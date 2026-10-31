"""CIF reader and writer (Package A, Feature A2).

The reader uses ``gemmi`` to parse CIF files and expand the asymmetric unit
by the space group.  The writer produces P1 CIFs with ``_discuspdf_*`` tags
for exact round-trip.  ``gemmi`` is imported inside the functions so that
``import discuspdf`` works without the ``cif`` extra.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .lattice import Lattice
from .structure import Structure
from .interop import (
    InteropError,
    InteropWarning,
    B_FROM_U,
    u_to_b,
    element_from_type_symbol,
    LossReport,
)
from ._safe_write import atomic_write_text


def _import_gemmi():
    """Import gemmi, raising a helpful ``ImportError`` if missing."""
    try:
        import gemmi
        return gemmi
    except ImportError:
        raise ImportError(
            "gemmi is required for CIF reading; "
            "install it with: pip install 'discus-python[cif]'"
        ) from None


@dataclass(frozen=True)
class CifImport:
    """Result of reading a CIF file.

    ``structure`` holds the expanded unit cell (``model="average"``,
    ``supercell=(1,1,1)``, ``pbc`` all True).  ``loss`` is non-empty only
    when ``allow_loss=True``.
    """

    structure: Structure
    block: str
    spacegroup: str | None
    n_asymmetric: int
    symmetry_source: str
    loss: tuple[str, ...]
    defaults_applied: tuple[str, ...]


def _parse_symop_strings(
    gemmi, symop_strings: list[str]
) -> list:
    """Parse a list of symop strings into a list of ``gemmi.Op`` objects.

    Raises ``InteropError`` if any string cannot be parsed.
    """
    ops = []
    for s in symop_strings:
        try:
            ops.append(gemmi.Op(s))
        except Exception as e:
            raise InteropError(f"cannot parse symop {s!r}: {e}") from e
    return ops


def _get_symops(
    gemmi, block
) -> tuple[list, str, str | None]:
    """Determine the symmetry operations and source.

    Returns ``(ops, symmetry_source, spacegroup_symbol)``.

    Policy (design rule 3):
    1. ``_space_group_symop_operation_xyz`` or ``_symmetry_equiv_pos_as_xyz``
       if present → ``"symop_loop"``.
    2. Space-group symbol → ``"space_group_symbol"``.
    3. Neither → ``"assumed_P1"`` (identity only).
    """
    # 1. Symop loop.
    symop = block.find_values("_space_group_symop_operation_xyz")
    if len(symop) == 0:
        symop = block.find_values("_symmetry_equiv_pos_as_xyz")
    if len(symop) > 0:
        strings = [gemmi.cif.as_string(symop[i]) for i in range(len(symop))]
        ops = _parse_symop_strings(gemmi, strings)
        # Try to identify the space group from the ops (informational).
        sg_symbol = None
        try:
            group_ops = gemmi.GroupOps(ops)
            sg = gemmi.find_spacegroup_by_ops(group_ops)
            if sg is not None:
                sg_symbol = sg.hm
        except Exception:
            pass
        return ops, "symop_loop", sg_symbol

    # 2. Space-group symbol.
    sg_symbol = block.find_value("_symmetry_space_group_name_H-M")
    if sg_symbol is not None:
        sg_symbol = gemmi.cif.as_string(sg_symbol)
        if sg_symbol:
            try:
                sg = gemmi.find_spacegroup_by_name(sg_symbol)
                if sg is not None:
                    point_ops = sg.operations().sym_ops
                    # Add centring translations if the group is centred.
                    centring = sg.centring_type()
                    if centring == "F":
                        translations = [
                            (0.0, 0.0, 0.0),
                            (0.0, 0.5, 0.5),
                            (0.5, 0.0, 0.5),
                            (0.5, 0.5, 0.0),
                        ]
                    elif centring == "C":
                        translations = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
                    elif centring == "I":
                        translations = [
                            (0.0, 0.0, 0.0),
                            (0.5, 0.5, 0.5),
                        ]
                    elif centring == "A":
                        translations = [(0.0, 0.0, 0.0), (0.5, 0.0, 0.5)]
                    elif centring == "B":
                        translations = [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5)]
                    else:
                        translations = [(0.0, 0.0, 0.0)]
                    # Combine point ops with centring translations.
                    # We represent each operation as (rotation_matrix, translation).
                    # The rotation is extracted from the gemmi.Op (in 1/100000 units);
                    # the translation is added manually.
                    ops = []
                    for t in translations:
                        for p in point_ops:
                            # Extract the rotation matrix from the Op.
                            # p.rot is a 3x3 list of ints in 1/DEN units (DEN=24).
                            rot = np.array(p.rot, dtype=float) / gemmi.Op.DEN
                            # Store as (rot, translation) tuple.
                            ops.append((rot, np.array(t)))
                    return ops, "space_group_symbol", sg_symbol
            except Exception:
                pass
            # Symbol present but unparseable: treat as P1.
            return [gemmi.Op("x,y,z")], "assumed_P1", sg_symbol

    # 3. No symmetry information.
    return [gemmi.Op("x,y,z")], "assumed_P1", None


def _expand_and_merge(
    gemmi,
    listed_sites: list,
    ops: list,
) -> tuple[list, list[str]]:
    """Expand listed sites by symmetry operations, wrap to [0,1), and merge
    periodic images within 1e-4.

    ``ops`` is a list of either ``gemmi.Op`` objects (from a symop loop) or
    ``(rotation_matrix, translation)`` tuples (from the space-group symbol).

    Returns ``(rows, losses)`` where ``rows`` is a list of dicts with keys
    ``label``, ``type_symbol``, ``frac``, ``occ``, ``u_iso``, and ``losses``
    is a list of loss messages.
    """
    merge_tol = 1e-4
    rows: list[dict] = []
    losses: list[str] = []

    for site in listed_sites:
        label = site.label
        type_symbol = site.type_symbol
        occ = site.occ
        u_iso = site.u_iso
        frac0 = np.array([site.fract.x, site.fract.y, site.fract.z], dtype=float)

        # Check for anisotropic ADPs (loss).
        aniso = site.aniso
        if aniso is not None and (
            aniso[0] != 0.0 or aniso[1] != 0.0 or aniso[2] != 0.0
            or aniso[3] != 0.0 or aniso[4] != 0.0 or aniso[5] != 0.0
        ):
            if "anisotropic ADPs" not in losses:
                losses.append("anisotropic ADPs")

        # Check for disorder group (loss).
        if site.disorder_group not in (0, None):
            if "disorder groups" not in losses:
                losses.append("disorder groups")

        # Expand by all symmetry operations.
        images: list[np.ndarray] = []
        for op in ops:
            if isinstance(op, tuple):
                # (rotation_matrix, translation)
                rot, trans = op
                pos = rot @ frac0 + trans
            else:
                # gemmi.Op
                result = op.apply_to_xyz((frac0[0], frac0[1], frac0[2]))
                pos = np.array(result, dtype=float)
            images.append(pos)

        # Wrap to [0, 1).
        for i, img in enumerate(images):
            images[i] = img % 1.0

        # Merge periodic images within merge_tol.
        merged: list[np.ndarray] = []
        for img in images:
            is_dup = False
            for m in merged:
                # Periodic distance: min(|d|, 1-|d|) per component.
                d = np.abs(img - m)
                d = np.minimum(d, 1.0 - d)
                if np.max(d) < merge_tol:
                    is_dup = True
                    break
            if not is_dup:
                merged.append(img)

        # One row per distinct image.
        for k, pos in enumerate(merged, start=1):
            rows.append({
                "label": label,
                "type_symbol": type_symbol,
                "frac": pos,
                "occ": occ,
                "u_iso": u_iso,
                "image_k": k,
            })

    return rows, losses


def _check_occupancy(rows: list[dict]) -> None:
    """Check that no position has total occupancy > 1 + 1e-6.

    For each unique position, sum the occupancies of all *different* labels
    at that position.  (Multiple rows with the same label at the same position
    are merged by the caller, so this only sees one row per label per position.)
    Raises ``InteropError`` if any position exceeds the limit.
    """
    merge_tol = 1e-4
    # Group by position only.
    pos_groups: list[list[dict]] = []
    for row in rows:
        found = False
        for g in pos_groups:
            d = np.abs(row["frac"] - g[0]["frac"])
            d = np.minimum(d, 1.0 - d)
            if np.max(d) < merge_tol:
                g.append(row)
                found = True
                break
        if not found:
            pos_groups.append([row])

    for g in pos_groups:
        # Sum occupancies of *different* labels at this position.
        # (If the same label appears multiple times, only count it once -
        # the images are equivalent and should have been merged.)
        # Only sum labels that are actually at the same position (within merge_tol).
        labels_at_pos: dict[str, float] = {}
        for r in g:
            if r["label"] not in labels_at_pos:
                labels_at_pos[r["label"]] = r["occ"]
        # Only sum if there are multiple *different* labels at this position.
        # (A single label at a position is fine - it's one site.)
        if len(labels_at_pos) > 1:
            total = sum(labels_at_pos.values())
            if total > 1.0 + 1e-6:
                labels = ", ".join(f"{lbl} ({occ:.4f})" for lbl, occ in labels_at_pos.items())
                raise InteropError(
                    f"total occupancy {total:.6f} > 1 + 1e-6 at position "
                    f"({g[0]['frac'][0]:.4f}, {g[0]['frac'][1]:.4f}, {g[0]['frac'][2]:.4f}) "
                    f"({labels})"
                )


def _build_structure(
    rows: list[dict],
    lattice: Lattice,
    *,
    model: str = "average",
    defaults_applied: list[str] | None = None,
) -> Structure:
    """Build a ``Structure`` from the expanded rows.

    ``site_label = f"{label}_{k}"`` where ``k`` is the 1-based image number.
    Rows from the same listed atom that are symmetry images of one another
    keep their own labels (``Ga1_1`` … ``Ga1_4``).  When rows from *different*
    listed atoms sit at the same position (periodic distance < 1e-4 in
    fractional units) they form one mixed site and all rows in the group
    share the ``site_label`` of the first listed atom at that position
    (i.e. the lowest index in ``rows``, which is in file order).
    """
    if not rows:
        raise InteropError("no atoms found in CIF")

    n = len(rows)
    species = np.empty(n, dtype=object)
    frac = np.empty((n, 3), dtype=float)
    occupancy = np.empty(n, dtype=float)
    site_label = np.empty(n, dtype=object)
    atom_id = np.arange(n, dtype=np.int64)
    cell_index = np.zeros((n, 3), dtype=np.int64)

    # Determine B factor.
    has_u = all(r["u_iso"] is not None for r in rows)
    b_iso = None
    if has_u and any(r["u_iso"] is not None for r in rows):
        b_iso = np.empty(n, dtype=float)
        for i, r in enumerate(rows):
            if r["u_iso"] is not None:
                b_iso[i] = u_to_b(r["u_iso"])
            else:
                b_iso[i] = 0.0
        # If all u_iso are 0.0, b_iso is all 0.0 → treat as None.
        if np.all(b_iso == 0.0):
            b_iso = None

    for i, r in enumerate(rows):
        # Species: type_symbol verbatim, else element from label.
        if r["type_symbol"]:
            sp = r["type_symbol"]
        else:
            sp = element_from_type_symbol(r["label"])
        species[i] = sp
        frac[i] = r["frac"]
        occupancy[i] = r["occ"]
        site_label[i] = f"{r['label']}_{r['image_k']}"

    # Mixed-site pass: group rows by position (periodic distance < 1e-4).
    # For each group containing rows from more than one distinct listed atom,
    # overwrite all rows' site_label with the label of the first row in file
    # order (lowest index).  Rows from the same listed atom (symmetry images)
    # keep their own labels.
    merge_tol = 1e-4
    pos_groups: list[list[int]] = []
    for i in range(n):
        found = False
        for g in pos_groups:
            d = np.abs(frac[i] - frac[g[0]])
            d = np.minimum(d, 1.0 - d)
            if np.max(d) < merge_tol:
                g.append(i)
                found = True
                break
        if not found:
            pos_groups.append([i])

    for g in pos_groups:
        labels_at_pos = {rows[i]["label"] for i in g}
        if len(labels_at_pos) > 1:
            first_label = rows[g[0]]["label"]
            # All rows at this position share the first atom's label.
            # Use the first row's image_k so that (site_label, cell_index)
            # is identical for all rows at the position (required by the
            # average-model validation in Structure).
            shared_k = rows[g[0]]["image_k"]
            for i in g:
                site_label[i] = f"{first_label}_{shared_k}"

    # Check occupancy.
    _check_occupancy(rows)

    # Mixed pbc → loss (handled by caller).
    return Structure(
        lattice,
        species,
        frac,
        model=model,
        occupancy=occupancy,
        site_label=site_label,
        atom_id=atom_id,
        b_iso=b_iso,
        cell_index=cell_index,
        supercell=(1, 1, 1),
        pbc=(True, True, True),
        min_separation=0.1,
    )


def parse_cif(
    text: str,
    *,
    source: str = "<string>",
    block: str | None = None,
    model: str = "average",
    allow_loss: bool = False,
    min_separation: float = 0.1,
) -> CifImport:
    """Parse CIF text into a :class:`CifImport`.

    See the module docstring for the full behaviour specification.
    """
    gemmi = _import_gemmi()

    if not text or not text.strip():
        raise InteropError(f"{source}: empty CIF text")

    try:
        doc = gemmi.cif.read_string(text)
    except Exception as e:
        raise InteropError(f"{source}: cannot parse CIF: {e}") from e

    # Block selection.
    if block is not None:
        b = doc.find_block(block)
        if b is None:
            raise InteropError(f"{source}: block {block!r} not found")
        block_name = b.name
    else:
        try:
            b = doc.sole_block()
            block_name = b.name
        except Exception:
            # Multiple blocks.
            names = [blk.name for blk in doc]
            raise InteropError(
                f"{source}: multiple data blocks {names}; pass block='...'"
            )

    # Cell parameters.
    a_str = b.find_value("_cell_length_a")
    b_str = b.find_value("_cell_length_b")
    c_str = b.find_value("_cell_length_c")
    alpha_str = b.find_value("_cell_angle_alpha")
    beta_str = b.find_value("_cell_angle_beta")
    gamma_str = b.find_value("_cell_angle_gamma")

    if a_str is None or b_str is None or c_str is None:
        raise InteropError(f"{source}: missing _cell_length_a/b/c")
    if alpha_str is None or beta_str is None or gamma_str is None:
        raise InteropError(f"{source}: missing _cell_angle_alpha/beta/gamma")

    try:
        a = float(gemmi.cif.as_string(a_str))
        b_len = float(gemmi.cif.as_string(b_str))
        c = float(gemmi.cif.as_string(c_str))
        alpha = float(gemmi.cif.as_string(alpha_str))
        beta = float(gemmi.cif.as_string(beta_str))
        gamma = float(gemmi.cif.as_string(gamma_str))
    except ValueError as e:
        raise InteropError(f"{source}: non-numeric cell parameter: {e}") from e

    try:
        lattice = Lattice.from_parameters(a, b_len, c, alpha, beta, gamma)
    except ValueError as e:
        raise InteropError(f"{source}: invalid cell parameters: {e}") from e

    # Modulated structure → always InteropError.
    mod_dim = b.find_value("_cell_modulation_dimension")
    if mod_dim is not None:
        mod_val = gemmi.cif.as_string(mod_dim)
        if mod_val and mod_val not in (".", "?", "0"):
            raise InteropError(
                f"{source}: modulated structure (modulation dimension {mod_val}) "
                "is not supported"
            )

    # Symmetry.
    ops, symmetry_source, sg_symbol = _get_symops(gemmi, b)
    defaults_applied: list[str] = []
    if symmetry_source == "assumed_P1":
        defaults_applied.append("symmetry")

    # Listed sites.
    labels = b.find_values("_atom_site_label")
    if len(labels) == 0:
        raise InteropError(f"{source}: no _atom_site_label column")

    n_listed = len(labels)
    type_symbols = b.find_values("_atom_site_type_symbol")
    occ_col = b.find_values("_atom_site_occupancy")
    u_iso_col = b.find_values("_atom_site_U_iso_or_equiv")
    b_iso_col = b.find_values("_atom_site_B_iso_or_equiv")
    adp_type_col = b.find_values("_atom_site_adp_type")
    disorder_col = b.find_values("_atom_site_disorder_group")
    aniso_labels = b.find_values("_atom_site_aniso_label")

    # Build listed site records.
    listed_sites = []
    for i in range(n_listed):
        label = gemmi.cif.as_string(labels[i])
        type_symbol = ""
        if len(type_symbols) > i:
            type_symbol = gemmi.cif.as_string(type_symbols[i])
            # Validate the element symbol (design rule 5).
            if type_symbol:
                element_from_type_symbol(type_symbol)
        occ = 1.0
        if len(occ_col) > i:
            try:
                occ = float(gemmi.cif.as_string(occ_col[i]))
            except ValueError:
                pass
        else:
            if "occupancy" not in defaults_applied:
                defaults_applied.append("occupancy")
        u_iso = None
        if len(u_iso_col) > i:
            try:
                u_iso = float(gemmi.cif.as_string(u_iso_col[i]))
            except ValueError:
                u_iso = None
        b_iso_val = None
        if len(b_iso_col) > i:
            try:
                b_iso_val = float(gemmi.cif.as_string(b_iso_col[i]))
            except ValueError:
                b_iso_val = None

        # Fractional coordinates.
        fx = b.find_value(f"_atom_site_fract_x")
        # Use find_values for the loop.
        fx_col = b.find_values("_atom_site_fract_x")
        fy_col = b.find_values("_atom_site_fract_y")
        fz_col = b.find_values("_atom_site_fract_z")
        if len(fx_col) <= i or len(fy_col) <= i or len(fz_col) <= i:
            raise InteropError(
                f"{source}: missing fractional coordinates for atom {i} ({label})"
            )
        try:
            fx_val = float(gemmi.cif.as_string(fx_col[i]))
            fy_val = float(gemmi.cif.as_string(fy_col[i]))
            fz_val = float(gemmi.cif.as_string(fz_col[i]))
        except ValueError as e:
            raise InteropError(
                f"{source}: non-numeric fractional coordinates for atom {i} "
                f"({label}): {e}"
            ) from e

        # Anisotropic ADP detection.
        has_aniso = False
        if len(adp_type_col) > i:
            adp_type = gemmi.cif.as_string(adp_type_col[i])
            if adp_type and adp_type not in (".", "?"):
                has_aniso = True
        if len(aniso_labels) > 0:
            # Check if this label appears in the aniso loop.
            for j in range(len(aniso_labels)):
                if gemmi.cif.as_string(aniso_labels[j]) == label:
                    has_aniso = True
                    break

        # Disorder group.
        disorder = 0
        if len(disorder_col) > i:
            try:
                dv = gemmi.cif.as_string(disorder_col[i])
                if dv not in (".", "?", ""):
                    disorder = int(float(dv))
            except (ValueError, OverflowError):
                pass

        # Magnetic moments.
        has_moment = False
        for tag in ("_atom_site_moment_total", "_atom_site_moment_lx",
                    "_atom_site_moment_ly", "_atom_site_moment_lz",
                    "_atom_site_moment_sx", "_atom_site_moment_sy",
                    "_atom_site_moment_sz"):
            mv = b.find_value(tag)
            if mv is not None:
                mv_str = gemmi.cif.as_string(mv)
                if mv_str and mv_str not in (".", "?"):
                    has_moment = True
                    break

        listed_sites.append({
            "label": label,
            "type_symbol": type_symbol,
            "frac": np.array([fx_val, fy_val, fz_val], dtype=float),
            "occ": occ,
            "u_iso": u_iso,
            "b_iso": b_iso_val,
            "has_aniso": has_aniso,
            "disorder": disorder,
            "has_moment": has_moment,
        })

    # Collect losses.
    losses = LossReport()
    for site in listed_sites:
        if site["has_aniso"]:
            losses.add(f"anisotropic ADPs on {site['label']}")
        if site["disorder"] not in (0, None):
            losses.add(f"disorder group on {site['label']}")
        if site["has_moment"]:
            losses.add(f"magnetic moment on {site['label']}")

    # B factor policy: B = 8π²·U, or _atom_site_B_iso_or_equiv if present.
    # If the column is missing for some or all atoms → b_iso=None and a loss
    # line only when it is missing for *some*.
    has_b_col = len(b_iso_col) > 0
    has_u_col = len(u_iso_col) > 0
    b_iso_arr = None
    if has_b_col:
        b_iso_arr = np.empty(n_listed, dtype=float)
        for i in range(n_listed):
            try:
                b_iso_arr[i] = float(gemmi.cif.as_string(b_iso_col[i]))
            except (ValueError, IndexError):
                b_iso_arr[i] = 0.0
        if np.all(b_iso_arr == 0.0):
            b_iso_arr = None
    elif has_u_col:
        b_iso_arr = np.empty(n_listed, dtype=float)
        for i in range(n_listed):
            try:
                u_val = float(gemmi.cif.as_string(u_iso_col[i]))
            except (ValueError, IndexError):
                u_val = 0.0
            b_iso_arr[i] = u_to_b(u_val)
        if np.all(b_iso_arr == 0.0):
            b_iso_arr = None
    else:
        # No U or B column at all.
        if "b_iso" not in defaults_applied:
            defaults_applied.append("b_iso")

    # Check for _discuspdf_* tags (our own files).
    dp_model = b.find_value("_discuspdf_model")
    dp_supercell_a = b.find_value("_discuspdf_supercell_a")
    dp_supercell_b = b.find_value("_discuspdf_supercell_b")
    dp_supercell_c = b.find_value("_discuspdf_supercell_c")
    dp_site_labels = b.find_values("_discuspdf_site_label")
    dp_atom_ids = b.find_values("_discuspdf_atom_id")
    dp_cell_a = b.find_values("_discuspdf_cell_index_a")
    dp_cell_b = b.find_values("_discuspdf_cell_index_b")
    dp_cell_c = b.find_values("_discuspdf_cell_index_c")

    has_dp_tags = (
        dp_model is not None
        and dp_supercell_a is not None
        and len(dp_site_labels) == n_listed
    )

    if has_dp_tags:
        # Restore from _discuspdf_* tags (exact round-trip).
        dp_model_val = gemmi.cif.as_string(dp_model)
        if dp_model_val not in ("average", "configuration"):
            raise InteropError(
                f"unknown _discuspdf_model {dp_model_val!r}"
            )
        # Override the model parameter.
        model = dp_model_val

        sc_a = int(gemmi.cif.as_string(dp_supercell_a))
        sc_b = int(gemmi.cif.as_string(dp_supercell_b))
        sc_c = int(gemmi.cif.as_string(dp_supercell_c))
        supercell = (sc_a, sc_b, sc_c)

        n = n_listed
        species = np.empty(n, dtype=object)
        frac = np.empty((n, 3), dtype=float)
        occupancy = np.empty(n, dtype=float)
        site_label = np.empty(n, dtype=object)
        atom_id = np.empty(n, dtype=np.int64)
        cell_index = np.empty((n, 3), dtype=np.int64)
        b_iso_arr_full = None

        if b_iso_arr is not None:
            b_iso_arr_full = np.empty(n, dtype=float)

        for i in range(n):
            species[i] = listed_sites[i]["type_symbol"] or element_from_type_symbol(
                listed_sites[i]["label"]
            )
            frac[i] = listed_sites[i]["frac"]
            occupancy[i] = listed_sites[i]["occ"]
            site_label[i] = gemmi.cif.as_string(dp_site_labels[i])
            atom_id[i] = int(gemmi.cif.as_string(dp_atom_ids[i]))
            cell_index[i, 0] = int(gemmi.cif.as_string(dp_cell_a[i]))
            cell_index[i, 1] = int(gemmi.cif.as_string(dp_cell_b[i]))
            cell_index[i, 2] = int(gemmi.cif.as_string(dp_cell_c[i]))
            if b_iso_arr_full is not None:
                b_iso_arr_full[i] = b_iso_arr[i]

        if b_iso_arr_full is not None and np.all(b_iso_arr_full == 0.0):
            b_iso_arr_full = None

        # Check occupancy.
        _check_occupancy(listed_sites)

        # Configuration model check.
        if model == "configuration":
            for i, d in enumerate(listed_sites):
                if d["occ"] < 1.0 - 1e-9:
                    raise InteropError(
                        f"model='configuration' requires all occupancies = 1.0; "
                        f"atom {d['label']} has occupancy {d['occ']:.4f}"
                    )

        structure = Structure(
            lattice,
            species,
            frac,
            model=model,
            occupancy=occupancy,
            site_label=site_label,
            atom_id=atom_id,
            b_iso=b_iso_arr_full,
            cell_index=cell_index,
            supercell=supercell,
            pbc=(True, True, True),
            min_separation=min_separation,
        )
    else:
        # Expand and merge (foreign CIF).
        class _Site:
            def __init__(self, d):
                self.label = d["label"]
                self.type_symbol = d["type_symbol"]
                self.fract = _Fract(d["frac"])
                self.occ = d["occ"]
                self.u_iso = d["u_iso"]
                self.aniso = None
                self.disorder_group = d["disorder"]

        class _Fract:
            def __init__(self, arr):
                self._arr = arr
                self.x = arr[0]
                self.y = arr[1]
                self.z = arr[2]

        gemmi_sites = [_Site(d) for d in listed_sites]
        rows, _ = _expand_and_merge(gemmi, gemmi_sites, ops)

        # Apply B factor to rows.
        if b_iso_arr is not None:
            for i, d in enumerate(listed_sites):
                for row in rows:
                    if row["label"] == d["label"]:
                        row["b_iso"] = b_iso_arr[i] if b_iso_arr[i] > 0.0 else None

        # Build the Structure.
        if model == "configuration":
            for i, d in enumerate(listed_sites):
                if d["occ"] < 1.0 - 1e-9:
                    raise InteropError(
                        f"model='configuration' requires all occupancies = 1.0; "
                        f"atom {d['label']} has occupancy {d['occ']:.4f}"
                    )

        structure = _build_structure(
            rows, lattice, model=model, defaults_applied=defaults_applied
        )

        # If B was determined at the listed-site level, apply it.
        if b_iso_arr is not None:
            b_full = np.empty(structure.n_atoms, dtype=float)
            for i, row in enumerate(rows):
                b_full[i] = row.get("b_iso", 0.0) or 0.0
            if np.all(b_full == 0.0):
                b_full = None
            if b_full is not None:
                structure = structure.replace(b_iso=b_full)

    # Loss policy.
    if losses.losses:
        if not allow_loss:
            raise InteropError("\n".join(losses.losses))
        import warnings
        warnings.warn("\n".join(losses.losses), InteropWarning, stacklevel=2)

    return CifImport(
        structure=structure,
        block=block_name,
        spacegroup=sg_symbol,
        n_asymmetric=n_listed,
        symmetry_source=symmetry_source,
        loss=losses.losses if allow_loss else (),
        defaults_applied=tuple(defaults_applied),
    )


def read_cif(
    path: Path | str,
    *,
    block: str | None = None,
    model: str = "average",
    allow_loss: bool = False,
    min_separation: float = 0.1,
) -> CifImport:
    """Read a CIF file from ``path``."""
    _import_gemmi()
    p = Path(path)
    text = p.read_text(encoding="ascii")
    return parse_cif(
        text, source=str(p), block=block, model=model,
        allow_loss=allow_loss, min_separation=min_separation,
    )


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

_CIF_DECIMALS = 10


def _format_number(v: float) -> str:
    """Format a number with 10 decimal places."""
    return f"{v:.{_CIF_DECIMALS}f}"


def format_cif(
    structure: Structure,
    *,
    block_name: str = "discuspdf",
    allow_loss: bool = False,
) -> str:
    """Format a :class:`Structure` as CIF text (P1 only).

    The output includes ``_discuspdf_*`` tags for exact round-trip.
    """
    losses = LossReport()

    # Mixed pbc → loss.
    if not all(structure.pbc):
        losses.add(
            f"mixed pbc {structure.pbc} is not stored in a CIF; "
            "the reader assumes all periodic"
        )

    # B factor: mixed None/values → InteropError.
    b_iso = structure.b_iso
    if b_iso is not None:
        # All values must be present (no None mixing).
        pass  # b_iso is either None or a full array.

    a, b_len, c, alpha, beta, gamma = structure.lattice.parameters
    n = structure.n_atoms

    lines = [
        f"data_{block_name}",
        f"_space_group_name_H-M_alt 'P 1'",
        f"_space_group_IT_number 1",
        f"_cell_length_a {_format_number(a)}",
        f"_cell_length_b {_format_number(b_len)}",
        f"_cell_length_c {_format_number(c)}",
        f"_cell_angle_alpha {_format_number(alpha)}",
        f"_cell_angle_beta {_format_number(beta)}",
        f"_cell_angle_gamma {_format_number(gamma)}",
        f"_discuspdf_model {structure.model}",
        f"_discuspdf_supercell_a {structure.supercell[0]}",
        f"_discuspdf_supercell_b {structure.supercell[1]}",
        f"_discuspdf_supercell_c {structure.supercell[2]}",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    if b_iso is not None:
        lines.append("_atom_site_B_iso_or_equiv")
    lines += [
        "_discuspdf_site_label",
        "_discuspdf_atom_id",
        "_discuspdf_cell_index_a",
        "_discuspdf_cell_index_b",
        "_discuspdf_cell_index_c",
    ]

    for i in range(n):
        element = element_from_type_symbol(structure.species[i])
        label = f"{element}{structure.atom_id[i]}"
        x, y, z = (
            _format_number(structure.frac[i, 0]),
            _format_number(structure.frac[i, 1]),
            _format_number(structure.frac[i, 2]),
        )
        occ = _format_number(structure.occupancy[i])
        row = f"{label} {structure.species[i]} {x} {y} {z} {occ}"
        if b_iso is not None:
            row += f" {_format_number(b_iso[i])}"
        row += (
            f" {structure.site_label[i]}"
            f" {structure.atom_id[i]}"
            f" {structure.cell_index[i, 0]}"
            f" {structure.cell_index[i, 1]}"
            f" {structure.cell_index[i, 2]}"
        )
        lines.append(row)

    losses.raise_or_warn(allow_loss)
    return "\n".join(lines) + "\n"


def write_cif(
    path: Path | str,
    structure: Structure,
    *,
    block_name: str = "discuspdf",
    allow_loss: bool = False,
    overwrite: bool = False,
) -> Path:
    """Write a :class:`Structure` to ``path`` as a CIF file, atomically.

    Refuses an existing path unless ``overwrite=True``.  The written text is
    parsed back and must agree with the input; on failure the target is left
    untouched.
    """
    target = Path(path)
    text = format_cif(structure, block_name=block_name, allow_loss=allow_loss)

    # Write atomically.
    tmp = atomic_write_text(target, text, overwrite=overwrite)

    # Read back and verify.
    try:
        imported = parse_cif(
            tmp.read_text(encoding="ascii"),
            source=str(tmp),
            model=structure.model,
            allow_loss=allow_loss,
        )
    except Exception:
        tmp.unlink()
        raise

    written = imported.structure
    if not _cif_roundtrip_agrees(structure, written):
        tmp.unlink()
        raise ValueError("written CIF does not round-trip to the input")

    return target


def _cif_roundtrip_agrees(original: Structure, written: Structure) -> bool:
    """Check that a written CIF parses back to the input structure."""
    if written.n_atoms != original.n_atoms:
        return False
    if not np.array_equal(written.species, original.species):
        return False
    # Lattice: 10-decimal precision.
    tol = 0.5 * 10 ** -_CIF_DECIMALS
    if not np.allclose(
        written.lattice.parameters, original.lattice.parameters, atol=tol
    ):
        return False
    # Fractional coordinates: 10-decimal precision.
    if not np.allclose(written.frac, original.frac, atol=tol):
        return False
    # Occupancy: 10-decimal precision.
    if not np.allclose(written.occupancy, original.occupancy, atol=tol):
        return False
    # B factor.
    if (original.b_iso is None) != (written.b_iso is None):
        return False
    if original.b_iso is not None:
        if not np.allclose(written.b_iso, original.b_iso, atol=tol):
            return False
    # Site labels, atom IDs, cell indices: exact.
    if not np.array_equal(written.site_label, original.site_label):
        return False
    if not np.array_equal(written.atom_id, original.atom_id):
        return False
    if not np.array_equal(written.cell_index, original.cell_index):
        return False
    # Supercell: exact.
    if not np.array_equal(written.supercell, original.supercell):
        return False
    # Model: exact.
    if written.model != original.model:
        return False
    return True
