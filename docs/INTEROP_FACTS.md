# INTEROP_FACTS.md — A2-1/A2-1b probe findings

Environment: Python 3.14.7, gemmi 0.7.5, ase 3.29.0, diffpy.structure 3.5.0.
All evidence files are in `tests/fixtures/interop/evidence/`.

## P1: gemmi small-structure API on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `gemmi.read_small_structure` reads the CIF; block name `9008845`; space group `F -4 3 m` (No. 216) | `P1.log` | `spacegroup_hm: F -4 3 m` / `spacegroup_number: 216` |
| 2 listed sites: Ga at (0,0,0), As at (0.25,0.25,0.25); both occ=1.0, u_iso=0.0, aniso=zero, disorder=0 | `P1.log` | `n_listed_sites: 2` / `listed[0]: label='Ga' … fract=<gemmi.Fractional(0, 0, 0)> occ=1.0` / `listed[1]: label='As' … fract=<gemmi.Fractional(0.25, 0.25, 0.25)> occ=1.0` |
| `get_all_unit_cell_sites()` returns 8 sites (2 × 4 F-centring images) | `P1.log` | `n_unit_cell_sites: 8` |
| The 8 unit-cell sites: Ga at (0,0,0), (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0); As at (0.25,0.25,0.25), (0.25,-0.25,-0.25), (-0.25,-0.25,0.25), (-0.25,0.25,-0.25) | `P1.log` | `uc[0]…uc[7]` (see log) |
| 4 of the 8 sites have coordinates outside [0,1) (negative values); gemmi does NOT wrap to [0,1) | `P1.log` | `uc[5]: … x=0.25 y=-0.25 z=-0.25 in_[0,1)=False` |
| `type_symbol` is present and equals the element (Ga, As) for both listed sites | `P1.log` | `type_symbol='Ga'` / `type_symbol='As'` |
| `_atom_site_label` values verbatim: `Ga`, `As` (no numeric suffix in the CIF) | `P1.log` | `label='Ga'` / `label='As'` |

## P2: gemmi small-structure API on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: P -1, 2 listed sites → 4 unit-cell sites; 2 sites outside [0,1) | `P2.log` | `n_listed_sites: 2` / `n_unit_cell_sites: 4` / `uc[1]: … in_[0,1)=False` |
| `ga07in03as_mixed.cif`: F -4 3 m, 3 listed sites (Ga1 occ=0.7, In1 occ=0.3, As1 occ=1.0) → 12 unit-cell sites; 4 outside [0,1) | `P2.log` | `n_listed_sites: 3` / `n_unit_cell_sites: 12` / `uc[9]: … in_[0,1)=False` |
| `p1_aniso_adp.cif`: P 1, 2 listed sites → 2 unit-cell sites; aniso tensor populated (Fe1: U11=0.013, U22=0.012, U33=0.011, U12=0.0001); u_iso=0.012/0.010 | `P2.log` | `aniso=<gemmi.SMat33d(0.013, 0.012, 0.011, 0.0001, 0, 0)>` / `u_iso=0.012` |
| `p1_disorder_groups.cif`: P 1, 3 listed sites → 3 unit-cell sites; disorder_group=1, 2, 0 | `P2.log` | `disorder=1` / `disorder=2` / `disorder=0` |
| `p1_occupancy_over_one.cif`: P 1, 2 listed sites at same position (Ga1 occ=0.7, In1 occ=0.5) → 2 unit-cell sites; gemmi does NOT reject total occupancy > 1 | `P2.log` | `uc[0]: … occ=0.7` / `uc[1]: … occ=0.5` |
| `spacegroup_number` is 0 for all 5 synthetic CIFs (gemmi does not resolve the number from the H-M symbol in `read_small_structure`) | `P2.log` | `spacegroup_number: 0` (all 5) |

## P3: ASE on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| ASE 3.29.0 expands the 2-atom asymmetric unit to 8 atoms | `P3.log` | `n_atoms: 8` / `species: ['Ga'×4, 'As'×4]` |
| Cell matrix matches the CIF: a=b=c=5.6537 | `P3.log` | `cell.array: [[5.6537 0 0] [0 5.6537 0] [0 0 5.6537]]` |
| `info` contains `spacegroup` (Spacegroup(216, setting=1)) and `unit_cell` ('conventional') but NO `occupancy` key (all occupancies are 1.0) | `P3.log` | `info: {'spacegroup': Spacegroup(216, setting=1), 'unit_cell': 'conventional'}` |
| `arrays` keys: `numbers`, `positions`, `spacegroup_kinds` (no U/B data) | `P3.log` | `arrays keys: ['numbers', 'positions', 'spacegroup_kinds']` |
| `get_scaled_positions()` returns fractional coordinates (unwrapped to [0,1)) | `P3.log` | `scaled positions: [[0 0 0] [0 0.5 0.5] …]` |

## P4: ASE on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: 4 atoms; `info['occupancy']` = `{'0': {'Cu': 1.0}, '1': {'O': 1.0}}` (keys are listed-site indices) | `P4.log` | `info: … 'occupancy': {'0': {'Cu': 1.0}, '1': {'O': 1.0}}` |
| `ga07in03as_mixed.cif`: 8 atoms (Ga×4, As×4; **In dropped from `species`**); `info['occupancy']` = `{'0': {'Ga': 0.7, 'In': 0.3}, '1': {'In': 0.3, 'Ga': 0.7}, '2': {'As': 1.0}}` — the In data IS in `info`, not silent | `P4.log` | `species: ['Ga'×4, 'As'×4]` / `info: … 'occupancy': {'0': {'Ga': 0.7, 'In': 0.3}, …}` |
| `p1_aniso_adp.cif`: 2 atoms; `info['occupancy']` = `{'0': {'Fe': 1.0}, '1': {'Fe': 1.0}}`; U/B not in `arrays` | `P4.log` | `info: … 'occupancy': {'0': {'Fe': 1.0}, '1': {'Fe': 1.0}}` |
| `p1_disorder_groups.cif`: 3 atoms; `info['occupancy']` = `{'0': {'Ca': 0.5}, '1': {'Ca': 0.5}, '2': {'O': 1.0}}` | `P4.log` | `info: … 'occupancy': {'0': {'Ca': 0.5}, '1': {'Ca': 0.5}, '2': {'O': 1.0}}` |
| `p1_occupancy_over_one.cif`: **1 atom only** (Ga); `info['occupancy']` = `{'0': {'Ga': 0.7, 'In': 0.5}, '1': {'In': 0.5, 'Ga': 0.7}}` — the In data IS in `info` | `P4.log` | `n_atoms: 1` / `species: ['Ga']` / `info: … 'occupancy': {'0': {'Ga': 0.7, 'In': 0.5}, …}` |
| **Atom-index ↔ key mapping:** the keys of `info['occupancy']` are indices of the *listed* CIF sites (0, 1, 2), NOT indices of the expanded atoms. For the mixed fixture, listed site 0 (Ga1) and listed site 1 (In1) both map to atom indices 0–3 (the 4 F-centring images); listed site 2 (As1) maps to atom indices 4–7. The key `'0'` covers atoms 0–3, key `'1'` covers atoms 0–3 (same positions), key `'2'` covers atoms 4–7. **The mapping from atom index to key is NOT recoverable from `Atoms` alone** — there is no array that maps atom index → listed-site index. | `P4.log` | (inferred from `species` + `scaled positions` + `info['occupancy']` keys) |

## P5: diffpy.structure on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| diffpy.structure 3.5.0 expands to 8 atoms; cell matches; space group No. 216 (F-43m), 96 sym matrices, 24 point-group matrices | `P5.log` | `n_atoms: 8` / `spacegroup: SpaceGroup #216 (F-43m, Cubic). Symmetry matrices: 96, point sym. matr.: 24` |
| `cif_sgname` is `F -4 2 3` (diffpy normalises the H-M symbol; this is just DiffPy's own string for group 216) | `P5.log` | `cif_sgname: F -4 2 3` |
| All 8 atoms: occ=1.0, U=0.0, B=0.0 (defaults) | `P5.log` | `atom 0: Ga … occ=1.0 U=0.0 B=0.0` |

## P6: diffpy.structure on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: 4 atoms; U=0.0, B=0.0 (defaults) | `P6.log` | `n_atoms: 4` / `atom 0: Cu … U=0.0 B=0.0` |
| `ga07in03as_mixed.cif`: 12 atoms (3×4); U and B populated from Uiso (B=8π²U) | `P6.log` | `n_atoms: 12` / `atom 0: Ga … U=0.01 B=0.7895683520871487` |
| `p1_aniso_adp.cif`: 2 atoms; diffpy computes isotropic-equivalent U from anisotropic tensor (U=(U11+U22+U33)/3); B from that | `P6.log` | `atom 0: Fe … U=0.012898155079923389 B=1.0183975051419596` |
| `p1_disorder_groups.cif`: 3 atoms; U=0.0, B=0.0 | `P6.log` | `n_atoms: 3` |
| `p1_occupancy_over_one.cif`: 2 atoms (Ga occ=0.7, In occ=0.5); diffpy does NOT reject total occupancy > 1 | `P6.log` | `n_atoms: 2` / `atom 0: Ga … occ=0.7` / `atom 1: In … occ=0.5` |

## P7: cross-library comparison on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| All four libraries agree on cell parameters to machine precision (diff=0.0000000000) | `P7.log` | `gemmi vs ASE cell: a diff=0.0000000000 …` |
| gemmi reads 2 atoms (asymmetric unit); ASE and diffpy expand to 8 | `P7.log` | `n_atoms: gemmi=2 ase=8 diffpy=8` |
| discuspdf `Lattice.from_parameters` round-trips exactly | `P7.log` | `discuspdf Lattice: a=5.6537 b=5.6537 c=5.6537 al=90.0 be=90.0 ga=90.0` |

## P8: B/U handling across libraries

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| gemmi reads anisotropic U11/U22/U33 and adp_type directly; does not compute isotropic equivalent | `P8.log` | `Fe1: adp=Uani U11=0.0130 U22=0.0120 U33=0.0110` |
| ASE does not store U or B in the `Atoms` object after CIF reading | `P8.log` | `ASE arrays: ['numbers', 'positions', 'spacegroup_kinds']` |
| diffpy computes U_iso=(U11+U22+U33)/3 and B=8π²U; defaults to U=0.0, B=0.0 when columns absent | `P8.log` | `Fe: U=0.012898… B=1.018398…` / `Ga: U=0.0 B=0.0` |

## P9: gemmi CIF writing with _discuspdf_* tags

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| A CIF string with `_discuspdf_model`, `_discuspdf_site_label`, `_discuspdf_atom_id` columns is written and re-read by gemmi with all columns and data items intact | `P9.log` | `model: average` / `dp_site_labels: ['Ga1_1', 'As1_1']` / `dp_atom_ids: ['0', '1']` |
| Labels like `Ga1_1` are NOT quoted in the written CIF and are read back as plain strings by gemmi | `P9.log` | `labels: ['Ga1_1', 'As1_1']` |
| `find_value` returns quoted strings with literal quotes (e.g. `"'P 1'"`); `gemmi.cif.as_string` DOES strip the quotes (see P13) | `P9.log` / `P13.log` | `sg raw: "'P 1'"` / `as_string: 'P 1'` |
| ASE reads the same CIF without error; `_discuspdf_*` columns are ignored (not in `arrays`); `info['occupancy']` is produced | `P9.log` | `n_atoms: 2` / `species: ['Ga', 'As']` / `arrays keys: ['numbers', 'positions', 'spacegroup_kinds']` |

## P10: ASE deep inspection

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `info['occupancy']` format: `{str(listed_site_index): {symbol: occupancy}}`. For the mixed fixture: `{'0': {'Ga': 0.7, 'In': 0.3}, '1': {'In': 0.3, 'Ga': 0.7}, '2': {'As': 1.0}}` | `P10.log` | `info: … 'occupancy': {'0': {'Ga': 0.7, 'In': 0.3}, '1': {'In': 0.3, 'Ga': 0.7}, '2': {'As': 1.0}}` |
| **Atom-index ↔ key mapping is NOT recoverable from `Atoms` alone.** The keys are listed-site indices, but `Atoms` has no array mapping atom index → listed-site index. For the mixed fixture, atoms 0–3 are Ga (from listed sites 0 AND 1), atoms 4–7 are As (from listed site 2). Both key `'0'` and key `'1'` refer to the same 4 atom positions. | `P10.log` | (inferred: `species` + `scaled positions` + `info['occupancy']` keys) |
| `atoms.cell.array` vs `Lattice.from_parameters(*params).matrix`: max abs diff = 0.0 for triclinic, 3.55e-16 for cubic (GaAs) | `P10.log` | `max abs diff (ASE cell.array vs Lattice.matrix): 0.00e+00` (triclinic) / `3.55e-16` (cubic) |
| `atoms.copy()` and `pickle` round-trip: custom arrays (`dp_test`) and `info` entries (`dp_model`) survive both | `P10.log` | `copy: dp_test in arrays: True, dp_model in info: True` / `pickle: dp_test in arrays: True, dp_model in info: True` |
| `Atoms(symbols, cell=<left-handed matrix>)`: constructor does NOT raise; `get_scaled_positions()` returns `[0, 0, 0]` (no error, no warning) | `P10.log` | `left-handed: constructor OK, det=-125.0` / `get_scaled_positions(): [[0. 0. 0.]]` |
| `Atoms(symbols, cell=zeros)`: constructor does NOT raise; `get_scaled_positions()` returns `[0, 0, 0]` (no error) | `P10.log` | `zero cell: constructor OK, det=0.0` / `get_scaled_positions(): [[0. 0. 0.]]` |

## P11: DiffPy deep inspection

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `Lattice(4.1,5.3,6.2,78,84,101).base` is in a DIFFERENT orientation from `Lattice.from_parameters(...).matrix`; max abs diff = 1.44 (crystallographic convention difference) | `P11.log` | `max abs diff (Lattice.base vs our matrix): 1.44e+00` |
| `Lattice.cell_parms` is the replacement for the deprecated `abcABG` (returns `(4.1, 5.3, 6.2, 78.0, 84.0, 101.0)`) | `P11.log` | `cell_parms: (4.1, 5.3, 6.2, 78.0, 84.0, 101.0)` |
| `Atom(atype="O2-")`: `.element` = `"O2-"` (charges preserved); there is no `.atype` attribute | `P11.log` | `element='O2-'` |
| `Atom(...).Uisoequiv` default is 0.0 | `P11.log` | `Uisoequiv default: 0.0` |
| `Uisoequiv` with `anisotropy=False` and U11/U22/U33 set: the last assignment wins (isotropic tensor overwrite), NOT a meaningful equivalent — see P14 for the correct anisotropic case | `P11.log` | `Uisoequiv with U11=0.013 U22=0.012 U33=0.011: 0.011` (anisotropy was False) |
| Iteration order preserves the order atoms were passed to the `Structure` constructor; `.xyz` is a `numpy.ndarray` | `P11.log` | `iteration order: ['O2-', 'Ga']` / `xyz type: <class 'numpy.ndarray'>` |

## P12: B/U conversion

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `gemmi` does NOT have a `u_to_b` function (`hasattr(gemmi, 'u_to_b')` is `False`) | `P12.log` | `gemmi has u_to_b: False` |
| The conversion is `B = 8π²·U`: `8π²·0.01 = 0.7895683520871487`, `8π²·0.008 = 0.6316546816697189` | `P12.log` | `8*pi^2*0.01 = 0.7895683520871487` / `8*pi^2*0.008 = 0.6316546816697189` |
| The only B-related function in gemmi is `binrfactor` (unrelated to U→B conversion) | `P12.log` | `gemmi B-related functions: ['binrfactor']` |

## P13: gemmi.cif.as_string quote-stripping (corrected)

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `gemmi.cif.as_string` DOES strip CIF single-quotes: `as_string(b.find_value('_symmetry_space_group_name_H-M'))` → `'F -4 3 m'` from `"'F -4 3 m'"` | `P13.log` | `find_value raw: "'F -4 3 m'"` / `as_string: 'F -4 3 m'` |
| Works on quoted values from `find_values` (loop items) too | `P13.log` | `find_values[0] raw: 'Ga'` / `as_string: 'Ga'` (identity on unquoted) |
| `as_string` is the identity on unquoted values; strips quotes on quoted values. Never use ad-hoc `.strip("'")`. | `P13.log` | `Summary: as_string strips CIF single-quotes; identity on unquoted values.` |

## P14: DiffPy Uisoequiv with anisotropy (corrected)

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| With `anisotropy=True` and U11=0.013, U22=0.012, U33=0.011 on a **cubic** lattice: `Uisoequiv` = 0.012 (the correct isotropic equivalent) | `P14.log` | `cubic, anisotropy=True, …: Uisoequiv=0.012` |
| Same tensor on a **triclinic** lattice (4.1/5.3/6.2/78/84/101): `Uisoequiv` = 0.012911887014119466 | `P14.log` | `triclinic, anisotropy=True, …: Uisoequiv=0.012911887014119466` |
| Isotropic case: `anisotropy=False`, `Uisoequiv=0.012` set directly → `anisotropy` stays False, and the diagonal U11=U22=U33=0.012 | `P14.log` | `isotropic, anisotropy=False, Uisoequiv=0.012, U11=0.012 U22=0.012 U33=0.012` |
| The earlier P11 finding ("Uisoequiv = U33 = 0.011") was a probe artefact: setting U11/U22/U33 on an atom whose `anisotropy` is False just overwrites an isotropic tensor, so the last assignment wins. With `anisotropy=True`, `Uisoequiv` is the correct lattice-weighted isotropic equivalent. | `P14.log` | (see the three numbers above) |

## Overseer rulings (verbatim)

Copied from `00_START_HERE.md`, "Amendments after A2-1b (overseer, 2026-10-31)":

- **ASE foreign atoms:** `info['occupancy']` is never interpreted (see rule 4 above).
- **gemmi** `get_all_unit_cell_sites()` returns UNWRAPPED fractions (e.g. -0.25); the reader wraps to [0,1) and merges periodic images within 1e-4 (as already specified).
- **Quotes:** `gemmi.cif.as_string(value)` removes CIF quotes (checked by the overseer: `as_string("'F -4 3 m'") == 'F -4 3 m'`, also for values from `find_values`). Never use ad-hoc `.strip("'")`.
- **DiffPy:** `from_diffpy` takes the lattice from `lattice.cell_parms()` → `Lattice.from_parameters`, never from `Lattice.base`; fractional `xyz` are used directly. Anisotropic atoms: `atom.anisotropy` is True; `Uisoequiv` is then correct (overseer: diag(0.013,0.012,0.011) → 0.012 cubic, 0.012912 triclinic) and is used for `b_iso`; the anisotropy is a loss. `to_diffpy` leaves `anisotropy` False and sets `Uisoequiv`.
- **B ↔ U:** gemmi has no `u_to_b`; `interop.u_to_b` is our own `8π²·U`.

## Conflicts with the fixed design

- **None.**
