"""DISCUS runs for a disorder-builder model, through discuspdf's own configuration and runner classes.

* 3D: :class:`~discuspdf.FourierConfig` (grid centred on 0, Bragg subtraction ``set aver``, LOTS, explicit
  8-value seed, radiation, wavelength, own neutron scattering lengths) staged with
  :func:`~discuspdf.prepare_fourier_run` and run by :class:`~discuspdf.Diffuse3DRunner`.
* 1D: :class:`~discuspdf.PDFConfig` with the supercell ``.stru`` as structure (``read stru``) run by
  :class:`~discuspdf.PDFRunner` (hash checks, r-grid check, ``status.json``, ``result.csv``).
* DISCUS Monte Carlo (atom-only models): :class:`~discuspdf.mc.MCConfig` from the page's ``mmc`` block.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

from .model import VOID

__all__ = ["seed8_from", "build_mc_config", "fourier_wavelength", "scattering_lengths_for", "fourier_config",
           "run_fourier", "pdf_config", "run_powder_pdf"]


def seed8_from(rng) -> tuple:
    return tuple(int(x) for x in rng.integers(1, 99999, size=8))


def build_mc_config(mc: dict, seed8):
    """The page's ``mmc`` block as an :class:`~discuspdf.mc.MCConfig`."""
    from ..fourier_macro import Seed
    from ..mc import (MCConfig, Neighborhood, Vector, OccupationalTarget, DisplacementTarget,
                      SpringTarget, LennardJonesTarget, Swap, SwapDisplacement, Shift)
    nbs = []
    for nbd in mc["neighborhoods"]:
        vecs = tuple(Vector(int(v["origin"]), int(v["target"]), tuple(int(x) for x in v["t"])) for v in nbd["vectors"])
        nbs.append(Neighborhood(vecs, symmetric=bool(nbd.get("symmetric", True))))
    targets = []
    for t in mc["targets"]:
        kind = t["type"]
        common = dict(neighborhood=int(t["neighborhood"]), a=t["a"], b=t["b"])
        val = {"correlation": float(t["value"])} if t.get("mode", "correlation") == "correlation" \
            else {"energy": float(t["value"])}
        if kind == "occupational":
            targets.append(OccupationalTarget(**common, **val))
        elif kind == "displacement":
            targets.append(DisplacementTarget(**common, direction=tuple(float(x) for x in t["direction"]), **val))
        elif kind == "spring":
            targets.append(SpringTarget(**common, distance=float(t["distance"]), k=float(t["k"])))
        elif kind == "lennard_jones":
            targets.append(LennardJonesTarget(**common, distance=float(t["distance"]), depth=float(t["depth"])))
        else:
            raise ValueError(f"unknown MC target {kind!r}")
    moves = []
    mv = mc.get("moves", {})
    if mv.get("swap", 0) > 0:
        moves.append(Swap(probability=float(mv["swap"]), pairing=mv.get("pairing", "all")))
    if mv.get("swap_displacement", 0) > 0:
        moves.append(SwapDisplacement(probability=float(mv["swap_displacement"]), pairing=mv.get("pairing", "all")))
    if mv.get("shift", 0) > 0:
        moves.append(Shift(sigma={k: float(v) for k, v in mv.get("shift_sigma", {}).items()},
                           probability=float(mv["shift"])))
    return MCConfig(tuple(nbs), tuple(targets), tuple(moves), Seed(tuple(seed8)),
                    moves_per_atom=float(mc.get("moves_per_atom", 10.0)),
                    feedback_per_atom=float(mc.get("feedback_per_atom", 1.0)),
                    temperature=float(mc.get("temperature", 1.0)))


def fourier_wavelength(w):
    """The page's λ field: a DISCUS symbol ('moa1', 'cua1', …) or a number in Å; empty means 'moa1'."""
    s = str(w if w is not None else "").strip()
    if not s:
        return "moa1"
    try:
        return float(s)
    except ValueError:
        return s.lower()


def scattering_lengths_for(cfg: dict, species, radiation: str) -> dict[str, float]:
    """The model's own neutron scattering lengths for the atom names that occur in the structure.

    X-ray runs get none (they always use the element's form factor)."""
    if radiation != "neutron":
        return {}
    custom = {str(k).upper(): float(v) for k, v in (cfg.get("scattering") or {}).items()}
    names = sorted({str(s).strip().upper() for s in species} - {VOID})
    return {nm: custom[nm] for nm in names if nm in custom}


def fourier_config(cfg: dict, calc3: dict, aver: float, outfile: str, seed8, lots_ok: bool, species=()):
    from ..fourier_macro import FourierConfig, Lots, Seed
    from ..lattice import Lattice
    from ..reciprocal_grid import ReciprocalGrid
    cell = cfg["cell"]
    lat = Lattice.from_parameters(cell["a"], cell["b"], cell["c"], cell["alpha"], cell["beta"], cell["gamma"])
    steps = [float(x) for x in calc3["steps"]]
    hmax = [float(x) for x in calc3["hmax"]]
    shape = tuple(2 * int(round(h / s)) + 1 for h, s in zip(hmax, steps))
    grid = ReciprocalGrid.centered(lat, shape, steps)
    lots = None
    if lots_ok and calc3.get("lots", {}).get("enabled"):
        L = calc3["lots"]
        lots = Lots(L.get("shape", "box"), tuple(int(x) for x in L["size"]),
                    "all" if str(L.get("count", "all")) == "all" else int(L["count"]), bool(L.get("periodic", True)))
    radiation = calc3.get("radiation") or "xray"
    # an explicit 8-value seed is always written: needed for 0 < aver < 100 and random lots, harmless otherwise
    return FourierConfig(grid, outputs=(("int", outfile),), aver_percent=aver, lots=lots, seed=Seed(tuple(seed8)),
                         temperature="use" if cfg.get("calc", {}).get("debye_waller") else "ignore",
                         symmetry=calc3.get("symmetry") or None, radiation=radiation,
                         wavelength=fourier_wavelength(calc3.get("wavelength")),
                         scattering_lengths=scattering_lengths_for(cfg, species, radiation))


def run_fourier(cfg, out: Path, stru_path: Path, n, launcher, seed8, *, dry_run: bool, log, aver, name,
                lots_ok=True, species=()):
    """Stage (and unless ``dry_run`` run) one DISCUS Fourier calculation; returns a ScatteringVolume or None."""
    from ..diffuse3d import Diffuse3DRunner
    from ..fourier_macro import StructureSource, prepare_fourier_run
    calc3 = cfg["calc"]["fourier"]
    outfile = f"{name}.h5"
    fc = fourier_config(cfg, calc3, aver, outfile, seed8, lots_ok, species)
    if fc.scattering_lengths:
        log("  neutron scattering lengths of their own: "
            + ", ".join(f"{nm} {b:g} fm" for nm, b in fc.scattering_lengths))
    runner = Diffuse3DRunner(out / "runs_3d", launcher=launcher, output_filename=outfile,
                             timeout_seconds=int(cfg["calc"].get("timeout_s", 86400)))
    src = StructureSource.stru(stru_path.name, crystal_cells=tuple(n))
    run = prepare_fourier_run(runner, fc, src, stru_path, title=f"{cfg.get('title', 'model')} {name}")
    log(f"  staged {name}: grid {fc.grid.shape}, set aver {aver:g} → {run.directory}")
    if dry_run:
        return None
    runner.execute(run, native_nexus=True)
    log(f"  DISCUS finished {name}")
    return run.load_scattering_volume()


def pdf_config(cfg: dict, stru_path: Path, species=(), log=None):
    """The model's 1D PDF settings as a :class:`~discuspdf.PDFConfig` on the supercell ``.stru``."""
    from ..pdf import PDFConfig, QmaxLimitWarning
    p = cfg["calc"]["pdf"]
    radiation = p.get("radiation", "xray")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", QmaxLimitWarning)
        pc = PDFConfig(structure=stru_path, supercell=(1, 1, 1), rmax=float(p["rmax"]), dr=float(p["dr"]),
                       qmax=float(p["qmax"]), qsig=float(p.get("qsig", 0.0)), radiation=radiation,
                       qmin=float(p.get("qmin", 0.5)), dq=float(p.get("dq", 0.01)),
                       wavelength=float(p["wavelength"]) if str(p.get("wavelength") or "").strip() else None,
                       scattering_lengths=scattering_lengths_for(cfg, species, radiation))
    for warning in w:
        if log:
            log(f"  WARNING: {warning.message}")
    return pc


def run_powder_pdf(cfg, out: Path, stru_path: Path, launcher, *, dry_run: bool, log, species=()):
    """1D PDF G(r) of the disordered supercell (Debye sum, ``read stru``); returns the PDFRun or None."""
    from ..pdf import PDFRunner
    pc = pdf_config(cfg, stru_path, species, log)
    if pc.scattering_lengths:
        log("  neutron scattering lengths of their own: "
            + ", ".join(f"{nm} {b:g} fm" for nm, b in pc.scattering_lengths))
    runner = PDFRunner(out / "runs_pdf", launcher=launcher, timeout_seconds=int(cfg["calc"].get("timeout_s", 86400)))
    run = runner.prepare(pc)
    log(f"  staged 1D PDF → {run.directory}")
    if dry_run:
        return None
    return runner.execute(run)


def species_of(structure) -> np.ndarray:
    return np.asarray(structure.species, dtype=object)
