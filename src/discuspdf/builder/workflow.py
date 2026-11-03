"""Run a disorder-builder model end to end: structure → optional DISCUS mmc → files → 3D → 1D PDF.

From Python::

    from discuspdf.builder import run
    run(model)                                   # model: dict, JSON text or path to model.json
    run(model, dry_run=True, out="results")      # stage every DISCUS input, start nothing

From the command line::

    python -m discuspdf.builder run model.json [--dry-run] [--numpy-preview] [--launcher "..."] ...

The script the web page exports embeds the model and calls :func:`main`, so it takes the same options.
"""
from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

import numpy as np

from .model import VOID, export_formats, load_model

ENGINE_VERSION = "2.0"   # 1.0 was the builder's self-contained engine.py


def run(model, *, out=None, seed=None, launcher=None, dry_run=False, structure_only=False,
        numpy_preview=False, overwrite=False, log=None) -> dict:
    """Build the model's structure and run every enabled calculation. Returns a summary dict
    (``out``, ``structure``, ``info`` and the result objects of the runs that were executed)."""
    from ..discus_io import write_stru
    from ..xyz_io import write_extxyz
    from . import outputs, runs
    from .realize import build_structure, clash_report

    cfg = load_model(model)
    formats = export_formats(cfg)
    seed = int(seed if seed is not None else cfg.get("seed", 1))
    out = Path(out or cfg.get("output_dir") or "discus_model_out").resolve()
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError(f"{out} exists and is not empty; use overwrite=True (--overwrite) or another --out")
    out.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def say(msg):
        log_lines.append(msg)
        (log or (lambda m: print(m, flush=True)))(msg)

    calc = cfg.get("calc", {})
    if launcher is None:
        launcher = calc.get("launcher") or "discus_suite -noautorun -macro"
    launcher = shlex.split(launcher) if isinstance(launcher, str) else list(launcher)
    rng_discus = np.random.default_rng([seed, 99])
    result: dict = {"out": out}

    say(f"== {cfg.get('title', 'model')} (seed {seed}) ==")
    say("Building the supercell …")
    structure, info = build_structure(cfg, seed, say)
    say(f"  {structure.n_atoms} rows ({info['n_void']} VOID padding), {info['atoms_per_cell']} per cell, "
        f"supercell {tuple(int(x) for x in structure.supercell)}, {info['n_molecules']} molecules")
    clashes = clash_report(structure, info, float(cfg.get("disorder", {}).get("clash_threshold", 0.8)), say)

    mc = cfg.get("mmc") or {}
    if mc.get("enabled"):
        if info["n_molecules"]:
            raise ValueError("DISCUS mmc is only offered for atom-only models (molecules would be torn apart)")
        from ..mc import MCRunner
        mcc = runs.build_mc_config(mc, runs.seed8_from(rng_discus))
        for c in mcc.cautions():
            say(f"  mmc note: {c}")
        runner = MCRunner(out / "runs_mc", launcher=launcher, timeout_seconds=int(calc.get("timeout_s", 86400)))
        mrun = runner.prepare(mcc, structure, title=f"{cfg.get('title', 'model')} mmc")
        say(f"  staged DISCUS mmc → {mrun.directory}")
        if not dry_run and not structure_only:
            res = runner.execute(mrun)
            structure = res.structure
            info["is_void"] = np.asarray(structure.species, dtype=object) == VOID
            info["molecule_id"] = np.full(structure.n_atoms, -1)
            say(f"  mmc finished: {res.status}")
            result["mmc"] = res

    stru_path = out / "supercell.stru"
    title = " ".join(str(cfg.get("title", "model")).split())[:60] or "model"
    # molecules in several orientations and VOID padding on top of real atoms: rows stay cell-major, but
    # the read-back can neither infer cells from positions nor require a minimum separation
    write_stru(stru_path, structure, title=title, overwrite=True,
               check_order=not info["n_molecules"], min_separation=None)
    write_extxyz(out / "supercell.xyz", structure, mask=~info["is_void"], overwrite=True)
    summary = {k: v for k, v in info.items() if k not in ("molecule_id", "is_void", "unit_lattice")}
    summary.update(seed=seed, clashes=clashes, engine=ENGINE_VERSION, n_rows=int(structure.n_atoms),
                   formula=structure.select(mask=~info["is_void"]).formula if (~info["is_void"]).any() else "",
                   outputs=formats)
    (out / "structure_summary.json").write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")
    (out / "model.json").write_text(json.dumps(cfg, indent=1), encoding="utf-8")
    say(f"  wrote {stru_path.name}, supercell.xyz, structure_summary.json")
    result.update(structure=structure, info=info, stru=stru_path)

    if numpy_preview:
        outputs.numpy_preview(cfg, structure, info, out, say)
    if structure_only:
        (out / "run_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        return result

    n = tuple(int(x) for x in structure.supercell)
    seed8 = runs.seed8_from(rng_discus)
    f3 = calc.get("fourier", {})
    custom_b = {str(k).upper() for k in (cfg.get("scattering") or {})}
    if custom_b and calc.get("debye_waller") and structure.b_iso is not None:
        names = np.char.upper(structure.species.astype(str))
        for nm in sorted(custom_b):
            if len(np.unique(structure.b_iso[names == nm])) > 1:
                say(f"  WARNING: {nm} has several B values, so DISCUS makes several atom types of that name; "
                    "`scat` sets the scattering length of the first one only. Give the atoms one B or switch "
                    "Debye–Waller off.")
    if f3.get("enabled"):
        say(f"3D diffuse scattering (DISCUS fourier, {f3.get('radiation') or 'xray'}, "
            f"λ {runs.fourier_wavelength(f3.get('wavelength'))}, Bragg subtracted) …")
        diffuse = runs.run_fourier(cfg, out, stru_path, n, launcher, seed8, dry_run=dry_run, log=say,
                                   species=structure.species, aver=float(f3.get("aver_percent", 100.0)),
                                   name="diffuse")
        if diffuse is not None:
            result["diffuse"] = diffuse
            if f3.get("delta_pdf", True):
                say("3D-ΔPDF (fractional FFT of the diffuse volume) …")
            outputs.write_volume_outputs(diffuse, out, "diffuse", formats, say,
                                         delta_pdf=bool(f3.get("delta_pdf", True)), real_name="3D-ΔPDF",
                                         real_quantity="delta_pdf")
        if f3.get("total"):
            say("Total 3D scattering (Bragg present) and 3D-PDF (Patterson) …")
            total = runs.run_fourier(cfg, out, stru_path, n, launcher, seed8, dry_run=dry_run, log=say,
                                     aver=0.0, name="total", species=structure.species, lots_ok=False)
            if total is not None:
                result["total"] = total
                outputs.write_volume_outputs(total, out, "total", formats, say, delta_pdf=True,
                                             real_name="3D-PDF", real_quantity="patterson",
                                             require_diffuse=False)

    pdf = calc.get("pdf", {})
    if pdf.get("enabled"):
        say("1D powder PDF (DISCUS powder, Debye sum over the supercell) …")
        if structure.n_atoms > 60000:
            say(f"  note: {structure.n_atoms} atoms; the Debye sum scales as N², this can take long")
        prun = runs.run_powder_pdf(cfg, out, stru_path, launcher, dry_run=dry_run, log=say,
                                   species=structure.species)
        if prun is not None:
            result["pdf"] = outputs.write_pdf_outputs(prun, out, formats, cfg.get("title", "model"), say)

    if dry_run:
        say("Dry run: DISCUS inputs are staged; run without --dry-run to execute them.")
    say(f"Done. Results in {out}")
    (out / "run_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return result


def _parser(prog=None, description=None) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog=prog, description=description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="read the model from a JSON file instead of the embedded one")
    ap.add_argument("--out", help="output directory (default: from the model)")
    ap.add_argument("--seed", type=int, help="override the model's random seed")
    ap.add_argument("--launcher", help='command that runs a DISCUS macro, e.g. "discus_suite -noautorun -macro"')
    ap.add_argument("--dry-run", action="store_true", help="build and stage everything; do not start DISCUS")
    ap.add_argument("--structure-only", action="store_true", help="stop after writing the structure")
    ap.add_argument("--numpy-preview", action="store_true", help="numpy diffuse plane (no DISCUS needed)")
    ap.add_argument("--overwrite", action="store_true", help="allow writing into an existing output directory")
    return ap


def main(argv=None, *, model=None, description=None) -> int:
    """Command line of an exported script: ``model`` is the embedded model (``--config`` overrides it)."""
    args = _parser(description=description).parse_args(argv)
    source = args.config or model
    if source is None:
        sys.exit("No model: pass --config model.json")
    try:
        run(source, out=args.out, seed=args.seed, launcher=args.launcher, dry_run=args.dry_run,
            structure_only=args.structure_only, numpy_preview=args.numpy_preview, overwrite=args.overwrite)
    except FileExistsError as e:
        sys.exit(str(e))
    return 0


def cli(argv=None) -> int:
    """``python -m discuspdf.builder run model.json [options]``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: python -m discuspdf.builder run MODEL.json [--out DIR] [--seed N] [--launcher CMD]\n"
              "                                   [--dry-run] [--structure-only] [--numpy-preview] [--overwrite]")
        return 0 if argv else 2
    if argv[0] != "run" or len(argv) < 2:
        print("usage: python -m discuspdf.builder run MODEL.json [options]", file=sys.stderr)
        return 2
    return main(["--config", argv[1], *argv[2:]])
