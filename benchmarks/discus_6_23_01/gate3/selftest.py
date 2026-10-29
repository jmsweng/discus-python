"""Gate 3 evaluator self-test: synthetic runs built from the Gate 2 v4 volumes and the source emulation must PASS when
clean and FAIL (or be INVALID) for each kind of corruption. Proves the evaluator can fail; says nothing about DISCUS."""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util  # noqa: E402

import analyze as AN  # noqa: E402
import native3dpdf as N  # noqa: E402

_spec = importlib.util.spec_from_file_location("gate3_generate", HERE / "generate.py")   # avoid gate2's generate.py
G = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(G)

LOG = """          ***********************************************************
          *         D I S C U S - S U I T E  Version 6.23.01        *
          ***********************************************************
   Temp. factors      : {temp}
                   Atom typ = CU         (#      2048 )
                   Atom typ = AU         (#      2048 )
"""


def write_native(path: Path, values, corner, inc):
    import h5py
    with h5py.File(path, "w") as f:
        g = f.create_group("entry").create_group("data")
        g["data_values"] = values; g["data_corner"] = corner; g["data_increment_vector"] = np.diag(inc)
        g["data_type_style"] = np.bytes_("single_pdf"); g["unit_cell_lengths"] = np.full(3, 4.0)


def make_run(work: Path, corrupt: str | None) -> Path:
    import h5py
    run = work / (corrupt or "clean")
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(HERE / "inputs", run)
    corner, inc = N.native_axes(G.SHAPE, G.VI)
    for kind in G.KINDS:
        p = G.v4_file(kind, run / f"gate3_block_{kind}_int.h5")
        with h5py.File(p, "r") as f:
            I = f["entry/data/data_values"][()]
        if corrupt == "int_changed":
            with h5py.File(p, "r+") as f:
                f["entry/data/data_values"][0, 0, 1] *= 1 + 1e-12
        if kind == "total":
            J, _ = N.bragg_reweight(I, G.CORNER, G.VI, 2048 if corrupt == "bragg_reweighted" else G.N_CELLS)
        else:
            J = I
        nat = N.native_3dpdf(J)
        c, d = corner, inc
        if corrupt == "unitary_missing":
            nat = nat * np.sqrt(np.prod(G.SHAPE))
        elif corrupt == "origin_shift":
            nat = np.roll(nat, 1, axis=0)
        elif corrupt == "true_axes_metadata":
            c, d = N.true_axes(G.SHAPE, G.VI)
        write_native(run / f"gate3_block_{kind}_3dpdf.h5", nat, c, d)
        (run / f"gate3_block_{kind}.log").write_text(LOG.format(temp="used" if corrupt == "temp_used" else "ignored"))
    if corrupt == "edited_macro":
        p = run / "gate3_block_diffuse.mac"; p.write_text(p.read_text().replace("set aver,100.0", "set aver,50.0"))
    return run


CASES = {None: ("PASS", None), "unitary_missing": ("FAIL", "E1"), "origin_shift": ("FAIL", "E1"),
         "true_axes_metadata": ("FAIL", "M1"), "bragg_reweighted": ("FAIL", "E1 total"), "int_changed": ("FAIL", "R0"),
         "temp_used": ("FAIL", "L "), "edited_macro": ("INVALID", None)}


def run_selftest(work: Path) -> list:
    rep = []
    for corrupt, (want, frag) in CASES.items():
        res = AN.analyze(make_run(work, corrupt))
        failed = [c["criterion"] for c in res.get("criteria", []) if not c["pass"]]
        ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
        rep.append((corrupt or "clean", want, res["status"], ok, failed[:4]))
    return rep


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="gate3_selftest_"))
    rep = run_selftest(work)
    for name, want, got, ok, failed in rep:
        print(f"{'ok ' if ok else 'BAD'} {name:20s} expected {want:7s} got {got:7s} failed: {failed}")
    sys.exit(0 if all(r[3] for r in rep) else 1)


if __name__ == "__main__":
    main()
