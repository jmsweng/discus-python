"""Evaluate the builder-backend probe against its pre-registered criteria (PREREG.md, criteria.json).

    python benchmarks/discus_6_23_01/e_builder_backend/analyze.py <evidence dir> [--json out.json]
    python benchmarks/discus_6_23_01/e_builder_backend/analyze.py --selftest

<evidence dir> holds work/<run>/ (as run by run_probe.sh), exit_codes.txt and input_sha256.txt.
Exit code 0 = every criterion passed, 1 = a criterion failed, 2 = selftest failed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_inputs as mk  # noqa: E402

from discuspdf.diffuse3d import _ENGINE_ERROR_PATTERN  # noqa: E402
from discuspdf.reference_scattering import ScatteringFactors, diffuse_intensity  # noqa: E402

CRITERIA = json.loads((HERE / "criteria.json").read_text(encoding="utf-8"))


def rel_l2(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / np.linalg.norm(b))


def max_rel(a, b) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(b)).max() / np.abs(b).max())


def oracle_mix(hkl: np.ndarray) -> np.ndarray:
    fac = ScatteringFactors.neutron(mk.B_DISCUS, "DISCUS 6.23.01 per_bcoh")
    return diffuse_intensity(mk.structure("mix"), fac, hkl=hkl)


def evaluate(fourier: dict, pdf: dict, hkl: np.ndarray, meta: dict) -> tuple[dict, list[str]]:
    """fourier: id -> 3D intensity array; pdf: id -> (r, G); meta: id -> dict(bragg, shape)."""
    res: dict = {}
    vals: dict = {}
    ref = fourier["f_mix"]
    vals["E1"] = rel_l2(ref, oracle_mix(hkl))
    vals["E2"] = rel_l2(fourier["f_relabel_scat"], ref)
    vals["E3"] = float(np.abs(fourier["f_relabel_noscat"]).max() / np.abs(ref).max())
    vals["E4"] = rel_l2(fourier["f_void"], ref)
    want_shape = list(mk.grid().shape)
    vals["E5"] = float(sum(m["bragg"] != "bragg_subtracted" or list(m["shape"]) != want_shape
                           for m in meta.values()))
    expected_r = np.arange(round(mk.PDF["rmax"] / mk.PDF["dr"]) + 1) * mk.PDF["dr"]
    bad_grid = 0
    for r, g in pdf.values():
        if len(r) != len(expected_r) or not np.allclose(r, expected_r, atol=1e-6) or not np.isfinite(g).all():
            bad_grid += 1
    vals["P0"] = float(bad_grid)
    gm = pdf["p_mix"][1]
    vals["P1"] = max_rel(pdf["p_relabel_scat"][1], gm)
    vals["P2"] = max_rel(pdf["p_relabel_noscat"][1], gm)
    vals["P3"] = max_rel(pdf["p_void"][1], gm)
    vals["P4"] = max_rel(pdf["p_ni_stru_n"][1], pdf["p_ni_cll_n"][1])
    vals["P5"] = max_rel(pdf["p_ni_stru_x"][1], pdf["p_ni_cll_x"][1])
    failed = []
    for key, c in CRITERIA.items():
        v = vals[key]
        ok = v <= c["max"] if "max" in c else v >= c["min"]
        res[key] = dict(c, value=v, passed=bool(ok))
        if not ok:
            failed.append(key)
    return res, failed


def check_integrity(ev: Path) -> tuple[list[str], list[str]]:
    """(problems, notes). Problems fail E0: a run missing or with a nonzero exit code, an input that differs
    from the committed file, a missing or extra input file, a missing log. Notes are informational (I0):
    log lines with an engine marker (``diffuse3d._ENGINE_ERROR_PATTERN``), e.g. a tagged DISCUS warning."""
    spec = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
    problems, notes = [], []
    codes = dict(l.split(" exit=") for l in (ev / "exit_codes.txt").read_text().splitlines() if " exit=" in l)
    ids = {r["id"] for r in spec["runs"]}
    if set(codes) != ids or any(v.strip() != "0" for v in codes.values()):
        problems.append(f"exit codes: {codes}")
    committed = {p.relative_to(HERE / "inputs").as_posix() for p in (HERE / "inputs").rglob("*") if p.is_file()}
    hashed = set()
    for line in (ev / "input_sha256.txt").read_text().splitlines():
        h, name = line.split(maxsplit=1)
        name = name.strip().removeprefix("./")
        hashed.add(name)
        f = HERE / "inputs" / name
        if not f.is_file():
            problems.append(f"input {name} is not a committed file")
        elif hashlib.sha256(f.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != h:
            problems.append(f"input {name} differs from the committed file")
    if committed - hashed:
        problems.append(f"committed inputs not used: {sorted(committed - hashed)}")
    for r in spec["runs"]:
        log = ev / f"{r['id']}.log"
        if not log.is_file():
            problems.append(f"{r['id']}: no engine log")
            continue
        for line in log.read_text(errors="replace").splitlines():
            if _ENGINE_ERROR_PATTERN.search(line):
                notes.append(f"{r['id']}: {line.strip()}")
    return problems, notes


def load_evidence(ev: Path):
    from discuspdf.discus_output import read_discus_volume
    from discuspdf.pdf import read_pdf
    spec = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
    problems, notes = check_integrity(ev)
    fourier, meta, pdf, hkl = {}, {}, {}, None
    for r in spec["runs"]:
        d = ev / "work" / r["id"]
        if r["type"] == "fourier":
            vol = read_discus_volume(d / f"{r['id']}.h5")
            fourier[r["id"]] = np.asarray(vol.values, float)
            meta[r["id"]] = {"bragg": vol.with_bragg, "shape": vol.values.shape}
            if hkl is None:
                hkl = np.einsum("ab,ijkb->ijka", vol.increments,
                                np.indices(vol.values.shape).transpose(1, 2, 3, 0)) + vol.corner
        else:
            pdf[r["id"]] = read_pdf(d / "result.pdf")
    return fourier, pdf, hkl, meta, problems, notes


def selftest() -> int:
    """The criteria accept oracle-consistent data and reject the failure modes they are meant to catch."""
    g = mk.grid()
    h, k, l = (np.asarray(x) for x in g.hkl_axes())
    hkl = np.stack(np.meshgrid(h, k, l, indexing="ij"), axis=-1)
    ref = oracle_mix(hkl)
    meta = {i: {"bragg": "bragg_subtracted", "shape": ref.shape} for i in ("f_mix", "f_relabel_scat", "f_relabel_noscat", "f_void")}
    r = np.arange(501) * 0.02
    gm = np.sin(3 * r) * np.exp(-r / 5)
    good_f = {"f_mix": ref, "f_relabel_scat": ref * (1 + 1e-12), "f_relabel_noscat": ref * 1e-9, "f_void": ref.copy()}
    good_p = {k: (r, gm.copy()) for k in ("p_mix", "p_relabel_scat", "p_void", "p_ni_stru_n", "p_ni_cll_n",
                                          "p_ni_stru_x", "p_ni_cll_x")}
    good_p["p_relabel_noscat"] = (r, gm * 1.1)
    ok = True
    res, failed = evaluate(good_f, good_p, hkl, meta)
    if failed:
        print("selftest: consistent data failed", failed)
        ok = False
    cases = {
        "E1": ({**good_f, "f_mix": ref * 1.01, "f_relabel_scat": ref * 1.01, "f_void": ref * 1.01,
                "f_relabel_noscat": ref * 1e-9}, good_p),
        "E2": ({**good_f, "f_relabel_scat": good_f["f_relabel_noscat"]}, good_p),   # scat ignored
        "E3": ({**good_f, "f_relabel_noscat": ref}, good_p),                         # NI63 not read as Ni
        "E4": ({**good_f, "f_void": ref * 1.25}, good_p),                            # VOID scatters
        "P1": (good_f, {**good_p, "p_relabel_scat": good_p["p_relabel_noscat"]}),
        "P2": (good_f, {**good_p, "p_relabel_noscat": good_p["p_mix"]}),
        "P3": (good_f, {**good_p, "p_void": (r, gm * 0.8)}),                         # VOID in normalisation
        "P4": (good_f, {**good_p, "p_ni_stru_n": (r, gm * 1.001)}),
        "P5": (good_f, {**good_p, "p_ni_stru_x": (r, gm * 1.001)}),
        "P0": (good_f, {**good_p, "p_mix": (r[:-1], gm[:-1]), "p_relabel_scat": (r[:-1], gm[:-1]),
                        "p_void": (r[:-1], gm[:-1]), "p_relabel_noscat": (r[:-1], 1.1 * gm[:-1])}),
    }
    for key, (f, p) in cases.items():
        _, failed = evaluate(f, p, hkl, meta)
        if key not in failed:
            print(f"selftest: case {key} was not caught (failed: {failed})")
            ok = False
    bad_meta = dict(meta, f_void={"bragg": "bragg_present", "shape": ref.shape})
    if "E5" not in evaluate(good_f, good_p, hkl, bad_meta)[1]:
        print("selftest: case E5 was not caught")
        ok = False
    ok = _selftest_integrity() and ok
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 2


def _selftest_integrity() -> bool:
    """check_integrity on a fake evidence directory laid out as run_probe.sh writes it."""
    import shutil
    import tempfile
    ok = True
    spec = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory() as tmp:
        ev = Path(tmp)
        shutil.copytree(HERE / "inputs", ev / "work")
        lines = []
        for p in sorted((ev / "work").rglob("*")):
            if p.is_file():
                lines.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  ./{p.relative_to(ev / 'work').as_posix()}")
        (ev / "input_sha256.txt").write_text("\n".join(lines) + "\n")
        (ev / "exit_codes.txt").write_text("".join(f"{r['id']} exit=0\n" for r in spec["runs"]))
        for r in spec["runs"]:
            (ev / f"{r['id']}.log").write_text("discus done\n")
        (ev / "f_mix.log").write_text(" ***FOUR*** something DISCUS printed\n")
        problems, notes = check_integrity(ev)
        if problems or len(notes) != 1:
            print("selftest integrity: clean evidence gave", problems, notes)
            ok = False
        bad = lines[0].split()[0].replace(lines[0][0], "0" if lines[0][0] != "0" else "1", 1)
        (ev / "input_sha256.txt").write_text("\n".join([bad + "  " + lines[0].split()[1]] + lines[2:]) + "\n")
        (ev / "exit_codes.txt").write_text("f_mix exit=1\n")
        (ev / "p_void.log").unlink()
        problems, _ = check_integrity(ev)
        want = ["exit codes", "differs from the committed file", "committed inputs not used", "p_void: no engine log"]
        missing = [w for w in want if not any(w in p for p in problems)]
        if missing:
            print("selftest integrity: not caught", missing, problems)
            ok = False
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence", type=Path, nargs="?")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.evidence is None:
        ap.error("give the evidence directory or --selftest")
    fourier, pdf, hkl, meta, problems, notes = load_evidence(a.evidence)
    res, failed = evaluate(fourier, pdf, hkl, meta)
    if problems:
        failed = ["E0"] + failed
        res["E0"] = {"what": "integrity (exit codes, input hashes, logs)", "problems": problems, "passed": False}
    res["I0"] = {"what": "engine-marker lines in the logs (informational)", "lines": notes}
    for key, row in res.items():
        if key == "E0":
            print(f"E0 FAIL integrity: {row['problems']}")
            continue
        if key == "I0":
            print(f"I0 informational: {len(row['lines'])} engine-marker line(s)" + "".join(f"\n   {l}" for l in row["lines"][:20]))
            continue
        bound = f"<= {row['max']:g}" if "max" in row else f">= {row['min']:g}"
        print(f"{key} {'PASS' if row['passed'] else 'FAIL'} {row['value']:.3e} {bound}  {row['what']}")
    print(f"criteria failed: {len(failed)} {failed}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=1, default=float) + "\n", encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
