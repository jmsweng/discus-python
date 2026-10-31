"""Load discuspdf-written MDHisto files in Mantid and record what Mantid sees (run with Mantid's Python).

    <mantid-python> benchmarks/mantid_7_0_0/check_in_mantid.py tests/fixtures/mantid_7_0_0 <written dir>

For each <name>.discuspdf.nxs:
* LoadMD must succeed; dimensions (names, frames, units, bins, extents), signal, errors² and the
  lattice are compared with Mantid's own load of the original fixture <name>.nxs;
* the workspace is saved again with SaveMD as <name>.resaved.nxs (to check that Mantid keeps the
  projection transform when it writes the file itself).
For skewed.discuspdf.nxs (no Mantid original) the values are checked against skewed.expected.json and
Mantid's basis vectors and origin are reported.
Writes mantid_check.json in the written dir and exits 1 if any check fails.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import mantid
from mantid.simpleapi import LoadMD, SaveMD, mtd


def describe(ws) -> dict:
    dims = []
    for i in range(ws.getNumDims()):
        d = ws.getDimension(i)
        dims.append({"name": d.name, "units": str(d.getUnits()), "frame": d.getMDFrame().name(),
                     "bins": d.getNBins(), "min": d.getMinimum(), "max": d.getMaximum()})
    lat = None
    if ws.getNumExperimentInfo() and ws.getExperimentInfo(0).sample().hasOrientedLattice():
        o = ws.getExperimentInfo(0).sample().getOrientedLattice()
        lat = [o.a(), o.b(), o.c(), o.alpha(), o.beta(), o.gamma()]
    return {"dims": dims, "lattice": lat}


def masked_array(ws, scratch: Path) -> np.ndarray:
    """Mantid's Python API does not expose the bin mask, so save the workspace with SaveMD and read
    the mask Mantid writes (file axes are slowest-first: reverse them)."""
    import h5py
    tmp = scratch / f"_mask_{ws.name()}.nxs"
    SaveMD(ws, Filename=str(tmp))
    with h5py.File(tmp, "r") as f:
        mask = np.asarray(f["MDHistoWorkspace/data/mask"][()]).transpose() != 0
    tmp.unlink()
    return mask


def main(src: Path, written: Path) -> int:
    report = {"mantid_version": mantid.__version__, "files": {}}
    failures = []
    for path in sorted(written.glob("*.discuspdf.nxs")):
        name = path.name.replace(".discuspdf.nxs", "")
        entry: dict = {}
        ours = LoadMD(str(path), OutputWorkspace=f"ours_{name}")
        entry["ours"] = describe(ours)
        s_ours = np.array(ours.getSignalArray())
        e_ours = np.array(ours.getErrorSquaredArray())
        try:
            m_ours = masked_array(ours, written)
        except Exception as exc:          # pragma: no cover - API differences between versions
            m_ours, entry["mask_error"] = None, str(exc)
        original = src / f"{name}.nxs"
        if original.is_file():
            ref = LoadMD(str(original), OutputWorkspace=f"ref_{name}")
            entry["reference"] = describe(ref)
            s_ref = np.array(ref.getSignalArray())
            e_ref = np.array(ref.getErrorSquaredArray())
            m_ref = masked_array(ref, written)
            # discuspdf masks NaN bins (no coverage) in addition to Mantid's own mask, by design.
            expected_mask = m_ref | ~np.isfinite(s_ref)
            same_mask = m_ours is not None and np.array_equal(m_ours, expected_mask)
            entry["mask_counts"] = {"mantid": int(m_ref.sum()), "nan": int((~np.isfinite(s_ref)).sum()),
                                   "ours": None if m_ours is None else int(m_ours.sum())}
            usable = ~m_ref & np.isfinite(s_ref)
            checks = {
                "dim_names": [d["name"] for d in entry["ours"]["dims"]] == [d["name"] for d in entry["reference"]["dims"]],
                "dim_frames": all(d["frame"] == "HKL" for d in entry["ours"]["dims"]),
                "dim_bins": [d["bins"] for d in entry["ours"]["dims"]] == [d["bins"] for d in entry["reference"]["dims"]],
                "dim_extents": bool(np.allclose([[d["min"], d["max"]] for d in entry["ours"]["dims"]],
                                                [[d["min"], d["max"]] for d in entry["reference"]["dims"]],
                                                rtol=0, atol=1e-5)),
                "signal_unmasked": bool(np.array_equal(s_ours[usable], s_ref[usable])),
                "errors_squared_unmasked": bool(np.array_equal(e_ours[usable], e_ref[usable])),
                "mask": bool(same_mask),
                "lattice": bool(np.allclose(entry["ours"]["lattice"], entry["reference"]["lattice"], atol=1e-9)),
            }
        else:
            exp = json.loads((written / f"{name}.expected.json").read_text())
            i0, i1, i2 = np.meshgrid(*[np.arange(n) for n in s_ours.shape], indexing="ij")
            f = 100.0 * i0 + 10.0 * i1 + i2
            masked_index = tuple(exp["masked_index"])
            ok = np.ones(f.shape, bool)
            ok[masked_index] = False
            checks = {
                "signal_formula": bool(np.array_equal(s_ours[ok], f[ok])),
                "errors_squared_formula": bool(np.allclose(e_ours[ok], 0.5 * f[ok])),
                "mask": bool(m_ours is not None and m_ours[masked_index] and m_ours.sum() == 1),
                "dim_frames": all(d["frame"] == "HKL" for d in entry["ours"]["dims"]),
                "lattice": bool(np.allclose(entry["ours"]["lattice"], exp["cell"], atol=1e-9)),
            }
        resaved = written / f"{name}.resaved.nxs"
        SaveMD(f"ours_{name}", Filename=str(resaved))
        entry["resaved"] = resaved.name
        entry["checks"] = checks
        failures += [f"{name}: {k}" for k, v in checks.items() if not v]
        report["files"][name] = entry
    report["failures"] = failures
    (written / "mantid_check.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v["checks"] for k, v in report["files"].items()}, indent=1))
    print("FAILURES:", failures or "none")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()))
