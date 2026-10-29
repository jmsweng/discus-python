"""Generate the pre-registered Gate 3 inputs and predictions (run BEFORE the engine run).

Question: what exactly is DISCUS 6.23.01's native `output -> value 3DPDF` (values AND axis metadata), and how does it
convert to the package's absolutely-normalized `delta_pdf_fractional`?

Design: the Gate 2 v4 block structure and grid, unchanged (`gate2_block.stru`, 128^3, H = -4 + k/16, B = 0 via
`temp ignore`). One Fourier run per macro writes `value int` and then `value 3DPDF` (the order matters: writing
3DPDF consumes dsi3d). Two macros: set aver,0 (total) and set aver,100 (diffuse).

Predictions are computed from the Gate 2 v4 engine volumes, recovered from git history (commit 8cd59ac) and checked
against EVIDENCE_CONTENT_HASHES.json, through native3dpdf.py (a parameter-free emulation written from the DISCUS
source). Nobody has seen a native 3DPDF output on this grid.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
REPO = BENCH.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(BENCH)); sys.path.insert(0, str(BENCH / "gate2_prereg"))
import native3dpdf as N  # noqa: E402
import oracle as O  # noqa: E402
import regenerate_evidence as R  # noqa: E402

INPUTS = HERE / "inputs"
V4_COMMIT = "8cd59ac"
V4_DIR = "benchmarks/discus_6_23_01/gate2_prereg/evidence/20261029T010305Z_LAPTOP-NGU3BT0N"
N_CELLS = O.N_CELL ** 3
VI = np.full(3, O.DH)
CORNER = np.full(3, -O.H_MAX)
SHAPE = (O.N_H,) * 3

MACRO = """# Gate 3 pre-registered run: block model, set aver,{aver}; value int then value 3DPDF (generated; do not edit)
discus
read
  stru gate2_block.stru
eval n[1]
eval n[2]
eval n[3]
fourier
  xray
  wvle moa1
  ll -4.0,-4.0,-4.0
  lr  3.9375,-4.0,-4.0
  ul -4.0, 3.9375,-4.0
  tl -4.0,-4.0, 3.9375
  na 128
  no 128
  nt 128
  abs h
  ord k
  top l
  temp ignore
  disp off
  lots off
  set technique:turbo
  set aver,{aver}
  show
  run
exit
output
  value int
  form nexus
  outf gate3_block_{kind}_int.h5
  run
  value 3DPDF
  form nexus
  outf gate3_block_{kind}_3dpdf.h5
  run
exit
exit
exit
"""

THRESHOLDS = {
    "R0_int_content_hash_equals_gate2_v4": True,
    "M1_native_metadata_abs_tol": 1e-9,
    "E1_native_vs_emulation_of_same_run_int_relL2_max": 1e-12,
    "E2_native_fingerprint_rtol": 1e-9,
    "C1_pipeline_vs_converted_native_relL2_max": 1e-12,
}
KINDS = {"total": ("0.0", "bragg_present"), "diffuse": ("100.0", "bragg_subtracted")}


def v4_file(kind: str, dest: Path) -> Path:
    """Write the Gate 2 v4 engine volume (from git history) to dest, verified against its recorded content hash."""
    rel = f"{V4_DIR}/gate2_block_{kind}.h5"
    raw = subprocess.run(["git", "-C", str(REPO), "show", f"{V4_COMMIT}:{rel}"], check=True,
                         capture_output=True).stdout
    dest.write_bytes(raw)
    want = R.load_registry()["evidence"][V4_DIR]["volumes"][f"gate2_block_{kind}.h5"]
    assert R.content_hash(dest) == want, "v4 volume from history does not match its recorded content hash"
    return dest


def v4_volume(kind: str) -> np.ndarray:
    import h5py
    tmp = v4_file(kind, HERE / f".tmp_v4_{kind}.h5")
    try:
        with h5py.File(tmp, "r") as f:
            return f["entry/data/data_values"][()]
    finally:
        tmp.unlink()


def conversion_factor() -> float:
    """rho_package(u) = k * native(u):  k = sqrt(n1 n2 n3) * |det(vi)| / V_cell  (valid when scalef == 1)."""
    return float(np.sqrt(np.prod(SHAPE)) * np.prod(VI) / O.A_CELL ** 3)


def main():
    INPUTS.mkdir(exist_ok=True)
    shutil.copy2(BENCH / "gate2_prereg" / "inputs" / "gate2_block.stru", INPUTS / "gate2_block.stru")
    files = {"gate2_block.stru": INPUTS / "gate2_block.stru"}
    for kind, (aver, _) in KINDS.items():
        p = INPUTS / f"gate3_block_{kind}.mac"
        p.write_text(MACRO.format(aver=aver, kind=kind))
        files[p.name] = p
    nat_corner, nat_inc = N.native_axes(SHAPE, VI)
    true_corner, true_inc = N.true_axes(SHAPE, VI)
    pred = {"schema": "gate3-prereg-v1",
            "status": "PRE-REGISTERED: written before any native 3DPDF engine output on this grid existed",
            "source": "DiffuseCode tag v.6.23.01 (19a214a); see native3dpdf.py docstring",
            "design": {"structure": "Gate 2 v4 block model, seed 20261101", "grid": "128^3, H = -4 + k/16",
                       "ncells": N_CELLS, "v4_volumes_from": f"{V4_COMMIT}:{V4_DIR}"},
            "native_metadata_predicted": {"data_corner": nat_corner.tolist(), "increment": nat_inc.tolist(),
                                          "data_type_style": "single_pdf"},
            "true_axes": {"corner": true_corner.tolist(), "increment": true_inc.tolist(),
                          "origin_index_0based": [int(N._rev_index(k)[0]) for k in SHAPE],
                          "native_label_at_true_origin": (nat_corner + nat_inc * 64).tolist()},
            "conversion_factor_k": conversion_factor(), "thresholds": THRESHOLDS, "kinds": {}}
    for kind, (_, bragg) in KINDS.items():
        I = v4_volume(kind)
        J, scalef = N.bragg_reweight(I, CORNER, VI, N_CELLS) if kind == "total" else (I, None)
        nat = N.native_3dpdf(J)
        pred["kinds"][kind] = {"with_bragg": bragg, "four_weight_scalef": scalef,
                               "int_content_hash_gate2_v4": R.load_registry()["evidence"][V4_DIR]["volumes"][
                                   f"gate2_block_{kind}.h5"],
                               "native_fingerprint": O.fingerprint(nat),
                               "native_argmax_index": [int(v) for v in np.unravel_index(np.argmax(nat), nat.shape)],
                               "native_max": float(nat.max())}
    pred["input_sha256"] = {n: R.sha256_file(p) for n, p in sorted(files.items())}
    pred["code_sha256"] = {n: R.sha256_file(HERE / n) for n in ("native3dpdf.py", "generate.py", "analyze.py",
                                                                 "selftest.py")}
    (HERE / "predictions.json").write_text(json.dumps(pred, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: pred[k] for k in ("native_metadata_predicted", "true_axes", "conversion_factor_k")}, indent=1))
    for kind, v in pred["kinds"].items():
        print(kind, "scalef", v["four_weight_scalef"], "argmax", v["native_argmax_index"], "max", v["native_max"])


if __name__ == "__main__":
    main()
