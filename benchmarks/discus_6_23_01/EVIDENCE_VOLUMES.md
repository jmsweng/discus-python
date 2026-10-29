# Engine HDF5 volumes of the Gate 1 / Gate 2 evidence runs

The 20 DISCUS output volumes (~270 MB) of these evidence directories are **not in the working tree**:

- `gate1_c1_dyadic/evidence/20261028T234015Z_LAPTOP-NGU3BT0N`
- `gate2_prereg/evidence/20261028T234015Z_LAPTOP-NGU3BT0N` (v2)
- `gate2_prereg/evidence/20261028T235829Z_LAPTOP-NGU3BT0N` (v3)
- `gate2_prereg/evidence/20261029T010305Z_LAPTOP-NGU3BT0N` (v4)
- `gate3/evidence/20261029T105111Z_LAPTOP-NGU3BT0N` (Gate 3; never committed, hashes reported by the operator)

What stays committed in each directory: the exact `.stru` and `.mac` that were run, the engine logs, `ENGINE_ENV.txt`
(engine binary SHA-256, registered commit), `RUN_STATUS.txt`, the evaluator reports and `MANIFEST.sha256`. The manifest
still lists the original HDF5 file hashes as a record. The originals also remain in the git history of the PR #8 branch,
up to commit `d78448e`.

## Recreating them

On a machine with DISCUS_SUITE 6.23.01 (the runs used the binary with SHA-256 `d4274f13…45706f3`), with numpy and h5py:

```sh
python benchmarks/discus_6_23_01/regenerate_evidence.py list
python benchmarks/discus_6_23_01/regenerate_evidence.py regenerate \
    benchmarks/discus_6_23_01/gate2_prereg/evidence/20261029T010305Z_LAPTOP-NGU3BT0N  ~/regen_gate2_v4  --analyze
```

The script:
1. checks every input against `MANIFEST.sha256`;
2. runs the macros in the recorded order into a new directory (never into the evidence directory);
3. compares each volume with the **content hash** in `EVIDENCE_CONTENT_HASHES.json`.

That hash covers the intensity array and the grid metadata, not the date stamp that DISCUS writes into each file, so a
regenerated file can never match the original file hash. With `--analyze`, it also reruns the evaluator of the registered
commit, from a temporary git worktree, on the new volumes.

- **MATCH:** the data are identical to the originals.
- **DIFFERS:** expected with another binary, compiler or thread count. The registered evaluator's verdict and numbers then show whether the result is scientifically the same.

The script was tested end-to-end with a stand-in engine. Date-stamped copies MATCH, a single-voxel change of 1×10⁻¹³ is
detected, and the registered evaluators reproduce PASS. **Real-engine confirmation:** the Gate 3 run recomputed the Gate 2 v4 block total and diffuse volumes from the same inputs on a different day. Both content hashes MATCH, so the engine is reproducible and the date-independent hash works as designed.
