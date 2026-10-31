# Mantid 7.0.0 MDHisto fixtures

Generated on 2026-10-31 in the cloud with Mantid 7.0.0 (conda channel `mantid`, Python 3.13, an
isolated environment removed afterwards). Nothing here was written by hand.

| File | Made by | What it pins down |
|---|---|---|
| `aligned_hkl.nxs` | `BinMD` (axis-aligned) of an HKL MDEventWorkspace | two events at known hkl bin centres |
| `projected_hh0.nxs` | `BinMD` with basis `[H,H,0]`, `[-H,H,0]`, `[0,0,L]` | the `transform_to_orig` convention |
| `w_matrix_before.nxs`, `w_matrix_after.nxs` | `ConvertToMD` (Uproj 1,1,0; Vproj 0,0,1; Wproj 1,-1,0) + `BinMD`, one extra event between them | the `W_MATRIX` convention and `in <L> A^-1` units |
| `histo_formula.nxs` | `CreateMDHistoWorkspace` with signal 100·i0 + 10·i1 + i2 + 1, errors 0.1·signal, `MaskMD` on the first H bin, one NaN | axis order, errors², mask |
| `expected.json` | `make_fixtures.py`, written before the files | all expected positions and values |
| `roundtrip/*.resaved.nxs` | discuspdf `write_mantid_mdhisto` → Mantid `LoadMD` → Mantid `SaveMD` | Mantid reads our files and keeps our projection |
| `roundtrip/mantid_check.json` | `check_in_mantid.py` inside Mantid | dimensions, signal, errors², mask and lattice as Mantid sees them (no failures) |

Scripts: `benchmarks/mantid_7_0_0/make_fixtures.py` and `check_in_mantid.py` (Mantid's Python),
`write_for_mantid.py` (discuspdf's Python).
