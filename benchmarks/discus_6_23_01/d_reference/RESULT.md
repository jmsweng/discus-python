# D reference probe — result: PASS (6/6 pre-registered criteria, v2)

**`discuspdf.reference_scattering` reproduces DISCUS 6.23.01 single-crystal intensities to 2.5 × 10⁻¹⁵ …
1.8 × 10⁻¹³ relative L2 on the commensurate grid (X-ray and neutron, total and `set aver,100`, isotropic B),
once DISCUS's own table approximations are emulated.**

| | |
|---|---|
| Registration | v1 `e74bf15` (inputs, macros, criteria); v2 `b7716f0` — amended after the engine run but before the overseer read any number (see `PREREG.md`, amendment log) |
| Engine run | owner's laptop, WSL Ubuntu-22.04, `discus_suite` SHA-256 `d4274f13…45706f3`, 2026-11-01 22:17 UTC, local agent task 001 |
| Evidence | `evidence/20261101T221751Z_LAPTOP-NGU3BT0N/` (logs, `engine.txt`, exit codes, input/output hashes, the 8 NeXus files, the v1 analysis as `analysis_v1.json`) |
| Integrity | 8/8 exit codes 0; input hashes = committed inputs; Bragg metadata and shapes as requested |
| Re-analysed | independently in the cloud sandbox (Python 3.13) from the committed evidence: identical numbers |

## Numbers (relative L2 of DISCUS − reference)

| Run | physics (exact s) | engine emulation | Criterion |
|---|---|---|---|
| `cuau_x_tot_c4` | 7.0e-4 | 2.5e-15 | ≤ 1e-9 PASS |
| `cuau_x_dif_c4` | 8.5e-4 | 5.1e-15 | ≤ 1e-9 PASS |
| `cuau_n_tot_c4` | 1.3e-4 | 4.1e-15 | ≤ 1e-9 PASS |
| `cuau_n_dif_c4` | 2.4e-3 | 1.8e-13 | ≤ 1e-9 PASS |
| `mono_x_tot_m` (monoclinic, displaced, step 0.2) | 1.1e-3 | 9.1e-6 | ≤ 2e-3 PASS |
| `mono_n_tot_m` | 2.1e-4 | 9.2e-6 | ≤ 2e-3 PASS |
| `cuau_x_tot_nc` (step 0.3, not commensurate) | 5.5e-4 | 1.9e-5 | informational |
| `cuau_x_dif_nc` | 7.2e-4 | 1.2e-4 | informational |

## Findings

1. **DISCUS's scattering factors are tabulated at s rounded to 0.001 Å⁻¹, Debye–Waller factor included.** For
   isotropic B, `four_formtab` stores `exp(−B s_r²)·f(s_r)`; the exact per-voxel `four_dbwtab` is used only for
   anisotropic ADPs. This costs 1e-4 … 2e-3 relative (largest for the neutron diffuse intensity, where the
   small, strongly s-dependent difference of two factors matters). v1 of the criteria assumed an exact DW and
   failed by exactly this amount (`analysis_v1.json`), as the review had predicted from the source.
2. **On a commensurate (dyadic) grid the agreement is at machine precision**, including `set aver,100`
   (amplitude subtraction of the cell-averaged structure, vacancies none) — the reference's `⟨F⟩` is DISCUS's.
3. **Non-orthogonal geometry is right.** The monoclinic, displaced crystal agrees to 9 × 10⁻⁶ on a
   non-commensurate 0.2 step; the residual is DISCUS's turbo phase-table approximation on non-dyadic steps
   (1.9e-5 for the cubic crystal on step 0.3, 1.2e-4 after `set aver`, in line with the ≈1e-3-of-max level
   seen before). An orthogonal-metric or Ga/As-swap error would have given 2e-2 / 6e-3 (measured before the run).
4. **Tables.** DISCUS uses b(Au) = 7.63 fm (Sears 1992); `periodictable` 2.1 gives 7.90 fm. The reference
   always states which table it uses.

Reproduce: `python benchmarks/discus_6_23_01/d_reference/analyze.py benchmarks/discus_6_23_01/d_reference/evidence/20261101T221751Z_LAPTOP-NGU3BT0N`
(also run by `tests/test_d_reference_probe.py`).
