# Builder-backend probe v2: VOID atoms in the 1D powder PDF (pre-registered)

Registered **before** the engine run, after v1 (`../RESULT.md`). `inputs/`, `runs.json`, `criteria.json`,
`make_inputs.py` and `analyze.py` are committed first and not changed afterwards. A failed criterion is
reported as a finding.

## Background
In v1, adding one VOID atom per cell (320 rows, 256 atoms) multiplied the powder G(r) by exactly 1.25 = 320/256
at every r (P3), while the 3D intensity was unchanged (E4). The disorder builder pads cells with VOID whenever
a site's occupants have different atom counts, so its 1D PDFs of such models are scaled by rows / atoms.

## Questions
1. **One-cell route.** Is the crystal written as a single cell of the whole supercell, with only the real atoms
   (`ncell 1,1,1,N`, no VOID rows), equivalent to the verified layout for the powder PDF? If so, the builder
   writes such a file for the 1D run.
2. **Scale law.** Does G_padded = (rows / atoms) × G_unpadded hold for another VOID fraction (vacancies:
   256 rows, 240 atoms) and for X-rays?

## Inputs (`make_inputs.py`; fcc, a = 3.52 Å, 4×4×4 cells, B = 1; PDF settings as v1)
| Structure | Rows / atoms | Layout |
|---|---|---|
| `mix` | 256 / 256 | v1's `mix` (64 Cu among Ni), `ncell 4,4,4,4` (byte-identical to v1's file) |
| `mix_big` | 256 / 256 | the same atoms as one 14.08 Å cell, `ncell 1,1,1,256` |
| `void` | 320 / 256 | v1's `void` (one VOID per cell), `ncell 4,4,4,5` |
| `vac` | 256 / 240 | `mix` with 16 rows renamed VOID (seed 20261104), `ncell 4,4,4,4` |
| `vac_big` | 240 / 240 | the 240 real atoms of `vac` as one cell, `ncell 1,1,1,240` |

Runs: `p_mix_n`, `p_mix_x`, `p_mix_big_n`, `p_mix_big_x`, `p_void_x`, `p_vac_void_n`, `p_vac_void_x`,
`p_vac_big_n`, `p_vac_big_x` (n = neutron, x = X-ray), run with `../run_probe.sh`.

## Criteria (`criteria.json`; metric max |ΔG| / max |G|)
| Key | Comparison | Pass if |
|---|---|---|
| E0 | every run exit 0, committed inputs used unchanged, logs present | no problem |
| V0 | r = 0…10 step 0.02, finite G, all runs | 0 failures |
| V1 | `p_mix_big_n` vs `p_mix_n` | ≤ 1e-6 |
| V2 | `p_mix_big_x` vs `p_mix_x` | ≤ 1e-6 |
| V3 | (240/256)·`p_vac_void_n` vs `p_vac_big_n` | ≤ 1e-6 |
| V4 | (240/256)·`p_vac_void_x` vs `p_vac_big_x` | ≤ 1e-6 |
| V5 | (256/320)·`p_void_x` vs `p_mix_x` | ≤ 1e-6 |
| V6 | `p_vac_big_n` vs `p_mix_n` (power: removing 16 atoms changes G; estimate ≈ 6e-2) | ≥ 1e-3 |
| I0 | engine-marker lines in the logs | informational |
| I1 | `p_mix_n` byte-identical to v1's `p_mix` (same files, same binary) | informational |

The limits follow v1: the pair sums are identical, so only summation order and the exact scale factor can
differ (v1 measured 3.5e-16 for the factor 1.25).

## What the outcomes mean for the builder
* **V1, V2 pass:** the builder writes the real atoms as one cell for the 1D PDF (no VOID rows) and keeps the
  padded file for the 3D runs. This is the planned fix.
* **V3 also passes:** the one-cell route is right with vacancies too, which is the builder's actual case (real
  atoms per cell vary). V3 failing while V1 passes would mean the one-cell file and the scale law disagree when
  atoms are missing; then neither is offered until that is understood.
* **V3–V5 pass:** the scale law holds generally, which explains v1 P3 and makes rescaling a fallback.
* **V1 or V2 fail:** the `ncell` layout itself changes the powder PDF; then the builder rescales by atoms / rows
  instead, if V3–V5 pass, and the 1D PDF of padded models is otherwise withdrawn.
