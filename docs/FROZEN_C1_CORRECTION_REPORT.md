# FROZEN C1 CORRECTION REPORT

**Date:** 2026-10-28
**Environment:** WSL2 Ubuntu 22.04, DISCUS_SUITE 6.23.01 (`/usr/local/bin/discus_suite`, SHA-256 `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3`)
**Run directory:** `~/discus_3d_frozen_corrected_20261028T022435Z/` (WSL2)

## 1. Objective

Correct the frozen-structure loading procedure (change `read cell` → `read stru`) and demonstrate that scattering from the frozen C1 model reproduces the original C1 data.

## 2. Root cause

The original frozen scatter macro used:

```
read
  cell C1_frozen.stru
```

Per the DISCUS 6.23.01 help file (`/usr/local/share/suite.hlp`):

- **`cell <name>`** — "Reads the contents of a **single asymmetric unit** from file <name> and then generates the unit cell."
- **`stru <name>`** — "Reads a **whole structure** from file <name>."

The 4096-atom saved `.stru` file was misinterpreted as a unit cell, causing DISCUS to expand it under F-43m symmetry (96 operations), producing 8192 "Atom typ" entries (each with 4 atoms) and an all-zero diffuse output.

## 3. Corrected macro

```
discus
read
  stru C1_frozen.stru
fourier
  xray
  wvle moa1
  ll -4.0,-4.0,-4.0
  lr  4.0,-4.0,-4.0
  ul -4.0, 4.0,-4.0
  tl -4.0,-4.0, 4.0
  na 65
  no 65
  nt 65
  abs h
  ord k
  top l
  temp use
  disp off
  lots off
  set technique:turbo
  set aver,0.0
  show
  run
exit
output
  value int
  form nexus
  outf C1_corrected_total.h5
  run
exit
fourier
  set aver,100.0
  show
  run
exit
output
  value int
  form nexus
  outf C1_corrected_diffuse.h5
  run
exit
exit
exit
```

The only change from the original frozen macro is `cell` → `stru`. All Fourier grid, X-ray form factors, temperature factors, and `set aver` settings are retained. No random substitution is performed.

## 4. Engine log verification

All 3 corrected runs report:

```
Lot #    1/   1 : complete crystal
           Atom typ = GA         (#      1723 )
           Atom typ = AS         (#      2373 )
```

| Metric | Original C1 | Incorrect frozen (`cell`) | Corrected frozen (`stru`) |
|---|---|---|---|
| `Atom typ` lines | 4 | 8192 | **4** |
| GA atoms | 1723 | 13784 (8× overcount) | **1723** |
| AS atoms | 2373 | 18984 (8× overcount) | **2373** |
| Total atoms | 4096 | 32768 (8× overcount) | **4096** |
| Lot description | — | — | **"complete crystal"** |

The corrected run reads the 4096-atom structure as a single lot ("complete crystal"), matching the original C1's atom counts exactly.

## 5. HDF5 numerical comparison

### 5.1 Corrected vs original C1

| Metric | Total | Diffuse |
|---|---|---|
| Shape | (65, 65, 65) | (65, 65, 65) |
| Max (corrected) | 1.733590×10¹⁰ | 7545.1438 |
| Max (original) | 1.733590×10¹⁰ | 7545.1438 |
| Max abs diff | **0.0** | **0.0** |
| Mean abs diff | **0.0** | **0.0** |
| Relative L2 error | **0.0** | **0.0** |
| Pearson correlation | **1.0** | **1.0** |
| Bit-identical | **True** | **True** |

Both the total and Bragg-subtracted diffuse arrays are **bit-identical** to the original C1. The maximum absolute difference is 0.0 across all 274,625 voxels.

### 5.2 Reproducibility: 3 corrected runs

| Comparison | Total | Diffuse |
|---|---|---|
| Run 1 vs Run 2 | bit-identical | bit-identical |
| Run 1 vs Run 3 | bit-identical | bit-identical |
| Run 2 vs Run 3 | bit-identical | bit-identical |

All 3 runs produce bit-identical numerical data. The HDF5 file checksums differ due to the `audit_creation_date` metadata field, but the `data_values` arrays are identical.

## 6. Checksums

| File | SHA-256 |
|---|---|
| `corrected/C1_frozen.stru` | `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9` |
| `corrected/C1_corrected_total_run1.h5` | `ac5788e4c3619a0fc8aca7625972a87ef83830bca2b1c83df9821af203371291` |
| `corrected/C1_corrected_diffuse_run1.h5` | `c2287bae7ab5e119f09ce07299e3d71475bc7f731592c55bdab59412855be75f` |
| `corrected/C1_corrected_total_run2.h5` | `bb15c3082634682d2b19719477aa38081b0db79de1f18b88e4aeba2a86c61aa1` |
| `corrected/C1_corrected_diffuse_run2.h5` | `28493c49f43c8d49069a4087ff47bbed6261b4c10c15e7c60f72a20f1f258ff5` |
| `corrected/C1_corrected_total_run3.h5` | `cefd4420bcd7a46b2be8beaef788cca97203e115c8c43752e2b91e8b1a9196db` |
| `corrected/C1_corrected_diffuse_run3.h5` | `4e46045400c969d6881413249d21eb7c41924b8372acda87c1d0d3092bc23297` |
| `original/C1_total.h5` | `da7e74553d4e03c90e56ad7c2e7bc059a044894b2be5bb03b598d14e55d06666` |
| `original/C1_diffuse.h5` | `ef3c131db731387a9a4653da6f28c53be0d027dd7638be1a14f5ec0ded927ee3` |

Note: The HDF5 file checksums differ between runs due to the `audit_creation_date` metadata field. The `data_values` arrays are bit-identical (verified with `np.array_equal`).

## 7. Pass/fail conclusions

| Check | Result |
|---|---|
| Exactly 4096 atoms represented | **PASS** (1723 GA + 2373 AS) |
| Composition matches frozen `.stru` model | **PASS** (identical SHA-256) |
| Atomic coordinates match frozen model | **PASS** (single lot, "complete crystal") |
| Total intensity nonzero and consistent with original C1 | **PASS** (bit-identical) |
| Diffuse intensity nonzero and consistent with original C1 | **PASS** (bit-identical, max 7545.14) |
| 3 repeated runs produce identical numerical data | **PASS** (bit-identical) |
| Corrected results match original C1 physical output | **PASS** (max abs diff = 0.0) |

**Overall: PASS.** The corrected frozen-structure workflow (`read stru`) reproduces the original C1 scattering data to numerical precision. The frozen C1 model is now a valid benchmark for independent scientific validation.

## 8. Next step: X-ray form-factor table benchmark

With the frozen C1 structure validated, the next step is to benchmark the X-ray form-factor tables using the same frozen atomic structure for both calculations:

1. Run the corrected frozen macro with the default form-factor table (Waasmaier–Kirfel).
2. Run the corrected frozen macro with `xray table:inter` (interpolated table).
3. Compare the total and diffuse intensities to quantify the isolated form-factor table effect.

This addresses P1 item 2 in the independent review's priority handoff.

## 9. Status

- **Frozen C1 correction:** **COMPLETE** (validated)
- **Form-factor table benchmark:** **PENDING** (next step)
- **Independent 3D-ΔPDF benchmark:** **NOT COMPLETE** (separate milestone)
- **PR #7:** **NOT MERGED** (per instructions)
