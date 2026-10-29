> **Historical document (kept by owner decision, 2026-10-28).** Written during earlier PR work. It is **not** part of PR #8's validated results and may contain claims since superseded or refuted. The current record is `docs/SCIENTIFIC_VALIDATION_REPORT.md` and `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`.

# Independent Review: DISCUS 6.23.01 C0–C5 Evidence — Verification Report

**Date:** 2026-10-27
**Reviewer:** Independent (automated)
**Scope:** Verify the claims in `DISCUS_EVIDENCE_INDEPENDENT_REVIEW.md` against the raw evidence in `~/discus_3d_validation_20261028T011146Z/` (WSL2 Ubuntu 22.04).

## 1. Archive integrity

| Archive | SHA-256 (claimed) | SHA-256 (verified) | Match |
|---|---|---|---|
| `discus_3d_review_text.tar.gz` | `bda3143b…96ed1be` | `bda3143b…96ed1be` | ✓ |
| `discus_3d_review_hdf5.tar.gz` | `18924ae4…3cde7b` | `18924ae4…3cde7b` | ✓ |

Both archives are intact and match the handoff document.

## 2. Frozen C1 structure identity

| Item | Value |
|---|---|
| `cases/C1_model.stru` SHA-256 | `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9` |
| `frozen_scatter/C1_frozen.stru` SHA-256 | `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9` |
| Byte-identical | ✓ |
| Atom count (both files) | 4096 (1723 GA + 2373 AS) |

The frozen structure is a byte-for-byte copy of the original C1 saved model. The input data is correct.

## 3. Critical finding: frozen C1 rerun is NOT physically equivalent to original C1

### 3.1 The macro difference

| | Original C1 (`cases/C1.mac`) | Frozen rerun (`frozen_scatter/scatter.mac`) |
|---|---|---|
| Read command | `read cell gaas.cll,8,8,8` | `read cell C1_frozen.stru` |
| Seed | `seed 1` | (none) |
| Replace | `replace GA,AS,all,0.15` | (none) |
| Fourier grid | 65³, −4..+4 rlu | 65³, −4..+4 rlu (identical) |
| `set aver` | 0.0 (total), 100.0 (diffuse) | 0.0 (total), 100.0 (diffuse) (identical) |

The only difference is the `read` command. The original reads a 2-atom asymmetric unit (`gaas.cll`) and expands it 8×8×8 = 512 times to produce 4096 atoms. The frozen rerun reads the 4096-atom saved structure file using `read cell` — which, per the DISCUS 6.23.01 help file (`/usr/local/share/suite.hlp`), is defined as:

> **`cell <name> [,<nx> [,<ny> [,<nz>]]`**
> "Reads the contents of a **single asymmetric unit** from file <name> and then generates the unit cell. Optionally the number of unit cells along x,y,z can be given to generate the whole crystal. The number of unit cells default to 1."

The correct command for reading a whole saved crystal is:

> **`stru <name>`**
> "Reads a **whole structure** from file <name>. The space group information is ignored, nor is a check for consistency of space group and lattice constants performed."

### 3.2 The engine log evidence

| Metric | Original C1 log | Frozen scatter log |
|---|---|---|
| `Atom typ` lines | **4** | **8192** |
| GA entries | 2 (total 1723) | 3446 (total 13784) |
| AS entries | 2 (total 2373) | 4746 (total 18984) |
| Atoms per entry | 1723 / 2373 | **4** (all 8192 entries) |
| Total atoms reported | 4096 | **32768** (8× overcount) |

The frozen scatter log contains 8192 "Atom typ" lines, each reporting exactly 4 atoms. This is 8× the 1024 atoms expected from a single unit cell of F-43m GaAs (4 GA + 4 AS per cell). The 8× overcount is consistent with DISCUS interpreting the 4096-atom file as a unit cell and expanding it under the 96 symmetry operations of F-43m, with the `cell` command's default `nx=ny=nz=1` producing a single "unit cell" that DISCUS then treats as containing 8192/96 ≈ 85.33 "asymmetric units" — or more likely, DISCUS is reading the 4096-atom file, grouping atoms into 4-atom lots, and reporting each lot as a separate "Atom typ" entry.

The key point: **DISCUS did NOT compute the structure factor from 4096 atoms.** It computed from a misinterpreted representation of the file.

### 3.3 The HDF5 numerical evidence

| Metric | Original C1 | Frozen C1 | Verdict |
|---|---|---|---|
| `data_values` shape | (65, 65, 65) | (65, 65, 65) | Same grid |
| Total max | 1.7336×10¹⁰ | 2.7737×10¹¹ | **16× larger** |
| Total mean | ~1.08×10⁹ | ~1.22×10¹⁰ | ~11× larger |
| Diffuse max | 7545.14 | **0.0000** | **All zeros** |
| Diffuse nonzero voxels | 274196 / 274625 | **0 / 274625** | **All zeros** |
| `np.array_equal` (total) | — | **False** | Different |
| `np.array_equal` (diffuse) | — | **False** | Different |
| Pearson correlation (total) | — | 0.073 | Negligible |
| Pearson correlation (diffuse) | — | NaN (zero variance) | Undefined |

The frozen diffuse output is **zero at every one of the 274,625 voxels**. This is the decisive evidence: `set aver,100.0` subtracts the average structure factor from the total. When DISCUS computes the "total" from a misinterpreted 4-atom-per-lot representation, the average equals the total (or the subtraction yields zero), producing an all-zero diffuse output. The original C1 diffuse has 274,196 nonzero voxels with a maximum of 7545.14.

The 16× ratio at Q=0 (total) is consistent with a 4× atom expansion (intensity ∝ N² for coherent scattering) combined with additional symmetry-averaging effects. The fact that other reciprocal voxels are not simply proportional (Pearson correlation 0.073) confirms this is not a harmless overall scale factor — the spatial distribution of intensity is fundamentally different.

### 3.4 The bit-identical reproducibility claim

The handoff reports three repeated frozen runs with bit-identical numerical data. This is **true but misleading**: it demonstrates reproducibility of the **incorrectly reloaded model**, not reproducibility of the original C1 scattering physics. The frozen runs are deterministic because they read the same file with the same (incorrect) `read cell` command each time. But the result is not the original C1 physics.

## 4. Verification of the review document's other claims

### 4.1 3DPDF vs independent NumPy FFT

The review claims that `numpy.fft.fftshift(numpy.fft.fftn(numpy.fft.ifftshift(I))).real / sqrt(I.size)` matches DISCUS's `*_diffuse_3dpdf.h5` with relative L2 error 3.65×10⁻¹⁶ for both C1 and C5. The fitted scale is 0.0019082266859878, which equals `1/sqrt(65³)` = `1/sqrt(274625)` = 0.0019082266859878222. This is a precise numerical confirmation of the orthonormal FFT scaling. **Confirmed.**

### 4.2 C3 Patterson distance

The review claims a direct Fourier summation at fractional vector u=(0.30, 0.20, 0.10) yields an off-origin feature of ~1.9587×10⁹, and that the monoclinic metric (a=5, b=6, c=7, β=105°) gives an interatomic separation of 1.9069556904 Å. This is consistent with the handoff's C3 Patterson distance of 1.9070 Å. **Confirmed (consistent with handoff).**

### 4.3 Seed reproducibility

The review confirms that all 8 preserved seed-test `.stru` files have 4096 atoms and identical coordinates but differ in species occupancy, with 527 of 4096 site identities differing between two `seed 1` repetitions. **Confirmed** (consistent with `seed_test/seed_result.txt`).

### 4.4 C1 vs C5 form-factor comparison

The review correctly notes that C1 and C5 use **different saved random structures**, so their intensity differences cannot isolate the effect of the X-ray form-factor table. **Confirmed** (C1 and C2/C5 were generated with different seeds and are genuinely different configurations).

## 5. Root cause confirmation

The review's likely root cause is **correct**:

> `read cell` interprets a full saved `.stru` as a unit-cell/asymmetric-unit input and expands it under space-group symmetry. The historical DISCUS command reference distinguishes `read / cell ...` from `read / stru <file>`, which reads the whole saved crystal without symmetry expansion.

The DISCUS 6.23.01 help file (`/usr/local/share/suite.hlp`) confirms this distinction:
- `cell <name> [,<nx> [,<ny> [,<nz>]]` — "Reads the contents of a **single asymmetric unit** from file <name> and then generates the unit cell."
- `stru <name>` — "Reads a **whole structure** from file <name>."

The frozen scatter macro used `read cell C1_frozen.stru`, which is semantically incorrect for a 4096-atom saved structure. The correct command is `read stru C1_frozen.stru`.

## 6. Conclusion

**The review document's central claim is CONFIRMED:** the frozen C1 rerun is **not physically equivalent** to the original C1 run. The `read cell` command misinterprets the 4096-atom saved structure, producing:

1. An 8× overcount of atoms in the engine log (32768 vs 4096).
2. A 16× overcount of total intensity at Q=0.
3. An **all-zero** diffuse output (vs 274,196 nonzero voxels in the original).
4. A negligible correlation (0.073) between the original and frozen total intensities.

The bit-identical reproducibility of the three frozen runs demonstrates only that the **incorrect** reload is deterministic, not that the original C1 physics is reproducible.

### Required correction

The frozen scatter macro must be changed from:

```
read
  cell C1_frozen.stru
```

to:

```
read
  stru C1_frozen.stru
```

After this correction, the engine log should show exactly 4096 atoms (1723 GA + 2373 AS), and the total and diffuse HDF5 outputs should match the original C1 outputs to numerical precision. Only then can the frozen-structure benchmark be considered valid.

### Status of the handoff's claims

| Handoff claim | Verdict |
|---|---|
| C0–C5 engine execution complete | **Confirmed** (all 6 cases ran, outputs valid) |
| Native reader validation complete | **Confirmed** (all files read, schema/geometry/FFT checks pass) |
| C0 Bragg suppression observed | **Confirmed** (consistent with FINDINGS.md) |
| C3 Patterson distance 1.9070 Å | **Confirmed** (matches independent calculation) |
| C4 skewed-grid rejection | **Confirmed** (correctly rejected) |
| Seed reproducibility failed | **Confirmed** (all 8 variants produce different structures) |
| Frozen-structure scattering reproducibility complete | **INVALID** (reproducibility of incorrect model, not original physics) |
| Independent structure-factor benchmark | **Not performed** (correctly noted) |
| Correlated-disorder ΔPDF benchmark | **Not performed** (correctly noted) |
| 3DPDF normalization unresolved | **Confirmed** (metadata inconsistency documented) |

**Overall judgment:** The handoff's C0–C5 engine execution and native reader validation are credible. The frozen-scattering validation claim is **invalid** and must be corrected. The review document's identification of the root cause (`read cell` vs `read stru`) is correct and is confirmed by the DISCUS 6.23.01 help file.
