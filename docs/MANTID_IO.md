# Mantid MDHisto files

`read_mantid_mdhisto` reads a binned single-crystal volume saved by Mantid, and
`write_mantid_mdhisto` writes one that Mantid can load. Mantid does not need to be installed for
either. Use them to:

- compare measured diffuse scattering with a DISCUS calculation on the same hkl points;
- send calculated or processed data back to Mantid for viewing or further reduction.

The files are what Mantid's `SaveMD` writes for an **MDHistoWorkspace** (the binned result of
`BinMD`, `MDNorm`, `CreateMDHistoWorkspace`, …) in the **HKL frame**.

## Quick start

```python
from discuspdf import read_mantid_mdhisto, write_mantid_mdhisto

sv = read_mantid_mdhisto("binned.nxs")     # a ScatteringVolume
sv.grid.hkl()                              # hkl of every bin centre, projection applied
sv.variance                                # Mantid's errors squared
sv.mask                                    # Mantid's mask plus bins with no data (NaN); None if nothing is masked
sv.grid.lattice.parameters                 # a, b, c, α, β, γ from the file

write_mantid_mdhisto(sv, "for_mantid.nxs")  # LoadMD in Mantid gives the same axes, values, mask and lattice
```

The result is an ordinary [`ScatteringVolume`](SCATTERING_DATA.md), so slicing, plotting, xarray and
HDF5 all work the same as for DISCUS data.

## What is read

- **Axes, including projections.** Workspaces binned along directions such as `[H,H,0]`,
  `[-H,H,0]` and `[0,0,L]` (`BinMD` projections, or `ConvertToMD` with `Uproj/Vproj/Wproj`) give a
  grid whose points have the correct hkl.
- **Values, errors and mask.** Errors² become `variance`. Mantid's mask is kept with the reason
  `Mantid mask`. Bins whose signal is NaN (no coverage, e.g. from `MDNorm`) are masked with the
  reason `NaN in Mantid signal (no coverage)`.
- **Lattice.** Built from a, b, c, α, β, γ in the standard orientation (a along x). Mantid's UB and
  U matrices are kept in `sv.provenance` but not applied.
- **Extra axes with one bin** (e.g. an energy-transfer window) are dropped and their range is
  recorded in the provenance.

As a safety check, each axis's direction is worked out three ways (from the stored transform,
from the axis name such as `[H,H,0]`, and from its length in Å⁻¹ in the units), and the file is
refused if they disagree.

## What is refused

The reader stops with a `MantidFormatError` that says why, rather than guessing, for:

| File | What to do |
|---|---|
| Q_lab or Q_sample frame | convert to HKL in Mantid first (applying UB and goniometer here is not supported) |
| not exactly three HKL axes | bin to a 3D hkl volume |
| an extra axis (e.g. `DeltaE`) with more than one bin | integrate over it in Mantid |
| uneven bins, unknown units, no lattice | re-bin, or set the lattice (`SetUB`) in Mantid |
| `SaveMD` version 1, or an MDEventWorkspace | bin to an MDHistoWorkspace and save with the current `SaveMD` |

## What is written

`write_mantid_mdhisto` writes each grid axis as an HKL axis named in Mantid's style (`[H,H,0]`), in
`r.l.u.`, with the projection stored so that Mantid places every bin at the right hkl. It also
writes errors² from `variance` (zeros if there is none), masked bins as 0 with the mask set, and
the lattice. Planes and lines are refused: a one-point axis has no direction to store.

## Reference: Mantid's file layout

For anyone reading these files by hand. These are observed on files written by Mantid 7.0.0.

| Item | Meaning |
|---|---|
| `data/D0 … D{n-1}` | bin **edges** of each dimension; attributes `frame`, `long_name` (e.g. `[H,H,0]`), `units` |
| `data/signal`, attribute `axes = "D2:D1:D0"` | the array's first axis is the **last** dimension |
| `errors_squared`, `mask` (1 = masked), `num_events` | same shape as `signal` |
| `transform_to_orig` (4×4, row-major) | its **columns** are the projection directions |
| `experiment0/logs/W_MATRIX` (3×3, row-major) | from `ConvertToMD`; its **columns** are Uproj, Vproj, Wproj |
| units `in 1.67 A^-1` | coordinates along that axis are in units of 1.67 Å⁻¹, which is 2π × \|direction\| (here `[H,H,0]` in a 5, 6, 7 Å, β = 105° cell) |
| `experiment0/sample/oriented_lattice/orientation_matrix` | the UB matrix |

So hkl = W_MATRIX · (transform · x + offset) for bin-centre coordinates x.

<details>
<summary>How this was checked</summary>

Test files in `tests/fixtures/mantid_7_0_0/` were written by **Mantid 7.0.0** from
`benchmarks/mantid_7_0_0/make_fixtures.py`. The expected position and value of every event were
written down in `expected.json` before Mantid ran.

- Events placed at known hkl (plain axes, `[H,H,0]`/`[-H,H,0]`/`[0,0,L]` projections, and a
  `ConvertToMD` projection) land on the grid point with that hkl, with the right weight and error.
- A workspace with values 100·i0 + 10·i1 + i2 + 1, a masked box and a NaN bin reads back with the
  right axis order, errors, mask and NaN handling.
- **Inside Mantid**: our written files load in `LoadMD` with the same axes, values, errors, mask and
  lattice as Mantid's own (report: `roundtrip/mantid_check.json`, no failures). A skewed monoclinic
  volume Mantid had never seen loads with the expected values. Mantid then re-saved each file, and
  reading those gives back the same grid and values.
- Two older third-party files (an inelastic file with 48 energy bins, and a file with two HKL axes)
  are refused with the reasons above.
- Deliberate bugs (transposed matrices, wrong axis order, disabled name/units checks, NaN masking,
  writer axis order and flip, written UB) are each caught.
</details>
