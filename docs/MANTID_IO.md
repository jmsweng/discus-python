# Mantid MDHisto files (Package B, Feature 18)

`read_mantid_mdhisto(path)` reads a Mantid `SaveMD` (version 2) **MDHistoWorkspace** in the **HKL
frame** into a `ScatteringVolume`; `write_mantid_mdhisto(sv, path)` writes one that Mantid's `LoadMD`
reads. No Mantid installation is needed for either.

```python
from discuspdf import read_mantid_mdhisto, write_mantid_mdhisto
sv = read_mantid_mdhisto("binned.nxs")       # grid in hkl, variance = errors², mask, provenance
sv.grid.axis_vectors                         # hkl step of each Mantid dimension (projection applied)
write_mantid_mdhisto(sv, "for_mantid.nxs")   # LoadMD in Mantid gives the same axes, values, mask, lattice
```

## File conventions (observed on files written by Mantid 7.0.0)

| Item | Meaning |
|---|---|
| `data/D0…D{n-1}` | bin **edges** of each dimension; attrs `frame` (`HKL`, `General Frame`, …), `long_name` (`[H,H,0]`), `units` (`r.l.u.` or `in <L> A^-1`) |
| `data/signal` attr `axes = "D{n-1}:…:D0"` | the first array axis is the **last** dimension (Mantid's linear index runs D0 fastest) |
| `errors_squared`, `mask` (int8, 1 = masked), `num_events` | same shape as `signal` |
| `transform_to_orig` (16 × float32, row-major 4×4) | original coordinates = M·[x, 1]; **columns** are the projection basis vectors (`BinMD` `[H,H,0]` → column (1,1,0)) |
| `experiment0/logs/W_MATRIX/value` (row-major 3×3, from `ConvertToMD`) | **columns** are Uproj, Vproj, Wproj; hkl = W·orig (identity when absent) |
| `units = "in <L> A^-1"` | L = 2π·\|basis vector\| in Å⁻¹ (e.g. 1.67 for [H,H,0] in a 5/6/7 Å, β = 105° cell) |
| `experiment0/sample/oriented_lattice/orientation_matrix` | **UB** (UBᵀ·UB = reciprocal metric; B in Busing–Levy form) |

So for an HKL-frame workspace, **hkl = W · (M · x + t)** for bin-centre coordinates x.

## Cross-checks and refusals

For every axis the reader checks that three independent records agree, and refuses the file
(`MantidFormatError`) if they do not:
1. the basis vector from `transform_to_orig` and `W_MATRIX`;
2. the axis name, when it has the `[aH,bK,cL]` form;
3. the Å⁻¹ length in `units`, when given (3 significant digits, as Mantid prints it).

Also refused, with the reason: `SaveMDVersion` other than 2; Q_lab / Q_sample frames (convert to HKL
in Mantid; this reader does not apply UB or goniometer); fewer or more than three HKL axes; extra
axes such as `DeltaE` with more than one bin (one-bin axes are dropped and their range recorded);
non-uniform bins; HKL axes coupled to non-HKL axes; units other than `r.l.u.` / `in <L> A^-1`; no
oriented lattice. NaN signal (no coverage, e.g. from `MDNorm`) is **masked** with reason
`NaN in Mantid signal (no coverage)`; Mantid's own mask keeps reason `Mantid mask`.

The lattice is rebuilt from a, b, c, α, β, γ in this package's canonical orientation (a ∥ x);
Mantid's UB and U are kept in `provenance` but not applied (Q in the sample frame needs the goniometer
too).

## Writer

Each grid axis becomes an HKL dimension whose basis vector is the grid step scaled so its largest
component is ±1 (a decreasing axis is written with a flipped basis vector and increasing edges);
names follow Mantid's `[H,H,0]` style and units are `r.l.u.`; `transform_to_orig` holds the basis.
Variance → `errors_squared` (zeros if absent), masked voxels are written as 0 with `mask = 1`. The
`experiment0` skeleton is the one Mantid writes (`LoadMD` needs it), with UB = Busing–Levy B
(U = identity). Planes and lines are refused (a one-point axis has no direction to store). No
`W_MATRIX` is written: the projection is entirely in `transform_to_orig`.

## Evidence

Fixtures in `tests/fixtures/mantid_7_0_0/` (see its README) were written by **Mantid 7.0.0** from
`benchmarks/mantid_7_0_0/make_fixtures.py`; every expected position and value was fixed in
`expected.json` before Mantid ran.

- Events placed at known hkl (aligned `BinMD`, projected `BinMD` with `[H,H,0]`, `[-H,H,0]`, `[0,0,L]`)
  land on the grid point with that hkl, with the right weight and variance.
- `ConvertToMD` with Uproj (1,1,0), Vproj (0,0,1), Wproj (1,−1,0): the single bin that changes when
  one event is added sits at hkl = W·x, as predicted.
- `CreateMDHistoWorkspace` with signal 100·i0 + 10·i1 + i2 + 1: axis order, errors², the masked box
  and the NaN bin come out as specified.
- **Writer, checked inside Mantid** (`check_in_mantid.py`, report `roundtrip/mantid_check.json`, no
  failures): for each fixture, the file we write loads in `LoadMD` with the same dimension names,
  frames, bins and extents, signal and errors², mask (plus the NaN bin, by design) and lattice as
  Mantid's own file; a skewed monoclinic volume Mantid never saw loads with the expected values and
  mask. Mantid then saved each of our workspaces again (`roundtrip/*.resaved.nxs`), and reading those
  gives back the same grid (to the float32 precision of the stored transform, ≤ 6e-8) and values.
- Two older real files: Shiver's 2023 file (inelastic, `DeltaE` with 48 bins) and Javelin's test file
  (two HKL axes, no transform) are refused with those reasons. They are not committed (GPL-3 /
  separate projects).
- Mutation checks: `W_MATRIX` transposed, transform rows/columns, axis order, name check, units check,
  NaN masking, writer axis order, decreasing-axis flip, written UB — each caught.

## Not included

MDEventWorkspace files, SaveMD version 1, Q_lab / Q_sample frames, applying UB/goniometer to get Q in
Mantid's sample frame, non-HKL `General Frame` volumes, `num_events`-based normalisation (the sum is
recorded in provenance).
