# Reciprocal-space grids (Package B, Feature 6)

`discuspdf.ReciprocalGrid` describes **where** reciprocal-space samples are: a regular (affine)
grid of points in hkl, tied to a `Lattice`. It holds no intensities; Feature 7 (xarray) attaches
data to it. It never resamples, and it never relabels an oblique grid as Cartesian Q.

## Conventions

| Quantity | Representation |
|---|---|
| Point at integer index `i = (i0, i1, i2)` | `hkl = corner + i @ axis_vectors` (r.l.u.) |
| `axis_vectors` | 3×3, **row** `k` = hkl step along array axis `k` |
| DISCUS 6.23.01 NeXus | `data_increment_vector[:, k]` (**column** `k`) = `axis_vectors[k]`; verified on the skewed real C4 fixture |
| Sample positions | voxel **centres**; edges at half-integer indices (`index_to_hkl(i - 0.5)`) |
| Cartesian Q | `Q = 2π · (hkl @ B)` Å⁻¹ for the grid's `Lattice` orientation, `|Q| = 2π/d` |
| Planes and lines | an axis with one point is *inactive*; its vector never enters a coordinate |
| Centred (FFT-ready) | hkl = 0 at index `n // 2` on every axis (numpy `fftshift` convention; even `n` has one extra negative point) |
| Real-space dual grid | rows `R[k]` with `axis_vectors @ R[k] = e_k / n_k` (fractional); `1/(n·dh)` for an hkl-aligned grid |

## Construction

```python
from discuspdf import Lattice, ReciprocalGrid

lat = Lattice.from_parameters(5, 6, 7, 90, 105, 90)

# Same commands as the DISCUS fourier menu (ll, lr, ul, tl, na, no, nt):
g = ReciprocalGrid.from_corners(lat, ll=(-4, -5, -6), lr=(4, -5, -6),
                                ul=(-3, 5, -6), tl=(-4, -5, 6), shape=(33, 41, 49))

g = ReciprocalGrid.from_ranges(lat, h=(-4, 4, 33), k=(-5, 5, 41), l=(0, 0, 1))   # hk0 plane
g = ReciprocalGrid.centered(lat, shape=(64, 64, 64), steps=(0.125, 0.125, 0.125))  # FFT-ready
g = ReciprocalGrid.from_discus_volume(read_discus_volume("C1_diffuse.h5"))
```

`from_discus_volume` rebuilds the lattice from the file's cell parameters in the canonical
orientation (a ∥ x, b in the xy plane). DISCUS files carry no orientation matrix; give a rotated
`Lattice` (`Lattice.from_matrix`) if an experimental orientation matters.

## What it answers

- Coordinates: `hkl()`, `q()`, `q_magnitude()`, `index_to_hkl(idx)`, `hkl_to_index(hkl)` (volumes only),
  `hkl_axes()` (hkl-aligned grids only), `q_axes()` (only when step k ∥ +Cartesian axis k).
- Geometry: `is_hkl_axis_aligned()`, `is_q_axis_aligned()`, `is_q_orthogonal()`, `handedness()`,
  `is_centered()`, `real_space_step_vectors()`.
- Derived grids: `plane(axis, index)`; `discus_corners()` gives the `ll/lr/ul/tl/na/no/nt` values for a macro.
- Persistence: `to_dict()` / `from_dict()` / JSON, schema `discuspdf.ReciprocalGrid` version 1
  (unknown schema or version is refused). The full lattice matrix is stored, so orientation survives.

### Which 3D transform is exact

`grid.fft_compatibility()` returns `FFTCompatibility(fractional, cartesian, reasons)`:

| Transform | Requires | Typical case |
|---|---|---|
| `delta_pdf_fractional` (fractional) | volume, hkl-aligned, increasing, centred | any cell, incl. monoclinic/triclinic |
| `delta_pdf_from_diffuse` / `ReciprocalVolume` (cartesian) | volume, Q step k ∥ +x/+y/+z, centred | orthogonal cell, hkl-aligned grid |
| neither | skewed grid, plane/line, uncentred, decreasing axis | `reasons` says which |

A skewed (oblique) volume is still a regular grid, and its DFT is exact on the dual real-space
grid given by `real_space_step_vectors()` (tested against a brute-force sum), but no current
transform function accepts such data. Feeding it into the Cartesian transform is refused by
`q_axes()`.

## Evidence

`tests/test_reciprocal_grid.py` (32 tests):

- Q of seven reflections in a monoclinic cell equals 2π/d from the textbook monoclinic formula
  (to 1e-12); the a*–c* angle equals 180° − β; hexagonal a*–b* angle is 60°.
- Real DISCUS 6.23.01 output: grids built from the `ll/lr/ul/tl/na/no/nt` lines of the committed
  C3 (aligned) and C4 (skewed) macros equal the grids read from the engine's HDF5 files.
- C3: `real_space_step_vectors()` equals the `u` spacing of `delta_pdf_fractional`.
- Mutation checks (transposed DISCUS increments, untransposed dual, missing 2π, centring, sign,
  inactive-axis and rank checks disabled): each makes at least one test fail.

## Not verified / not included

- What DISCUS writes in `data_increment_vector` for a **plane or line** (an axis with one point):
  no such engine file exists in the repository yet (probe in the Package B handoff).
- Transforms for skewed volumes, resampling, and masks/uncertainties (Feature 7 onward).
- `DiscusVolume` still has its own geometry code; consolidating it onto `Lattice`/`ReciprocalGrid`
  is a separate PR with Gate 2/3 evidence parity.
