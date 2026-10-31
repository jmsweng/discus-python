# Reciprocal-space grids

A `ReciprocalGrid` describes **where** reciprocal-space data points are: which h, k, l each point
has, and its Q in Å⁻¹. It holds no intensities; [`ScatteringVolume`](SCATTERING_DATA.md) attaches
the data. You use a grid to:

- set up the points DISCUS should calculate (the `ll/lr/ul/tl/na/no/nt` lines of a Fourier macro);
- look up the hkl or Q of any data point;
- find out whether a 3D-ΔPDF can be calculated on that grid, and if not, why.

## Quick start

```python
from discuspdf import Lattice, ReciprocalGrid

lat = Lattice.from_parameters(5, 6, 7, 90, 105, 90)          # a, b, c (Å), α, β, γ (°)

g = ReciprocalGrid.from_ranges(lat, h=(-4, 4, 33), k=(-5, 5, 41), l=(-6, 6, 49))
g.shape                       # (33, 41, 49)
g.index_to_hkl((16, 20, 24))  # array([0., 0., 0.])
g.q_magnitude()               # |Q| in Å⁻¹ for every point, shape (33, 41, 49)
```

## Making a grid

| You have | Use |
|---|---|
| an h, k, l range for each axis | `ReciprocalGrid.from_ranges(lat, h=(start, stop, n), k=..., l=...)` |
| the corners you would type in DISCUS's `fourier` menu | `ReciprocalGrid.from_corners(lat, ll=..., lr=..., ul=..., tl=..., shape=(na, no, nt))` |
| a grid you want to Fourier-transform to a 3D-ΔPDF | `ReciprocalGrid.centered(lat, shape=(64, 64, 64), steps=(0.125, 0.125, 0.125))` |
| a DISCUS output file | `ReciprocalGrid.from_discus_volume(read_discus_volume("file.h5"))`, or load it as a `ScatteringVolume` (also works for planes and lines) |

```python
# A plane: give one point along the axis you want to keep fixed
hk0 = ReciprocalGrid.from_ranges(lat, h=(-4, 4, 33), k=(-5, 5, 41), l=(0, 0, 1))

# A skewed grid, exactly as in a DISCUS macro (the second axis also steps in h)
c4 = ReciprocalGrid.from_corners(lat, ll=(-4, -5, -6), lr=(4, -5, -6),
                                 ul=(-3, 5, -6), tl=(-4, -5, 6), shape=(33, 41, 49))

# Going the other way: the values to put in a DISCUS macro
hk0.discus_corners()   # {'ll': [...], 'lr': [...], 'ul': [...], 'tl': [...], 'na': 33, 'no': 41, 'nt': 1}
```

A grid read from a DISCUS file gets its lattice from the file's a, b, c, α, β, γ, in the standard
orientation (a along x, b in the xy plane). DISCUS files do not record a crystal orientation; if
the orientation matters (e.g. to compare with an experiment), build the `Lattice` yourself with
`Lattice.from_matrix`.

## Asking a grid questions

| Question | Method |
|---|---|
| hkl of every point / of one point | `g.hkl()` · `g.index_to_hkl((i, j, k))` |
| which point has this hkl (3D grids) | `g.hkl_to_index((h, k, l))` (can be fractional) |
| Q vector / \|Q\| of every point, Å⁻¹ | `g.q()` · `g.q_magnitude()` |
| 1D h, k, l axes (when each axis runs along one of h, k, l) | `g.hkl_axes()` |
| one layer of a 3D grid | `g.plane(axis, index)` |
| are two grids the same | `g.isclose(other)` |
| can I compute a 3D-ΔPDF here | `g.fft_compatibility()` (below) |
| save / load the grid | `g.to_dict()` / `ReciprocalGrid.from_dict(d)` (JSON-friendly; keeps the full orientation) |

## Can this grid be Fourier-transformed?

`g.fft_compatibility()` answers this before you try, and lists the reasons when the answer is no:

```python
c4.fft_compatibility()
# FFTCompatibility(fractional=False, cartesian=False, reasons=(
#   'hkl = 0 is not at index n//2 on every axis (centred FFT would be mis-phased)',
#   'fractional: grid is skewed in hkl (non-diagonal axis_vectors)', ...))
```

| Result | Meaning | Typical grid |
|---|---|---|
| `fractional=True` | `delta_pdf_fractional` / `ScatteringVolume.delta_pdf()` work | any cell (incl. monoclinic, triclinic), axes along h, k, l, hkl = 0 in the centre |
| `cartesian=True` | the older Cartesian-Q transform (`delta_pdf_from_diffuse`) also works | orthogonal cell, axes along h, k, l, centred |
| both False | no transform in this package accepts the grid | skewed grids, planes, lines, grids not centred on hkl = 0 |

To get a grid that passes, use `ReciprocalGrid.centered(...)`: it puts hkl = 0 at index `n // 2`
on every axis, as NumPy's `fftshift` does.

## Reference

| Item | Convention |
|---|---|
| point at index `(i0, i1, i2)` | hkl = `corner + i @ axis_vectors` (r.l.u.) |
| `axis_vectors` | 3×3; **row** k is the hkl step along array axis k |
| DISCUS 6.23.01 NeXus files | `data_increment_vector[:, k]` (**column** k) is the step of axis k |
| Q | `Q = 2π · hkl · B` in Å⁻¹, so \|Q\| = 2π/d |
| points | voxel **centres**; voxel edges are half a step either side |
| planes and lines | an axis with one point is ignored in all coordinates (DISCUS writes a zero step, or `tl − ll`, for it) |
| centred | hkl = 0 at index `n // 2` on every axis (for even n there is one more negative point) |
| real-space step of the matching 3D-ΔPDF | `g.real_space_step_vectors()` (fractional); 1/(n·Δh) for a grid along h, k, l |

The grid never interpolates or resamples, and it never treats an oblique hkl grid as if it were a
Cartesian Q grid.

<details>
<summary>How this was checked</summary>

`tests/test_reciprocal_grid.py`:
- |Q| for seven reflections in a monoclinic cell equals 2π/d from the textbook formula (to 1e-12);
  the a*–c* angle is 180° − β and the hexagonal a*–b* angle is 60°.
- Grids built from the corner lines of the committed DISCUS 6.23.01 C3 (aligned) and C4 (skewed)
  macros equal the grids in the files DISCUS wrote.
- `real_space_step_vectors()` matches the real-space spacing of `delta_pdf_fractional` on C3.
- Deliberate bugs (transposed DISCUS steps, missing 2π, wrong centring, sign, plane handling) each
  make at least one test fail.

`DiscusVolume` (the older reader) still has its own geometry code; moving it onto
`Lattice`/`ReciprocalGrid` is planned separately.
</details>
