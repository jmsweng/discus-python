# Plotting and VTK export

These tools draw reciprocal-space data and 3D-ΔPDF maps **in their real geometry**, with axes in
Å⁻¹ or Å, and write them for ParaView:

- `plot_slice`: a 2D layer through the data, with true lengths and angles (a monoclinic layer is
  drawn as a parallelogram, not a rectangle);
- `plot_isosurfaces`: 3D surfaces at chosen values (needs `pip install -e ".[viz]"`, i.e. scikit-image);
- `write_vtk`: files that ParaView (or any VTK reader) opens with the correct coordinates.

All of them take a `SampledVolume`: the values plus the Cartesian position of every point, taken
from the data's own lattice and grid, never guessed from array indices.

## Quick start

```python
from discuspdf import ScatteringVolume, SampledVolume, plot_slice, plot_isosurfaces, write_vtk

sv = ScatteringVolume.from_discus_file("C1_diffuse.h5")
q = SampledVolume.from_scattering_volume(sv)          # positions in Q (Å⁻¹)
plot_slice(q, axis=2)                                 # the middle layer along the third axis

rho = SampledVolume.from_fractional_fft(sv.delta_pdf())   # positions in r (Å)
plot_slice(rho, axis=2)                               # red positive, blue negative, white zero
plot_isosurfaces(rho, levels=(+2e-3, -2e-3))
write_vtk(rho, "c1_delta_pdf.vts")                    # open in ParaView
```

## Making a SampledVolume

| You have | Use | Positions |
|---|---|---|
| a `ScatteringVolume` | `SampledVolume.from_scattering_volume(sv)` | Q, Å⁻¹ (keeps mask, variance, Bragg label) |
| a 3D-ΔPDF from `delta_pdf_fractional` or `sv.delta_pdf()` | `SampledVolume.from_fractional_fft(result)` | r, Å |
| a result of the older `delta_pdf_from_diffuse` | `SampledVolume.from_fft_result(result)` | r, Å |

The orientation is that of the lattice: by default a along x and b in the xy plane.

## Slices

```python
plot_slice(vol, axis=2, index=10)       # layer 10 perpendicular to array axis 2
plot_slice(vol, axis=0, log=True)       # log colour scale (unsigned data)
plot_slice(vol, axis=2, ax=my_axes, title="...", vmax=1e5)
```

- Axis labels name the crystal direction, e.g. `Q ∥ [1 0 0]* (Å⁻¹)`; the title gives the layer's
  offset from the origin.
- Data that can be negative (Bragg-subtracted intensity, ΔPDF, Patterson) get a red–blue scale
  that is symmetric about zero. Other data use `viridis`.
- Masked points are grey; they are never filled in.

## Isosurfaces

```python
plot_isosurfaces(rho, levels=(+2e-3, -2e-3))           # one surface per level, equal axis scaling
verts, faces = isosurface(rho, 2e-3)                   # the mesh itself, in Å
```

Masked points need an explicit `fill_value=`.

## VTK / ParaView

```python
write_vtk(vol, "out.vts")       # works for every grid
write_vtk(vol, "out.vti")       # smaller, but only for grids whose axes run along x, y, z
data = read_vtk("out.vts")      # read back without VTK installed
```

| Extension | Use it for |
|---|---|
| `.vts` (StructuredGrid) | anything: non-orthogonal cells, skewed grids, oblique planes |
| `.vti` (ImageData) | only grids whose axes run along +x, +y, +z; anything else is refused, because ParaView would draw it in the wrong place |

Each file holds the values (masked points as NaN), the mask (1 = excluded) and any extra arrays
such as variance. A small `<file>.json` next to it records units, origin, steps, crystal directions
and provenance. VTK does not need to be installed to write the files.

## Not included

An interactive viewer (ParaView or napari can open the files). The older `plot_delta_pdf_slice`
remains for existing code, but labels its axes by index only; use `plot_slice` for new work.

<details>
<summary>How this was checked</summary>

`tests/test_viz.py`:
- The position of every point equals `ReciprocalGrid.q()` (reciprocal space) or `u @ A` (real space).
- On the real DISCUS C3 file (monoclinic), the strongest Patterson peak lies within one voxel of the
  interatomic vector u = ±(0.3, 0.2, 0.1), converted to Å with the cell metric.
- VTK 9.7 reads our `.vts` and `.vti` files with matching points, values, mask and variance (this
  test needs `vtk` installed and is skipped in CI; the `read_vtk` round trip always runs).
- Slices: drawn cell edges have the true lengths and angle; masks and the symmetric scale are kept.
- An isosurface of |r|² = 1 on an orthogonal and on a skewed grid lies within 0.02 of radius 1.
- Deliberate bugs (point order, plane basis, FFT origin, axis-alignment test, isosurface mapping,
  colour scale, masking) are each caught.
</details>
