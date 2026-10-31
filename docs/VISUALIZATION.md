# Visualisation and VTK export (Package B, Feature 20)

`discuspdf.SampledVolume` holds values on a regular grid together with the **Cartesian** position
of every point: point `(i0, i1, i2)` sits at `origin + i @ steps`. Its constructors take the
calibrated objects of this package, so coordinates come from source metadata, never from indices:

| Constructor | Source | Coordinates |
|---|---|---|
| `SampledVolume.from_scattering_volume(sv)` | `ScatteringVolume` (grid, mask, variance, Bragg label) | Q = 2π·hkl·B, Å⁻¹ |
| `SampledVolume.from_fractional_fft(res)` | `delta_pdf_fractional` result | r = u·A, Å |
| `SampledVolume.from_fft_result(res)` | `delta_pdf_from_diffuse` result | Cartesian r, Å |

The orientation is the `Lattice` orientation (canonical: a ∥ x, b in the xy plane, unless the lattice
was built from a rotated matrix).

```python
from discuspdf import read_discus_volume, delta_pdf_fractional, ScatteringVolume
from discuspdf.viz import SampledVolume, plot_slice, plot_isosurfaces, write_vtk

vol = read_discus_volume("C1_diffuse.h5")
q = SampledVolume.from_scattering_volume(ScatteringVolume.from_discus_volume(vol))
plot_slice(q, axis=2)                          # hk-plane through the centre, Å⁻¹, true angles

rho = SampledVolume.from_fractional_fft(delta_pdf_fractional(vol))
plot_slice(rho, axis=2)                        # signed: symmetric RdBu_r colour scale
plot_isosurfaces(rho, levels=(+2e-3, -2e-3))   # needs pip install -e ".[viz]" (scikit-image)
write_vtk(rho, "c1_delta_pdf.vts")             # ParaView; non-orthogonal or oblique -> .vts
```

## Slices

`plot_slice(vol, axis, index)` draws one layer perpendicular to array `axis` with `pcolormesh` on
the projected voxel corners:
- in-plane coordinates on an orthonormal basis of the plane (first axis along x), so cell edges have
  their real lengths and the angle between the axes is the real angle (`aspect='equal'`);
- labels name the direction in crystal coordinates, e.g. `Q ∥ [1 0.2 0]* (Å⁻¹)`, and the title gives
  the plane's offset from the origin;
- signed data (Bragg-subtracted intensity, ΔPDF, Patterson) use `RdBu_r` with a symmetric range;
  unsigned data use `viridis`, optionally `log=True`;
- masked voxels are grey (they are never filled with a number).

## Isosurfaces

`isosurface(vol, level)` runs marching cubes (scikit-image, optional extra `viz`) in index space and
maps the vertices with the grid's affine geometry. Masked voxels need an explicit `fill_value`.
`plot_isosurfaces(vol, levels)` draws them with equal axis scaling.

## VTK / ParaView

`write_vtk(vol, path)` writes VTK XML with raw appended binary data (UInt64 headers, little endian);
no VTK installation is needed to write.

| Extension | VTK type | When |
|---|---|---|
| `.vti` | ImageData (origin + spacing) | **only** if every active axis runs along +x, +y, +z; otherwise refused |
| `.vts` | StructuredGrid (explicit point coordinates) | any regular grid: oblique planes, non-orthogonal cells, skewed hkl grids |

Point data: the values (named after `quantity`, masked voxels written as NaN), `mask` (UInt8,
1 = excluded) and extra arrays such as `variance`. A sidecar `<file>.json` records space, unit,
origin, steps, the crystal index vectors and the provenance metadata. `read_vtk` reads these files back
without VTK.

## Evidence

`tests/test_viz.py` (17 tests):
- Q of every point equals `ReciprocalGrid.q()`; real-space points equal `u @ A`.
- Real DISCUS C3 file (monoclinic, Bragg present): the strongest interatomic Patterson peak inside one
  cell lies within one voxel of u = ±(0.3, 0.2, 0.1) converted to Å with the cell metric.
- VTK 9.7 (`vtkXMLStructuredGridReader`, `vtkXMLImageDataReader`) reads our files: points, values
  (incl. NaN), mask, variance, ImageData origin/spacing and point lookup all match. (Skipped when
  `vtk` is not installed; CI does not install it. The package's own `read_vtk` round trip always runs.)
- Slices: plotted cell edges have the true step lengths and angle; masked voxels stay masked; signed
  scale is symmetric.
- Isosurface of |r|² = 1 on an orthogonal and on a skewed grid: vertex radii within 0.02 of 1.
- Mutation checks: VTK point order, StructuredGrid point order, plane basis, FFT origin, axis-alignment
  test, isosurface mapping, symmetric scale and masking — each caught.

## Not included

An interactive viewer (ParaView or napari can open the `.vts`/`.vti` files), DISCUS plane files
(read them into `ScatteringVolume` first), and the older `plot_delta_pdf_slice`, which stays for
compatibility but labels index-derived axes only.
