# Package B — Reciprocal-space grids, scattering datasets and data exchange

**Status:** TODO · **Priority:** Very high · **Suggested PR family:** `feat/reciprocal-*`

**Reference in Javelin:** [grid.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/grid.py), [fourier.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/fourier.py), [io.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/io.py), [DISCUS tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/discus.rst).

## Objective

Work with measured or simulated reciprocal-space intensities as calibrated, labeled multidimensional data. Add flexible Q-space selections and export without silently changing unit conventions or sacrificing scientific provenance.

## Planned features

### 6. Flexible reciprocal-space grids
- [ ] `ReciprocalGrid` defined by corners / orientation vectors / ranges / number of bins; support arbitrary reciprocal-space planes and volumes.
- [ ] Distinguish h,k,l (r.l.u.), Cartesian `Q` (rad/Å), voxel centers vs edges, and direct ↔ reciprocal lattice basis.
- [ ] Validate invertibility, axis ordering, handedness, shape and non-orthogonal coordinate transformations.
- [ ] For **3D FFT**, Cartesian orthogonal/regular-grid prerequisites must be enforced; arbitrary oblique planes should not be fed into the existing Cartesian transform without documented resampling/metric or an appropriate nonuniform algorithm.

### 7. xarray-based, metadata-rich scattering data
- [ ] Add `xarray.DataArray` / `Dataset` adapter or integration to the existing `ReciprocalVolume`, retaining numeric type, labeled axes, errors, masks, units, lattice and provenance.
- [ ] Make missing data, masked Bragg regions and background corrections explicit. Preserve validation rules already being developed by OpenCode/Claude.
- [ ] Support memory-efficient access/chunking only after size and numerical tests justify it.

### 8. Average-structure/Bragg subtraction
- [ ] Expose **verified DISCUS** choices for average structure subtraction / Bragg exclusion; save method and parameters in provenance.
- [ ] Confirm whether a given method subtracts coherent **amplitudes** or **intensities**, and specify mathematically what output contains.
- [ ] Test a perfect crystal (diffuse → nearly zero away from Bragg features) and a known disordered supercell; preserve negative-valued corrected data when physically warranted.

### 9. LOTS averaging
- [ ] Expose LOTS (lots of tiny scatters) settings through Python **using validated DISCUS features**, such as subvolume size, number of lots, seeds and averaging choices.
- [ ] Test model-size and random-seed dependence; record sampling variance and normalization. Do not conflate finite-domain average and Bragg removal automatically.

### 18. Mantid MDHisto import/export
- [ ] Explore Mantid MDHisto/NeXus HDF5 readers and exporters with real fixtures; read signal, error/variance, masks, bin edges/centers, orientation, projection matrix and lattice.
- [ ] Reject or explicitly convert unsupported axis frames and ambiguous metadata.
- [ ] Round-trip fixtures with verified coordinate positions and intensity values; do not assume Javelin's older file layout covers current Mantid releases.

### 19. Labeled persistence
- [ ] Save/load a versioned HDF5 or Zarr/xarray schema with intensity, axes, uncertainty, masks, cell information, input hashes, method and engine version.
- [ ] Guard against silent unit/schema mismatch and unsafe deserialization; write schema migration tests.

### 20. 3D visualization and VTK exports
- [ ] Add Q-volume and 3D-ΔPDF slice / isosurface exploration with physical labels and sensible color conventions.
- [ ] Export VTI only when grid axes are physically axis-aligned; otherwise support an appropriate VTK structured grid or clearly state restrictions.
- [ ] Ensure coordinate calibration and orientation survive export.

## Suggested PRs

1. ReciprocalGrid + lattice/coordinate tests (coordinate with Package A).
2. xarray adapter and validated native DISCUS outputs, once Claude's format verification has landed.
3. Validated Bragg/LOTS controls, tied to known DISCUS v6 commands and physical benchmarks.
4. Mantid and persistence adapters.
5. Visualization and VTK exporters.

## Acceptance criteria

- [ ] 3D datasets carry sufficient source-derived coordinate metadata to distinguish r.l.u. from Å⁻¹ and understand 2π factors.
- [ ] Reference peaks at known h,k,l map to independently checked Cartesian Q and real-space distances.
- [ ] Masks/errors/provenance survive round trips and clearly report unsupported data.
- [ ] New data paths do not break existing `ReciprocalVolume` or `FFTResult` usage without migration notes.
- [ ] Real DISCUS/experimental files and synthetic numerical fixtures have separate clearly labeled test coverage.

## Dependencies / cautions

**Requires Claude's real 3D DISCUS macro/HDF5 investigation.** The existing `Diffuse3DRunner` is an adapter, not evidence that any HDF5 layout or physics is verified. Coordinate calibration and FFT normalization must be benchmarked before quantitative interpretations. Prefer calling DISCUS for core scattering, not porting Javelin's computation verbatim.
