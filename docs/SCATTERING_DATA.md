# Labelled scattering data (Package B, Features 7 and 19)

`discuspdf.ScatteringVolume` is reciprocal-space data with everything needed to interpret it:
the `ReciprocalGrid` (hkl, lattice, Cartesian Q), the values, their variance, an explicit mask,
units, how Bragg scattering was treated, the corrections applied and provenance. It is a NumPy
container; xarray is an optional view on it, and it has its own versioned HDF5 format.

Install the optional pieces with `pip install -e ".[data]"` (h5py, xarray).

## Rules the class enforces

| Topic | Rule |
|---|---|
| Mask | boolean, **True = excluded**; each `with_mask` call records a reason in `mask_reasons` |
| Missing data | NaN/Inf allowed **only on masked voxels**; unmasked missing data is refused |
| Uncertainty | `variance` = σ² in `units²`, ≥ 0 on unmasked voxels |
| Numeric type | stored as float64 (integer and float32 input is converted; nothing else is accepted) |
| Units | `units` is required (e.g. `"DISCUS intensity (whole crystal)"`, `"e^2 per crystal"`) |
| Bragg | `bragg_present`, `bragg_subtracted` (e.g. DISCUS `set aver`), `bragg_masked` (needs a mask) or `unknown`; a source's declaration, not a proof |
| Negative values | kept; nothing is clipped |
| Corrections | `with_correction(values, "background: …")` records each step in `corrections` |
| Real-space transforms | refuse masked voxels unless `fill_masked=<value>` is given; the fill is recorded |

## Typical use

```python
from discuspdf import read_discus_volume, ScatteringVolume

sv = ScatteringVolume.from_discus_volume(read_discus_volume("C1_diffuse.h5"))
sv.bragg                     # 'bragg_subtracted' (from data_type_with_bragg)
sv.provenance["sha256"]      # hash of the source file; DISCUS metadata under 'discus_metadata'

rho = sv.delta_pdf()         # same result as delta_pdf_fractional(vol)

# Total scattering: mask Bragg regions instead of subtracting them
tot = ScatteringVolume.from_discus_volume(read_discus_volume("C1_total.h5"))
masked = tot.mask_bragg(radius=0.15)                 # Å⁻¹ around every integer hkl -> 'bragg_masked'
rho2 = masked.delta_pdf(fill_masked=0.0)             # explicit; recorded in rho2.provenance["corrections"]

sv.save_hdf5("c1_diffuse_labelled.h5")
again = ScatteringVolume.load_hdf5("c1_diffuse_labelled.h5")
```

`delta_pdf()` uses `delta_pdf_fractional` (any cell metric; the grid must be hkl-aligned, increasing
and centred). `to_reciprocal_volume()` gives the existing Cartesian `ReciprocalVolume`, only when the
grid's Q steps are parallel to x, y, z. Both refuse what `ReciprocalGrid.fft_compatibility()` refuses.

`mask_near_bragg(grid, radius)` masks by Euclidean distance in Q to the nearest reflection; the
radius must be below half the shortest reciprocal-lattice vector.

## Reading DISCUS output directly (volumes, planes, lines)

```python
sv = ScatteringVolume.from_discus_file("hk0_plane.h5")              # any DISCUS 6.23.01 form-nexus hkl file
sv = ScatteringVolume.from_discus_file("job_int.h5", config=cfg)    # also checked against its FourierConfig
sv = job.load_scattering_volume()                                   # Diffuse3DRun: uses the staged fourier_config.json
```

`from_discus_file` builds the grid with `ReciprocalGrid`, so planes and lines work (DISCUS writes a
zero step, or `tl − ll`, for a one-point axis; `read_discus_volume` refuses those files). Accepted
contents are `intensity`, `real part` and `imaginary part` (`value int/real/imag`); real-space files
(`3DPDF`, axes `uvw`) are refused. With a `FourierConfig` the reader refuses a file whose name, content
(`value`), Bragg label (`set aver`), shape, corner or steps (to the 6 digits written in the macro)
disagree with the configuration, and stores the configuration (seed, lots, labels) in
`provenance['fourier_config']`. `Diffuse3DRunner.execute(..., native_nexus=True)` now validates with
this reader, so plane runs end `complete` and a contradicting output ends `failed`.

## xarray

```python
ds = sv.to_xarray()
# dims (i0, i1, i2) with integer index coords; h, k, l on hkl-aligned grids; qx, qy, qz on
# Cartesian-aligned grids; intensity / variance / mask variables; metadata as JSON strings in attrs.
ds.swap_dims(i0="h", i1="k", i2="l")["intensity"].sel(l=0.0)   # label-based selection

sub = ScatteringVolume.from_xarray(ds.isel(i0=slice(10, 50, 2)))  # grid rebuilt for the selection
```

The dims stay `i0, i1, i2` because a skewed grid has no independent h, k, l axes. `from_xarray`
rebuilds the grid from the index coordinates, so `isel` with any constant step (including
reversal) is safe. It refuses: missing schema attributes (`attrs` dropped), transposed variables,
uneven selections, h/k/l or Q coordinates that no longer match the grid, a non-boolean mask (e.g.
after `ds * 2`), and unknown schema versions. Arithmetic in xarray does not update `units`,
`bragg` or `corrections`; use `with_correction` for corrections you want recorded.

## HDF5 schema `discuspdf.ScatteringVolume`, version 1

| Item | Content |
|---|---|
| root attrs | `schema`, `schema_version` (1), `written_by`, `written_utc`, `grid` (ReciprocalGrid JSON with full lattice matrix), `intensity_units`, `bragg`, `mask_reasons`, `corrections`, `provenance` (JSON), `intensity_sha256`, `variance_sha256` / `mask_sha256` when present |
| `intensity` | float64, grid shape |
| `variance` | float64, optional |
| `mask` | uint8 0/1, optional (1 = excluded) |

Loading checks the schema name and version (readers are registered per version, so a future
version adds a reader or migration rather than changing v1), refuses non-numeric datasets and mask
values other than 0/1, and verifies the SHA-256 of intensity, variance and mask. Nothing is unpickled. `save_hdf5`
refuses to overwrite unless `overwrite=True`; `compression="gzip"` is optional. Chunking and Zarr
are not offered: no size test has shown a need yet.

## Evidence

`tests/test_scattering_data.py`:

- On the real DISCUS C3 file, `ScatteringVolume.delta_pdf` reproduces `delta_pdf_fractional` to 1e-13
  (β round-trips through the lattice matrix as 105.00000000000001).
- On an orthorhombic grid, the fractional path and the Cartesian path (`to_reciprocal_volume` →
  `delta_pdf_from_diffuse`) give the same complex volume.
- `mask_near_bragg` equals a brute-force enumeration of all reflections, also on a triclinic, skewed grid.
- xarray round trips (aligned and the skewed C4 file), four `isel` selections incl. reversal and
  single layers, and five refusal cases; HDF5 round trips (with and without gzip, real C3 provenance)
  and nine refusal cases (version, schema, edited intensity, edited variance, edited mask, mask
  value, missing attribute, string dataset, not HDF5).
- Mutation checks: 8 of 9 caught. Not caught: restricting the Bragg-mask search to `round(hkl)`
  only; within the allowed radius no tested lattice (up to γ = 165°) makes rounding pick the
  wrong reflection, so the 27-point search is a safety net.

## Not included

Chunked/lazy access and Zarr (no size test justifies them yet), Mantid MDHisto (Feature 18),
`set aver`/LOTS controls (Features 8 and 9, waiting for the engine probes), visualisation (Feature 20).
