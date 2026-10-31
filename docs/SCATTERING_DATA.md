# Scattering data with labels

A `ScatteringVolume` is a block of reciprocal-space data that carries everything you need to
interpret it later:

- **where** each value is: a [`ReciprocalGrid`](RECIPROCAL_GRID.md) (hkl, lattice, Q);
- the **values**, their **variance** (σ²) and **units**;
- a **mask** of excluded points, with the reason for each masking step;
- how **Bragg scattering** was handled (present, subtracted, masked or unknown);
- the **corrections** applied, and where the data came from (**provenance**: file, SHA-256, DISCUS settings).

It is plain NumPy underneath. It can be turned into an xarray Dataset, and saved to and loaded from
its own HDF5 format. Install the optional pieces with `pip install -e ".[data]"` (h5py, xarray).

## Quick start

```python
from discuspdf import ScatteringVolume

sv = ScatteringVolume.from_discus_file("C1_diffuse.h5")   # any DISCUS 6.23.01 `form nexus` hkl file
sv.intensity.shape        # the values, as a NumPy array
sv.grid.hkl()             # hkl of every point
sv.bragg                  # 'bragg_subtracted' (DISCUS ran with `set aver`)
sv.provenance["sha256"]   # hash of the file it came from

rho = sv.delta_pdf()      # the 3D-ΔPDF (needs Bragg-subtracted or Bragg-masked data)
```

## Loading data

| Source | Use |
|---|---|
| a DISCUS 6.23.01 output file (volume, plane or line) | `ScatteringVolume.from_discus_file(path)` |
| … and check it matches the `FourierConfig` that produced it | `ScatteringVolume.from_discus_file(path, config=cfg)` |
| a run made with `Diffuse3DRunner` | `job.load_scattering_volume()` (uses the configuration saved with the run) |
| a `DiscusVolume` you already read | `ScatteringVolume.from_discus_volume(vol)` |
| a Mantid MDHisto file | `read_mantid_mdhisto(path)` (see [MANTID_IO.md](MANTID_IO.md)) |
| your own arrays | `ScatteringVolume(grid, intensity, units="...", bragg="...", variance=..., mask=...)` |
| a saved file / an xarray Dataset | `ScatteringVolume.load_hdf5(path)` / `ScatteringVolume.from_xarray(ds)` |

`from_discus_file` reads `value int`, `real` and `imag` files. It refuses real-space files
(`value 3DPDF`). When you give it a `config`, it refuses a file whose name, content, Bragg setting,
shape, corner or steps disagree with that configuration, and it stores the configuration in
`sv.provenance["fourier_config"]`.

## Dealing with Bragg peaks

A 3D-ΔPDF needs the Bragg peaks removed. There are two ways:

1. **Subtract them in DISCUS** with `set aver,100` (see [FOURIER_CONTROLS.md](FOURIER_CONTROLS.md)).
   The file is then labelled `bragg_subtracted` and `sv.delta_pdf()` works directly.
2. **Mask them** in total-scattering data, then say what value the masked points should take:

```python
total = ScatteringVolume.from_discus_file("C1_total.h5")      # bragg_present
masked = total.mask_bragg(radius=0.15)                         # masks within 0.15 Å⁻¹ of every integer hkl
rho = masked.delta_pdf(fill_masked=0.0)                        # the fill value is recorded in rho's provenance
```

`delta_pdf` on data with Bragg peaks present is refused. For a Patterson function of total
scattering, pass `require_diffuse=False`.

## Masks, corrections and missing data

```python
sv2 = sv.with_mask(my_mask, "detector gap")                   # True = excluded; adds to the existing mask
sv3 = sv2.with_correction(sv2.intensity - bkg, "background: empty-can run 1234")
sv3.mask_reasons      # ('detector gap',)
sv3.corrections       # ('background: empty-can run 1234',)
filled = sv3.filled(0.0)                                      # a copy of the values with masked points set to 0
```

- NaN or infinite values are allowed **only on masked points**.
- Variance must be ≥ 0 on unmasked points, in units².
- Negative values are kept (Bragg-subtracted data can be negative); nothing is clipped.
- Every method returns a new object; the original is not changed.
- Transforms refuse masked points unless you give `fill_masked=`, so a fill value is never chosen silently.

## Working in xarray

```python
ds = sv.to_xarray()
ds.swap_dims(i0="h", i1="k", i2="l")["intensity"].sel(l=0.0)     # the l = 0 layer, by label

sub = ScatteringVolume.from_xarray(ds.isel(i0=slice(10, 50, 2)))  # back to a ScatteringVolume
```

The dimensions are `i0, i1, i2` (point indices), with `h, k, l` added as coordinates when the grid
runs along h, k, l, and `qx, qy, qz` when it runs along x, y, z. A skewed grid has no separate
h, k, l axes, so the index dimensions are what always works.

Going back with `from_xarray` works after any evenly spaced selection (`isel` with a constant step,
including reversing). It refuses data that no longer matches the grid: uneven selections,
transposed variables, a mask that is no longer True/False (e.g. after `ds * 2`), or missing
attributes. Arithmetic in xarray does not update `units`, `bragg` or `corrections`; use
`with_correction` for steps you want recorded.

## Saving and loading

```python
sv.save_hdf5("c1_diffuse_labelled.h5")                 # refuses to overwrite unless overwrite=True
again = ScatteringVolume.load_hdf5("c1_diffuse_labelled.h5")
```

Loading checks the file's format name and version and a SHA-256 checksum of every array, so an
edited or damaged file is refused rather than read wrongly.

## Reference

**Bragg labels** (`sv.bragg`): `bragg_present`, `bragg_subtracted` (e.g. DISCUS `set aver`),
`bragg_masked` (needs a mask) or `unknown`. A label records what the source says was done; it is
not a test of the data.

**HDF5 format** `discuspdf.ScatteringVolume`, version 1:

| Item | Content |
|---|---|
| `intensity` | float64, grid shape |
| `variance` | float64, optional |
| `mask` | uint8, 1 = excluded, optional |
| root attributes | `schema`, `schema_version`, `written_by`, `written_utc`, `grid` (JSON, full lattice matrix), `intensity_units`, `bragg`, `mask_reasons`, `corrections`, `provenance` (JSON), SHA-256 of each array |

Optional `compression="gzip"`. Nothing is unpickled on loading. A future format version will add
a reader rather than change version 1.

**Other details**
- Values are stored as float64 (integer and float32 input is converted).
- `mask_near_bragg(grid, radius)` is the function behind `mask_bragg`: it masks by distance in Q to
  the nearest reflection. The radius must be less than half the shortest reciprocal-lattice vector.
- `to_reciprocal_volume()` converts to the older Cartesian `ReciprocalVolume`, only for grids whose
  steps run along x, y, z.

<details>
<summary>How this was checked</summary>

`tests/test_scattering_data.py` and `tests/test_discus_scattering.py`:
- On the real DISCUS C3 file, `delta_pdf()` reproduces `delta_pdf_fractional` to 1e-13; on an
  orthorhombic grid the fractional and Cartesian routes give the same result.
- `mask_near_bragg` equals a brute-force search over all reflections, also on a skewed triclinic grid.
- xarray round trips (aligned, and the skewed C4 file), four kinds of selection, and five refusals;
  HDF5 round trips (with and without gzip) and nine refusals (wrong version, edited arrays, bad mask
  values, missing attributes, non-HDF5 files).
- DISCUS planes and lines load with the right grid; files that contradict their configuration are refused.
- Deliberate bugs in the new code are caught by the tests (8 of 9 for the data class; the one not
  caught is a redundant safety check in the Bragg-mask search).
</details>
