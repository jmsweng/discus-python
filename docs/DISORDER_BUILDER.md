# Disorder builder back end (`discuspdf.builder`)

The [discus-disorder-builder](https://github.com/jmsweng/discus-disorder-builder) web page is a user interface:
it edits a disordered crystal model and exports it. Building the structure and running DISCUS happen here.

## Running a model

```bash
python -m discuspdf.builder run model.json            # build, then every enabled DISCUS calculation
python -m discuspdf.builder run model.json --dry-run  # stage every DISCUS input, start nothing
python model.py --dry-run                             # the script the page exports: same options
```

```python
from discuspdf.builder import run
res = run("model.json", out="results", dry_run=True)
res["structure"]                       # the configuration Structure (with VOID padding rows)
```

Options: `--out`, `--seed` (overrides the model's), `--launcher` (e.g.
`"wsl.exe -d Ubuntu-22.04 -- discus_suite -noautorun -macro"`), `--structure-only`, `--numpy-preview`,
`--overwrite`.

## What a run does

1. **Structure** (`builder.realize`): an occupant (atom, molecule orientation or vacancy) is drawn for every
   site of every cell, at random or with the exact composition. If the model asks for it, the occupants are
   ordered by Monte Carlo (`builder.ordering`: neighbour energies J and forbidden-contact penalties; works for
   molecules and vacancies). Then the atoms are placed with rigid-body translation and libration of molecules
   and Gaussian atomic displacements. Every cell is padded with `VOID` atoms to the same row count, in DISCUS
   cell-major order. A clash report lists atoms of different molecules or sites that are closer than the
   model's threshold.
2. **DISCUS mmc** (atom-only models, optional): `discuspdf.mc`.
3. **Files**: `supercell.stru` (`write_stru`; molecule configurations use `check_order=False,
   min_separation=None`, see [DISCUS_IO](DISCUS_IO.md)), `supercell.xyz` (`write_extxyz`, without VOIDs),
   `structure_summary.json`, `model.json`.
4. **3D** (`FourierConfig` → `Diffuse3DRunner`): diffuse scattering with `set aver`, optional LOTS, then the
   3D-ΔPDF. Optionally the total scattering and its 3D-PDF (Patterson).
5. **1D PDF** (`PDFConfig` with the `.stru` → `PDFRunner`): the Debye sum over the whole supercell.

Own neutron scattering lengths (model key `scattering`, e.g. `{"NI62": -8.7}`) reach DISCUS as `scat` lines
through `FourierConfig(scattering_lengths=…)` and `PDFConfig(scattering_lengths=…)`, only for names that occur
in the structure and only for neutron runs.

## Output formats (model key `outputs`, the page's Script step)

| Key | Default | Files |
|---|---|---|
| `hdf5` | on | `diffuse.h5`, `total.h5` (`ScatteringVolume.save_hdf5`), `delta_pdf.h5`, `pdf3d.h5` |
| `csv` | on | `pdf_1d.csv` (`r (Angstrom),G(r)`) |
| `png` | on | `*_slices.png`, `pdf_1d.png` |
| `vtk` | off | `delta_pdf.vts`, `pdf3d_patterson.vts` |
| `npz` | off | `diffuse.npz`, `total.npz`, `delta_pdf.npz` (the builder engine 1.0 files) |
| `txt` | off | `pdf_1d.txt` |

Each DISCUS run also keeps its own directory under `runs_3d/`, `runs_pdf/` and `runs_mc/`: the macro, the
staged structure, the configuration, the log and the raw output.

## The model format

The JSON the page exports: `schema: "discuspdf-disorder-builder"`, `version: 1`. Orientation groups arrive as
explicit rotation matrices with weights; the page expands them from Euler angles and symmetry operations, so no
symmetry parsing happens here. `builder_state` (the page's own state, for re-import) is carried along and
ignored. A model with a higher version is refused with a request to update discus-python.

## Atom names (`discuspdf.atom_names`)

DISCUS identifies atoms by name. `discus_atom_name(name)` returns the table entry DISCUS uses: a name it knows
(225 names, including ions, `D` and `VOID`), or else its guess from the name up to a `+`/`-`, then the first two
letters, then the first letter. The table is generated from DiffuseCode `lib_f90/element_data_mod.f90`
(`_discus_names.py` records the source's SHA-256). `element_of`, `atomic_number`, `mass_number` and
`is_isotope_label` follow from it. The page's atom-name table is generated from this module.

## Compatibility with the builder's engine 1.0

The code was moved from the builder's self-contained `engine.py` (builder @ 92fe94e) without changing what it
computes. For the same model and seed, `supercell.stru`, `supercell.xyz`, the summary, every staged 3D and mmc
macro and every staged structure file are byte-identical to engine 1.0. This was checked for all 19 parity
models: every builder example in three calculation variants, and a coverage model with VOID padding,
displacements, site symmetry and a centre-of-mass anchor. `tests/test_builder.py` pins five of them. The 1D PDF
macro differs only in the staged file name (`stru model.stru` instead of `stru supercell.stru`).

Changes: the default result files (above); the `--numpy-preview` uses `reference_scattering` with real
scattering factors (neutron b from DISCUS's table, X-ray f₀ from `xraydb` when installed, else f = Z) instead
of f = Z; and the run log now ends with its "Done" line.
