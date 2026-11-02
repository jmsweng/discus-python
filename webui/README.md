# DISCUS Disorder Builder

A static web page for building disordered crystal models and exporting a Python script that runs them through
[discus-python](https://github.com/jmsweng/discus-python) and DISCUS_SUITE 6.23.01.

## What it does

- **Cell and crystal**: lattice parameters (discus-python convention, a ∥ x, b in xy) and supercell size.
- **Molecule builder**: atoms in Cartesian Å, templates (H2O, CO2, NH4+, CH4, C6H6, SF6, NO3-, …), pasted XYZ, or
  click-to-draw in a rotatable 3D view; anchor at centroid, centre of mass or one atom.
- **Orientation library**: a group is a base orientation (ZYZ / ZXZ / XYZ Euler angles) plus symmetry operations
  (`C4 [0 0 1]`, `S4 [001]`, `m [1 -1 0]`, `-1`, `R 30 [1 1 1]`, Jones symbols `-y,x,z`, explicit matrices), with
  optional group closure. Placements that coincide because of the molecule's own symmetry are merged. Each
  orientation has a relative weight.
- **Sites**: each site holds one occupant per cell (atom, molecule in one orientation of a group, or vacancy) with
  given probabilities, drawn independently or with exact composition.
- **Disorder**: correlated occupants (composition-conserving swap Monte Carlo with pair energies along neighbour
  vectors; works for atoms, vacancies and molecule orientations), rigid-body translation and libration of
  molecules, Gaussian atomic displacements, and a contact check.
- **DISCUS Monte Carlo** (atom-only models): neighbourhoods, occupational / displacement / spring / Lennard-Jones
  targets and moves, through `discuspdf.mc.MCConfig`.
- **Calculations**: 3D diffuse scattering with Bragg subtraction (`FourierConfig`, `set aver`, LOTS), 3D-ΔPDF by the
  fractional FFT, optional total scattering and 3D-PDF (Patterson), and the 1D powder PDF of the supercell.
- **Preview**: a rotatable view of a random draw of the crystal (coloured by orientation, element or site) and an
  in-browser direct-sum diffuse plane (f = Z, averaged over several random crystals).

## The generated script

One file, the model embedded as JSON (re-import the `.py` or the exported `.json` into the page).

```
pip install "discuspdf[data] @ git+https://github.com/jmsweng/discus-python"
python model.py --dry-run          # build the crystal and stage every DISCUS run, no DISCUS needed
python model.py --numpy-preview    # also a numpy diffuse plane
python model.py                    # full run; --launcher "…" to change how DISCUS is started
```

Outputs: `supercell.stru`, `supercell.xyz`, `structure_summary.json`, `diffuse.npz` + slice plots,
`delta_pdf.npz` / `delta_pdf.vts`, optionally `total.npz` / `pdf3d_patterson.vts`, `pdf_1d.txt` / `.png`, and one
run folder per DISCUS job with its macro, inputs and log.

How the crystal reaches DISCUS: every unit cell lists the same number of atoms (DISCUS `ncell` layout). When a site's
occupants need different atom counts (molecules of different size, vacancies), the unused slots are filled with
DISCUS's non-scattering `VOID` atoms placed where the site's largest occupant would put its atoms.

## Files

- `index.html` – the whole app in one file, built by `build.py` (committed so it also opens straight from a checkout).
- `src/` – sources: `core.js` (maths, orientations, realisation, preview, export), `viewer.js`, `app.js`,
  `examples.js`, `style.css`, `shell.html`.
- `engine.py` – the Python template the page fills in.
- `build.py` – inlines `src/` and `engine.py` into `index.html`. Run it after editing anything in `src/` or `engine.py`.
- `tests/node_test.js` – Node checks of the core; writes the generated example scripts to `$OUT` (default: the temp
  folder), which you can then run with `--dry-run`.

## Hosting on GitHub Pages

The workflow `.github/workflows/pages.yml` rebuilds `index.html` and deploys it on every push to `main` that touches
`webui/`. One-time setup: repository Settings → Pages → Build and deployment → Source: **GitHub Actions**. The site
is then at `https://<owner>.github.io/discus-python/`. The page is fully static; it only loads Google Fonts.

## Not yet verified against DISCUS

The structure building, staged macros and file formats were tested with dry runs; the DISCUS runs themselves have not
been executed yet. Points to check on the first live run: `VOID` padding atoms in `read stru`, the 1D PDF macro with
`read stru` instead of `read cell`, and the Windows → WSL launcher.
