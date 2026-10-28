# S1 probe protocol — real DISCUS_SUITE 3D output (for the WSL2 operator)

Historical probe protocol: this text was **written before the experiments**. Parts were subsequently executed on WSL2; see [FINDINGS.md](FINDINGS.md) for what was actually observed and [the reproduction plan](../docs/SCIENTIFIC_BENCHMARK_REPRODUCTION.md) for what is still needed. The original baseline was
`main` = `b5f98d350ca8eb7aeb7775b979cfd2164b2fc2ed` (local working tree was on
`fix/post-merge-review` = `abd52f0`). Claude's workspace had no `discus_suite`
(`which discus_suite` returned nothing on the Linux VM), at the time, so this protocol alone makes no 3D execution claim. Later real-engine observations are documented in FINDINGS.md. The steps below produce the observations
needed to write it correctly.

## Rules
- Run in Ubuntu 22.04 / WSL2 with the project venv; do not mix Windows Python.
- No `sudo`, installs, or large sweeps. Cap each calculation at roughly 64^3 voxels
  and a few minutes until sizes are known; record runtime and memory.
- Keep every log, macro, seed and hash. Failed attempts are evidence: keep them.

## Step 0 — environment and documentation (`probe_environment.sh`)
`bash probe_environment.sh` records OS/WSL, binary path/hash/linked libraries
(is HDF5 linked?), the version banner, the manual PDFs present, keyword hits for
HDF5 / LOTS / Bragg / 3D output, and the in-suite `help` text of the `fourier`
menu. Return the output directory.

## Step 1 — read the version-matched sources (operator or Claude, from Step 0 output)
Answer and cite section/page for each, for **6.23.01** only:
1. Which command computes a 3D reciprocal-space intensity on a grid, and how are
   the grid corners/steps (`ll`, `lr`, `lu`, `lv`-style settings or equivalents)
   and number of points set? Is the grid in hkl (rlu) or in Å⁻¹?
2. Which output formats exist for that grid (check for HDF5), what file name does
   the macro choose (never assume `diffuse.h5`), and which dataset paths/attributes
   does the file contain?
3. How are Bragg peaks selected or excluded (LOTS / ensemble averaging /
   subtracting the average structure / separate Bragg-only calculation)?
4. How is disorder generated (e.g. occupational or displacive via the `mmc`,
   `chem`, `surface`/`thermal`-type tools) and how is the random seed set?

## Step 2 — three small calculations
Candidate macros are in `macros/case.mac.in`; `bash run_cases.sh` generates and runs C0/C1/C2 (GaAs 8x8x8, random Ga->As antisites p=0.15, seeds 1 and 2; 65^3 hkl grid, -4..4), writes total (Bragg present, `set aver,0.0`), Bragg-subtracted (`set aver,100.0`) and DISCUS 3D-PDF outputs in both `nexus` and `hdf5` formats, and inspects every .h5. The macro syntax is taken from the manual and is UNVERIFIED until it runs.
| id | Model | Purpose |
|---|---|---|
| C0 | perfect crystal, non-cubic or low-symmetry cell if feasible | Bragg-only control: diffuse intensity should be ~0 away from Bragg positions |
| C1 | the same cell with one simple, interpretable disorder mechanism (e.g. random substitution at one site, known short-range correlation) | genuine diffuse scattering; pair correlation known by construction |
| C2 | C1 with a different seed, or a doubled supercell | variance / finite-size check |

For each run record: macro, input structure(s) and SHA-256, seed, cell and
supercell, grid definition, wall time, peak memory, `discus_suite` stdout/stderr.

## Step 3 — inspect outputs without interpreting them
`python inspect_output.py <outputfile> --json out_C0.json` (works for HDF5 and for
other formats). Return the JSON for C0, C1, C2. Needed facts: file name and
format, dataset paths, shapes, dimension order, coordinate arrays or attributes,
units attributes (if truly present), min/max/finite counts, where the maximum
sits, and whether the grid center is index 0 or `n//2`.

## Step 4 — return for Claude
Directory tarball from Step 0, the three macros, the inspector JSON, and the
answers to Step 1 with citations. Claude then builds the adapter, Q-geometry and
benchmark from these observations. Until they exist, S1 is **not complete**.

## Acceptance for S1
A real macro and real output reproduced by the local agent, with paths, logs,
array shape and sampling documented (handoff §1).
