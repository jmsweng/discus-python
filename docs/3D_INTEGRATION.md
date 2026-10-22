# Real 3D diffuse-scattering integration (not yet scientifically verified)

`Diffuse3DRunner` stages an **independently validated DISCUS_SUITE macro**
and all its referenced input files into a unique run directory. It expects
the macro to produce `diffuse.h5` and saves logs and input hashes. This is
an execution adapter, not a claim that a valid DISCUS 3D macro has been
established. HDF5 data paths must be inspected from **actual** installed-version
output; no dataset paths are invented.

Example **only after** independently validating the macro and file format:

```python
from discuspdf import Diffuse3DRunner, OutputSpec
runner = Diffuse3DRunner()
job = runner.prepare("verified_3d.mac", assets=("model.cll",))
# Without a spec the run is recorded as "produced_unverified":
job = runner.execute(job)
# With a spec validated against the actual file, the run is "complete":
job = runner.execute(job, output_spec=OutputSpec(
    intensity_dataset="/ACTUAL/PATH",
    q_datasets=("/ACTUAL/QX", "/ACTUAL/QY", "/ACTUAL/QZ"),
    expected_shape=(128, 128, 128),
))
vol = job.load_volume(
    intensity_dataset="/ACTUAL/PATH",
    q_datasets=("/ACTUAL/QX", "/ACTUAL/QY", "/ACTUAL/QZ"),
    coordinate_system="cartesian_q",
    units="1/angstrom",
    preprocessing=("bragg_excluded",),  # only if independently confirmed
)
```

The sample paths are **placeholders**, and this example must not be called
a validated workflow. Confirm Bragg removal/exclusion, reciprocal-lattice
metric, symmetry, masks, zero-centering, and normalization with real reference
outputs before quantitative interpretations.

## Output validation and status semantics

`Diffuse3DRunner.execute(run, output_spec=None)` validates the produced
`diffuse.h5` in two layers:

* **Structural (always):** the file must open as HDF5.
* **Declarative (only with an `OutputSpec`):** the caller-supplied
  `OutputSpec(intensity_dataset, q_datasets, expected_shape=None)` names the
  actual dataset paths (inspect the real file; no DISCUS layout is assumed).
  The intensity dataset must be numeric, 3D, and finite; each Q dataset must
  be a finite 1D axis whose length matches the corresponding intensity
  dimension; `expected_shape`, when given, must match exactly.

Status transitions are `prepared` → `running` → one of:

* `complete` — a user-supplied `output_spec` was fully validated;
* `produced_unverified` — the file is structurally valid HDF5 but no spec was
  validated; this is **not** a claim of scientific completeness;
* `failed` — preflight (staged-file hashes, macro presence, existing output),
  engine (nonzero exit, definitive error in the log, timeout, missing binary),
  or output-validation failure. `status.json` records the return code when an
  engine process ran, the failure reason, and `engine_started` for preflight
  failures. PGPLOT/X11 noise such as `xdpyinfo`/`pgxwin_server` messages is
  not treated as fatal.

The runner does not infer units, Bragg subtraction, or diffraction
geometry. Source metadata must be verified by the operator. Scientific
benchmarking is still required.
