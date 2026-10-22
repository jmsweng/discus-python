# Real 3D diffuse-scattering integration (not yet scientifically verified)

`Diffuse3DRunner` stages an **independently validated DISCUS_SUITE macro**
and all its referenced input files into a unique run directory. It expects
the macro to produce `diffuse.h5` and saves logs and input hashes. This is
an execution adapter, not a claim that a valid DISCUS 3D macro has been
established. HDF5 data paths must be inspected from **actual** installed-version
output; no dataset paths are invented.

Example **only after** independently validating the macro and file format:

```python
from discuspdf import Diffuse3DRunner
runner = Diffuse3DRunner()
job = runner.prepare("verified_3d.mac", assets=("model.cll",))
job = runner.execute(job)
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

The runner does not infer units, Bragg subtraction, or diffraction
geometry. Source metadata must be verified by the operator. Scientific
benchmarking is still required.
