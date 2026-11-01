# Monte Carlo disorder models with DISCUS MMC (Package D, Features 15–17)

`discuspdf.mc` writes, runs and reads DISCUS 6.23.01 `mmc` (multiple-energy Monte Carlo) simulations from
Python. You describe the run with an `MCConfig`, read the generated macro, run it with `MCRunner`, and get an
`MCResult` with the final structure, DISCUS's feedback history, and an independent measurement of every target
with Package C. Only behaviour that was read in the DISCUS source **and** checked by the engine probe
([`d_mc_probes/RESULT.md`](../benchmarks/discus_6_23_01/d_mc_probes/RESULT.md)) is offered.

```python
from discuspdf import read_stru
from discuspdf.fourier_macro import Seed
from discuspdf.mc import MCConfig, MCRunner, Neighborhood, OccupationalTarget, Swap

nn = Neighborhood.from_translations([(1, 0, 0), (0, 1, 0), (0, 0, 1)])     # ±a, ±b, ±c (reverse added)
cfg = MCConfig(neighborhoods=(nn,),
               targets=(OccupationalTarget(0, "Cu", "Au", correlation=-0.3),),
               moves=(Swap(),), seed=Seed((11, 22, 33, 44, 55, 66, 77, 88)),
               moves_per_atom=10, feedback_per_atom=1)
print(cfg.render("mc_in.stru", n_atoms=4096))     # read the macro first
print(cfg.cautions())                             # known engine behaviour that affects this run

runner = MCRunner("runs_mc")                       # launcher: discus_suite (WSL on Windows, see below)
result = runner.execute(runner.prepare(cfg, structure))
result.target_report(tol=0.02)                     # DISCUS's final target vs achieved, 'reached'
result.measure()                                   # the same quantities measured with Package C
result.history(0), result.acceptance()             # feedback trajectory, accepted good/neutral/bad
```

## What DISCUS's MMC does (and what the API makes of it)

| Fact (source + probe) | Consequence in `discuspdf.mc` |
|---|---|
| Acceptance: downhill always; uphill **and ΔE = 0** with `1/(1+e^{ΔE/kT})` (heat bath, not Metropolis; ΔE rounded to 1e-4). Probe: J = 0 → exactly ½; kT → ∞ → ½ of non-downhill moves | stated in every run's `cautions`; choose energies with typical ΔE ≫ 1e-4 |
| `set cyc` counts **valid trial moves**; same-species swap draws are redrawn, not counted | `moves_per_atom × N`; the result must report exactly that many (`Gen`) or the run fails |
| Same 8-value seed + serial → byte-identical; OpenMP → irreproducible and more moves than requested | every macro has `seed` (8 values, required) and `set parallel, useomp:off` |
| `set vec` pairs are directed (origin → target only) | `Neighborhood(symmetric=True)` adds the reverse vectors |
| Default shift σ is 0 | `Shift` requires σ per species |
| `CORR`: energy ignored, PID feedback with damping; scaled ×0.005 for \|c\|<0.1, J = 0 for \|c\|<0.001 | at most one feedback target per run; cautions for small targets |
| Lennard-Jones needs four numbers | always written |
| `set mode` with an atom list sets the wrong probability (bug) | moves never use atom lists |
| `set valid` is a no-op; `set finish` thresholds unverified | not offered; convergence is judged from the history and `target_report` |
| `save` drops `ncell`, upper-cases names | `MCRunner.load` restores supercell and spelling |

## Targets and moves

| Class | DISCUS | Energy (units of kT) |
|---|---|---|
| `OccupationalTarget(nb, a, b, correlation=c \| energy=J)` | `corr … CORR/ENER` | J·(N_like − N_unlike) per centre; J < 0 clusters, J > 0 orders |
| `DisplacementTarget(nb, a, b, direction, correlation= \| energy=)` | `cd … dir:[u,v,w]` | J·(u_c·d̂)(u_n·d̂) (Å²), counted only when the moved atom is the vector origin |
| `SpringTarget(nb, a, b, distance, k)` | `spring` | k (d − d0)², averaged over the centre's neighbours |
| `LennardJonesTarget(nb, a, b, distance, depth, n_rep=12, n_att=6)` | `lennard` | minimum −\|depth\| at d0, averaged over the neighbours |
| `Swap(p, pairing)` / `SwapDisplacement(p, pairing)` / `Shift(sigma, p)` | `swchem` / `swdisp` / `shift` + `set move` | `pairing`: all, local, slocal, site; σ in unit-cell fractions |

## Results and diagnostics

* `status` is `complete` only when DISCUS exited 0, printed no error, did all requested moves and wrote the
  structure; otherwise `MCRunner` raises and `status.json` says why. **`complete` does not mean the targets were
  reached**: use `target_report(tol)` (feedback targets get `reached` True/False; fixed-energy targets `None`).
* `measure()` recomputes each target on the final structure with Package C — occupational α (pooled formula,
  equal to DISCUS's printed value for binary systems with symmetric pair counts), uncentred displacement
  correlation from the site averages, mean and spread of spring/LJ pair distances. On the probe runs it matches
  DISCUS's printed values to their 3 decimals.
* `replica_summary(results)` gives the seed-to-seed mean and spread of achieved and measured values for runs
  that differ only in the seed.

## Running on Windows
DISCUS runs in WSL: `MCRunner(runs_dir, launcher=("wsl.exe", "-d", "Ubuntu-22.04", "--", "discus_suite",
"-noautorun", "-macro"))`; the run directory path must not contain spaces.

## Not covered
Molecules, angle/environment/connectivity neighbourhoods, `unid`, `cn`, `group`, `heli`, `bucking`, `repulsive`,
more than one feedback target, `set finish` convergence, parallel runs, partial occupancies (vacancies as a species
are untested).
