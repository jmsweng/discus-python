# Running DISCUS_SUITE from Python

Three pieces make every DISCUS_SUITE command usable from Python, including those without a
typed discuspdf API:

| Module | What it does |
|---|---|
| `discuspdf.engine` | `Engine`: one DISCUS_SUITE installation (launcher, `stock` or `fork`). `Engine.run` is the single way the package starts the engine. `fork_reasons` / `choose_engine` decide which build a macro needs. |
| `discuspdf.macro` | `Macro`: builds a macro line by line, with menus as `with` blocks; `render()` checks it against the [command catalogue](SUITE_CATALOGUE.md). |
| `discuspdf.suite_runner` | `SuiteRunner`: prepares a run directory (macro, inputs, hashes), runs it, and checks the declared outputs. |

```python
from discuspdf.engine import Engine
from discuspdf.macro import Macro
from discuspdf.suite_runner import SuiteRunner

m = Macro("GaAs 2x2x2 supercell")
with m.program("discus"):
    with m.menu("read"):                 # read takes one command and returns by itself
        m.command("cell", "gaas.cll", 2, 2, 2)
    m.command("save", "supercell.stru")  # 'save <file>' writes at once (no menu)

runner = SuiteRunner("runs_suite", Engine(("wsl.exe", "-d", "Ubuntu-22.04", "--",
                                           "discus_suite", "-noautorun", "-macro")))
run = runner.prepare(m, {"gaas.cll": "tests/fixtures/discus_6_23_01/gaas.cll"}, ["supercell.stru"])
print(run.macro.read_text())             # inspect before running
runner.execute(run)                      # status.json: complete, with output hashes
```

## The engine

`Engine.run(directory, "run.mac")` starts DISCUS the way every validated run in this
repository did:

- in the run directory, with `exit` on stdin;
- with no `DISPLAY` and `PGPLOT=none`, because otherwise `pgxwin_server` can keep the output
  pipe open after DISCUS has exited;
- with the log written to `run.log` as stdout, a `--- STDERR ---` line, then stderr.

A launcher that cannot be found (`FileNotFoundError`) or a timeout
(`subprocess.TimeoutExpired`) is raised to the caller.

`PDFRunner`, `Diffuse3DRunner` and `MCRunner` take `engine=` and use `Engine.run`. Their
`launcher=` and `timeout_seconds=` arguments still work and build a stock engine. Their
macros, staged files, statuses and error checks are unchanged.

`Engine.from_env("stock")` reads the launcher from `DISCUSPDF_ENGINE`, defaulting to
`discus_suite -noautorun -macro`. `Engine.from_env("fork")` reads `DISCUSPDF_ENGINE_FORK`;
the fork has no default.

## Stock or fork

The rule, agreed with the owner, is: stock DISCUS 6.23.01 unless the macro uses something
only the fork computes correctly. `fork_reasons(macro)` lists those uses:

| Use | Why the fork |
|---|---|
| REFINE `run` | MPI derivative steps, `MAX(dP/sig)` criterion, seed restore with Gaussian noise, MPI result transfer |
| REFINE `newparam …,points:2` | exists only on the fork |
| RMC `set lots,…` (not `off`) | lot moves, updates and averages |
| Fourier `set technique:nufft` | coordinate scaling on non-integer grids |
| a second Fourier `run` after `temp use` in one session | Debye–Waller cache key; a single run is correct on stock (Gate 0 v2.1) |

`SuiteRunner` uses the fork when the list is not empty. If the runner has no fork, it refuses
to prepare the run (`EngineChoiceError`). `engine="stock"` or `engine="fork"` forces a
choice. The chosen engine and the reasons are recorded in `run.json`.

The specialised runners do not route. Their macros contain none of these features except
`temp use`, which they run once per session.

## The macro builder

- `command(word, *params, **options)` writes `word p1,p2,name:value`. Integers are written as
  they are, floats with up to 12 significant digits, and strings unchanged.
- The command word must exist in the current menu. Abbreviations are resolved the way DISCUS
  resolves them.
- `menu(word, …)` is a `with` block that enters the sub-menu and writes `exit` at its end. The
  `read` menu takes exactly one command and needs no `exit`. A command that would open a menu
  must be written as a `menu(...)` block.
- `raw(line)` writes a line unchanged. It is the escape hatch for assignments, `do`/`if`
  blocks and `@macro` calls. Raw lines are still checked for unknown commands, but menu
  changes they cause are not tracked.
- `render()` walks the finished macro through the catalogue. It refuses unknown commands, lines
  that the engine would read in a different menu than the code intends, and a macro that ends
  inside a menu.

The builder checks names and menus, not what the arguments mean. Whether a macro computes
the right thing is only known from an engine run.

## The general runner

`SuiteRunner.prepare(macro, inputs, outputs)`:

- writes `run.mac` (LF line endings);
- copies the inputs under the names the macro uses;
- records SHA-256 of the macro and the inputs, the engine and the fork reasons in `run.json`;
- sets the status to `prepared`.

It refuses unknown commands, file names with paths or spaces, reserved names and duplicate
outputs.

`execute(run)`:

- re-checks every hash and refuses outputs that already exist;
- runs the engine;
- fails on a nonzero exit, an engine error in the log (ANSI colours removed first) or a
  missing or empty output.

`complete` means the engine ran cleanly and wrote the files, not that the results are
physically right.

## Live checks

`tests/test_suite_runner_live.py` (opt-in, `DISCUS_LIVE=1`) has four small cases:

- a DISCUS read/save on the 2-atom GaAs cell;
- a KUPLOT function written with `ksave`;
- a run that must fail on a missing input;
- the read/save case on the fork, when `DISCUSPDF_ENGINE_FORK` is set.
