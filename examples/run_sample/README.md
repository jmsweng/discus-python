# Sample DISCUS 1D PDF run

This directory contains a **real** DISCUS_SUITE 1D powder-PDF run produced by
`discuspdf.PDFRunner` (GaAs, 5×5×5 supercell, r = 0–12 Å in 0.02 Å steps, X-ray).

| File | Description |
|---|---|
| `run.mac` | The DISCUS macro that was executed (v6 `powder`→`output` PDF workflow) |
| `model.cll` | The input unit-cell file (copied into the run directory) |
| `config.json` | Recorded parameters + SHA-256 of the input + launcher used |
| `run.log` | DISCUS stdout/stderr for the run |
| `status.json` | Runner status (`complete`, return code, output presence) |
| `result.pdf` | The numeric 1D PDF output (two columns: `r` in Å, `G(r)`) — **not** an Adobe PDF |
| `MANIFEST.sha256` | SHA-256 of the six run files above (checked in CI) |

Regenerate an equivalent run with:

```python
from discuspdf import PDFConfig, PDFRunner
cfg = PDFConfig("examples/gaas.cll", supercell=(5, 5, 5),
                rmax=12.0, dr=0.02, qmax=30.0, qsig=0.0, radiation="xray")
job = PDFRunner(runs_directory="runs").run(cfg)
r, gr = job.load()
```

See `docs/VERIFIED_DISCUS_MACRO.md` for the verified macro syntax.

**Note (2026-10-31):** although the macro says `set qmax,30`, DISCUS 6.23.01 cannot use Q above about 8.0 Å⁻¹ at λ = 1.5406 Å, so this G(r) is the same as for Qmax ≈ 8 Å⁻¹. Current `discuspdf` produces the same macro and result but issues a `QmaxLimitWarning`; see `docs/VERIFIED_DISCUS_MACRO.md`, "Qmax limit and wavelength".

## Provenance (re-checked 2026-11-02)

- **One coherent run.** `run.mac`, `run.log` and `result.pdf` belong to the same DISCUS_SUITE 6.23.01 run of 2026-10-22, made with `qmax=30`. `run.log` echoes every command of `run.mac`, including `outf result.pdf`.
- **`run.mac` was restored.** Commit `6ff7fac` (switch to CSV output) edited the archived `run.mac` by hand to `outf result.csv` without rerunning. The file now has its original content again (blob `a8453bf`). It is byte-identical to what `render_pdf_macro` writes today for the parameters in `config.json`.
- **`config.json` is kept exactly as the runner wrote it.** Its `structure` path is the absolute path on the machine that made the run. The structure is identified by `input_sha256`, which equals both `model.cll` and the repository's `examples/gaas.cll`. Older runners did not record `qmin`, `dq`, `wavelength` or `macro_sha256`. The macro shows the values used: qmin 0.5, dq 0.01, λ 1.5406 Å.
- **Engine binary not recorded.** The run predates binary hashing, so the exact DISCUS binary is unknown. The gate evidence of 2026-10-28/29 used SHA-256 `d4274f13…45706f3`.
- **Checks.** `tests/test_run_sample_fixture.py` (standard CI, no DISCUS) checks the manifest, the macro against the generator and against the log, the structure hash and the r grid. `tests/test_integration_discus.py` (opt-in, `DISCUS_INTEGRATION=1`) reruns the recorded settings and compares G(r) with `result.pdf` (relative L2 < 1e-6).
- **qmax 30 vs 25.** Both are above the DISCUS limit at this wavelength and give identical G(r) (see the Qmax note above), so the test's earlier `qmax=25` would not have changed the numbers. The test now uses the recorded 30.

## Engine re-check (2026-11-02)

`DISCUS_INTEGRATION=1 python -m pytest -q -rP tests/test_run_sample_fixture.py tests/test_integration_discus.py` on the owner's WSL machine, fresh LF clone of `fix/gaas-fixture-provenance` at `d770b74`: **7 passed**. DISCUS_SUITE 6.23.01 rerun of the recorded settings (qmax 30) against `result.pdf`: same 601-point r grid, **relative L2 = 2.454×10⁻¹⁶** (limit 1×10⁻⁶), i.e. identical to floating-point precision. Engine `/usr/local/bin/discus_suite`, SHA-256 `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3`, the same binary as the Gate 1–3 evidence runs.
