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
| `result.csv` | The numeric 1D PDF output (two columns: `r` in Å, `G(r)`) |

Regenerate an equivalent run with:

```python
from discuspdf import PDFConfig, PDFRunner
cfg = PDFConfig("examples/gaas.cll", supercell=(5, 5, 5),
                rmax=12.0, dr=0.02, qmax=30.0, qsig=0.0, radiation="xray")
job = PDFRunner(runs_directory="runs").run(cfg)
r, gr = job.load()
```

See `SETUP_REPORT.md` and `docs/VERIFIED_DISCUS_MACRO.md` for the full setup and the
verified macro syntax.
