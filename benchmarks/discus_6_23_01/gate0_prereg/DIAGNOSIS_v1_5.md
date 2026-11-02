# Gate 0 v1.5 — diagnosis of the FAIL (run 20261102T063824_LAPTOP-NGU3BT0N)

## What happened

v1.5 run `20261102T063824_LAPTOP-NGU3BT0N`: verdict **FAIL** only on **E4 B** ("log shows B values", 0/1 for total and diffuse). Every other criterion passed in both arms (E1 ≤ 8.5e-15, E2 ≤ 8.5e-15, E3 = 0, S1 7.1e-3, S2 1.5e-2, S3 0/15 561). Cause: the analyzer's E4 regex looks for lines like `Cu B = 0.5`; DISCUS 6.23.01 prints no per-atom B values — its log says `Temp. factors      : used` (arm B) / `ignored` (arm A). The selftest wrote **synthetic** logs containing `Cu B = 0.5`, so the selftest could never reveal that the criterion is unsatisfiable on real output. (That the engine used exactly B = 0.5/0.8 is already shown by E1/E2 for arm B at 7e-15 against an oracle computed with those B — a wrong B would give percent-level errors.)

## Verbatim E4 lines of the analyzer (`gate0_analyze.py`)

```python
    temp = re.search(r"Temp\. factors\s*:\s*(\w+)", txt)
    # Look for B values in the log
    b_cu = re.search(r"B\s*=\s*([\d.]+).*Cu|Cu.*B\s*=\s*([\d.]+)", txt)
    b_au = re.search(r"B\s*=\s*([\d.]+).*Au|Au.*B\s*=\s*([\d.]+)", txt)
```

and the E4 decision:

```python
        if arm == "A":
            ok = lg["temp"] is not None and lg["temp"].lower() == "ignored"
            crit.append(_crit(f"E4 {arm} {kind}: log shows 'temp ignore'", float(ok), 1.0, ">="))
        else:
            # Check B values are present in the log
            ok = lg["b_cu"] is not None and lg["b_au"] is not None
            crit.append(_crit(f"E4 {arm} {kind}: log shows B values", float(ok), 1.0, ">="))
```

## Verbatim `Temp. factors` lines from the committed v1.5 logs

- `gate0_B_total_00.log`: `   Temp. factors      : used`
- `gate0_B_diffuse_00.log`: `   Temp. factors      : used`
- `gate0_A_total_00.log`: `   Temp. factors      : ignored`
- `gate0_A_diffuse_00.log`: `   Temp. factors      : ignored`

(No line of the form `B = <number>` appears anywhere in any arm-B log; the only `temp`-related lines are `discus/fourier> temp use` and `Temp. factors      : used`.)

## Criterion table (from `GATE0_RESULTS_PREREG.md`)

| criterion | value | limit | pass |
|---|---|---|---|
| E1 A: max total vs realization engine oracle (relL2) | 8.502e-15 | <= 1.000e-09 | PASS |
| E2 A: max diffuse vs realization engine oracle (relL2) | 7.833e-15 | <= 1.000e-06 | PASS |
| E3 A: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 A total: log shows 'temp ignore' | 1.000e+00 | >= 1.000e+00 | PASS |
| E4 A diffuse: log shows 'temp ignore' | 1.000e+00 | >= 1.000e+00 | PASS |
| S1 A: |mean x - 1| | 7.115e-03 | <= 3.937e-02 | PASS |
| S2 A: |var x - 1| | 1.501e-02 | <= 1.114e-01 | PASS |
| S3 A: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 A: max |E_phys - E_eng| / E_eng | 1.038e-03 | <= 5.000e-03 | PASS |
| E1 B: max total vs realization engine oracle (relL2) | 6.662e-15 | <= 1.000e-09 | PASS |
| E2 B: max diffuse vs realization engine oracle (relL2) | 8.471e-15 | <= 1.000e-06 | PASS |
| E3 B: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 B total: log shows B values | 0.000e+00 | >= 1.000e+00 | **FAIL** |
| E4 B diffuse: log shows B values | 0.000e+00 | >= 1.000e+00 | **FAIL** |
| S1 B: |mean x - 1| | 7.115e-03 | <= 3.937e-02 | PASS |
| S2 B: |var x - 1| | 1.501e-02 | <= 1.114e-01 | PASS |
| S3 B: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 B: max |E_phys - E_eng| / E_eng | 1.899e-03 | <= 5.000e-03 | PASS |
