# Gate 3 pre-registered experiment: PASS

predictions.json sha256 69b7f4751dd2c0a4f297529194ccefa6f70460cdd2c7f39538e9438fa8ecadfe

0/15 criteria failed.

| criterion | value | limit | pass |
|---|---|---|---|
| R0 diffuse: int volume content hash equals Gate 2 v4 (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| M1 diffuse: native corner/increment equal source prediction (max abs dev) | 0.000e+00 | <= 1.000e-09 | PASS |
| M1 diffuse: native data_type_style == single_pdf | 1.000e+00 | >= 1.000e+00 | PASS |
| E1 diffuse: native vs emulation of same-run int (relL2) | 3.041e-16 | <= 1.000e-12 | PASS |
| E2 diffuse: native matches registered fingerprint (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| P1 diffuse: native argmax at predicted index [64, 64, 64] (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| C1 diffuse: package delta_pdf_fractional == k * native on true axes (relL2) | 5.083e-16 | <= 1.000e-12 | PASS |
| L diffuse: log atom counts, no errors, Temp. factors not 'used' | 1.000e+00 | >= 1.000e+00 | PASS |
| R0 total: int volume content hash equals Gate 2 v4 (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| M1 total: native corner/increment equal source prediction (max abs dev) | 0.000e+00 | <= 1.000e-09 | PASS |
| M1 total: native data_type_style == single_pdf | 1.000e+00 | >= 1.000e+00 | PASS |
| E1 total: native vs emulation of same-run int (relL2) | 2.933e-16 | <= 1.000e-12 | PASS |
| E2 total: native matches registered fingerprint (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| P1 total: native argmax at predicted index [64, 64, 64] (fraction) | 1.000e+00 | >= 1.000e+00 | PASS |
| L total: log atom counts, no errors, Temp. factors not 'used' | 1.000e+00 | >= 1.000e+00 | PASS |

```json
{
 "diffuse": {
  "int_content_hash": "756d781ee47bafe9024efa8f30a54e589548fbd6dc221f3ddd8636fe97033fa1",
  "four_weight_scalef_used": null,
  "native_argmax_index": [
   64,
   64,
   64
  ],
  "native_label_at_argmax": [
   0.1259842519685037,
   0.1259842519685037,
   0.1259842519685037
  ],
  "log": {
   "ok": true,
   "counts": {
    "AU": 2048,
    "CU": 2048
   },
   "errors": [],
   "temp": "ignored"
  }
 },
 "total": {
  "int_content_hash": "e76cea5fa13e857b96f0f27ca118c51d35347391007c360a26d7bc5e60eea964",
  "four_weight_scalef_used": 1.0,
  "native_argmax_index": [
   64,
   64,
   64
  ],
  "native_label_at_argmax": [
   0.1259842519685037,
   0.1259842519685037,
   0.1259842519685037
  ],
  "log": {
   "ok": true,
   "counts": {
    "AU": 2048,
    "CU": 2048
   },
   "errors": [],
   "temp": "ignored"
  }
 }
}
```
