# Powder settings probe (powder-settings-v1): FAIL (27/30 failed)

| criterion | value | limit | pass |
|---|---|---|---|
| R delta002: engine produced G(r) | no output | – | FAIL |
| R delta005: engine produced G(r) | no output | – | FAIL |
| R delta010: engine produced G(r) | no output | – | FAIL |
| T1 base: G vs sine transform of own F(Q) (shift1_zero) | 3.666e-01 | 2e-03 | FAIL |
| T1 base: |scale - 1| | 5.701e-01 | 2e-03 | FAIL |
| T1 qmin02: G vs sine transform of own F(Q) (shift1_zero) | 4.208e-01 | 2e-03 | FAIL |
| T1 qmin02: |scale - 1| | 5.884e-01 | 2e-03 | FAIL |
| T1 qmin10: G vs sine transform of own F(Q) (shift1_zero) | 3.639e-01 | 2e-03 | FAIL |
| T1 qmin10: |scale - 1| | 5.696e-01 | 2e-03 | FAIL |
| T1 qmin20: G vs sine transform of own F(Q) (shift1_zero) | 5.629e-01 | 2e-03 | FAIL |
| T1 qmin20: |scale - 1| | 3.899e-01 | 2e-03 | FAIL |
| T1 dq005: G vs sine transform of own F(Q) (shift0_zero) | 3.636e-01 | 2e-03 | FAIL |
| T1 dq005: |scale - 1| | 5.687e-01 | 2e-03 | FAIL |
| T1 dq020: G vs sine transform of own F(Q) (shift1_zero) | 3.665e-01 | 2e-03 | FAIL |
| T1 dq020: |scale - 1| | 5.701e-01 | 2e-03 | FAIL |
| T1 dq050: G vs sine transform of own F(Q) (shift1_zero) | 3.663e-01 | 2e-03 | FAIL |
| T1 dq050: |scale - 1| | 5.701e-01 | 2e-03 | FAIL |
| T1 qmax50: G vs sine transform of own F(Q) (shift0_zero) | 3.747e-01 | 2e-03 | FAIL |
| T1 qmax50: |scale - 1| | 5.662e-01 | 2e-03 | FAIL |
| T1 qmax60: G vs sine transform of own F(Q) (shift1_zero) | 3.606e-01 | 2e-03 | FAIL |
| T1 qmax60: |scale - 1| | 5.701e-01 | 2e-03 | FAIL |
| T3 qmin02: F(Q) at 701 common Q vs base | 2.796e-02 | 1e-06 | FAIL |
| T3 qmin10: F(Q) at 651 common Q vs base | 1.005e-02 | 1e-06 | FAIL |
| T3 qmin20: F(Q) at 551 common Q vs base | 4.738e-01 | 1e-06 | FAIL |
| T3 dq005: F(Q) at 701 common Q vs base | 5.340e-02 | 1e-06 | FAIL |
| T3 dq020: F(Q) at 351 common Q vs base | 1.069e-01 | 1e-06 | FAIL |
| T3 dq050: F(Q) at 141 common Q vs base | 4.175e-01 | 1e-06 | FAIL |
| T3 qmax50: F(Q) at 451 common Q vs base | 4.120e-16 | 1e-06 | PASS |
| T3 qmax60: F(Q) at 551 common Q vs base | 7.416e-16 | 1e-06 | PASS |
| T4 base_plain vs base: G(r) | 4.402e-16 | 1e-12 | PASS |

## Effect of each setting on G(r): max|G - G_base| / max|G_base|

| case | settings | effect |
|---|---|---|
| qmin02 | qmin 0.2, qmax 7.5, dq 0.01, delta 0.0 | 1.423e-01 |
| qmin10 | qmin 1.0, qmax 7.5, dq 0.01, delta 0.0 | 2.041e-02 |
| qmin20 | qmin 2.0, qmax 7.5, dq 0.01, delta 0.0 | 4.356e-01 |
| dq005 | qmin 0.5, qmax 7.5, dq 0.005, delta 0.0 | 1.283e-03 |
| dq020 | qmin 0.5, qmax 7.5, dq 0.02, delta 0.0 | 2.571e-03 |
| dq050 | qmin 0.5, qmax 7.5, dq 0.05, delta 0.0 | 1.046e-02 |
| qmax50 | qmin 0.5, qmax 5.0, dq 0.01, delta 0.0 | 5.847e-01 |
| qmax60 | qmin 0.5, qmax 6.0, dq 0.01, delta 0.0 | 3.773e-01 |
