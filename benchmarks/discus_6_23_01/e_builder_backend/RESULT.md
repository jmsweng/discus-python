# Builder-backend probe: result of v1 (DISCUS_SUITE 6.23.01)

**10 of 11 criteria pass; P3 fails.** Run `20261103T113126Z_LAPTOP-NGU3BT0N` (evidence in
`evidence/20261103T113126Z_LAPTOP-NGU3BT0N/`, binary SHA-256 `d4274f13…45706f3`, WSL Ubuntu-22.04). The kit was
registered at `d2159a3` (PR #69). All 12 runs exited 0 and no engine markers were found (I0). The analyzer was
re-run in the cloud on the committed evidence and gave the same verdicts; values agree except E1 at the 1e-17
level (oracle rounding).

| Key | Value | Limit | |
|---|---|---|---|
| E1 `f_mix` vs `reference_scattering` | 2.5e-14 | ≤ 1e-9 | PASS |
| E2 `scat NI63` = Cu reproduces `f_mix` | 0 | ≤ 1e-6 | PASS: bit-identical, so the custom b is not stored in lower precision |
| E3 NI63 without `scat` is read as Ni | 1.5e-28 | ≤ 1e-6 | PASS |
| E4 VOID padding, 3D | 0 | ≤ 1e-9 | PASS: bit-identical |
| E5 metadata and shape | 0 | 0 | PASS |
| P0 r grid, finite G | 0 | 0 | PASS |
| P1 `scat` kept for the powder run | 3.0e-16 | ≤ 1e-6 | PASS |
| P2 powder G distinguishes the lengths | 3.2e-3 | ≥ 1e-4 | PASS (estimate was 3.5e-3) |
| **P3 VOID padding, powder PDF** | **0.25** | ≤ 1e-6 | **FAIL** |
| P4 neutron `read stru` = `read cell 4,4,4` | 4.0e-16 | ≤ 1e-6 | PASS |
| P5 X-ray `read stru` = `read cell 4,4,4` | 4.0e-16 | ≤ 1e-6 | PASS |

## What this establishes
- **`scat NAME, 0,0,0,0,0,0,0,0, b`** sets the neutron scattering length of an atom name in Fourier runs and
  is still in effect for a following `powder` run (E2, P1, P2). `FourierConfig(scattering_lengths=…)` and
  `PDFConfig(scattering_lengths=…)` can be offered.
- **DISCUS reads an unknown name such as NI63 as its element** (Ni) from the first two letters (E3).
- **VOID atoms do not scatter in single-crystal Fourier runs** (E4).
- **The powder PDF from `read stru` equals the verified `read cell` route** on the same crystal, neutron and
  X-ray (P4, P5). `PDFConfig` with a `.stru` can be offered.

## Finding: VOID atoms enter the powder PDF normalisation (P3)
`p_void` (the `mix` crystal plus one VOID per cell, 320 rows, 256 atoms) gives exactly **1.25 × G(r) of
`p_mix`** at every r: the ratio is 1.250000000 with a spread of 3.5e-16, and a fit G_void = a·G_mix + c·r gives
a = 1.25 and c = −1e-17. 1.25 = 320/256 = rows / scattering atoms. So DISCUS 6.23.01 counts the VOID rows in
the normalisation of the Debye-sum PDF while they add nothing to the pair sum. The effect is a pure scale factor:
the low-r slope −4πρ₀ is scaled by the same 1.25, so it is not a change of the density term alone.

Consequence (PREREG "What a failure would mean"): a supercell padded with VOID atoms must not be used as it
is for the 1D PDF. This concerns the disorder builder's models with vacancies or with molecules of different
sizes on one site, and it affected the builder's engine 1.0 as well. 3D results are not affected (E4).

Follow-up: v2 (`v2/`, registered separately) tests two ways out: the same crystal written as **one cell without
VOID rows**, and the **scale law** with another VOID fraction and with X-rays.
