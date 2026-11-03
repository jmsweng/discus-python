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

# Result of v2 (VOID atoms in the powder PDF)

**All 7 criteria pass.** Run `20261103T114220Z_LAPTOP-NGU3BT0N` (evidence `v2/evidence/…`, same binary). The kit
was registered at `ac0ff79` before the run. 9 runs, exit 0, no engine markers; re-analysed in the cloud with
identical values.

| Key | Value | Limit |
|---|---|---|
| V1 one cell without VOID rows = verified layout, neutron | 4.0e-16 | ≤ 1e-6 |
| V2 same, X-ray | 4.0e-16 | ≤ 1e-6 |
| V3 vacancies, neutron: (240/256)·padded = one cell | 7.3e-16 | ≤ 1e-6 |
| V4 vacancies, X-ray | 5.2e-16 | ≤ 1e-6 |
| V5 VOID padding, X-ray: (256/320)·padded = unpadded | 4.0e-16 | ≤ 1e-6 |
| V6 power: vacancies change G | 4.9e-2 | ≥ 1e-3 |
| V0 r grid, finite G | 0 | 0 |

- **Fix confirmed:** writing the real atoms as one cell of the whole crystal (`ncell 1,1,1,N`, no VOID rows)
  gives the same powder G(r) as the verified layout, for neutrons and X-rays, and with vacancies. The disorder
  builder uses it for its 1D PDF (`supercell_pdf.stru`).
- **Scale law:** G_padded = (rows / atoms) · G_true holds exactly (to 1e-15) for VOID fractions 1/5 and 1/16,
  neutron and X-ray. This explains v1 P3.
- **I1 (informational):** `p_mix_n` repeats v1's `p_mix` with byte-identical inputs and binary but is not
  byte-identical: the values differ in the last bits (max 3e-16 relative). The powder Debye sum of 6.23.01 is
  not bit-reproducible between runs; comparisons of powder runs need a tolerance, as all criteria here have.
