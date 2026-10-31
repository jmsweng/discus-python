# DISCUS 6.23.01 fixtures

Small `.cll` / `.stru` files and the probe logs that support the verified facts in
`docs/DISCUS_FORMATS.md`. All files were produced by, or verified against, the real
DISCUS 6.23.01 (`/usr/local/bin/discus_suite`, WSL `Ubuntu-22.04`).

DISCUS version: **6.23.01** (banner `D I S C U S - S U I T E  Version 6.23.01`).
Date: 2026-10-30.

## Fixtures

| File | How produced |
|---|---|
| `gaas.cll` | Byte-for-byte copy of the tracked `examples/gaas.cll` (2-atom F-43m GaAs unit cell, 4-column atom lines, no `ncell`/`format`). |
| `gaas_111.stru` | Written by DISCUS: `read cell gaas.cll,1,1,1` then `save gaas_111.stru` (probe Q3, macro `q3.mac`). 8 atoms, centred 0…1 coords, no `ncell`. |
| `gaas_222.stru` | Written by DISCUS: `read cell gaas.cll,2,2,2` then `save gaas_222.stru` (probe Q3, macro `q3.mac`). 64 atoms, centred −1…+0.75 coords, no `ncell`. |
| `p1_2x2x2_a4.stru` | DISCUS-style file with `ncell 2,2,2,1`, P1, a = 4 Å, 8 atoms of two species (Cu/Ag), B = 0. Built from a DISCUS-written 2×2×2 P1 cell (macro `p1_222.mac`) and completed to the spec; **verified to read back** in DISCUS (`p1_verify.mac`: n[1]=8, n[2]=2, n[3]=1, n[4]=0). |

## Evidence logs

`evidence/` holds the DISCUS console logs for the probes (Q1–Q13, `00_env.txt`,
`smoke.log`). Each log is the raw `discus_suite -noautorun -macro <name>.mac`
output. The macros/inputs that produced them are in
`benchmarks/discus_6_23_01/a3_format_probes/`.

- `00_env.txt` — environment summary (Python, WSL distro, DISCUS banner, 87-passed).
- `smoke.log` — smoke test: `read cell gaas.cll,2,2,2` → n[1]=64.
- `Q1.log` … `Q13.log`, `Q7_v1.log` … `Q7_v9.log` — probe logs (see `docs/DISCUS_FORMATS.md`).

## Integrity

`MANIFEST.sha256` lists the SHA-256 of every file in this directory (format
`<hash>  <name>`). Recompute with `Get-FileHash -Algorithm SHA256` and compare.
