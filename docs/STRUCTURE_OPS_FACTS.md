# Structure Ops — Facts (A4-1)

Evidence: `tests/fixtures/structure_ops/evidence/P1.log`

## Facts

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| P1: NumPy 2.5.4; `default_rng(1234).choice([0..9], 4, replace=False) = [7, 3, 8, 6]`; `normal(0,1,(2,3)) = [[-1.604, 0.064, 0.741], [0.153, 0.864, 2.913]]`; `random(3) = [0.977, 0.380, 0.923]` | P1.log | `numpy.__version__ = 2.5.4` … `default_rng(1234).choice(np.arange(10), size=4, replace=False) = [7 3 8 6]` |
| P2: `with_frac` accepts positions outside [0,1) (no range validation); `wrapped()` wraps to [0,1); `distance_matrix` works; no overlap check failure for shifted atoms (min distance = 0.0) | P1.log | `with_frac accepted: n_atoms=8, frac range=[-0.3000, 1.9500]` … `wrapped() frac range=[0.0500, 0.9500]` … `distance_matrix shape=(8, 8), min=0.0000` |
| P3: 2×2×2 GaAs supercell: lattice (11.3074, 11.3074, 11.3074, 90, 90, 90); supercell=[2,2,2]; 64 atoms; cell_index starts [0,0,0]×8 then [0,0,1]×8 (iz fastest); frac range [0, 0.875]; atom_id 0–19 | P1.log | `lattice.parameters = (11.3074, 11.3074, 11.3074, 90.0, 90.0, 90.0)` … `cell_index[:20] = [[0, 0, 0], …, [0, 0, 1], …]` … `frac min = [0. 0. 0.]` … `frac max = [0.875 0.875 0.875]` |
| P4: 4096-atom (8×8×8) construction with default overlap check: 0.061 s; 64000-atom (20×20×20): 0.800 s | P1.log | `4096 atoms: 0.061 s` … `64000 atoms: 0.800 s` |
| P5: The `_gaas_supercell` helper creates rows NOT in DISCUS order (iz fastest vs DISCUS ix fastest); `format_stru` raises `DiscusFormatError`; `Structure.select(rows=…)` does NOT reorder rows (returns ascending); a manually DISCUS-ordered structure round-trips with `frac` identical but `site_label` lost (.stru does not store site labels); `format_cll` raises `ValueError: a .cll is a unit cell; use write_stru for supercells` | P1.log | `format_stru raised: DiscusFormatError: rows are not in DISCUS order; the writer never reorders (row 8 is in cell 0,0,1, expected cell 1)` … `Round trip equal: False` … `frac diff max: 0.00e+00` … `site_label equal: False` … `format_cll raised: ValueError: a .cll is a unit cell; use write_stru for supercells` |
| P5b: DISCUS-ordered GaAs 4×4×4 with all atoms displaced by Gaussian noise: σ=0.05 and σ=0.4 round-trip ok (max |Δfrac| = 1.25e-07); σ=1.0 and σ=2.0 raise `DiscusFormatError: atom order does not match ncell cell-major order`; single atom +3.0 Å raises `DiscusFormatError` | P1.log | `sigma=0.05: round trip ok, max \|df\| = 1.25e-07` … `sigma=0.4: round trip ok, max \|df\| = 1.25e-07` … `sigma=1.0: DiscusFormatError: <string>:0: atom order does not match ncell cell-major order (row 91, slot 3, axis 0): ''` … `sigma=2.0: DiscusFormatError: <string>:0: atom order does not match ncell cell-major order (row 1, slot 1, axis 0): ''` … `single atom +3.0 A: DiscusFormatError: <string>:0: atom order does not match ncell cell-major order (row 0, slot 0, axis 0): ''` |
| P6: `Structure.sites()` returns sites in **sorted** order (by site_label), NOT in order of first appearance. For the GaAs CIF fixture: sites are [As_1, As_2, As_3, As_4, Ga_1, Ga_2, Ga_3, Ga_4] (sorted), but first appearance is [Ga_1, …, As_4] | P1.log | `sites() order: [('As_1', (0, 0, 0)), …, ('Ga_4', (0, 0, 0))]` … `First appearance order: [('Ga_1', …), …, ('As_4', …)]` … `Same order: False` |
| P7: `Lattice.from_matrix(A * [2,3,4][:, None])` builds a scaled lattice; parameters (10, 15, 20, 90, 90, 90); angles preserved | P1.log | `Scaled lattice.parameters = (10.0, 15.0, 20.0, 90.0, 90.0, 90.0)` |

## Conflicts

none (see Amendments in 00_START_HERE.md)

## Overseer rulings (verbatim)

- `sample_configuration` site order: `Structure.sites()` groups, ordered by `min(members)` (probe P6: `sites()` is sorted by label, not by first appearance).
- `repeat` lattice: `Lattice.from_matrix(matrix * reps[:, None])`.
- Facts used (probes): `with_frac`/`Structure` accept positions outside [0,1) (P2); the `.stru` writer checks row order from `cell_index` only and never reorders (P5), the existing test helper `_gaas_supercell` is NOT in DISCUS order (iz fastest), so in tests build DISCUS-ordered inputs with `repeat` itself; `format_cll` refuses `supercell != (1,1,1)`; the `.stru` reader returns `DiscusStru` (use `.structure`), assigns `site_label = site{slot+1}` and `atom_id = 0..N-1`, so round-trip comparisons ignore `site_label`/`atom_id`. The reader's order check tolerates an atom offset of up to 0.5 unit cell from the per-slot median on each axis (source: `discus_io.py` lines ~800).
