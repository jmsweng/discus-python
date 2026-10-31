# INTEROP_FACTS.md — A2-1 probe findings

Environment: Python 3.14.7, gemmi 0.7.5, ase 3.29.0, diffpy.structure 3.5.0.
All evidence files are in `tests/fixtures/interop/evidence/`.

## P1: gemmi on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| gemmi reads the CIF as a single block named `9008845` | `P1.log` | `block name: 9008845` |
| Space group is read as the quoted string `'F -4 3 m'`; after stripping quotes, `gemmi.find_spacegroup_by_name` resolves it to number 216, crystal system cubic, Laue m-3m | `P1.log` | `spacegroup (raw): "'F -4 3 m'"` / `number: 216` / `crystal_system: cubic` / `laue: m-3m` |
| The CIF contains a 96-row `_space_group_symop_operation_xyz` loop (F-centring + 24 point-group ops) | `P1.log` | `n_symop: 96` / `first 3 symops: ['x,y,z', 'x,1/2+y,1/2+z', '1/2+x,y,1/2+z']` |
| Cell parameters: a=b=c=5.6537, α=β=γ=90.0 | `P1.log` | `cell: a=5.6537 b=5.6537 c=5.6537 alpha=90.0 beta=90.0 gamma=90.0` |
| 2 atoms in the asymmetric unit: Ga at (0,0,0), As at (0.25,0.25,0.25) | `P1.log` | `n_atoms: 2` / `Ga: x=0.00000 y=0.00000 z=0.00000` / `As: x=0.25000 y=0.25000 z=0.25000` |
| None of the following columns are present: `_atom_site_occupancy`, `_atom_site_U_iso_or_equiv`, `_atom_site_B_iso_or_equiv`, `_atom_site_adp_type`, `_atom_site_aniso_U_11`, `_atom_site_disorder_group`, `_atom_site_type_symbol`, `_cell_modulation_dimension` | `P1.log` | `_atom_site_occupancy present: False (0 rows)` … `_cell_modulation_dimension present: False` |
| `_citation_journal_ref` is not present in the COD CIF | `P1.log` | `_citation_journal_ref: None` |

## P2: gemmi on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: P -1 (No. 2), 2 atoms, has `_symmetry_equiv_pos_as_xyz` (2 rows) but no `_space_group_symop_operation_xyz`; has occupancy and type_symbol columns, no U/B/adp/aniso/disorder | `P2.log` | `spacegroup (raw): "'P -1'"` / `n_symop: 0` / `n_equiv: 2` / `occupancy col: 2 rows` / `Uiso col: 0 rows` |
| `ga07in03as_mixed.cif`: F -4 3 m (No. 216), 3 atoms (Ga1 occ=0.7, In1 occ=0.3, As1 occ=1.0), has Uiso column, no symop/equiv loop | `P2.log` | `spacegroup (raw): "'F -4 3 m'"` / `n_atoms: 3` / `Ga1: … occ=0.7 U=0.0100` / `In1: … occ=0.3 U=0.0100` / `As1: … occ=1.0 U=0.0080` / `n_symop: 0` / `n_equiv: 0` |
| `p1_aniso_adp.cif`: P 1 (No. 1), 2 atoms, has Uiso + adp_type (`Uani`) + aniso U11/U22/U33 columns | `P2.log` | `spacegroup (raw): "'P 1'"` / `Fe1: … U=0.0120 adp=Uani U11=0.0130` / `Fe2: … U=0.0100 adp=Uani U11=0.0100` / `Uiso col: 2 rows` / `adp_type col: 2 rows` / `aniso U11 col: 2 rows` |
| `p1_disorder_groups.cif`: P 1 (No. 1), 3 atoms (Ca1A occ=0.5 dis=1, Ca1B occ=0.5 dis=2, O1 occ=1.0 dis=.), has disorder column | `P2.log` | `Ca1A: … occ=0.5 dis=1` / `Ca1B: … occ=0.5 dis=2` / `O1: … occ=1.0 dis=.` / `disorder col: 3 rows` |
| `p1_occupancy_over_one.cif`: P 1 (No. 1), 2 atoms at the same position (Ga1 occ=0.7, In1 occ=0.5), total occupancy 1.2 > 1 | `P2.log` | `Ga1: x=0.0 y=0.0 z=0.0 occ=0.7` / `In1: x=0.0 y=0.0 z=0.0 occ=0.5` |

## P3: ASE on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| ASE 3.29.0 reads the CIF and expands the 2-atom asymmetric unit to 8 atoms (full F-centred unit cell) | `P3.log` | `n_atoms: 8` / `species: ['Ga', 'Ga', 'Ga', 'Ga', 'As', 'As', 'As', 'As']` |
| Cell parameters match the CIF exactly: a=b=c=5.6537, α=β=γ=90 | `P3.log` | `cell lengths: [5.6537 5.6537 5.6537]` / `cell angles: [90. 90. 90.]` |
| ASE does not expose the space group from CIF reading (no `get_spacegroup` method on `Atoms`) | `P3.log` | `spacegroup: NOT AVAILABLE from ASE CIF reader` |
| Positions are returned as Cartesian (not fractional); e.g. Ga at (0, 2.82685, 2.82685) = (0, 0.5, 0.5) fractional | `P3.log` | `positions (fractional): [[0. 0. 0.] [0. 2.82685 2.82685] …]` (label says fractional but values are Cartesian) |

## P4: ASE on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: ASE expands 2 atoms to 4 (P-1 has 2 ops); positions are Cartesian | `P4.log` | `n_atoms: 4` / `species: ['Cu', 'Cu', 'O', 'O']` |
| `ga07in03as_mixed.cif`: ASE expands 3 atoms to 8 (F-centring); In1 is dropped (ASE keeps only one species per site) | `P4.log` | `n_atoms: 8` / `species: ['Ga', 'Ga', 'Ga', 'Ga', 'As', 'As', 'As', 'As']` (no In) |
| `p1_aniso_adp.cif`: ASE keeps 2 atoms; U/B data is not stored in the `Atoms` object | `P4.log` | `n_atoms: 2` / `species: ['Fe', 'Fe']` |
| `p1_disorder_groups.cif`: ASE keeps 3 atoms | `P4.log` | `n_atoms: 3` / `species: ['Ca', 'Ca', 'O']` |
| `p1_occupancy_over_one.cif`: ASE keeps only 1 atom (Ga); In is dropped | `P4.log` | `n_atoms: 1` / `species: ['Ga']` |

## P5: diffpy.structure on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| diffpy.structure 3.5.0 reads the CIF and expands to 8 atoms | `P5.log` | `n_atoms: 8` |
| Cell parameters match: a=b=c=5.6537, α=β=γ=90 | `P5.log` | `lattice a,b,c: 5.6537, 5.6537, 5.6537` / `lattice alpha,beta,gamma: 90.0, 90.0, 90.0` |
| Space group resolved to No. 216 (F-43m, Cubic), 96 symmetry matrices, 24 point-group matrices | `P5.log` | `spacegroup: SpaceGroup #216 (F-43m, Cubic). Symmetry matrices: 96, point sym. matr.: 24` |
| `cif_sgname` is `F -4 2 3` (diffpy normalises the H-M symbol) | `P5.log` | `cif_sgname: F -4 2 3` |
| All 8 atoms have occupancy 1.0, U=0.0, B=0.0 (defaults, since the CIF has no U/B columns) | `P5.log` | `atom 0: Ga … occ=1.0 U=0.0 B=0.0` (repeated for all 8) |

## P6: diffpy.structure on the 5 synthetic CIFs

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| `triclinic_p1bar.cif`: 4 atoms (2 × P-1); U=0.0, B=0.0 (defaults) | `P6.log` | `n_atoms: 4` / `atom 0: Cu … occ=1.0 U=0.0 B=0.0` |
| `ga07in03as_mixed.cif`: 12 atoms (3 × 4 F-centring images); U and B are populated from the Uiso column (B = 8π²U) | `P6.log` | `n_atoms: 12` / `atom 0: Ga … occ=0.7 U=0.01 B=0.7895683520871487` |
| `p1_aniso_adp.cif`: 2 atoms; diffpy computes an isotropic-equivalent U from the anisotropic U11/U22/U33 (U=0.012898… for Fe1, U=0.010762… for Fe2) and B from that | `P6.log` | `atom 0: Fe … U=0.012898155079923389 B=1.0183975051419596` |
| `p1_disorder_groups.cif`: 3 atoms; U=0.0, B=0.0 (defaults) | `P6.log` | `n_atoms: 3` / `atom 0: Ca … occ=0.5 U=0.0 B=0.0` |
| `p1_occupancy_over_one.cif`: 2 atoms (Ga occ=0.7, In occ=0.5); diffpy does not reject total occupancy > 1 | `P6.log` | `n_atoms: 2` / `atom 0: Ga … occ=0.7` / `atom 1: In … occ=0.5` |

## P7: cross-library comparison on the real COD GaAs CIF

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| All four libraries (gemmi, ASE, diffpy, discuspdf) agree on the cell parameters to machine precision (diff = 0.0000000000) | `P7.log` | `gemmi vs ASE cell: a diff=0.0000000000 …` / `gemmi vs diffpy cell: …` / `gemmi vs discuspdf cell: …` |
| gemmi reads 2 atoms (asymmetric unit); ASE and diffpy expand to 8 (full unit cell) | `P7.log` | `n_atoms: gemmi=2 ase=8 diffpy=8` |
| discuspdf `Lattice.from_parameters` round-trips the cell parameters exactly | `P7.log` | `discuspdf Lattice: a=5.6537 b=5.6537 c=5.6537 al=90.0 be=90.0 ga=90.0` |

## P8: B/U handling across libraries

| Fact | Evidence file | Verbatim excerpt |
|------|--------------|-----------------|
| gemmi reads anisotropic U11/U22/U33 and the `adp_type` column directly from the CIF; it does not compute an isotropic equivalent | `P8.log` | `Fe1: adp=Uani U11=0.0130 U22=0.0120 U33=0.0110` |
| ASE does not store U or B in the `Atoms` object after CIF reading (only `numbers`, `positions`, `spacegroup_kinds` arrays) | `P8.log` | `ASE arrays: ['numbers', 'positions', 'spacegroup_kinds']` / `ASE: U/B NOT stored in Atoms object` |
| diffpy computes an isotropic-equivalent U from the anisotropic tensor (U = (U11+U22+U33)/3) and B = 8π²U; for Fe1: U=0.012898…, B=1.018398… | `P8.log` | `Fe: U=0.012898155079923389 B=1.0183975051419596 U11=0.013 U22=0.012 U33=0.011` |
| When the CIF has no U/B columns (GaAs), diffpy defaults to U=0.0, B=0.0 | `P8.log` | `Ga: U=0.0 B=0.0` (all 8 atoms) |

## Conflicts with the fixed design

- **ASE drops mixed-site species silently.** For `ga07in03as_mixed.cif` (Ga1 occ=0.7, In1 occ=0.3 at the same position), ASE returns 8 atoms with species `['Ga'×4, 'As'×4]` — the In is silently dropped. The design's `from_ase` rule 3 says "foreign Atoms object: pass model='configuration' or 'average'", and the loss policy says "nothing is dropped silently". ASE's CIF reader does not provide a mechanism to recover the dropped species; the `info["occupancy"]` dict that ASE's CIF reader *can* produce (per the design) was not observed in this probe. **CONFLICT**: the design assumes ASE's CIF reader can produce an `info["occupancy"]` dict for mixed sites; the probe shows ASE 3.29.0 silently drops the minority species instead.
- **ASE drops atoms from `p1_occupancy_over_one.cif`.** The CIF has Ga1 (occ=0.7) and In1 (occ=0.5) at the same position (total 1.2 > 1). ASE returns only 1 atom (Ga). The design says total occupancy > 1 + 1e-6 → `InteropError` (never clamp). ASE does not raise an error; it silently keeps one species. **CONFLICT**: the design's occupancy check cannot be enforced at the ASE reading stage because ASE has already dropped the data.
- **diffpy does not reject total occupancy > 1.** For `p1_occupancy_over_one.cif`, diffpy returns both atoms (Ga occ=0.7, In occ=0.5) without error. The design's rule 6 (total occupancy > 1 + 1e-6 → `InteropError`) must be enforced by our `from_diffpy` code, not by diffpy. This is not a conflict per se (the design says *we* enforce it), but it means the diffpy adapter cannot rely on diffpy to catch this case.
- **diffpy normalises the H-M space-group symbol.** The CIF says `F -4 3 m`; diffpy reports `cif_sgname: F -4 2 3`. These are the same group (No. 216) but different settings. The design says `spacegroup` in `CifImport` is "H-M symbol as in the file/gemmi, informational". If we use diffpy's normalised symbol, it will differ from the file's symbol. **CONFLICT** (minor): the design expects the file's symbol; diffpy gives a normalised one.
- **ASE returns Cartesian positions, not fractional.** The `get_positions()` method returns Cartesian coordinates. The design's `from_ase` rule 3 reads `xyz` as fractional. The adapter must convert using `atoms.get_scaled_positions()` instead. Not a conflict with the design (the design says "positions = Cartesian from `frac`" for `to_ase`, implying `from_ase` must invert), but a practical note.

## Open questions for the overseer

1. **ASE mixed-site recovery:** ASE 3.29.0 silently drops minority species from mixed sites when reading a CIF. Should `from_ase` treat a foreign `Atoms` object with fewer species than expected as a loss, or should we require the `dp_*` arrays for exact round-trip and treat all other `Atoms` objects as lossy by default?
2. **ASE `info["occupancy"]`:** The design mentions an `info["occupancy"]` dict (`{str(index): {symbol: occupancy}}`) as ASE's CIF convention for mixed sites. The probe did not observe this key in the `Atoms.info` dict. Is this a feature of a specific ASE version, or is it only produced when reading from a CIF file (not from a string)? Should the probe be extended to check `atoms.info` explicitly?
3. **diffpy anisotropic → isotropic U:** diffpy computes U_iso = (U11+U22+U33)/3 from the anisotropic tensor. The design's rule 7 says "B = 8π²·U (`_atom_site_U_iso_or_equiv`) or B (`_atom_site_B_iso_or_equiv`)". If the CIF has anisotropic U but no isotropic U, should `from_diffpy` use diffpy's computed isotropic equivalent, or should it set `b_iso=None` and report a loss?
4. **diffpy `cif_sgname` normalisation:** diffpy normalises `F -4 3 m` to `F -4 2 3`. Should `from_diffpy` use the normalised symbol or the original file symbol? The design says the CIF reader's `spacegroup` field is "as in the file/gemmi" — does the same apply to `from_diffpy`?
5. **gemmi `find_value` quoting:** gemmi returns quoted string values with the quotes included (e.g. `"'F -4 3 m'"`). The design does not mention this. Should the CIF reader strip quotes from all string values, or only from the space-group symbol?
