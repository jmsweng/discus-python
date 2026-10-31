# Interop test fixtures

CIF files used by the interop adapter tests.

## ga07in03as_mixed.cif

- Provenance: authored for this repository (public domain, synthetic test fixture).
- Content: mixed-occupancy Ga/In/As site (Ga 0.7 / In 0.3) in `F -4 3 m`, with isotropic U.

## triclinic_p1bar.cif

- Provenance: authored for this repository (public domain, synthetic test fixture).
- Content: triclinic `P -1` cell with explicit `_symmetry_equiv_pos_as_xyz` list.

## p1_aniso_adp.cif

- Provenance: authored for this repository (public domain, synthetic test fixture).
- Content: `P 1` cell with anisotropic ADPs (`Uani`) in a separate
  `_atom_site_aniso_U_*` loop.

## p1_disorder_groups.cif

- Provenance: authored for this repository (public domain, synthetic test fixture).
- Content: `P 1` cell with `_atom_site_disorder_group` labels (Ca1A/Ca1B split
  occupancy, O1 with `.`).

## p1_occupancy_over_one.cif

- Provenance: authored for this repository (public domain, synthetic test fixture).
- Content: `P 1` cell whose summed site occupancies exceed 1.0 (0.7 + 0.5) —
  an intentionally invalid fixture for validation tests.

## gaas_cod_9008845.cif

- Provenance: downloaded from the Cambridge Crystallographic Data Centre (COD),
  COD id 9008845.
- Journal reference: Wyckoff, R. W. G., Crystal Structures, 1963, vol. 1, pp. 85-237.
- DOI: none present.
- Licence: public domain (CC0, per the COD site).
- Download date: 2026-10-31.
- SHA-256: 2B011A3F401B638163CF663211726C2B72BF5C8A3BA32511DF0D867382FBA943
