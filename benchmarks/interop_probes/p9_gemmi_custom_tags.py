"""P9: gemmi CIF writing with _discuspdf_* tags and re-reading."""
import io
import os

def main():
    print("=== P9: gemmi CIF with _discuspdf_* tags ===")

    text = '''data_discuspdf_test
_discuspdf_model    average
_cell_length_a      5.6537
_cell_length_b      5.6537
_cell_length_c      5.6537
_cell_angle_alpha   90
_cell_angle_beta    90
_cell_angle_gamma   90
_symmetry_space_group_name_H-M   'P 1'

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_discuspdf_site_label
_discuspdf_atom_id
Ga1_1   0.0   0.0   0.0   1.0   Ga1_1   0
As1_1   0.25  0.25  0.25  1.0   As1_1   1
'''

    # Re-read with gemmi
    import gemmi
    doc = gemmi.cif.read_string(text)
    block = doc.sole_block()
    print("--- re-read with gemmi ---")
    print(f"model: {block.find_value('_discuspdf_model')}")
    labels = block.find_loop("_atom_site_label")
    print(f"labels: {[labels[i] for i in range(len(labels))]}")
    dp_labels = block.find_loop("_discuspdf_site_label")
    print(f"dp_site_labels: {[dp_labels[i] for i in range(len(dp_labels))]}")
    dp_ids = block.find_loop("_discuspdf_atom_id")
    print(f"dp_atom_ids: {[dp_ids[i] for i in range(len(dp_ids))]}")
    sg = block.find_value("_symmetry_space_group_name_H-M")
    print(f"sg raw: {sg!r}")
    print(f"sg as_string: {gemmi.cif.as_string(sg)!r}")

    # Re-read with ASE
    from ase.io import read
    atoms = read(io.StringIO(text), format="cif")
    print("--- re-read with ASE ---")
    print(f"n_atoms: {len(atoms)}")
    print(f"species: {atoms.get_chemical_symbols()}")
    print(f"info: {atoms.info}")
    print(f"arrays keys: {list(atoms.arrays.keys())}")
    print(f"scaled positions:\n{atoms.get_scaled_positions()}")

    print("=== P9 done ===")

if __name__ == "__main__":
    main()
