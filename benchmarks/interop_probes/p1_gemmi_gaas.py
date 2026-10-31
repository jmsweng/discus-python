"""P1: gemmi on the real COD GaAs CIF."""
import gemmi
import sys
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print(f"=== P1: gemmi on gaas_cod_9008845.cif ===")
    doc = gemmi.cif.read(CIF)
    block = doc.sole_block()
    print(f"block name: {block.name}")

    sg = block.find_value("_symmetry_space_group_name_H-M")
    print(f"spacegroup (raw): {sg!r}")
    if sg:
        sg_clean = sg.strip("'")
        print(f"spacegroup (clean): {sg_clean!r}")
        sg_obj = gemmi.find_spacegroup_by_name(sg_clean)
        if sg_obj:
            print(f"  number: {sg_obj.number}")
            print(f"  crystal_system: {sg_obj.crystal_system_str()}")
            print(f"  laue: {sg_obj.laue_str()}")
            ops = sg_obj.operations()
            print(f"  n_cen_ops: {len(ops.cen_ops)}")
            print(f"  n_sym_ops: {len(ops.sym_ops)}")
        else:
            print("  gemmi.find_spacegroup_by_name: NOT FOUND")

    a = float(block.find_value("_cell_length_a"))
    b = float(block.find_value("_cell_length_b"))
    c = float(block.find_value("_cell_length_c"))
    al = float(block.find_value("_cell_angle_alpha"))
    be = float(block.find_value("_cell_angle_beta"))
    ga = float(block.find_value("_cell_angle_gamma"))
    print(f"cell: a={a} b={b} c={c} alpha={al} beta={be} gamma={ga}")

    symops = block.find_loop("_space_group_symop_operation_xyz")
    print(f"n_symop: {len(symops)}")
    print(f"first 3 symops: {[symops[i] for i in range(min(3, len(symops)))]}")

    labels = block.find_loop("_atom_site_label")
    xs = block.find_loop("_atom_site_fract_x")
    ys = block.find_loop("_atom_site_fract_y")
    zs = block.find_loop("_atom_site_fract_z")
    print(f"n_atoms: {len(labels)}")
    for i in range(len(labels)):
        print(f"  {labels[i]}: x={xs[i]} y={ys[i]} z={zs[i]}")

    # Check optional columns
    occ = block.find_loop("_atom_site_occupancy")
    print(f"_atom_site_occupancy present: {len(occ) > 0} ({len(occ)} rows)")
    uiso = block.find_loop("_atom_site_U_iso_or_equiv")
    print(f"_atom_site_U_iso_or_equiv present: {len(uiso) > 0} ({len(uiso)} rows)")
    biso = block.find_loop("_atom_site_B_iso_or_equiv")
    print(f"_atom_site_B_iso_or_equiv present: {len(biso) > 0} ({len(biso)} rows)")
    adp = block.find_loop("_atom_site_adp_type")
    print(f"_atom_site_adp_type present: {len(adp) > 0} ({len(adp)} rows)")
    aniso = block.find_loop("_atom_site_aniso_U_11")
    print(f"_atom_site_aniso_U_11 present: {len(aniso) > 0} ({len(aniso)} rows)")
    dis = block.find_loop("_atom_site_disorder_group")
    print(f"_atom_site_disorder_group present: {len(dis) > 0} ({len(dis)} rows)")
    tsym = block.find_loop("_atom_site_type_symbol")
    print(f"_atom_site_type_symbol present: {len(tsym) > 0} ({len(tsym)} rows)")
    mod = block.find_value("_cell_modulation_dimension")
    print(f"_cell_modulation_dimension present: {mod is not None}")

    jr = block.find_value("_citation_journal_ref")
    print(f"_citation_journal_ref: {jr!r}")

    print("=== P1 done ===")

if __name__ == "__main__":
    main()
