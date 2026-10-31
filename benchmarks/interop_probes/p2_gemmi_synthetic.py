"""P2: gemmi on the 5 synthetic CIFs."""
import gemmi
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def inspect(fname):
    cif = os.path.join(FIXTURES, fname)
    doc = gemmi.cif.read(cif)
    block = doc.sole_block()
    print(f"\n--- {fname} ---")
    print(f"block name: {block.name}")

    sg = block.find_value("_symmetry_space_group_name_H-M")
    print(f"spacegroup (raw): {sg!r}")
    if sg:
        sg_clean = sg.strip("'")
        sg_obj = gemmi.find_spacegroup_by_name(sg_clean)
        if sg_obj:
            print(f"  number: {sg_obj.number}, crystal_system: {sg_obj.crystal_system_str()}")
            ops = sg_obj.operations()
            print(f"  n_cen_ops: {len(ops.cen_ops)}, n_sym_ops: {len(ops.sym_ops)}")
        else:
            print(f"  gemmi.find_spacegroup_by_name({sg_clean!r}): NOT FOUND")

    a = block.find_value("_cell_length_a")
    b = block.find_value("_cell_length_b")
    c = block.find_value("_cell_length_c")
    al = block.find_value("_cell_angle_alpha")
    be = block.find_value("_cell_angle_beta")
    ga = block.find_value("_cell_angle_gamma")
    print(f"cell: a={a} b={b} c={c} alpha={al} beta={be} gamma={ga}")

    symops = block.find_loop("_space_group_symop_operation_xyz")
    equiv = block.find_loop("_symmetry_equiv_pos_as_xyz")
    print(f"n_symop: {len(symops)}")
    print(f"n_equiv: {len(equiv)}")
    if len(equiv) > 0:
        print(f"  equiv[0]: {equiv[0]!r}")

    labels = block.find_loop("_atom_site_label")
    xs = block.find_loop("_atom_site_fract_x")
    ys = block.find_loop("_atom_site_fract_y")
    zs = block.find_loop("_atom_site_fract_z")
    occ = block.find_loop("_atom_site_occupancy")
    uiso = block.find_loop("_atom_site_U_iso_or_equiv")
    biso = block.find_loop("_atom_site_B_iso_or_equiv")
    adp = block.find_loop("_atom_site_adp_type")
    aniso = block.find_loop("_atom_site_aniso_U_11")
    dis = block.find_loop("_atom_site_disorder_group")
    tsym = block.find_loop("_atom_site_type_symbol")
    mod = block.find_value("_cell_modulation_dimension")

    print(f"n_atoms: {len(labels)}")
    for i in range(len(labels)):
        line = f"  {labels[i]}: x={xs[i]} y={ys[i]} z={zs[i]}"
        if len(occ) > 0:
            line += f" occ={occ[i]}"
        if len(uiso) > 0:
            line += f" U={uiso[i]}"
        if len(biso) > 0:
            line += f" B={biso[i]}"
        if len(adp) > 0:
            line += f" adp={adp[i]}"
        if len(aniso) > 0:
            line += f" U11={aniso[i]}"
        if len(dis) > 0:
            line += f" dis={dis[i]}"
        if len(tsym) > 0:
            line += f" sym={tsym[i]}"
        print(line)
    print(f"  occupancy col: {len(occ)} rows")
    print(f"  Uiso col: {len(uiso)} rows")
    print(f"  Biso col: {len(biso)} rows")
    print(f"  adp_type col: {len(adp)} rows")
    print(f"  aniso U11 col: {len(aniso)} rows")
    print(f"  disorder col: {len(dis)} rows")
    print(f"  type_symbol col: {len(tsym)} rows")
    print(f"  modulation: {mod!r}")

def main():
    print("=== P2: gemmi on 5 synthetic CIFs ===")
    for f in ["triclinic_p1bar.cif", "ga07in03as_mixed.cif", "p1_aniso_adp.cif",
              "p1_disorder_groups.cif", "p1_occupancy_over_one.cif"]:
        inspect(f)
    print("=== P2 done ===")

if __name__ == "__main__":
    main()
