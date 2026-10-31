"""P8: B/U handling across libraries."""
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def main():
    print("=== P8: B/U handling ===")

    # P1 aniso file
    cif = os.path.join(FIXTURES, "p1_aniso_adp.cif")
    print(f"\n--- p1_aniso_adp.cif ---")

    # gemmi
    import gemmi
    doc = gemmi.cif.read(cif)
    block = doc.sole_block()
    labels = block.find_loop("_atom_site_label")
    adp = block.find_loop("_atom_site_adp_type")
    u11 = block.find_loop("_atom_site_aniso_U_11")
    u22 = block.find_loop("_atom_site_aniso_U_22")
    u33 = block.find_loop("_atom_site_aniso_U_33")
    uiso = block.find_loop("_atom_site_U_iso_or_equiv")
    biso = block.find_loop("_atom_site_B_iso_or_equiv")
    print(f"gemmi: n_atoms={len(labels)}")
    for i in range(len(labels)):
        u11_v = u11[i] if len(u11) > i else "N/A"
        u22_v = u22[i] if len(u22) > i else "N/A"
        u33_v = u33[i] if len(u33) > i else "N/A"
        adp_v = adp[i] if len(adp) > i else "N/A"
        print(f"  {labels[i]}: adp={adp_v} U11={u11_v} U22={u22_v} U33={u33_v}")
    print(f"  Uiso col: {len(uiso)} rows")
    print(f"  Biso col: {len(biso)} rows")

    # ASE
    from ase.io import read
    atoms = read(cif)
    print(f"ASE: n_atoms={len(atoms)}")
    print(f"  ASE arrays: {list(atoms.arrays.keys())}")
    # ASE does not store U/B in the Atoms object from CIF reading
    print(f"  ASE: U/B NOT stored in Atoms object (only numbers, positions, spacegroup_kinds)")

    # diffpy
    from diffpy.structure import Structure
    s = Structure().read(cif)
    st = s.stru
    print(f"diffpy: n_atoms={len(st)}")
    for i, atom in enumerate(st):
        print(f"  {atom.element}: U={atom.Uisoequiv} B={atom.Bisoequiv} U11={atom.U11} U22={atom.U22} U33={atom.U33}")

    # GaAs (no B/U columns)
    cif2 = os.path.join(FIXTURES, "gaas_cod_9008845.cif")
    print(f"\n--- gaas_cod_9008845.cif (no B/U columns) ---")
    doc2 = gemmi.cif.read(cif2)
    block2 = doc2.sole_block()
    uiso2 = block2.find_loop("_atom_site_U_iso_or_equiv")
    biso2 = block2.find_loop("_atom_site_B_iso_or_equiv")
    adp2 = block2.find_loop("_atom_site_adp_type")
    print(f"gemmi: Uiso={len(uiso2)} rows, Biso={len(biso2)} rows, adp={len(adp2)} rows")

    atoms2 = read(cif2)
    print(f"ASE: n_atoms={len(atoms2)}")

    s2 = Structure().read(cif2)
    st2 = s2.stru
    print(f"diffpy: n_atoms={len(st2)}")
    for i, atom in enumerate(st2):
        print(f"  {atom.element}: U={atom.Uisoequiv} B={atom.Bisoequiv}")

    print("=== P8 done ===")

if __name__ == "__main__":
    main()
