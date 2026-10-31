"""P7: cross-library comparison on the real COD GaAs CIF."""
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print("=== P7: cross-library comparison on gaas_cod_9008845.cif ===")

    # gemmi
    import gemmi
    doc = gemmi.cif.read(CIF)
    block = doc.sole_block()
    g_a = float(block.find_value("_cell_length_a"))
    g_b = float(block.find_value("_cell_length_b"))
    g_c = float(block.find_value("_cell_length_c"))
    g_al = float(block.find_value("_cell_angle_alpha"))
    g_be = float(block.find_value("_cell_angle_beta"))
    g_ga = float(block.find_value("_cell_angle_gamma"))
    g_n = len(block.find_loop("_atom_site_label"))
    print(f"gemmi: a={g_a} b={g_b} c={g_c} al={g_al} be={g_be} ga={g_ga} n={g_n}")

    # ASE
    from ase.io import read
    atoms = read(CIF)
    a_l = atoms.cell.lengths()
    a_ang = atoms.cell.angles()
    a_n = len(atoms)
    print(f"ASE:   a={a_l[0]} b={a_l[1]} c={a_l[2]} al={a_ang[0]} be={a_ang[1]} ga={a_ang[2]} n={a_n}")

    # diffpy
    from diffpy.structure import Structure
    s = Structure().read(CIF)
    st = s.stru
    d_a = st.lattice.a
    d_b = st.lattice.b
    d_c = st.lattice.c
    d_al = st.lattice.alpha
    d_be = st.lattice.beta
    d_ga = st.lattice.gamma
    d_n = len(st)
    print(f"diffpy: a={d_a} b={d_b} c={d_c} al={d_al} be={d_be} ga={d_ga} n={d_n}")

    # discuspdf
    from discuspdf.lattice import Lattice
    lat = Lattice.from_parameters(g_a, g_b, g_c, g_al, g_be, g_ga)
    p = lat.parameters
    print(f"discuspdf Lattice: a={p[0]} b={p[1]} c={p[2]} al={p[3]} be={p[4]} ga={p[5]}")

    # Compare
    print("\n--- agreement check ---")
    print(f"gemmi vs ASE cell: a diff={abs(g_a - a_l[0]):.10f} b diff={abs(g_b - a_l[1]):.10f} c diff={abs(g_c - a_l[2]):.10f}")
    print(f"gemmi vs diffpy cell: a diff={abs(g_a - d_a):.10f} b diff={abs(g_b - d_b):.10f} c diff={abs(g_c - d_c):.10f}")
    print(f"gemmi vs discuspdf cell: a diff={abs(g_a - p[0]):.10f} b diff={abs(g_b - p[1]):.10f} c diff={abs(g_c - p[2]):.10f}")
    print(f"n_atoms: gemmi={g_n} ase={a_n} diffpy={d_n}")
    print(f"NOTE: gemmi reads 2 atoms (asymmetric unit), ASE and diffpy expand to 8 (full unit cell)")

    print("=== P7 done ===")

if __name__ == "__main__":
    main()
