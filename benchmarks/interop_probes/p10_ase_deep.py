"""P10: ASE deep inspection (info, arrays, cell matrix, copy, pickle, degenerate cells)."""
import os
import numpy as np

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def inspect_ase(fname):
    cif = os.path.join(FIXTURES, fname)
    print(f"\n--- {fname} ---")
    from ase.io import read
    atoms = read(cif)
    print(f"n_atoms: {len(atoms)}")
    print(f"species: {atoms.get_chemical_symbols()}")
    print(f"info: {atoms.info}")
    print(f"arrays keys: {list(atoms.arrays.keys())}")
    print(f"cell.array:\n{atoms.cell.array}")
    print(f"pbc: {atoms.pbc}")
    print(f"scaled positions:\n{atoms.get_scaled_positions()}")

    # Compare cell matrix with discuspdf Lattice
    from discuspdf.lattice import Lattice
    a, b, c = atoms.cell.lengths()
    al, be, ga = atoms.cell.angles()
    lat = Lattice.from_parameters(a, b, c, al, be, ga)
    our_matrix = lat.matrix
    ase_matrix = atoms.cell.array
    max_diff = np.max(np.abs(ase_matrix - our_matrix))
    print(f"max abs diff (ASE cell.array vs Lattice.matrix): {max_diff:.2e}")

def test_copy_pickle():
    print("\n--- copy and pickle round trip ---")
    from ase.io import read
    import pickle
    atoms = read(os.path.join(FIXTURES, "ga07in03as_mixed.cif"))
    atoms.set_array("dp_test", np.array([1, 2, 3, 4, 5, 6, 7, 8]))
    atoms.info["dp_model"] = "average"

    c = atoms.copy()
    print(f"copy: dp_test in arrays: {'dp_test' in c.arrays}, dp_model in info: {'dp_model' in c.info}")

    p = pickle.loads(pickle.dumps(atoms))
    print(f"pickle: dp_test in arrays: {'dp_test' in p.arrays}, dp_model in info: {'dp_model' in p.info}")

def test_degenerate_cells():
    print("\n--- degenerate cells ---")
    from ase import Atoms

    # Left-handed
    left_matrix = np.array([[-5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]])
    try:
        a = Atoms("Ga", positions=[[0, 0, 0]], cell=left_matrix, pbc=True)
        print(f"left-handed: constructor OK, det={np.linalg.det(a.cell.array):.1f}")
        print(f"  get_scaled_positions(): {a.get_scaled_positions()}")
    except Exception as e:
        print(f"left-handed: {type(e).__name__}: {e}")

    # Zero
    zero_matrix = np.zeros((3, 3))
    try:
        a = Atoms("Ga", positions=[[0, 0, 0]], cell=zero_matrix, pbc=True)
        print(f"zero cell: constructor OK, det={np.linalg.det(a.cell.array):.1f}")
        print(f"  get_scaled_positions(): {a.get_scaled_positions()}")
    except Exception as e:
        print(f"zero cell: {type(e).__name__}: {e}")

def main():
    print("=== P10: ASE deep inspection ===")
    for f in ["ga07in03as_mixed.cif", "triclinic_p1bar.cif", "gaas_cod_9008845.cif"]:
        inspect_ase(f)
    test_copy_pickle()
    test_degenerate_cells()
    print("=== P10 done ===")

if __name__ == "__main__":
    main()
