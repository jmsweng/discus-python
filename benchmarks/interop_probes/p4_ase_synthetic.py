"""P4: ASE on the 5 synthetic CIFs (with info)."""
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def inspect(fname):
    cif = os.path.join(FIXTURES, fname)
    print(f"\n--- {fname} ---")
    try:
        from ase.io import read
        atoms = read(cif)
        print(f"n_atoms: {len(atoms)}")
        print(f"species: {atoms.get_chemical_symbols()}")
        print(f"cell.array:\n{atoms.cell.array}")
        print(f"pbc: {atoms.pbc}")
        print(f"info: {atoms.info}")
        print(f"arrays keys: {list(atoms.arrays.keys())}")
        print(f"scaled positions:\n{atoms.get_scaled_positions()}")
    except Exception as e:
        print(f"ERROR reading: {type(e).__name__}: {e}")

def main():
    print("=== P4: ASE on 5 synthetic CIFs ===")
    for f in ["triclinic_p1bar.cif", "ga07in03as_mixed.cif", "p1_aniso_adp.cif",
              "p1_disorder_groups.cif", "p1_occupancy_over_one.cif"]:
        inspect(f)
    print("=== P4 done ===")

if __name__ == "__main__":
    main()
