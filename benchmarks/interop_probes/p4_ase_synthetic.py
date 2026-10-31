"""P4: ASE on the 5 synthetic CIFs."""
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
        print(f"cell: {atoms.cell}")
        print(f"cell lengths: {atoms.cell.lengths()}")
        print(f"cell angles: {atoms.cell.angles()}")
        print(f"positions (fractional): {atoms.get_positions()}")
        print(f"spacegroup: NOT AVAILABLE from ASE CIF reader")
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
