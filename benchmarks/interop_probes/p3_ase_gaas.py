"""P3: ASE on the real COD GaAs CIF."""
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print("=== P3: ASE on gaas_cod_9008845.cif ===")
    try:
        import ase
        from ase.io import read
        print(f"ase version: {ase.__version__}")
        atoms = read(CIF)
        print(f"n_atoms: {len(atoms)}")
        print(f"species: {atoms.get_chemical_symbols()}")
        print(f"cell: {atoms.cell}")
        print(f"cell lengths: {atoms.cell.lengths()}")
        print(f"cell angles: {atoms.cell.angles()}")
        print(f"positions (fractional): {atoms.get_positions()}")
        # ASE does not expose spacegroup from CIF reading
        print(f"spacegroup: NOT AVAILABLE from ASE CIF reader")
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")
    print("=== P3 done ===")

if __name__ == "__main__":
    main()
