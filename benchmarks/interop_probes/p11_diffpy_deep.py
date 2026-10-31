"""P11: DiffPy deep inspection (Lattice.base, Atom, Uisoequiv, iteration)."""
import os
import numpy as np

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def main():
    print("=== P11: DiffPy deep inspection ===")
    from diffpy.structure import Lattice, Atom, Structure

    # Lattice.base vs our matrix (triclinic)
    dl = Lattice(4.1, 5.3, 6.2, 78, 84, 101)
    print(f"diffpy Lattice.base:\n{dl.base}")
    print(f"type: {type(dl.base)}")

    from discuspdf.lattice import Lattice as OurLattice
    ol = OurLattice.from_parameters(4.1, 5.3, 6.2, 78.0, 84.0, 101.0)
    print(f"our matrix:\n{ol.matrix}")
    max_diff = np.max(np.abs(np.array(dl.base) - ol.matrix))
    print(f"max abs diff (Lattice.base vs our matrix): {max_diff:.2e}")

    # Lattice.cell_parms (replacement for deprecated abcABG)
    print(f"cell_parms: {dl.cell_parms}")

    # Atom(atype="O2-")
    a = Atom(atype="O2-", xyz=[0, 0, 0])
    print(f"\nAtom(atype='O2-'): element={a.element!r}")

    # Uisoequiv default
    a2 = Atom(atype="Ga", xyz=[0, 0, 0])
    print(f"Uisoequiv default: {a2.Uisoequiv}")

    # Uisoequiv with U tensor
    a3 = Atom(atype="Fe", xyz=[0, 0, 0])
    a3.U11 = 0.013
    a3.U22 = 0.012
    a3.U33 = 0.011
    print(f"Uisoequiv with U11=0.013 U22=0.012 U33=0.011: {a3.Uisoequiv}")

    # Iteration order and xyz type
    s = Structure(lattice=dl, atoms=[a, a2])
    print(f"iteration order: {[atom.element for atom in s]}")
    print(f"xyz type: {type(a.xyz)}")

    print("=== P11 done ===")

if __name__ == "__main__":
    main()
