"""P5: diffpy.structure on the real COD GaAs CIF."""
import os
import importlib.metadata

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print("=== P5: diffpy.structure on gaas_cod_9008845.cif ===")
    try:
        print(f"diffpy.structure version: {importlib.metadata.version('diffpy.structure')}")
        from diffpy.structure import Structure
        s = Structure().read(CIF)
        st = s.stru
        print(f"n_atoms: {len(st)}")
        print(f"lattice: {st.lattice}")
        print(f"lattice a,b,c: {st.lattice.a}, {st.lattice.b}, {st.lattice.c}")
        print(f"lattice alpha,beta,gamma: {st.lattice.alpha}, {st.lattice.beta}, {st.lattice.gamma}")
        print(f"spacegroup: {s.spacegroup}")
        print(f"cif_sgname: {s.cif_sgname}")
        for i, atom in enumerate(st):
            print(f"  atom {i}: {atom.element} label={atom.label} pos={atom.xyz} occ={atom.occupancy} U={atom.Uisoequiv} B={atom.Bisoequiv}")
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")
    print("=== P5 done ===")

if __name__ == "__main__":
    main()
