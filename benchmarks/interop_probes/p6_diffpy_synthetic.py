"""P6: diffpy.structure on the 5 synthetic CIFs."""
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def inspect(fname):
    cif = os.path.join(FIXTURES, fname)
    print(f"\n--- {fname} ---")
    try:
        from diffpy.structure import Structure
        s = Structure().read(cif)
        st = s.stru
        print(f"n_atoms: {len(st)}")
        print(f"lattice a,b,c: {st.lattice.a}, {st.lattice.b}, {st.lattice.c}")
        print(f"lattice alpha,beta,gamma: {st.lattice.alpha}, {st.lattice.beta}, {st.lattice.gamma}")
        print(f"spacegroup: {s.spacegroup}")
        print(f"cif_sgname: {s.cif_sgname}")
        for i, atom in enumerate(st):
            print(f"  atom {i}: {atom.element} label={atom.label} pos={atom.xyz} occ={atom.occupancy} U={atom.Uisoequiv} B={atom.Bisoequiv}")
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")

def main():
    print("=== P6: diffpy.structure on 5 synthetic CIFs ===")
    for f in ["triclinic_p1bar.cif", "ga07in03as_mixed.cif", "p1_aniso_adp.cif",
              "p1_disorder_groups.cif", "p1_occupancy_over_one.cif"]:
        inspect(f)
    print("=== P6 done ===")

if __name__ == "__main__":
    main()
