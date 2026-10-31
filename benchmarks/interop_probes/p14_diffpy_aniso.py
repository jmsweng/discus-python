"""P14: DiffPy Uisoequiv with anisotropy (corrected)."""
from diffpy.structure import Lattice, Atom

def main():
    print("=== P14: DiffPy Uisoequiv with anisotropy ===")

    # Cubic lattice, anisotropy=True
    lat_cub = Lattice(5.0, 5.0, 5.0, 90, 90, 90)
    a = Atom("Fe", [0, 0, 0], lattice=lat_cub)
    a.anisotropy = True
    a.U11, a.U22, a.U33 = 0.013, 0.012, 0.011
    print(f"cubic, anisotropy=True, U11=0.013 U22=0.012 U33=0.011: Uisoequiv={a.Uisoequiv}")

    # Triclinic lattice, anisotropy=True
    lat_tri = Lattice(4.1, 5.3, 6.2, 78, 84, 101)
    a2 = Atom("Fe", [0, 0, 0], lattice=lat_tri)
    a2.anisotropy = True
    a2.U11, a2.U22, a2.U33 = 0.013, 0.012, 0.011
    print(f"triclinic, anisotropy=True, U11=0.013 U22=0.012 U33=0.011: Uisoequiv={a2.Uisoequiv}")

    # Isotropic case: anisotropy=False, Uisoequiv set directly
    a3 = Atom("Fe", [0, 0, 0], lattice=lat_cub)
    a3.Uisoequiv = 0.012
    print(f"isotropic, anisotropy={a3.anisotropy}, Uisoequiv={a3.Uisoequiv}, "
          f"U11={a3.U11} U22={a3.U22} U33={a3.U33} (diag all 0.012)")

    print("=== P14 done ===")

if __name__ == "__main__":
    main()
