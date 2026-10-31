"""P2: gemmi small-structure API on the 5 synthetic CIFs."""
import gemmi
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")

def inspect(fname):
    cif = os.path.join(FIXTURES, fname)
    s = gemmi.read_small_structure(cif)
    print(f"\n--- {fname} ---")
    print(f"block name: {s.name}")
    print(f"spacegroup_hm: {s.spacegroup_hm}")
    print(f"spacegroup_number: {s.spacegroup_number}")
    print(f"cell: {s.cell}")
    print(f"n_listed_sites: {len(s.sites)}")

    for i, site in enumerate(s.sites):
        print(f"  listed[{i}]: label={site.label!r} type_symbol={site.type_symbol!r} "
              f"fract={site.fract} occ={site.occ} u_iso={site.u_iso} "
              f"aniso={site.aniso} disorder={site.disorder_group!r}")

    all_sites = s.get_all_unit_cell_sites()
    print(f"n_unit_cell_sites: {len(all_sites)}")
    for i, site in enumerate(all_sites):
        x, y, z = site.fract.x, site.fract.y, site.fract.z
        in_range = 0 <= x < 1 and 0 <= y < 1 and 0 <= z < 1
        print(f"  uc[{i}]: label={site.label!r} type_symbol={site.type_symbol!r} "
              f"x={x} y={y} z={z} in_[0,1)={in_range} "
              f"occ={site.occ} u_iso={site.u_iso} "
              f"aniso={site.aniso} disorder={site.disorder_group!r}")

def main():
    print("=== P2: gemmi small-structure on 5 synthetic CIFs ===")
    for f in ["triclinic_p1bar.cif", "ga07in03as_mixed.cif", "p1_aniso_adp.cif",
              "p1_disorder_groups.cif", "p1_occupancy_over_one.cif"]:
        inspect(f)
    print("=== P2 done ===")

if __name__ == "__main__":
    main()
