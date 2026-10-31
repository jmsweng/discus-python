"""P1: gemmi small-structure API on the real COD GaAs CIF."""
import gemmi
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print(f"=== P1: gemmi small-structure on gaas_cod_9008845.cif ===")
    s = gemmi.read_small_structure(CIF)
    print(f"block name: {s.name}")
    print(f"spacegroup_hm: {s.spacegroup_hm}")
    print(f"spacegroup_number: {s.spacegroup_number}")
    print(f"cell: {s.cell}")
    print(f"n_listed_sites: {len(s.sites)}")

    # Listed sites
    for i, site in enumerate(s.sites):
        print(f"  listed[{i}]: label={site.label!r} type_symbol={site.type_symbol!r} "
              f"fract={site.fract} occ={site.occ} u_iso={site.u_iso} "
              f"aniso={site.aniso} disorder={site.disorder_group!r}")

    # All unit-cell sites
    all_sites = s.get_all_unit_cell_sites()
    print(f"n_unit_cell_sites: {len(all_sites)}")
    for i, site in enumerate(all_sites):
        x, y, z = site.fract.x, site.fract.y, site.fract.z
        in_range = 0 <= x < 1 and 0 <= y < 1 and 0 <= z < 1
        print(f"  uc[{i}]: label={site.label!r} type_symbol={site.type_symbol!r} "
              f"x={x} y={y} z={z} in_[0,1)={in_range} "
              f"occ={site.occ} u_iso={site.u_iso} "
              f"aniso={site.aniso} disorder={site.disorder_group!r}")

    print("=== P1 done ===")

if __name__ == "__main__":
    main()
