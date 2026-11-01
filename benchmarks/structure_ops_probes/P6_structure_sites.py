"""P6: Structure.sites() for the mixed fixture read through read_cif."""
from discuspdf.cif_io import read_cif

cif = read_cif("tests/fixtures/interop/gaas_cod_9008845.cif")
s = cif.structure
print(f"n_atoms={s.n_atoms}, model={s.model}")

sites = s.sites()
print(f"n_sites={len(sites)}")
print("Order of sites (is it order of first appearance, or sorted?):")
for i, site in enumerate(sites):
    print(f"  [{i}] site_label={site['site_label']!r}, cell_index={site['cell_index']}, members={site['members']}, total_occupancy={site['total_occupancy']:.4f}")

# Check if the order matches first appearance in the structure
first_appearance = []
seen = set()
for i in range(s.n_atoms):
    key = (s.site_label[i], tuple(s.cell_index[i]))
    if key not in seen:
        seen.add(key)
        first_appearance.append(key)

site_order = [(site["site_label"], site["cell_index"]) for site in sites]
print(f"\nFirst appearance order: {first_appearance}")
print(f"sites() order:         {site_order}")
print(f"Same order: {first_appearance == site_order}")
