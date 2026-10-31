"""P13: gemmi.cif.as_string quote-stripping (corrected)."""
import gemmi
import os

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "fixtures", "interop")
CIF = os.path.join(FIXTURES, "gaas_cod_9008845.cif")

def main():
    print("=== P13: gemmi.cif.as_string quote-stripping ===")
    doc = gemmi.cif.read(CIF)
    b = doc.sole_block()

    # Quoted space-group value
    v = b.find_value("_symmetry_space_group_name_H-M")
    print(f"find_value('_symmetry_space_group_name_H-M') raw: {v!r}")
    print(f"as_string: {gemmi.cif.as_string(v)!r}")

    # Quoted loop value (use a file with a quoted loop value)
    # The GaAs CIF labels are unquoted; use the mixed CIF which has quoted values
    doc2 = gemmi.cif.read(os.path.join(FIXTURES, "ga07in03as_mixed.cif"))
    b2 = doc2.sole_block()
    sg2 = b2.find_value("_symmetry_space_group_name_H-M")
    print(f"\nmixed CIF find_value('_symmetry_space_group_name_H-M') raw: {sg2!r}")
    print(f"as_string: {gemmi.cif.as_string(sg2)!r}")

    # find_values on a loop (labels are unquoted in these CIFs)
    col = b.find_values("_atom_site_label")
    print(f"\nfind_values('_atom_site_label')[0] raw: {col[0]!r}")
    print(f"as_string: {gemmi.cif.as_string(col[0])!r}")

    # Confirm: as_string strips quotes from quoted strings, is identity on unquoted
    print(f"\nSummary: as_string strips CIF single-quotes; identity on unquoted values.")
    print("=== P13 done ===")

if __name__ == "__main__":
    main()
