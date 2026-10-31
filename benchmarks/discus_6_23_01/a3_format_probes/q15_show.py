"""Q15 show: print a table of B written vs B saved for the two saved files.

For each saved .stru file, parse the 4th number after the species name (Biso column)
and print one row per atom.
"""
import os

PROBE_DIR = r"C:\tmp\a3probe"

def parse_stru_biso(path):
    """Return list of (species, x, b_iso) for each atom line in a .stru file."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Skip header lines
            if line.startswith(("title", "spcgr", "cell", "format", "atoms", "ncell")):
                continue
            # Atom line: species followed by 12 comma-separated numbers
            parts = line.split()
            if len(parts) < 2:
                continue
            species = parts[0]
            # Parse the numeric fields (comma-separated)
            nums = []
            for p in parts[1:]:
                # Remove trailing comma
                p = p.rstrip(",")
                try:
                    nums.append(float(p))
                except ValueError:
                    pass
            if len(nums) >= 4:
                x, y, z, b = nums[0], nums[1], nums[2], nums[3]
                rows.append((species, x, b))
    return rows

def main():
    for label, fname in [("read-stru", "q15_stru_out.stru"), ("read-cell", "q15_cll_out.stru")]:
        path = os.path.join(PROBE_DIR, fname)
        if not os.path.exists(path):
            print(f"[{label}] {fname} NOT FOUND")
            continue
        rows = parse_stru_biso(path)
        print(f"\n=== {label} ({fname}) ===")
        print(f"{'idx':>3}  {'species':>6}  {'x':>10}  {'Biso (saved)':>14}")
        for i, (sp, x, b) in enumerate(rows):
            print(f"{i:>3}  {sp:>6}  {x:>10.6f}  {b:>14.6f}")

if __name__ == "__main__":
    main()
