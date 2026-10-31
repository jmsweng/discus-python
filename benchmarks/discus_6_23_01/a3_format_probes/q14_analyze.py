"""Q14: count distinct values per column for the XYZBPMMOS atom lines.

Columns (after the atom name): x, y, z, Biso, Property, MoleNo, MoleAt, Occ, St, Sh, Sk, Sl
We report distinct values for the non-coordinate columns:
  Biso (idx 3), Property (4), MoleNo (5), MoleAt (6), Occ (7), St (8), Sh (9), Sk (10), Sl (11)
"""
import sys
from collections import Counter


def analyze(path):
    with open(path, "r") as f:
        lines = f.read().splitlines()
    # find the atoms header line
    start = None
    for i, l in enumerate(lines):
        if l.strip().lower().startswith("atoms"):
            start = i + 1
            break
    if start is None:
        print(f"{path}: no atoms line")
        return
    cols = {
        "Biso": Counter(),
        "Property": Counter(),
        "MoleNo": Counter(),
        "MoleAt": Counter(),
        "Occ": Counter(),
        "St": Counter(),
        "Sh": Counter(),
        "Sk": Counter(),
    }
    n_atoms = 0
    for l in lines[start:]:
        s = l.strip()
        if not s:
            continue
        # atom line: name then 11 comma-separated fields
        # (x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk) -> 12 parts total
        parts = [p.strip() for p in s.split(",")]
        if len(parts) < 12:
            continue
        n_atoms += 1
        # parts[0]=name, parts[1..11] = x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk
        cols["Biso"][parts[3]] += 1
        cols["Property"][parts[4]] += 1
        cols["MoleNo"][parts[5]] += 1
        cols["MoleAt"][parts[6]] += 1
        cols["Occ"][parts[7]] += 1
        cols["St"][parts[8]] += 1
        cols["Sh"][parts[9]] += 1
        cols["Sk"][parts[10]] += 1
    print(f"=== {path} ({n_atoms} atom lines) ===")
    for name in ["Biso", "Property", "MoleNo", "MoleAt", "Occ", "St", "Sh", "Sk"]:
        c = cols[name]
        distinct = sorted(c.keys())
        print(f"  {name:9s}: {len(distinct)} distinct -> {distinct[:6]}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        analyze(p)
