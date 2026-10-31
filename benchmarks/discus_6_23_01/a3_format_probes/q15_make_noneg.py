"""Q15 probe (retry 4): write q15_noneg.stru by copying the exact gaas_111.stru atom line format.

The gaas_111.stru atom line is 123 chars:
  pos 1-11:   species (right-aligned in 11 chars)
  pos 12-20:  x (9 chars, right-aligned)
  pos 21:     comma
  pos 22-27:  6 spaces
  pos 28-35:  y (8 chars, right-aligned)
  pos 36:     comma
  pos 37-43:  7 spaces
  pos 44-51:  z (8 chars, right-aligned)
  pos 52:     comma
  pos 53-58:  6 spaces
  pos 59-66:  B (8 chars, right-aligned)
  pos 67:     comma
  pos 68-123: "       1,       0,       0,   1.000000, _,   0,   0,   0"
"""
import os

OUT_DIR = r"C:\tmp\a3probe"
os.makedirs(OUT_DIR, exist_ok=True)

atoms = [
    (0.00, 0.00, 0.00, 0.000000),
    (0.14, 0.00, 0.00, 0.000001),
    (0.28, 0.00, 0.00, 0.000010),
    (0.42, 0.00, 0.00, 0.001000),
    (0.56, 0.00, 0.00, 0.350000),
    (0.70, 0.00, 0.00, 2.000000),
]

TAIL = "       1,       0,       0,   1.000000, _,   0,   0,   0"

def make_atom_line(sp, x, y, z, b):
    # gaas layout: sp(10, left-aligned) x(9) , 7sp y(8) , 7sp z(8) , 6sp B(8) , TAIL
    sp_field = sp.ljust(10)
    x_field = f"{x:.6f}".rjust(9)
    y_field = f"{y:.6f}".rjust(8)
    z_field = f"{z:.6f}".rjust(8)
    b_field = f"{b:.6f}".rjust(8)
    return f"{sp_field}{x_field},       {y_field},       {z_field},      {b_field},{TAIL}"

# Verify against gaas
test = make_atom_line("GA", 0.0, 0.0, 0.0, 1.0)
print(f"Test line ({len(test)} chars): {test}")
gaas = "GA         0.000000,       0.000000,       0.000000,      1.000000,       1,       0,       0,   1.000000, _,   0,   0,   0"
print(f"Gaas line ({len(gaas)} chars): {gaas}")
print(f"Match: {test == gaas}")

stru_lines = [
    "title Q15 B-probe (no negative B)",
    "spcgr P 1             , setting:abc",
    "cell    5.000000,  5.000000,  5.000000, 90.000000, 90.000000, 90.000000",
    "format numbers,XYZBPMMOS",
    "atoms      x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl",
]
for x, y, z, b in atoms:
    stru_lines.append(make_atom_line("CU", x, y, z, b))

with open(os.path.join(OUT_DIR, "q15_noneg.stru"), "w") as f:
    f.write("\n".join(stru_lines) + "\n")

print(f"\nWrote q15_noneg.stru ({len(stru_lines)} lines)")
print(f"First atom line ({len(stru_lines[5])} chars): {stru_lines[5]}")
