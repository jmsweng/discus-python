"""Q15 probe: write q15.stru and q15.cll into C:\\tmp\\a3probe.

q15.stru: 7 CU atoms with Biso values 0.0, 1e-6, 1e-5, 1e-3, 0.35, 2.0, -0.1
q15.cll:  same 7 atoms as "CU x y z B" with the same B values.

Atom line format matches gaas_111.stru exactly (123 chars):
  sp(10, left) x(9) , 7sp y(8) , 7sp z(8) , 6sp B(8) , TAIL
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
    (0.84, 0.00, 0.00, -0.100000),
]

TAIL = "       1,       0,       0,   1.000000, _,   0,   0,   0"

def make_atom_line(sp, x, y, z, b):
    sp_field = sp.ljust(10)
    x_field = f"{x:.6f}".rjust(9)
    y_field = f"{y:.6f}".rjust(8)
    z_field = f"{z:.6f}".rjust(8)
    b_field = f"{b:.6f}".rjust(8)
    return f"{sp_field}{x_field},       {y_field},       {z_field},      {b_field},{TAIL}"

stru_lines = [
    "title Q15 B-probe",
    "spcgr P 1             , setting:abc",
    "cell    5.000000,  5.000000,  5.000000, 90.000000, 90.000000, 90.000000",
    "format numbers,XYZBPMMOS",
    "atoms      x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl",
]
for x, y, z, b in atoms:
    stru_lines.append(make_atom_line("CU", x, y, z, b))

with open(os.path.join(OUT_DIR, "q15.stru"), "w") as f:
    f.write("\n".join(stru_lines) + "\n")

# --- q15.cll ---
cll_lines = [
    "title Q15 B-probe",
    "spcgr P1",
    "cell 5.0 5.0 5.0 90.0 90.0 90.0",
    "atoms",
]
for x, y, z, b in atoms:
    cll_lines.append(f"CU    {x:7.4f}  {y:7.4f}  {z:7.4f}  {b:8.4f}")

with open(os.path.join(OUT_DIR, "q15.cll"), "w") as f:
    f.write("\n".join(cll_lines) + "\n")

print("Wrote q15.stru and q15.cll to", OUT_DIR)
print(f"First atom line ({len(stru_lines[5])} chars): {stru_lines[5]}")
