# Atomic-structure pictures and viewer export (Package D, Feature 21)

`discuspdf.structure_viz` draws a `Structure` in a notebook with matplotlib and writes files for external
viewers through ASE. It never changes coordinates: atoms are drawn and written at `frac @ lattice.matrix`
(Å), the cell exactly as `lattice.matrix` (rows a, b, c), nothing wrapped or rotated. Importing
`discuspdf` does not import matplotlib or ASE (tested), so headless runs are unaffected.

```python
from discuspdf.structure_viz import plot_structure, plot_displacements, plot_bonds, export_structure
from discuspdf import neighbor_list

plot_structure(configuration, color_by="species")                  # or color_by="site", select=..., cell="parent"
plot_displacements(configuration, reference, scale=10)             # arrows × 10, stated in the title
plot_bonds(neighbor_list(reference, 4.0).select(shell=1))          # first-shell pairs
export_structure(configuration, "model.xyz", format="extxyz")      # VESTA, OVITO, Jmol, ase gui …
```

| Function | What it shows |
|---|---|
| `plot_structure` | atoms coloured by species or site label (fixed Okabe–Ito colours per label), marker area × occupancy, supercell or parent-cell frame |
| `plot_displacements` | arrows from each reference site to its atom; atoms matched by `(site_label, cell_index)` (so the reference needs one atom per site), true minimum image in the supercell (`Lattice.minimum_image`) with all pbc True, plain difference with none; mixed pbc refused |
| `displacement_vectors` | the numbers behind those arrows: reference positions and Cartesian displacements |
| `plot_bonds` | one segment per directed pair of a `NeighborList` selection, from the origin atom to the periodic image that forms the pair; distances grouped into colours where they jump by more than `shell_tol` (0.05 Å), legend (distance range, pair count) for up to 8 groups |
| `export_structure` | `ase.io.write`, then reads the file back and refuses it if the cell, positions (1e-5 Å) or elements changed |

## Formats

| Format | Keeps a rotated cell? | Notes |
|---|---|---|
| `extxyz` | yes (full 3×3 lattice) | best for OVITO / ASE; also keeps no DISCUS-specific columns |
| `vasp` (POSCAR) | yes | |
| `xsf` | yes | XCrySDen, VESTA |
| `cif` | **no** — stores cell parameters only, and ASE wraps atoms into the cell | fine for a `Lattice.from_parameters` cell with atoms inside it; anything else is refused by the read-back check. Slow for big supercells (> 2 min for 31 000 atoms): use extxyz |

`xsf` and some other formats drop the cell of a non-periodic structure (pbc all False); the check says so.

The ASE extra (`pip install 'discuspdf-starter[ase]'`) is needed only for `export_structure`.

## Not covered
Interactive 3D viewers (use `ase.visualize.view(to_ase(structure))` or open an exported file), bonds by
chemistry rules, polyhedra, molecule drawing, magnetic moments.
