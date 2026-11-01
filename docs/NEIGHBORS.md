# Neighbor lists (Package C, Feature 12)

`discuspdf.neighbors` lists every pair of crystal sites closer than a cutoff, on any lattice
(including triclinic) and with periodic images handled exactly. It is the basis of the
occupational and displacement correlations of Package C.

```python
from discuspdf import neighbor_list

nl = neighbor_list(structure, cutoff=6.0)          # or structure.neighbors(6.0)
nl.shells()                                        # distance shells with coordination numbers
first = nl.select(shell=1)                         # nearest neighbors of every site
ab = nl.select(origin_site="Ga1", target_site="As1", translation=(0, 0, 0), symmetric=True)
nl.summary()                                       # coordination table per site pair and shell
nl.to_csv("pairs.csv")                             # one row per pair
```

## What a pair is

The list is built on one `Structure` (normally a supercell from `repeat`, a DISCUS `.stru`, or a
configuration from `sample_configuration`). Its **nodes** are the crystal sites, i.e. the
distinct (`site_label`, `cell_index`) pairs:

* in a `"configuration"` structure every site holds one atom, so node = atom
  (`node_atom_id`, `node_species` give it);
* in an `"average"` structure a mixed site (several species on one position) is **one** node;
  its `node_atom_id` is −1 and `node_species` is `None`. Species selection then needs a
  configuration (the correlation functions read species from a configuration).

A **directed pair** `(i → j, shift)` means: node *j*, moved by the integer `shift` of the
structure's own cell (the supercell), seen from node *i*:

| Quantity | Definition | Units |
|---|---|---|
| `vector` | `(frac[j] + shift − frac[i]) @ A`, `A = structure.lattice.matrix` | Å, Cartesian |
| `distance` | `‖vector‖` | Å |
| `translation` | `cell_index[j] + shift · supercell − cell_index[i]` | parent unit cells (integers) |
| `cell_vector` | `(frac[j] + shift − frac[i]) · supercell` | fractions of the parent cell |

`translation` is the DISCUS way of naming a neighbor: `(origin site, target site, translation)`
is what `set vec, n, is1, is2, dx, dy, dz` defines. `cell_vector = translation + (position of the
target site − position of the origin site)` in the parent cell.

`shift` refers to the coordinates **as stored**; coordinates outside [0, 1) (DISCUS writes
centred crystals with negative fractions) are allowed and change only `shift`, not the pairs.

## Which pairs are listed

* Every pair with `distance ≤ cutoff` (relative tolerance 1e-12), in **both directions**:
  `(i → j, s)` and `(j → i, −s)`. `unique()` keeps one direction of each unordered pair
  (`i < j`, or `i = j` with `shift` lexicographically positive).
* **All periodic images** within the cutoff, not only the minimum image. If the cutoff is longer
  than the cell, a node can be its own neighbor (`i = j`, `shift ≠ 0`), and two nodes can be
  neighbors through several images. Only `i = j` with `shift = 0` is excluded.
* `pbc` all True: periodic; all False: no images (`shift = 0`); mixed PBC raises
  `NotImplementedError`, as elsewhere in the package.
* The full metric is used (no per-axis rounding), so skewed cells are exact.

## Completeness of the search (why no pair is missed)

Positions are wrapped into [0, 1) and binned into `nb_k` slabs per lattice direction, each at
least one cutoff thick where the cell allows it (`nb_k = floor(d_k / cutoff)`, `d_k` = spacing of
the (100)/(010)/(001) planes of the cell). Any vector of length ≤ cutoff has fractional
components `|Δf_k| ≤ cutoff / d_k`, which crosses at most `floor(cutoff · nb_k / d_k) + 1`
slab boundaries, so searching that many slabs on each side (with wrap-around and the
corresponding image shift) reaches every image. The search is checked against an independent
brute-force enumeration in `tests/test_neighbors.py` for cubic, tetragonal, two triclinic cells,
cells smaller than the cutoff, unwrapped coordinates and open boundaries.

Speed: an fcc 25×25×25 supercell (62 500 atoms) with cutoff 3 Å (750 000 pairs) takes about 1.5 s.

## Selecting pairs — `select(...)`

All criteria given must hold (logical AND). The result keeps the node table, so node indices
stay comparable between selections.

| Argument | Keeps pairs whose … |
|---|---|
| `origin_site`, `target_site` | site label is in the given label(s) |
| `origin_species`, `target_species` | species is in the given species (needs one atom per site) |
| `distance=value` / `distance=(min, max)` | distance equals the value within `distance_tol` (1e-6 Å) / lies in the range |
| `translation=(da, db, dc)` | parent-cell translation equals the integers exactly |
| `cell_vector=(x, y, z)` | parent-cell fractional vector matches within `vector_tol` |
| `direction=…` | direction is within `angle_tol` degrees (default 1e-3°) of the given one; `direction_frame="fractional"` reads it as [uvw] of the **parent** cell, `"cartesian"` (default) as Cartesian |
| `shell=k` | lies in the k-th distance shell **of the selection after the other criteria** |
| `symmetric=True` | also keeps the reversed pairs (origin/target swapped, translation, vector and direction negated) |

`select(origin_site=A, target_site=B, translation=t, symmetric=True)` gives exactly the pairs
DISCUS `chem` uses for `set vec, n, A, B, t`: an atom on site A in cell c is paired with site B in
cell c + t, and an atom on B with site A in c − t (checked against DISCUS index arithmetic in the
tests). When A = B, DISCUS takes only the `+t` direction for each atom (one pair per atom),
whereas `symmetric=True` keeps both `+t` and `−t`; correlations that pool both directions are
unaffected, pair counts differ by a factor 2.

## Shells and coordination

`shells(tol=1e-3)` sorts the pairs by distance and starts a new shell wherever the gap to the
previous distance exceeds `tol` Å. Each `Shell` gives the mean/min/max distance, the number of
directed pairs, and `coordination = n_pairs / n_origin_nodes` (over all origin nodes of the
list, so an origin with no neighbor in that shell counts as zero). Shells of a **displaced**
configuration are blurred; build the list on the undisplaced reference (same sites and
`cell_index`) and evaluate correlations on the configuration (Features 13, 14).

The last shell can be incomplete if the cutoff cuts through it. Choose the cutoff a little
above the last distance you need.

`coordination()` gives the per-node count of the current selection; `summary()` gives a table
per (origin site, target site, shell) with the distance, pair count and coordination number.

## Export

`to_table()` returns plain arrays: node indices, atom ids (−1 for mixed sites), site labels,
species, distance, Cartesian vector (`rx, ry, rz`), translation (`ta, tb, tc`), parent-cell
vector (`fa, fb, fc`) and shift (`shift_a, shift_b, shift_c`) — ready for `pandas.DataFrame`.
`to_csv(path)` writes the same columns, full float precision.

## Not covered

Mixed PBC; molecules; neighbor definitions by bond angle (DISCUS `set ang`); symmetry-generated
vector sets (DISCUS `apply_symmetry`).
