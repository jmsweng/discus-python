# Lattice core (Package A, Feature A5)

`discuspdf.Lattice` is an immutable crystal-lattice value type. It is the single shared
lattice implementation for Package A; `Structure` (A1) will build on it.

## Conventions

These match the existing, validated `DiscusVolume.direct_basis()` / `reciprocal_basis()`:

| Quantity | Representation |
|---|---|
| Direct basis `A` | 3×3, **rows** are **a**, **b**, **c** in Å |
| Fractional → Cartesian | row vectors: `r = f @ A` |
| Reciprocal basis `B` | `B = inv(A).T`, rows **a\***, **b\***, **c\*** in Å⁻¹, **no 2π** (`A @ B.T = I`) |
| Cartesian Q | `Q = 2π · (hkl @ B)` in Å⁻¹ (2π included); hkl in r.l.u. |
| d-spacing | `d_hkl = 1 / |hkl @ B|` (Å); `hkl = (0,0,0)` raises |
| Metric | `G = A @ A.T`; reciprocal metric `G* = inv(G) = B @ B.T` |
| Volume | `det(A)` (Å³), must be > 0 |
| Angles | degrees in/out; lengths Å |

## Construction

- `Lattice.from_parameters(a, b, c, alpha, beta, gamma)` — canonical orientation
  (**a** ∥ x, **b** in the xy plane, right-handed). Same rejection rules as
  `discus_output._validated_direct_basis`: lengths finite > 0; angles finite in (0°, 180°);
  metric eigenvalue ratio ≤ 1e-12 → degenerate; c_z² ≤ 0 → invalid; det ≤ 0 → invalid.
  It does **not** reconstruct an experimental orientation matrix.
- `Lattice.from_matrix(A)` — preserves any supplied right-handed orientation.
  Rejects `det(A) ≤ 0` (singular or left-handed) with a `ValueError`; there is no
  left-handed opt-in.

## API

- properties: `matrix` (A), `reciprocal_matrix` (B, no 2π), `parameters` (a, b, c, α, β, γ),
  `volume`, `metric`, `reciprocal_metric`
- `frac_to_cart(f)`, `cart_to_frac(r)` — accept shape `(3,)` or `(..., 3)`, return the same leading shape
- `hkl_to_q(hkl)`, `q_to_hkl(q)` — 2π convention; units stated in docstrings
- `d_spacing(hkl)`
- `minimum_image(df) -> (df_min, n_shift)` for fractional displacements; `distance_pbc(f1, f2)` in Å
- `__eq__` / `isclose(other, rtol, atol)`; readable `__repr__`

All stored arrays are read-only; getters return copies or read-only views. Non-finite inputs
raise `ValueError` and never return NaN.

## Minimum image — algorithm

Rounding each fractional component is **not** correct in skewed cells. `minimum_image` uses
an exact bounded search:

1. `n0 = round(df)`, `r0 = (df − n0) @ A`, `R = |r0|`.
2. Any better image satisfies `|r| ≤ R`. Since `f_i = r · b_i` (rows of `B`), each component
   of the offset `k = n − n0` lies in `[−R|b_i|, R|b_i|]`. Enumerate integer offsets in that
   range per axis (vectorised, in chunks so peak memory is independent of the number of
   displacements), keep the shortest.
3. **Ties:** among images within `atol = 1e-12·R + 1e-15` Å of the shortest, return the one
   whose **reduced fractional vector `df − n` is lexicographically smallest** (component 0,
   then 1, then 2; components equal if they differ by < 1e-12). In a cubic cell this maps
   half-cell ties to −½, i.e. the half-open interval [−½, ½). The result is deterministic
   and **translation invariant**: `minimum_image(df + n)[0] == minimum_image(df)[0]` for
   integer `n` (a function of `df mod 1`).

## Not yet (later PRs)

- `Structure` (A1), CIF/ASE/DiffPy interop (A2), DISCUS `.cll`/`.stru` I/O (A3),
  supercells and atom manipulation (A4).
- Consolidating `DiscusVolume` / `discus_output` onto `Lattice` (separate PR, with Gate 2/3
  evidence parity).

## Example

`examples/lattice_basics.py` (run with `python examples/lattice_basics.py`); covered by
`tests/test_lattice.py::test_example_script_runs`.

Notebook: [`notebooks/03_lattice_basics.ipynb`](../notebooks/03_lattice_basics.ipynb) —
worked examples of the API (lattice construction, validation, fractional/Cartesian
round trips, Q and d-spacings, minimum image, immutability).
