# Magnetic neutron scattering — investigation and experimental reference (Package D, Feature 11)

**Status: experimental.** DISCUS 6.23.01 calculates no magnetic scattering, so discus-python has no DISCUS
magnetic workflow to wrap. `discuspdf.magnetic_reference` is a small independent reference with analytical
tests; it is not exported from `discuspdf` and is not a validated production feature.

## What DISCUS 6.23.01 supports (source read at tag `v.6.23.01`)

| Item | Finding |
|---|---|
| Data model | per-atom moment `cr_magn(0:3, N)` (`crystal_mod.f90:135`); component 0 is a flag, 1–3 the vector; settable via `magn[i,j]` (`discus_upd_par.f90:865-870`) |
| Symmetry | space-group expansion rotates the moment as a **polar** vector (`spcgr_apply.f90:2519-2524`): no det(R) factor (an axial vector needs one), no time reversal — magnetic space groups are not supported |
| `.stru` read / save | read sets moments to 0 (`structur.f90:1483, 1494, 2887`, marked `MAGNETIC_WORK`); `save write magnetic` fills `wr_magn` but never writes it (`discus_save.f90:1002-1007`) |
| HDF5 "unified" | read copies `/entry/data/magnetic_spins` (`structur.f90:7272`); write allocates zeros and never fills them (`discus_save.f90:1345-1347`) |
| Scattering | no Fourier, powder, 3D-PDF or PDF routine uses `cr_magn` (`discus_exp2pdf_run.f90:339` only zeroes it); `element_data_mod.f90` has no ⟨j0⟩/⟨j2⟩ coefficients; no magnetic command in the help |

**Verdict:** DISCUS can carry a moment per atom but cannot compute magnetic diffuse scattering, and its
`.stru` files do not keep moments. discus-python therefore does not add moments to `Structure` or the
DISCUS writers (Package A open item stays open).

## What Javelin does (`javelin/fourier.py:394-450`, HEAD a3538193)

Unpolarised `Σ_α |F⊥,α|²` with ⟨j0⟩ only (ion charge hard-coded to +3, `:413`), no prefactor (arbitrary
units), atoms without a table silently skipped (`:414-417`), NaN at Q = 0 (expected by its test). Its
projection `M − (M·Q)Q/|Q|²` mixes reciprocal-lattice components `2π(h,k,l)` with the Cartesian |Q|², so it
is right only for a Cartesian cell with a = 1 Å; for cubic a = 5 Å a moment parallel to Q gives 576× the
perpendicular value instead of 0. Javelin's test cannot see this (moments ∥ z, Q in the a*b* plane). It also
uses `np.complex`, removed in NumPy 1.24. Design ideas only; no Javelin code is used.

## Specification of the reference

Elastic, unpolarised neutrons, spin-only dipole approximation (Squires, *Thermal Neutron Scattering*, ch. 7):

    dσ/dΩ = p² |F⊥(Q)|²,  F(Q) = Σ_j f_j(|Q|) μ_j e^{iQ·r_j},  F⊥ = F − (F·Q̂)Q̂,  p = γr₀/2 = 2.6954 fm

| Convention | Choice |
|---|---|
| Moments | Cartesian, in μB, in the **same frame as the positions** (Å). No fractional/crystal-axis moments are accepted, so the frame cannot be mixed up |
| Q | Cartesian Å⁻¹ with 2π (`Lattice.hkl_to_q`) |
| Form factor | `f = ⟨j0⟩ + (1 − 2/g)⟨j2⟩`, `⟨j0⟩ = A e^{−as²} + B e^{−bs²} + C e^{−cs²} + D`, `⟨j2⟩ = s² × (same form)`, `s = |Q|/4π`; coefficients given explicitly with a source, or from `periodictable` (Brown, Int. Tables C 4.4.5) via `MagneticFormFactor.from_periodictable(element, charge, g)`. `g ≠ 2` requires ⟨j2⟩. `g` is the spectroscopic g of a spin-dominated (3d) ion, L = (g − 2)S; the rare-earth Landé-g_J form ⟨j0⟩ + (2/g_J − 1)⟨j2⟩ has the opposite sign and is not provided (refused for lanthanides/actinides with g ≠ 2) |
| Units | intensity in fm² for the whole sample (`p²` included); `magnetic_structure_factor` returns F⊥ in μB |
| Q = 0 | Q̂ is undefined: intensity NaN and `valid` mask False for `|Q| < q_min` (default 1e-9 Å⁻¹) |
| Refused | missing form factor for a species, non-(N,3) or non-finite inputs, `g ≠ 2` without ⟨j2⟩ |
| Not provided | polarisation analysis, inelastic scattering, orbital terms beyond the dipole (1 − 2/g)⟨j2⟩ term, magnetic space groups, Debye–Waller factor, nuclear–magnetic interference (zero for unpolarised neutrons) |

```python
from discuspdf.magnetic_reference import MagneticFormFactor, magnetic_intensity
ff = MagneticFormFactor.from_periodictable("Mn", 2)                     # pip install periodictable
I, valid = magnetic_intensity(positions_cartesian, moments_cartesian, lattice.hkl_to_q(hkl), {"Mn": ff},
                              species=species)                        # fm², NaN where Q = 0
```

## Analytical benchmarks (`tests/test_magnetic_reference.py`)

| Case | Expected |
|---|---|
| Single moment at angle θ to Q | `p² f² μ² sin²θ` (θ = 0 … 135°) |
| μ ∥ Q in a triclinic cell (a = 5.3 Å) | exactly 0, and `p² μ²` for μ ⊥ Q — the check Javelin's projection fails |
| Ferromagnetic / antiferromagnetic pair | `× 4cos²(Q·d/2)` / `× 4sin²(Q·d/2)` |
| AF chain, N = 8 | peaks `p² N²` at h = ½, 3/2, zero at integer h; ferromagnetic chain the opposite |
| Powder average | ⟨sin²θ⟩ = 2/3 (Fibonacci sphere, 2·10⁴ directions, 1e-4) |
| Q = 0 | NaN + mask |
| Form factor | `periodictable` j0/j2 values (V³⁺, Fe³⁺ with g = 1.9 — the same formula, so a consistency check, not an independent literature value); j0(V³⁺, 2π/5 Å⁻¹) = 0.88829, consistent with Javelin's expected intensity 493.1621 = 625·j0² (`tests/fourier_test.py`) |

## Before this could become a production feature

An independent numerical comparison (e.g. Spinvert, SCATTY or a published magnetic diffuse pattern), moments in
`Structure` with I/O that keeps them (DISCUS cannot, see above), and an agreed treatment of magnetic
symmetry.
