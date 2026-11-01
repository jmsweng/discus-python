# Independent structure-factor reference (Package D, Feature 10)

`discuspdf.reference_scattering` computes single-crystal scattering amplitudes by a plain sum over atoms.
It exists to **check** DISCUS and to build analytical fixtures. It is deliberately small: no lots, no
anomalous dispersion, no anisotropic ADPs, no form-factor tables of its own, no instrument effects.
For production calculations use DISCUS (`FourierConfig`).

```python
import numpy as np
from discuspdf import ScatteringFactors, reference_scattering as ref

fac = ScatteringFactors.neutron({"Cu": 7.718, "Au": 7.63}, source="Sears 1992")   # fm
# or ScatteringFactors.from_xraydb(["Cu", "Au"])        X-ray, Waasmaier–Kirfel (pip install xraydb)
# or ScatteringFactors.from_periodictable(["Cu", "Au"]) neutron b_c (pip install periodictable)

hkl = np.stack(np.meshgrid(np.linspace(-2, 2, 41), np.linspace(-2, 2, 41), [0.0], indexing="ij"), -1)
I_total   = ref.intensity(structure, fac, hkl=hkl, debye_waller=True)            # |F|², whole crystal
I_diffuse = ref.diffuse_intensity(structure, fac, hkl=hkl, debye_waller=True)    # |F − ⟨F⟩|²
```

## Conventions

| Quantity | Definition |
|---|---|
| Scattering vector | `Q = 2π · hkl @ B` (Å⁻¹), `B = inv(A)ᵀ`; `hkl` in r.l.u. of the **parent cell** (default, DISCUS convention) or of the supercell (`cell="supercell"`); or pass Cartesian `q=` |
| s | `sin θ / λ = |Q| / 4π` (Å⁻¹) |
| Amplitude | `F(Q) = Σ_j occ_j f_j(s) exp(−B_j s²) exp(+i Q·r_j)`; the Debye–Waller factor only with `debye_waller=True` (needs `b_iso`) |
| Average structure | `⟨F⟩ = L(Q) Σ_sites ⟨f e^{−Bs²} e^{iQ·(r − R_c)}⟩_cells`, `R_c` = origin of the atom's parent cell (`cell_index`), `L(Q) = Σ_c e^{iQ·R_c}` over all cells; a vacancy contributes zero to its site's average. Needs a periodic configuration model |
| Diffuse | `|F − ⟨F⟩|²` — amplitude subtraction, as DISCUS `set aver,100` |
| Units | X-ray `f` in electrons → electrons²; neutron `b` in fm → fm². `normalize="per_atom"` divides by Σ occ, `"per_cell"` by the number of parent cells |

Scattering factors always carry a stated source. Keys match `Structure.species` exactly or
case-insensitively; a species without a factor raises `KeyError`. The `xraydb`/`periodictable`
constructors accept plain element labels only: ions (`O2-`) and elements whose coherent length is
complex must be given explicitly. Note that tables differ: DISCUS 6.23.01 uses b(Au) = 7.63 fm (Sears 1992),
`periodictable` 2.1 gives 7.90 fm.

## Benchmarks (`tests/test_reference_scattering.py`)

| Case | Independent expectation | Tolerance |
|---|---|---|
| Phase sign | one atom at x = ¼: F(100) = +i f | 1e-12 |
| Single atom | `f(s)² e^{−2Bs²}` from the DABAX W-K transcription (PR #8); neutron `b²` flat | 1e-13 rel |
| Two atoms | `f₁² + f₂² + 2 f₁ f₂ cos(Q·d)` | 1e-12 rel |
| Perfect crystal | Laue function `Π sin²(πN_k h_k)/sin²(πh_k)` at 300 random hkl; fcc Bragg and extinctions | 1e-9 rel |
| Average of a perfect crystal | ⟨F⟩ = F, diffuse = 0 | 1e-10 |
| Correlated disorder (random and Ising-ordered binary alloy) | Warren–Cowley sum `Δf² Σ_R C(R) e^{−iQ·R}` with C(R) from rolled occupation grids, on the commensurate grid | 1e-9 rel |
| Off-grid disorder | brute-force pair sum `Δf² |Σ σ_j e^{iQ·r_j}|²` | 1e-9 rel |
| Package C cross-check | C(100)/c(1−c) from the rolled grid = Package C Pearson r_AA | 1e-12 |
| Optional tables | `xraydb` f₀ vs DABAX for Cu, Au (≤ 2e-6); `periodictable` b(Cu) = 7.718 fm | |

### Against DISCUS 6.23.01

`benchmarks/discus_6_23_01/d_reference/` compares the reference with DISCUS intensities (X-ray and neutron,
total and `set aver,100`, cubic and monoclinic cells, B ≠ 0). Criteria were fixed before the run in
[`PREREG.md`](../benchmarks/discus_6_23_01/d_reference/PREREG.md); the result is recorded in that folder's
`RESULT.md`.

## Performance

The sum is vectorised over blocks of Q points (at most 4·10⁶ atom×Q phase factors in memory at a time).
Measured in the cloud sandbox (one core, NumPy 2): about **30 ns per atom per Q point** for ≥ 4 000 atoms,
e.g. 13 824 atoms × 10 000 Q points in 4.2 s. Cost is O(N_atoms × N_Q); a 64³-point volume of a 10⁴-atom
crystal takes ~80 s. No Numba/Cython path: correctness first (Feature 10's rule), DISCUS for big jobs.

## Not covered

Anomalous dispersion (f′, f″), anisotropic ADPs, electron scattering, magnetic scattering (Feature 11),
partial `set aver`, lots, instrument resolution, absorption, absolute cross-section units.
