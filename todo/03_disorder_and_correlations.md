# Package C — Neighbor analysis and disorder correlations

**Status:** Implemented (PRs #30, #31, #32; see [TASKS.md](TASKS.md)) · **Priority:** High · **Suggested PR family:** `feat/correlations-*`

**Reference in Javelin:** [neighborlist.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/neighborlist.py), [structure.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/structure.py#L570-L694), [neighbors tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/neighbors.rst).

## Objective

Give researchers transparent Python measurements of local occupational and displacement disorder, independent of scattering visualization. These quantitative descriptors can later serve as simulation controls and optimization targets.

## Planned features

### 12. Neighbor lists by sites, species, shells and geometry
- [x] Represent neighbor pairs as atom/site IDs plus lattice translations; preserve relation to a supercell and boundary conditions.
- [x] Select by source/target site, elemental species, exact lattice vector, Cartesian distance tolerance, first/second/etc. neighbor shell, and direction.
- [x] Use the **full lattice metric** for non-orthogonal crystals. Implement periodic minimum-image and boundary handling explicitly; avoid counting the same pair twice unless intended.
- [x] Add neighbor count/coordination number summaries and exportable tables.

### 13. Occupational correlation measurements
- [x] Implement and document at least one rigorously defined measure: e.g. pair probabilities and Warren–Cowley short-range-order coefficients, with clear species and sublattice conventions.
- [x] Support binary and, when mathematically defined, multicomponent occupancy; distinguish site-specific occupancies from global concentrations.
- [x] Handle pure end members and zero-denominator cases predictably (error or NaN with reason).
- [x] Provide both signed neighbor-vector and shell-averaged summaries, including sample counts and uncertainty where appropriate.

### 14. Displacement correlation measurements
- [x] Calculate displacements relative to a **specified average/reference structure**, with fractional/Cartesian conventions explicit.
- [x] Support directional projections, vector covariances, normalized cross-correlations and cross-species/site analysis.
- [x] Define behavior for zero variance, unoccupied sites, supercell wrapping and non-orthogonal cells.
- [x] Where DISCUS reports analogous correlation metrics, compare reference outputs while documenting differing definitions.

## Suggested API (illustrative only)

```python
neighbors = structure.neighbors(shell=1, origin_site="A", target_site="B")
sro = structure.occupational_correlation(neighbors, species=("A", "B"))
disp = structure.displacement_correlation(neighbors, direction=(1, 0, 0))
```

## Suggested PRs

1. Neighbor-list model, periodic cell metric and selection logic.
2. Occupational pair probabilities and correlation metrics.
3. Directional displacement covariances and reference-structure handling.
4. Comparison notebook with tabular and chart views; links to DISCUS metrics.

## Acceptance criteria

- [x] Known cubic, tetragonal and triclinic toy cells produce independently checked neighbor distances and multiplicities.
- [x] Random occupational disorder approaches zero SRO within estimated uncertainty; ordered examples reproduce known sign/magnitude.
- [x] Constructed displacement fields reproduce analytic correlations along specified directions.
- [x] Correlation definitions, expected domains and zero-variance handling are documented and tested.
- [x] Tests pass across supported site/occupancy/PBC cases.

**Evidence:** [`tests/test_neighbors.py`](../tests/test_neighbors.py) (brute-force oracle, textbook shells, triclinic cells), [`tests/test_correlations.py`](../tests/test_correlations.py) (B2/L1₂ values, random-occupation standard errors, analytic displacement waves), [`benchmarks/discus_6_23_01/c_probes/RESULT.md`](../benchmarks/discus_6_23_01/c_probes/RESULT.md) (29/29 DISCUS 6.23.01 `chem corr` values to 5 × 10⁻¹³). Mixed PBC is not supported and raises `NotImplementedError`, as elsewhere in the package. Remaining scope gaps are listed under Package C in [TASKS.md](TASKS.md).

## Dependencies / cautions

Prefer Package A's tested `Structure` and lattice basis. Do not treat one Javelin correlation formula as interchangeable with all DISCUS definitions; explicitly map conventions before comparing or using as a fit target.
