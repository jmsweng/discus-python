# Package C — Neighbor analysis and disorder correlations

**Status:** TODO · **Priority:** High · **Suggested PR family:** `feat/correlations-*`

**Reference in Javelin:** [neighborlist.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/neighborlist.py), [structure.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/structure.py#L570-L694), [neighbors tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/neighbors.rst).

## Objective

Give researchers transparent Python measurements of local occupational and displacement disorder, independent of scattering visualization. These quantitative descriptors can later serve as simulation controls and optimization targets.

## Planned features

### 12. Neighbor lists by sites, species, shells and geometry
- [ ] Represent neighbor pairs as atom/site IDs plus lattice translations; preserve relation to a supercell and boundary conditions.
- [ ] Select by source/target site, elemental species, exact lattice vector, Cartesian distance tolerance, first/second/etc. neighbor shell, and direction.
- [ ] Use the **full lattice metric** for non-orthogonal crystals. Implement periodic minimum-image and boundary handling explicitly; avoid counting the same pair twice unless intended.
- [ ] Add neighbor count/coordination number summaries and exportable tables.

### 13. Occupational correlation measurements
- [ ] Implement and document at least one rigorously defined measure: e.g. pair probabilities and Warren–Cowley short-range-order coefficients, with clear species and sublattice conventions.
- [ ] Support binary and, when mathematically defined, multicomponent occupancy; distinguish site-specific occupancies from global concentrations.
- [ ] Handle pure end members and zero-denominator cases predictably (error or NaN with reason).
- [ ] Provide both signed neighbor-vector and shell-averaged summaries, including sample counts and uncertainty where appropriate.

### 14. Displacement correlation measurements
- [ ] Calculate displacements relative to a **specified average/reference structure**, with fractional/Cartesian conventions explicit.
- [ ] Support directional projections, vector covariances, normalized cross-correlations and cross-species/site analysis.
- [ ] Define behavior for zero variance, unoccupied sites, supercell wrapping and non-orthogonal cells.
- [ ] Where DISCUS reports analogous correlation metrics, compare reference outputs while documenting differing definitions.

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

- [ ] Known cubic, tetragonal and triclinic toy cells produce independently checked neighbor distances and multiplicities.
- [ ] Random occupational disorder approaches zero SRO within estimated uncertainty; ordered examples reproduce known sign/magnitude.
- [ ] Constructed displacement fields reproduce analytic correlations along specified directions.
- [ ] Correlation definitions, expected domains and zero-variance handling are documented and tested.
- [ ] Tests pass across supported site/occupancy/PBC cases.

## Dependencies / cautions

Prefer Package A's tested `Structure` and lattice basis. Do not treat one Javelin correlation formula as interchangeable with all DISCUS definitions; explicitly map conventions before comparing or using as a fit target.
