# Package A — Python structure and interoperability

**Status:** TODO · **Priority:** Very high · **Suggested PR family:** `feat/structure-*`

**Reference in Javelin:** [structure.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/structure.py), [unitcell.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/unitcell.py), [io.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/io.py), [ASE tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/ase.rst), [DiffPy tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/diffpy.rst).

## Objective

Make structures first-class, editable Python values while DISCUS remains the primary scattering/disorder engine. Scientists should no longer need to hand-edit every `.cll` file to explore atomic arrangements. Preserve the semantics required by DISCUS and 3D reciprocal-space calculations.

## Planned features

### 1. Editable Python `Structure` model
- [ ] Represent cell lengths/angles or lattice matrix, species/atomic numbers, fractional and Cartesian positions, occupancy, site IDs, supercell shape, periodic boundaries and optional displacement/magnetic metadata.
- [ ] Provide shape/finite/element/occupancy validation and clear copying semantics.
- [ ] Avoid silently treating a periodic crystal as a non-periodic molecule, or an average CIF model as a fully occupied disordered supercell.

### 2. CIF, ASE and DiffPy interoperability
- [ ] Provide `from_cif`, `from_ase`, `to_ase` and optional DiffPy adaptors using maintained installed libraries rather than custom full CIF parsing.
- [ ] Preserve fractional coordinates, lattice conventions, site occupancies (including mixed sites) and PBC. Surface unavoidable loss of information explicitly.
- [ ] Install ASE/DiffPy as optional extras, not mandatory for existing 1D PDF calculations.

### 3. DISCUS structure import/export
- [ ] Add version/formats-aware readers and writers for relevant DISCUS `.cll` and `.stru` files; distinguish unit-cell files from full expanded supercells.
- [ ] Confirm format syntax and ordering from DISCUS documentation and real fixture files; Javelin's `read_stru` can inform format tests but must not be assumed authoritative for all DISCUS releases.
- [ ] Preserve site indexing, coordinate order and atom counts. Refuse unsupported information with a helpful error.

### 4. Supercell and atom manipulation
- [ ] `repeat((nx,ny,nz))`, species substitution, site selection, random and directional displacements, composition checks.
- [ ] Reproducible random seeds; atom ID mapping maintained through changes; do not silently shuffle crystallographic sites.
- [ ] Support round-trip `Structure → DISCUS input → Structure` on representative model files.

### 5. Coordinate and unit-cell utilities
- [ ] Fractional ↔ Cartesian transforms, direct/reciprocal metric tensors, `Q=2π B hkl` convention, d-spacings, distances under PBC.
- [ ] Validate general triclinic cells, handedness, singular/inconsistent lattices and units (Å and Å⁻¹).
- [ ] Provide one shared lattice implementation for later Package B, not two inconsistent transform routines.

## Proposed interface (illustrative only)

```python
structure = Structure.from_cif("model.cif")
supercell = structure.repeat((8, 8, 8))
supercell.substitute(site="A1", old="A", new="B", fraction=0.10, seed=1234)
supercell.to_discus("model.stru")
```

## Suggested PRs

1. Core lattice/structure dataclasses and transform/property tests.
2. ASE + CIF adapters and edge-case/occupancy tests.
3. DISCUS file I/O with versioned fixtures, supercell transforms and round trips.

## Acceptance criteria

- [ ] Standard cubic and non-orthogonal cells give correct Cartesian distances and reciprocal Q.
- [ ] CIF/ASE and supported DISCUS format round trips preserve cell, atom identity, occupancy, site mapping and positions within documented tolerance.
- [ ] Known unsupported occupancy/symmetry cases fail loudly rather than silently corrupting models.
- [ ] Existing `PDFConfig` and `PDFRunner` keep working.
- [ ] Tests, docs and runnable Jupyter example included.

## Dependencies / cautions

Do not redefine the validated DISCUS 1D macro. Agree on common `Structure` and reciprocal-basis conventions with the Claude 3D work before finalizing the public API. License any reused Javelin source appropriately.
