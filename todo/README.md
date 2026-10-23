# Javelin-inspired development roadmap

**Status: proposed / not started.** This directory is a planning backlog, not evidence that its features already work.

## Goal

Extend [discus-python](https://github.com/jmsweng/discus-python) with the most useful Python-native ideas from [Javelin](https://github.com/rosswhitfield/javelin) while keeping the installed DISCUS_SUITE as the primary, validated physics engine. Prioritize scientific correctness, provenance, inspectable macros, real-data validation, and notebook usability.

Javelin source reference: [rosswhitfield/javelin](https://github.com/rosswhitfield/javelin) (MIT license, older codebase). DISCUS already offers many of these **underlying physics** capabilities; the new value here is clean Python-facing access, interoperability, and validated workflows. Do not assume an item is absent from DISCUS simply because it is absent from the current Python wrapper.

## Work packages

| Package | Document | Javelin-inspired feature IDs | Focus | Suggested order |
|---|---|---|---|---|
| A | [Structure and interoperability](01_structure_and_interoperability.md) | 1–5 | Editable atomic structures, CIF/ASE/DiffPy, DISCUS I/O, coordinates | 1 |
| B | [Reciprocal-space and 3D data](02_reciprocal_space_and_data.md) | 6–9, 18–20 | Flexible Q grids, xarray, LOTS/Bragg controls, Mantid, persistence and VTK | 2, after verified 3D bridge |
| C | [Disorder and correlations](03_disorder_and_correlations.md) | 12–14 | Neighbor vectors, occupational and displacement correlations | 3 |
| D | [Monte Carlo and advanced physics](04_monte_carlo_and_advanced_physics.md) | 10–11, 15–17, 21 | Python-configured MC, target correlations, optional reference scattering and magnetism, visualization | 4 |

These four packages cover all 21 features in the initial Javelin comparison.

## Suggested development order

1. Establish validated DISCUS 3D inputs, outputs, and coordinate conventions before adding dependent analysis features.
2. Package A foundations: structure, coordinate/lattice operations, interchange and round-trip tests.
3. Package B foundations: xarray/grid models and validated native/experimental data adapters. Add 3D visualization/export after coordinates and masking are reliable.
4. Package C: neighbor and correlation measurements validated on analytical toy models.
5. Package D: small Python API for verified DISCUS Monte Carlo commands, then optional independent physics engines and additional viewers.

## Cross-cutting definition of done

- [ ] Public, documented Python API with units and conventions explicit
- [ ] Deterministic toy-model unit tests, edge cases, and round-trip tests where relevant
- [ ] Real DISCUS integration tests where DISCUS behavior is claimed
- [ ] Comparison to authoritative/independent numerical benchmarks for scientific algorithms
- [ ] Version-aware format/command handling; refuse unknown schemas rather than guessing
- [ ] Preserve and report cell metric, reciprocal basis, supercell indexing, atomic identity, occupancy, PBC, masks, and provenance where relevant
- [ ] Python 3.10+ and installed dependency compatibility checked
- [ ] Notebook example using the implemented API; no unverified features shown as working
- [ ] One focused PR, GitHub CI green, written agent test report, user approval before merge

## Reuse and license policy

Prefer wrapping DISCUS and using maintained third-party libraries (e.g. ASE, xarray, DiffPy where useful) rather than porting Javelin's aging Cython stack wholesale. Javelin is [MIT-licensed](https://github.com/rosswhitfield/javelin/blob/master/LICENSE); retain required license and attribution notices if actual source code is copied. Design ideas and separate clean implementations need clear scientific verification either way.
