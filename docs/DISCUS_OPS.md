# DISCUS structure operations (`discuspdf.discus_ops`)

Typed wrappers for DISCUS commands that change a crystal: random substitution, thermal
displacements, modulation waves, atom edits. A `StructureJob`:

1. reads a crystal into DISCUS: a `.cll` unit cell repeated, a `.stru` file, or a discuspdf
   `Structure`;
2. writes the 8-value seed. A job with a step that draws random numbers (`ReplaceRandom`,
   `Thermal`, a density wave or a random phase, `Raw(..., random=True)`) refuses to run
   without one, because DISCUS would seed from the clock. `allow_unseeded=True` overrides this;
3. applies the steps in order;
4. saves the crystal and reads it back as a `Structure`.

It runs through the general [`SuiteRunner`](SUITE_ENGINE.md), so the macro, inputs, engine and
hashes of every run are recorded.

```python
from discuspdf.discus_ops import StructureJob, ReplaceRandom, Thermal, Wave
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

job = StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)),
                   [ReplaceRandom("CU", "AU", 0.25),
                    Wave("longitudinal", (1, 0, 0), 21.69, amplitude=0.05),
                    Thermal()],
                   seed=Seed((11, 22, 33, 44, 55, 66, 77, 88)))
print(job.macro().render())
result = job.run(SuiteRunner("runs_ops", engine), {"cu_b05.cll": "tests/fixtures/discus_ops/cu_b05.cll"})
result.structure     # the modified crystal; result.run is the run record
```

## Steps

| Step | DISCUS | Notes (from the 6.23.01 source) |
|---|---|---|
| `ReplaceRandom(old, new, p, site="all")` | `replace <old>,<new>,{all\|<site>},<p>` | Each atom of type `old` (on that site) becomes `new` with probability `p`, using `ran1`. `p` must be in (0, 1]. A new type takes the old type's B (`modify_mod.f90`). |
| `ReplaceAtom(index, atom, b=None)` | `replace <index>,<atom>[,<b>]` | One atom, 1-based DISCUS order. `b` is needed if the type is new. |
| `Switch(i, j)` | `switch i,j` | Swaps the chemistry of two atoms; their positions stay. |
| `Remove(start, end=None)` | `remove …` | Turns atoms into vacancies. They stay in the list. |
| `Purge(types=False)` | `purge [type:yes]` | Drops vacancies from the list (see below). |
| `Thermal(two_d=False, molecules=False)` | `thermal [mol][,2d]` | Gaussian displacements from each type's ADPs: sigma = sqrt(B/8π²) per Cartesian direction, for isotropic B (`thermal.f90`). B < 0.001 is read as 1.0. |
| `Wave(kind, vector, length, …)` | `waves` menu … `run` | Longitudinal, transverse or density waves, with sinus, box or triangle profiles. Atom modes only (no molecules or rotation). |
| `Raw(lines)` | the lines unchanged | Escape hatch. Lines are still checked against the catalogue. |

`Wave` details, from `waves.f90`:

- **Reset and select:** the step starts with `reset`, which deselects every atom. It therefore
  always writes `select` (`all` by default).
- **Vectors:** the wave vector and the oscillation direction are normalised, so only their
  directions matter.
- **Profile:** `sinus` is `amplitude·cos(2π·arg) + shift`.
- **Density waves:** the profile f runs between `plow` and `phigh`. An atom of type `from` is
  replaced by `to` when `ran1 > f`, so with probability 1 − f, i.e. between 1 − `phigh` and
  1 − `plow`. With `viceversa`, atoms of type `to` become `from` with probability f.
  Amplitude, shift and branch do not apply.
- **Quotes:** DISCUS keeps quotes in parameters literally, so the step writes `box`, `random` and
  `viceversa` without quotes, unlike the help text.
- **Optical branch:** it reverses the displacement only of atoms whose name contains `-`.
- **`site:` option:** `select … site:` is parsed but never used by `run`, so it is not offered.

## Reading the result back

DISCUS `save` writes no `ncell` line, so the crystal size comes from the source: the cell
count for `read cell`, the `ncell` line of an input `.stru`, or `Structure.supercell`.

- **Vacancies** from `remove` stay in the saved file as `VOID` atoms, so the crystal is still a
  whole number of cells. DISCUS writes their Property flag as 0 (`do_remove` clears
  `PROP_NORMAL`), and that column is not kept: it is listed in `result.stru.dropped`.
- **Order check:** the saved atoms are checked for DISCUS's cell-major order, which `read cell`
  and `read stru` keep. A purge that happens to remove a whole number of cells' worth of atoms
  would otherwise be mislabelled silently. `run(..., check_order=False)` skips the check, e.g.
  after displacements of more than half a cell.
- **After `purge`** the atom count is usually no longer a multiple of the cell count. The job
  then raises `ResultNotAStructure`, whose `.run` still gives the saved file
  (`run.output("result.stru")`).

## Verification

These wrappers are at the *engine-ran* level of [the catalogue's verification levels](SUITE_CATALOGUE.md#verification-levels).
`tests/test_discus_ops_live.py` (opt-in, `DISCUS_LIVE=1`) checks each step on Cu fcc 6×6×6
(864 atoms):

- the ideal crystal is read and saved unchanged;
- the substituted fraction is within 4σ of `p`, and is the same with the same seed;
- the thermal ⟨u²⟩ is within 12 % of B/8π²;
- `replace`, `switch` and `remove` give the expected result;
- `purge` leaves 861 atoms;
- a longitudinal wave moves atoms only along its vector, with maximum displacement close to
  the amplitude.

No pre-registered physics check is planned for these steps: they generate models, they don't
compute scattering.

## Phase 2b: stacking faults, domains and the engine's correlation analysis

| API | DISCUS menu | Notes (source facts in the linked pages) |
|---|---|---|
| `discus_stack.Stacking(layers, correlation, translations, modulus, number, first, average, …)` | `stack` … `create`, `run` | A `CrystalBuilder`: the job builds the crystal instead of reading one. Layer types are `.stru` files. Row i of `correlation` is P(next layer type \| previous i). `create` makes the layer origins and `run` the atoms, so the result is read back as one box. For close packing, give the straight average `(0, 0, 1/3)`. [`DISCUS_STACK_FACTS.md`](DISCUS_STACK_FACTS.md) |
| `discus_domain.Domains(content, inputfile, name, character, shape, orient, fuzzy=2.55)` | `domain` (pseudo-atom mode) | Inserts content atoms inside a shape around each pseudo atom. Host atoms within `fuzzy` Å of an inserted atom become vacancies. The result is not whole unit cells, so use `StructureJob(..., readback="box")`. [`DISCUS_DOMAIN_FACTS.md`](DISCUS_DOMAIN_FACTS.md) |
| `discus_chem.ChemAnalysis(neighbourhoods, calculations)`, `run_chem(job, runner, inputs)` | `chem` … `corr occ` / `corr disp` | A neighbourhood is an `mc.Neighborhood` (the same vectors as `mmc`, numbered the same way) or a `DistanceNeighbourhood`. Returns DISCUS's own correlation values, parsed from the log. The parser is checked on the committed Package C probe logs. DISCUS 6.23.01 aborts at exit (code 134) after `corr occ`; `run_chem` keeps the printed results and sets `engine_crashed_at_exit`. |

`StructureJob(readback=…)`:

- `"cells"` (default): a `Structure` of whole unit cells.
- `"box"`: one periodic box, through `discus_ops.read_box`. For a `CrystalBuilder` the box size is DISCUS's own estimate from the atom extents.
- `"none"`: only the run.

Live checks: `tests/test_discus_stack_domain_chem_live.py` (opt-in, `DISCUS_LIVE=1`). It covers:

- an ABAB stack;
- fcc ABC without faults;
- a Markov fcc sequence with reversed steps, and that the same seed gives the same file;
- a sphere domain: inserted Ag count equal to the vacancy count;
- no removal with `fuzzy = 0`;
- seed independence of a domain with no size spread;
- `chem corr occ` on a random alloy.

## Phase 2b part 2: RMC fits to diffuse planes

| API | DISCUS | Notes ([`DISCUS_RMC_FACTS.md`](DISCUS_RMC_FACTS.md)) |
|---|---|---|
| `discus_rmc.PlaneGrid(ll, lr, ul, nx, ny)` | `fourier` `ll/lr/ul/na/no`, `rmc data` corners | Point (i, j) = ll + i·(lr−ll)/(nx−1) + j·(ul−ll)/(ny−1); `hkl()` gives every point. nx ≤ ny (stock RMC loops over nx rows). `from_reciprocal_grid` / `to_reciprocal_grid(lattice)` convert to and from a `ReciprocalGrid` of shape (nx, ny, 1); `FourierPlane` and `RmcPlane` take either. |
| `discus_rmc.FourierPlane(grid, output, radiation, wavelength, debye_waller)` | `fourier` … `run`; `output` `value int`, `form stan` | Writes a NIPL plane in the orientation `rmc data` reads. Use with `StructureJob(..., readback="none")`. |
| `discus_rmc.RmcFit(planes, moves, select, tries, sigma, …)` | `rmc` … `run` | One plane (DISCUS 6.23.01 mixes up the tables of several). Moves `RmcSwap`, `RmcShift` (a sigma for every selected type). Starts with `reset` and always writes `set mode`, `dese all`/`sele`. `set aver` > 0 and lots run on the fork. |
| `discus_rmc.run_rmc(job, runner, inputs)` | — | Checks each plane file against its grid, runs the job, parses every log block (`RmcLog.progress`, `.final`). |
| `discus_rmc.read_nipl`, `write_nipl` | NIPL files | z has shape (ny, nx). |

Acceptance is exp(−Δχ²/(σ²/2)) for uphill moves; Δχ² is in intensity² units, so σ has to be set
relative to the data (σ = 0: downhill only). `set cycl` counts valid tries in total, not per atom.

Live checks: `tests/test_discus_rmc_live.py` (opt-in, `DISCUS_LIVE=1`), Cu fcc 4×4×4 (256 atoms),
one 21 × 31 plane hk0:

- R1: the plane file's size and orientation (L1_2 superstructure peak at h = 1, fundamental at h = 2);
- R2, pre-registered: a swap-only fit with σ = 0 from a random Cu0.75Au0.25 start to the L1_2 plane
  accepts no uphill move, never raises χ², keeps the Au count and ends with a first-neighbour
  Warren–Cowley α < −0.15 (perfect L1_2: −1/3);
- R3: the same seed gives the same fitted file.

## Phase 2c part B: powder patterns and the real-space `pdf` menu

| API | DISCUS | Notes |
|---|---|---|
| `discus_powder.PowderPattern(outputs, qmin, qmax, dq, radiation, wavelength \| energy, stepr, temperature, polarisation)` | `powder` (`set calc, debye`) then `output` `form powd,q,…`, `value inte\|S(Q)\|F(Q)\|I(Q)` | Debye sum over the crystal as one cluster (no periodic images). Always writes `set profile, off` (reset restores a pseudo-Voigt), `set stepr`, `set temp`, `set disp`, `set period` (not restored by reset). Values at label Q are computed at Q − dq (`read_powder_file(..., shift_labels=True)`). qmax ≤ 18 (DISCUS rescales S above 18 Å⁻¹). Facts: [`DISCUS_POWDER_FACTS.md`](DISCUS_POWDER_FACTS.md). |
| `discus_pdfmenu.PdfCalc(output, rmax, dr, radiation, weight_s, qmax, periodic, exact, sphere, density, thermal, delta, gamma, qalpha)` | `pdf` … `calc`, `save pdf` | G(r) = H/(N dr_i r) − 4πrρ0 from the weighted pair histogram on DISCUS's internal 0.001-Å step; `save pdf` samples it every dr, so dr > 0.001 needs Gaussian peaks, `qalpha` or a Qmax. Weights use the Waasmaier table at s = sinθ/λ = `weight_s`. Always writes `set radiation` (without it DISCUS weights X-rays by the constant c of the form factor only), a three-parameter `set boundary` and `set delta/gamma/qalpha` (kept in a session). `set qsig` is not offered (without Gaussian peaks it zeroes G). Facts: [`DISCUS_PDFMENU_EXP2PDF_FACTS.md`](DISCUS_PDFMENU_EXP2PDF_FACTS.md). |
| `discus_pdfmenu.PdfPartials(model, pairs)`, `read_pdf_partials(step, run)` | `pdf` … `isele`/`jsele`, `calc`, `save pdf` | Partial G_ab(r) with `model`'s settings. A same-type partial is DISCUS's own; a cross partial selected directly counts half of the a–b pairs (probed, 1e-8), so the step selects {a, b} on both sides (file `…_A_B_union.dat`) and returns G_{a,b} − G_aa − G_bb. Selections are reset to all at the end. |
| `discus_pdfmenu.PdfData`, `PdfFit(model, data, moves, select, tries, …)`, `run_pdf_fit(job, runner, inputs)` | `pdf` … `data`, `run`, `save pdf` | RMC fit to G(r) data (r, G, 0, σ on the grid dr…rmax). Moves and settings as `RmcFit` (`RmcShift`/`RmcSwap`). One `calc` over rmax + 1 before `data` (after `reset all` or in a fresh session `data` fails: "Too many points in PDF"). `set frange` cannot lower the upper end of the window (pdf.f90:393-394 takes the larger rmax), so `PdfData.write(fit_rmax=…)` weights the top of the range out (σ 1e6); use rmax − 0.5 Å with Gaussian peaks (they lack partners beyond rmax). The step resets the RMC settings first and sets the scale back to fixed 1 after saving (a fitted scale would carry into later `calc`s). The saved fitted G(r) and `calc` of the fitted crystal differ by ~8e-4 of max |G| (two code paths). Log parsed by `parse_rmc_log`. |
| `discus_patterson.PattersonMap(intensity, grid, plane, output)`, `InverseMap(real, imag, grid, plane, output, imag_output)`, `RealPlane`, `run_map`, `read_map` | `patterson` / `inverse` … `run`, then `output` | Maps rho = (n_cell/N)/V_cell Σ A exp(+2πi h·r) of a NIPL plane from `FourierPlane` (`value="int"`, `"real"`, `"imag"`). `run_map` strips the `#` header lines the menus cannot read; the steps ask for one extra ordinate row, which DISCUS computes wrongly (phase not reset), and `read_map` drops it. `diff-four`, SHELXL/HKLF4 input and the `set` options are not offered. Facts: [`DISCUS_PATTERSON_FACTS.md`](DISCUS_PATTERSON_FACTS.md). |

Live checks (`tests/test_powder_pdfmenu_live.py`, pre-registered): Debye S(Q) and F(Q) of a 108-atom Cu
cluster (stepr 0.001) against a Python emulation (DISCUS's binning and sine table) within 1e-5 / 1e-4; `pdf` menu G(r)
for Cu (neutron) and Cu3Au (X-ray, f(0) weights) against a periodic Python pair histogram within 1e-6.
Probes for the two `pdf`-menu defects: `tools/probes/pdf_menu_probes.py`.

Not yet wrapped: the `pdf` menu fit (`data`, `run`), `exp2pdf`, partial PDFs.
