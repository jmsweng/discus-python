/* Example models. Each is a complete builder state. */
(function (root) {
"use strict";
const h = 1.03 / Math.sqrt(3); // N–H 1.03 Å along <111>
const r3 = x => +x.toFixed(4);

function base() {
  return {
    title: "", seed: 1234, output_dir: "discus_model_out",
    cell: { a: 4, b: 4, c: 4, alpha: 90, beta: 90, gamma: 90 }, supercell: [10, 10, 10],
    molecules: [], orientations: [], sites: [],
    disorder: { atomic_displacements: {}, molecule_rigid: { translation_sigma: 0, libration_sigma_deg: 0 }, clash_threshold: 0.8,
      correlations: { enabled: false, sweeps: 50, kT: 1, match: "orientation", interactions: [] } },
    mmc: { enabled: false, neighborhoods: [], targets: [], moves: { swap: 1, shift: 0, swap_displacement: 0, pairing: "all", shift_sigma: {} },
      moves_per_atom: 10, feedback_per_atom: 1, temperature: 1 },
    calc: { launcher: "discus_suite -noautorun -macro", debye_waller: false, timeout_s: 86400,
      fourier: { enabled: true, steps: [0.1, 0.1, 0.1], hmax: [3, 3, 3], aver_percent: 100, delta_pdf: true, total: false, symmetry: "",
        lots: { enabled: false, shape: "box", size: [5, 5, 5], count: "all", periodic: true } },
      pdf: { enabled: true, radiation: "xray", qmin: 0.5, qmax: 7.5, dq: 0.01, qsig: 0, rmax: 20, dr: 0.02, wavelength: "" },
      preview: { u: [1, 0, 0], v: [0, 1, 0], range: 4 } },
  };
}

function nh4cl() {
  const s = base();
  s.title = "NH4Cl orientational disorder";
  s.cell = { a: 3.87, b: 3.87, c: 3.87, alpha: 90, beta: 90, gamma: 90 };
  s.supercell = [10, 10, 10];
  s.molecules = [{ id: "m1", name: "NH4+", anchor: 0, atoms: [
    { el: "N", xyz: [0, 0, 0] }, { el: "H", xyz: [r3(h), r3(h), r3(h)] }, { el: "H", xyz: [r3(h), r3(-h), r3(-h)] },
    { el: "H", xyz: [r3(-h), r3(h), r3(-h)] }, { el: "H", xyz: [r3(-h), r3(-h), r3(h)] }] }];
  s.orientations = [{ id: "g1", name: "Two tetrahedra", molecule: "m1", convention: "ZYZ", euler: [0, 0, 0],
    ops: "C4 [0 0 1]", frame: "lab", closure: true, weights: [1, 1] }];
  s.sites = [
    { label: "Cl", frac: [0, 0, 0], b_iso: 1, draw: "random", occupants: [{ type: "atom", species: "Cl", p: 1, b: "" }] },
    { label: "NH4", frac: [0.5, 0.5, 0.5], b_iso: 1, draw: "exact", occupants: [{ type: "molecule", group: "g1", p: 1, b: "" }] },
  ];
  s.disorder.molecule_rigid = { translation_sigma: 0.02, libration_sigma_deg: 3 };
  s.disorder.atomic_displacements = { Cl: 0.04 };
  s.disorder.correlations = { enabled: true, sweeps: 40, kT: 1, match: "orientation", interactions: [
    { site_a: 1, site_b: 1, t: [1, 0, 0], J: -0.6 }, { site_a: 1, site_b: 1, t: [0, 1, 0], J: -0.6 }, { site_a: 1, site_b: 1, t: [0, 0, 1], J: -0.6 }] };
  s.calc.fourier.steps = [0.1, 0.1, 0.1];
  s.calc.pdf.rmax = 18;
  s.calc.preview = { u: [1, 0, 0], v: [0, 1, 0], range: 4 };
  return s;
}

function cuau() {
  const s = base();
  s.title = "Cu3Au short-range order";
  s.cell = { a: 3.75, b: 3.75, c: 3.75, alpha: 90, beta: 90, gamma: 90 };
  s.supercell = [8, 8, 8];
  const occ = () => [{ type: "atom", species: "Cu", p: 0.75, b: "" }, { type: "atom", species: "Au", p: 0.25, b: "" }];
  s.sites = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]].map((f, i) => ({ label: `M${i + 1}`, frac: f, b_iso: 1, draw: "exact", occupants: occ() }));
  s.disorder.atomic_displacements = { Cu: 0.03, Au: 0.03 };
  // the 12 nearest neighbours of each fcc site, as DISCUS vectors between the four sites
  const nn = [];
  const pos = s.sites.map(x => x.frac);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let tx = -1; tx <= 1; tx++) for (let ty = -1; ty <= 1; ty++) for (let tz = -1; tz <= 1; tz++) {
    const d = [pos[j][0] + tx - pos[i][0], pos[j][1] + ty - pos[i][1], pos[j][2] + tz - pos[i][2]];
    if (Math.abs(Math.hypot(...d) - Math.SQRT1_2) < 1e-6) nn.push({ origin: i + 1, target: j + 1, t: [tx, ty, tz] });
  }
  s.mmc = { enabled: true, neighborhoods: [{ vectors: nn, symmetric: true }],
    targets: [{ type: "occupational", neighborhood: 0, a: "Cu", b: "Au", mode: "correlation", value: -0.15, direction: [1, 0, 0], distance: 2.65, k: 1, depth: 1 }],
    moves: { swap: 1, shift: 0, swap_displacement: 0, pairing: "all", shift_sigma: {} }, moves_per_atom: 20, feedback_per_atom: 1, temperature: 1 };
  s.calc.fourier.steps = [0.125, 0.125, 0.125];
  s.calc.fourier.hmax = [3, 3, 3];
  s.calc.pdf.rmax = 14;
  return s;
}

function rotor() {
  const s = base();
  s.title = "Linear rotor with vacancies";
  s.cell = { a: 5.6, b: 5.6, c: 5.6, alpha: 90, beta: 90, gamma: 90 };
  s.supercell = [8, 8, 8];
  s.molecules = [{ id: "m1", name: "CO2", anchor: "centroid", atoms: [{ el: "C", xyz: [0, 0, 0] }, { el: "O", xyz: [1.16, 0, 0] }, { el: "O", xyz: [-1.16, 0, 0] }] }];
  s.orientations = [{ id: "g1", name: "Along <111>", molecule: "m1", convention: "ZYZ", euler: [45, -35.2644, 0], ops: "C4 [0 0 1]\nC2 [1 1 0]", frame: "lab", closure: true, weights: [] },
    { id: "g2", name: "Along <100>", molecule: "m1", convention: "ZYZ", euler: [0, 0, 0], ops: "C3 [1 1 1]", frame: "lab", closure: true, weights: [] }];
  s.sites = [{ label: "Mol", frac: [0, 0, 0], b_iso: 1, draw: "random", occupants: [
    { type: "molecule", group: "g1", p: 0.6, b: "" }, { type: "molecule", group: "g2", p: 0.3, b: "" }, { type: "vacancy", p: 0.1 }] }];
  s.disorder.molecule_rigid = { translation_sigma: 0.03, libration_sigma_deg: 5 };
  s.disorder.correlations = { enabled: true, sweeps: 40, kT: 1, match: "orientation", interactions: [
    { site_a: 0, site_b: 0, t: [1, 0, 0], J: -0.8 }, { site_a: 0, site_b: 0, t: [0, 1, 0], J: -0.8 }, { site_a: 0, site_b: 0, t: [0, 0, 1], J: -0.8 }] };
  s.calc.pdf.rmax = 16;
  s.calc.fourier.steps = [0.125, 0.125, 0.125];
  s.calc.preview = { u: [1, 0, 0], v: [0, 1, 0], range: 3 };
  return s;
}

const EXAMPLES = { nh4cl: { name: "NH4Cl orientational order", make: nh4cl }, rotor: { name: "Linear rotor, two orientation families", make: rotor },
  cuau: { name: "Cu3Au short-range order (DISCUS mmc)", make: cuau }, blank: { name: "Empty cubic cell", make: () => { const s = base(); s.title = "New model";
    s.sites = [{ label: "A", frac: [0, 0, 0], b_iso: 1, draw: "random", occupants: [{ type: "atom", species: "Fe", p: 1, b: "" }] }]; return s; } } };
if (typeof module !== "undefined" && module.exports) module.exports = EXAMPLES; else root.DB_EXAMPLES = EXAMPLES;
})(typeof window !== "undefined" ? window : globalThis);
