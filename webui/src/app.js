/* Disorder Builder UI. State lives in S; every edit goes through commit(). */
(function () {
"use strict";
const DB = window.DB, EX = window.DB_EXAMPLES;
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => Array.from(el.querySelectorAll(s));
const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = (x, d = 3) => (isFinite(x) ? (+x).toFixed(d).replace(/\.?0+$/, "") : "—");
const STORE = "discuspdf-disorder-builder-v1";
const TEMPLATE = (document.getElementById("engine-template") || {}).textContent || "";

// ------------------------------------------------------------------------------------------ state
let S;
try { const raw = localStorage.getItem(STORE); if (raw) S = upgrade(JSON.parse(raw)); } catch (e) { S = null; }
if (!S) S = EX.rotor.make();
let section = "cell", molSel = S.molecules[0] && S.molecules[0].id, grpSel = S.orientations[0] && S.orientations[0].id;
let previewSeed = 0, placing = false, placeEl = "C", atomSel = -1, colorMode = "orientation", cellsShown = 3;
let lastReal = null;

function upgrade(s) {
  const b = EX.blank.make();
  const merge = (d, x) => { if (x == null) return d; if (typeof d !== "object" || Array.isArray(d)) return x; const o = { ...d }; for (const k in x) o[k] = merge(d[k], x[k]); return o; };
  return merge(b, s);
}
function save() { try { localStorage.setItem(STORE, JSON.stringify(S)); } catch (e) { /* storage unavailable: work in memory */ } }
function getPath(path) { return path.split(".").reduce((o, k) => (o == null ? o : o[k]), S); }
function setPath(path, v) {
  const ks = path.split("."); let o = S;
  for (let i = 0; i < ks.length - 1; i++) o = o[ks[i]];
  o[ks[ks.length - 1]] = v;
}
const uid = p => p + Math.random().toString(36).slice(2, 7);

// ------------------------------------------------------------------------------------------ commit / validation
let validation = { errs: [], warns: [] };
function commit({ render = true, preview = true } = {}) {
  save();
  try { validation = DB.validate(S); } catch (e) { validation = { errs: [e.message], warns: [] }; }
  renderStatus();
  if (render) renderSection();
  if (preview) schedulePreview();
  markDiffuseStale();
}
function renderStatus() {
  const e = validation.errs.length, w = validation.warns.length;
  $("#status").innerHTML = e ? `<span class="chip err"><span class="dot" style="background:currentColor"></span>${e} problem${e > 1 ? "s" : ""}</span>`
    : `<span class="chip ok"><span class="dot" style="background:currentColor"></span>Ready to export</span>`;
  if (w) $("#status").innerHTML += ` <span class="chip warn">${w} note${w > 1 ? "s" : ""}</span>`;
  const b = $("#badge-script"); if (b) b.textContent = e ? `${e}!` : "";
}

// ------------------------------------------------------------------------------------------ sections
const SECTIONS = [
  ["cell", "Cell"], ["molecules", "Molecules"], ["orientations", "Orientations"], ["sites", "Sites"],
  ["disorder", "Disorder"], ["mc", "Monte Carlo"], ["calc", "Calculations"], ["script", "Script"],
];
function renderRail() {
  $("#rail").innerHTML = SECTIONS.map(([k, t], i) => `<button type="button" data-sec="${k}" aria-current="${k === section}">
    <span class="k">${i + 1}</span><span class="t">${t}</span>${k === "script" ? '<span class="badge" id="badge-script"></span>' : ""}</button>`).join("");
  renderStatus();
}
function renderSection() {
  const fn = RENDER[section];
  const ed = $("#editor");
  const y = window.scrollY;
  ed.innerHTML = fn();
  if (AFTER[section]) AFTER[section]();
  window.scrollTo(0, y);
}
const num = (path, attrs = "") => `<input type="number" step="any" data-b="${path}" data-t="num" value="${esc(getPath(path))}" ${attrs}>`;
const txt = (path, attrs = "") => `<input type="text" data-b="${path}" data-t="str" value="${esc(getPath(path))}" ${attrs}>`;
const chk = (path, label, attrs = "") => `<label class="chk"><input type="checkbox" data-b="${path}" data-t="bool" ${getPath(path) ? "checked" : ""} ${attrs}> ${label}</label>`;
const sel = (path, opts, attrs = "") => `<select data-b="${path}" data-t="${typeof getPath(path) === "number" ? "num" : "str"}" ${attrs}>${opts.map(([v, l]) =>
  `<option value="${esc(v)}" ${String(getPath(path)) === String(v) ? "selected" : ""}>${esc(l)}</option>`).join("")}</select>`;
const trio = (path, attrs = "") => `<span class="trio">${[0, 1, 2].map(k => num(`${path}.${k}`, attrs)).join("")}</span>`;
const field = (label, inner) => `<label class="f"><span>${label}</span>${inner}</label>`;
const head = (t, p) => `<div class="sec-head"><h2>${t}</h2><p>${p}</p></div>`;
const msgs = (list, cls) => list.map(m => `<div class="msg ${cls}">${esc(m)}</div>`).join("");

function allSpecies() {
  const set = new Set();
  S.sites.forEach(s => s.occupants.forEach(o => { if (o.type === "atom" && o.species) set.add(o.species); }));
  S.sites.forEach(s => s.occupants.forEach(o => {
    if (o.type === "molecule") { const g = S.orientations.find(x => x.id === o.group); const m = g && S.molecules.find(x => x.id === g.molecule); if (m) m.atoms.forEach(a => set.add(a.el)); }
  }));
  return [...set];
}
function lattice() { try { return DB.latticeMatrix(S.cell); } catch (e) { return null; } }

const RENDER = {};
const AFTER = {};

// ---- 1 cell
RENDER.cell = () => {
  const A = lattice(), n = S.supercell.map(Number);
  const V = A ? DB.det(A) : NaN;
  return head("Unit cell and crystal", "The cell uses the discus-python convention: a along x, b in the xy plane. All positions are fractions of this cell; the crystal is the cell repeated n₁ × n₂ × n₃ times.") + `
  <div class="block"><h3>Model</h3><div class="fields">
    ${field("Title", txt("title", 'class="w-xl"'))}${field("Random seed", num("seed", 'step="1"'))}${field("Output folder", txt("output_dir", 'class="w-l"'))}</div></div>
  <div class="block"><h3>Lattice parameters</h3><div class="fields">
    ${["a", "b", "c"].map(k => field(`${k} (Å)`, num("cell." + k))).join("")}
    ${[["alpha", "α"], ["beta", "β"], ["gamma", "γ"]].map(([k, l]) => field(`${l} (°)`, num("cell." + k))).join("")}</div>
    <div class="fields">${field("Supercell n₁ n₂ n₃", trio("supercell", 'step="1" min="1" class="w-s"'))}
      <span class="btn-row">${[6, 8, 10, 12, 16].map(k => `<button type="button" class="btn small" data-a="supercell" data-v="${k}">${k}³</button>`).join(" ")}</span></div>
    ${A ? `<dl class="facts"><dt>Cell volume</dt><dd>${fmt(V, 3)} Å³</dd><dt>Crystal</dt><dd>${n.join(" × ")} = ${n[0] * n[1] * n[2]} cells, ${fmt(n[0] * DB.norm(A[0]), 1)} × ${fmt(n[1] * DB.norm(A[1]), 1)} × ${fmt(n[2] * DB.norm(A[2]), 1)} Å</dd>
      <dt>a, b, c (Cartesian)</dt><dd>${A.map(r => "(" + r.map(x => fmt(x, 4)).join(", ") + ")").join("  ")}</dd></dl>` : `<div class="msg err">These parameters do not form a valid cell.</div>`}
  </div>
  <div class="block"><h3>Start from an example</h3><p class="hint">Replaces the current model. Export it first if you want to keep it.</p>
    <div class="fields">${Object.entries(EX).map(([k, e]) => `<button type="button" class="btn" data-a="example" data-v="${k}">${esc(e.name)}</button>`).join("")}</div></div>`;
};

// ---- 2 molecules
const TEMPLATES = {
  "H2O": [["O", 0, 0, 0], ["H", .7572, .5865, 0], ["H", -.7572, .5865, 0]],
  "CO2": [["C", 0, 0, 0], ["O", 1.16, 0, 0], ["O", -1.16, 0, 0]],
  "N2": [["N", .549, 0, 0], ["N", -.549, 0, 0]],
  "C2H2": [["C", .601, 0, 0], ["C", -.601, 0, 0], ["H", 1.663, 0, 0], ["H", -1.663, 0, 0]],
  "NH4+": tetra("N", "H", 1.03), "CH4": tetra("C", "H", 1.09), "ClO4-": tetra("Cl", "O", 1.44), "SO4 2-": tetra("S", "O", 1.47),
  "NO3-": [["N", 0, 0, 0], ...[0, 120, 240].map(t => ["O", 1.24 * Math.cos(t * Math.PI / 180), 1.24 * Math.sin(t * Math.PI / 180), 0])],
  "SF6": [["S", 0, 0, 0], ["F", 1.56, 0, 0], ["F", -1.56, 0, 0], ["F", 0, 1.56, 0], ["F", 0, -1.56, 0], ["F", 0, 0, 1.56], ["F", 0, 0, -1.56]],
  "C6H6": [...[0, 60, 120, 180, 240, 300].map(t => ["C", 1.39 * Math.cos(t * Math.PI / 180), 1.39 * Math.sin(t * Math.PI / 180), 0]),
    ...[0, 60, 120, 180, 240, 300].map(t => ["H", 2.48 * Math.cos(t * Math.PI / 180), 2.48 * Math.sin(t * Math.PI / 180), 0])],
  "CH3CN": [["C", 0, 0, -1.47], ["C", 0, 0, 0], ["N", 0, 0, 1.157], ...[0, 120, 240].map(t => ["H", 1.03 * Math.cos(t * Math.PI / 180), 1.03 * Math.sin(t * Math.PI / 180), -1.83])],
};
function tetra(c, x, d) { const h = d / Math.sqrt(3); return [[c, 0, 0, 0], [x, h, h, h], [x, h, -h, -h], [x, -h, h, -h], [x, -h, -h, h]]; }

RENDER.molecules = () => {
  const m = S.molecules.find(x => x.id === molSel);
  const list = S.molecules.map(x => `<button type="button" class="btn" aria-pressed="${x.id === molSel}" data-a="mol-select" data-v="${x.id}">${esc(x.name)} <span class="tag">${x.atoms.length}</span></button>`).join(" ");
  let body = `<p class="hint">No molecule yet. Start from a template or an empty molecule.</p>`;
  if (m) {
    const mi = S.molecules.indexOf(m), loc = DB.moleculeLocal(m);
    const bonds = DB.bondsOf(loc.species, loc.xyz);
    body = `<div class="fields">${field("Name", txt(`molecules.${mi}.name`, 'class="w-l"'))}
      ${field("Placed at the site by", sel(`molecules.${mi}.anchor`, [["centroid", "centroid"], ["com", "centre of mass"], ...m.atoms.map((a, i) => [i, `atom ${i + 1} (${a.el})`])]))}
      ${field("Template", `<select id="tpl"><option value="">Replace atoms with…</option>${Object.keys(TEMPLATES).map(k => `<option>${k}</option>`).join("")}</select>`)}
      <button type="button" class="btn" data-a="mol-dup">Duplicate</button><button type="button" class="btn ghost" data-a="mol-del">Delete molecule</button></div>
    <div class="mol-layout">
      <div class="tbl-wrap"><table class="tbl"><thead><tr><th>#</th><th>Element</th><th>x (Å)</th><th>y (Å)</th><th>z (Å)</th><th></th></tr></thead><tbody>
      ${m.atoms.map((a, i) => `<tr class="${i === atomSel ? "sel" : ""}"><td class="num"><button type="button" class="btn ghost small" data-a="atom-pick" data-v="${i}">${i + 1}</button></td>
        <td><input type="text" class="w-s" data-b="molecules.${mi}.atoms.${i}.el" data-t="str" value="${esc(a.el)}"></td>
        ${[0, 1, 2].map(k => `<td>${num(`molecules.${mi}.atoms.${i}.xyz.${k}`)}</td>`).join("")}
        <td><button type="button" class="btn ghost" data-a="atom-del" data-v="${i}" aria-label="Remove atom ${i + 1}">✕</button></td></tr>`).join("")}
      </tbody></table>
      <div class="fields" style="margin-top:10px"><button type="button" class="btn" data-a="atom-add">Add atom</button>
        <button type="button" class="btn" data-a="mol-centre">Move anchor to origin</button>
        <span class="hint">${bonds.length} bonds drawn from covalent radii. Coordinates are Cartesian, in the cell's frame (a ∥ x).</span></div>
      <details style="margin-top:10px"><summary>Paste XYZ coordinates</summary>
        <textarea id="xyz-in" rows="5" placeholder="C 0.0 0.0 0.0&#10;O 1.16 0 0"></textarea>
        <button type="button" class="btn small" data-a="xyz-import">Replace atoms</button></details></div>
      <div><div class="canvas-box ${placing ? "placing" : ""}" id="molbox"><canvas id="molcv"></canvas>
        <div class="overlay">${placing ? "click to place " + esc(placeEl) : "drag to rotate · shift-drag to pan · wheel to zoom"}</div></div>
        <div class="fields" style="margin-top:8px"><button type="button" class="btn small" aria-pressed="${placing}" data-a="place-toggle">Draw atoms</button>
          ${field("Element", `<input type="text" id="place-el" class="w-s" value="${esc(placeEl)}">`)}
          <button type="button" class="btn small" data-a="mol-view" data-v="z">View ⟂ xy</button><button type="button" class="btn small" data-a="mol-view" data-v="x">View ⟂ yz</button></div>
        <p class="hint">Draw atoms places a new atom where you click, in the plane of the screen through the anchor.</p></div>
    </div>`;
  }
  return head("Molecules", "Build rigid molecules from atoms. A molecule is placed on a site by its anchor point, then turned into an orientation from its library.") +
    `<div class="block"><h3>Molecules <span class="grow"></span><button type="button" class="btn small" data-a="mol-add">New molecule</button></h3>
      <div class="fields">${list || '<span class="hint">None</span>'}</div>${body}</div>`;
};
let molViewer = null;
AFTER.molecules = () => {
  const cv = $("#molcv"); if (!cv) return;
  molViewer = new Viewer(cv, { aspect: 1, placing: () => placing, onClick: (p) => {
    const m = S.molecules.find(x => x.id === molSel); if (!m) return;
    const el = ($("#place-el").value || "C").trim(); placeEl = el;
    const loc = DB.moleculeLocal(m), o = DB.sub(m.atoms.length ? m.atoms[0].xyz.map(Number) : [0, 0, 0], loc.xyz[0] || [0, 0, 0]);
    const xyz = DB.add(p, o).map(x => +x.toFixed(4));
    m.atoms.push({ el, xyz }); atomSel = m.atoms.length - 1; commit();
  } });
  if (window._molR) molViewer.R = window._molR;
  drawMolecule();
  cv.addEventListener("pointerup", () => { window._molR = molViewer.R; });
  const tpl = $("#tpl"); if (tpl) tpl.addEventListener("change", () => {
    const m = S.molecules.find(x => x.id === molSel); const t = TEMPLATES[tpl.value]; if (!m || !t) return;
    m.atoms = t.map(([el, x, y, z]) => ({ el, xyz: [x, y, z].map(v => +v.toFixed(4)) })); m.anchor = "centroid";
    if (/^Mol|^New/.test(m.name)) m.name = tpl.value;
    commit();
  });
};
function drawMolecule() {
  const m = S.molecules.find(x => x.id === molSel); if (!m || !molViewer) return;
  const loc = DB.moleculeLocal(m);
  molViewer.set({ atoms: loc.xyz.map((p, i) => ({ p, color: DB.colorOf(loc.species[i]), r: 0.18 + 0.2 * DB.covOf(loc.species[i]), label: String(i + 1), ring: i === atomSel ? "#0d6a80" : null })),
    bonds: DB.bondsOf(loc.species, loc.xyz), lines: [], center: [0, 0, 0], axes: lattice() || undefined });
}

// ---- 3 orientations
RENDER.orientations = () => {
  const g = S.orientations.find(x => x.id === grpSel);
  const list = S.orientations.map(x => `<button type="button" class="btn" aria-pressed="${x.id === grpSel}" data-a="grp-select" data-v="${x.id}">${esc(x.name)}</button>`).join(" ");
  let body = `<p class="hint">${S.molecules.length ? "Create a group to start a library." : "Build a molecule first."}</p>`;
  if (g) {
    const gi = S.orientations.indexOf(g), A = lattice();
    const mol = S.molecules.find(x => x.id === g.molecule);
    let ex = { orientations: [], errors: ["choose a molecule"], nOps: 0 };
    if (mol && A) ex = DB.expandGroup(g, mol, A);
    const w = k => (g.weights && g.weights[k] != null ? g.weights[k] : 1);
    const wsum = ex.orientations.reduce((s, o, k) => s + (+w(k) || 0), 0) || 1;
    body = `<div class="fields">${field("Name", txt(`orientations.${gi}.name`, 'class="w-l"'))}
        ${field("Molecule", sel(`orientations.${gi}.molecule`, S.molecules.map(m => [m.id, m.name])))}
        <button type="button" class="btn" data-a="grp-dup">Duplicate</button><button type="button" class="btn ghost" data-a="grp-del">Delete group</button></div>
      <div class="fields">${field("Euler convention", sel(`orientations.${gi}.convention`, [["ZYZ", "ZYZ"], ["ZXZ", "ZXZ"], ["XYZ", "XYZ"]]))}
        ${field("α β γ (°) — base orientation", trio(`orientations.${gi}.euler`))}
        ${field("Operations act in", sel(`orientations.${gi}.frame`, [["lab", "crystal frame (after the base)"], ["molecule", "molecule frame (before the base)"]]))}
        ${chk(`orientations.${gi}.closure`, "Generate the full group")}</div>
      <div class="fields" style="align-items:start"><label class="f" style="flex:1;min-width:240px"><span>Symmetry operations, one per line</span>
        <textarea rows="5" data-b="orientations.${gi}.ops" data-t="str" spellcheck="false">${esc(g.ops)}</textarea></label>
        <div class="doc hint" style="max-width:44ch"><p><code>C4 [0 0 1]</code> n-fold axis along a lattice direction, <code>&lt;x y z&gt;</code> for Cartesian; <code>S4 [001]</code>; <code>m [1 -1 0]</code> mirror normal; <code>-1</code> inversion; <code>R 30 [1 1 1]</code> any angle; <code>-y,x,z</code> Jones symbol (translations ignored); <code>M</code> + nine numbers for a Cartesian matrix.</p>
        <div class="fields">${["C2 [0 0 1]", "C3 [1 1 1]", "C4 [0 0 1]", "m [0 0 1]", "-1"].map(o => `<button type="button" class="btn small" data-a="op-add" data-v="${o}">${o}</button>`).join("")}</div></div></div>
      ${msgs(ex.errors, "err")}
      <p class="hint">${ex.orientations.length} distinct orientation${ex.orientations.length === 1 ? "" : "s"} from ${ex.nOps} operation${ex.nOps === 1 ? "" : "s"}${ex.nOps > ex.orientations.length ? " (operations that map the molecule onto itself are merged)" : ""}. Relative weights set how often each orientation is drawn when a site picks this group. All thumbnails share one tilted view of the cell axes.</p>
      <div class="orient-grid">${ex.orientations.map((o, k) => `<div class="orient ${+w(k) > 0 ? "" : "off"}"><canvas data-orient="${k}"></canvas>
        <div class="row"><span class="tag">#${k + 1}</span>${o.improper ? '<span class="tag improper" title="contains a mirror: the molecule is inverted">mirror</span>' : ""}${o.multiplicity > 1 ? `<span class="tag" title="operations giving this placement">×${o.multiplicity}</span>` : ""}</div>
        <div class="row"><span>weight</span><input type="number" step="any" min="0" data-b="orientations.${gi}.weights.${k}" data-t="num" value="${esc(w(k))}"></div>
        <div class="row"><span class="tag">${fmt(100 * (+w(k) || 0) / wsum, 1)}%</span><button type="button" class="btn ghost small" data-a="mat-show" data-v="${k}">matrix</button></div></div>`).join("")}</div>
      <pre class="code" id="matbox" hidden></pre>`;
  }
  return head("Orientation library", "A group is a base orientation (Euler angles) plus symmetry operations. Every operation applied to the base gives one preferred orientation; sites then choose among a group's orientations with the weights below.") +
    `<div class="block"><h3>Groups <span class="grow"></span><button type="button" class="btn small" data-a="grp-add">New group</button></h3><div class="fields">${list || '<span class="hint">None</span>'}</div>${body}</div>`;
};
AFTER.orientations = () => {
  const g = S.orientations.find(x => x.id === grpSel), A = lattice();
  const mol = g && S.molecules.find(x => x.id === g.molecule);
  if (!g || !mol || !A) return;
  const ex = DB.expandGroup(g, mol, A), loc = DB.moleculeLocal(mol), bonds = DB.bondsOf(loc.species, loc.xyz);
  $$("canvas[data-orient]").forEach(cv => {
    const o = ex.orientations[+cv.dataset.orient];
    const v = new Viewer(cv, { static: true, aspect: 1 });
    v.viewAlong(DB.add(DB.add(DB.scale(A[0], 0.42), DB.scale(A[1], 0.28)), A[2]), A[2]);
    v.set({ atoms: loc.xyz.map((p, i) => ({ p: DB.mv(o.R, p), color: DB.colorOf(loc.species[i]), r: 0.16 + 0.18 * DB.covOf(loc.species[i]) })), bonds, lines: [], center: [0, 0, 0], axes: A });
  });
};

// ---- 4 sites
const OCC_COLORS = ["#0d6a80", "#d9822b", "#6b8e23", "#8e44ad", "#c0392b", "#2c7fb8", "#b8860b", "#7f8c8d"];
RENDER.sites = () => {
  const groups = S.orientations.map(g => [g.id, g.name]);
  const cards = S.sites.map((s, si) => {
    const tot = s.occupants.reduce((a, o) => a + (+o.p || 0), 0);
    const occ = s.occupants.map((o, oi) => {
      const p = `sites.${si}.occupants.${oi}`;
      const what = o.type === "atom" ? `<input type="text" class="w-s" data-b="${p}.species" data-t="str" value="${esc(o.species)}" aria-label="species">`
        : o.type === "molecule" ? (groups.length ? sel(`${p}.group`, groups) : '<span class="hint">no orientation groups yet</span>') : '<span class="hint">empty site</span>';
      return `<tr><td><span class="dot" style="background:${OCC_COLORS[oi % 8]}"></span></td><td>${sel(`${p}.type`, [["atom", "atom"], ["molecule", "molecule"], ["vacancy", "vacancy"]])}</td><td>${what}</td>
        <td>${num(`${p}.p`, 'min="0" max="1"')}</td><td>${o.type !== "vacancy" ? `<input type="number" step="any" data-b="${p}.b" data-t="numornull" value="${esc(o.b)}" placeholder="site" class="w-s">` : ""}</td>
        <td><button type="button" class="btn ghost" data-a="occ-del" data-v="${si}.${oi}" aria-label="Remove occupant">✕</button></td></tr>`;
    }).join("");
    return `<div class="site"><div class="site-head">${field("Label", `<input type="text" class="w-s" data-b="sites.${si}.label" data-t="str" value="${esc(s.label)}">`)}
      ${field("x y z (fractional)", trio(`sites.${si}.frac`))}${field("B (Å²)", num(`sites.${si}.b_iso`, 'class="w-s"'))}
      ${field("Draw", sel(`sites.${si}.draw`, [["random", "independent per cell"], ["exact", "exact composition"]]))}
      <span class="grow" style="flex:1"></span><button type="button" class="btn ghost" data-a="site-del" data-v="${si}">Delete site</button></div>
      <div class="probbar">${s.occupants.map((o, oi) => `<span style="width:${100 * Math.max(0, +o.p || 0) / Math.max(tot, 1)}%;background:${OCC_COLORS[oi % 8]}"></span>`).join("")}</div>
      <div class="tbl-wrap"><table class="tbl"><thead><tr><th></th><th>Occupant</th><th>Species / orientation group</th><th>Probability</th><th>B override</th><th></th></tr></thead><tbody>${occ}</tbody></table></div>
      <div class="fields"><button type="button" class="btn small" data-a="occ-add" data-v="${si}.atom">+ Atom</button><button type="button" class="btn small" data-a="occ-add" data-v="${si}.molecule">+ Molecule</button>
        <button type="button" class="btn small" data-a="occ-add" data-v="${si}.vacancy">+ Vacancy</button><button type="button" class="btn small" data-a="occ-norm" data-v="${si}">Scale to 1</button>
        <span class="sum ${Math.abs(tot - 1) > 1e-6 ? "bad" : ""}">Σp = ${fmt(tot, 4)}</span></div></div>`;
  }).join("");
  return head("Sites", "Each site sits at one fractional position in every unit cell. In each cell it holds one occupant, drawn with these probabilities: an atom, a molecule in one of its group's orientations, or nothing.") +
    `<div class="block"><h3>Sites in the unit cell <span class="grow"></span><button type="button" class="btn small" data-a="site-add">Add site</button></h3>
    <div class="cards">${cards}</div>
    <p class="hint">“Exact composition” fixes the number of each occupant in the crystal (shuffled); “independent” draws every cell separately. B is used only with Debye–Waller switched on (Calculations).</p></div>`;
};

// ---- 5 disorder
RENDER.disorder = () => {
  const sp = allSpecies(), d = S.disorder, c = d.correlations;
  const siteOpts = S.sites.map((s, i) => [i, s.label]);
  const rep = lastReal && lastReal.report;
  const inter = c.interactions.map((it, k) => {
    const r = rep && rep.find(x => x.it === undefined ? false : x.it.site_a === +it.site_a && x.it.site_b === +it.site_b && String(x.it.t) === String(it.t.map(Number)) && +x.it.J === +it.J);
    return `<tr><td>${sel(`disorder.correlations.interactions.${k}.site_a`, siteOpts)}</td><td>${sel(`disorder.correlations.interactions.${k}.site_b`, siteOpts)}</td>
      <td>${trio(`disorder.correlations.interactions.${k}.t`, 'step="1" class="w-s"')}</td><td>${num(`disorder.correlations.interactions.${k}.J`, 'class="w-s"')}</td>
      <td class="num">${r ? `${fmt(r.before, 2)} → ${fmt(r.after, 2)}` : "—"}</td>
      <td><button type="button" class="btn ghost" data-a="int-del" data-v="${k}" aria-label="Remove">✕</button></td></tr>`;
  }).join("");
  return head("Disorder", "Occupancies come from the Sites step. Here you add correlations between neighbouring occupants and the displacive disorder applied after the occupants are placed.") + `
  <div class="block"><h3>Correlated occupants ${chk("disorder.correlations.enabled", "on")}</h3>
    <p class="hint">A swap Monte Carlo that keeps the composition. Each listed pair of neighbours adds an energy J (in kT) when both hold the same occupant: J &lt; 0 makes like neighbours (domains), J &gt; 0 alternation. It works for atoms, vacancies and molecule orientations alike and runs in the generated script before DISCUS sees the crystal.</p>
    <div class="fields">${field("Sweeps", num("disorder.correlations.sweeps", 'step="1" class="w-s"'))}${field("kT", num("disorder.correlations.kT", 'class="w-s"'))}
      ${field("“Same occupant” means", sel("disorder.correlations.match", [["orientation", "same atom / same orientation"], ["group", "same atom / same molecule group"]]))}</div>
    <div class="tbl-wrap"><table class="tbl"><thead><tr><th>Site</th><th>Neighbour site</th><th>Cell offset t</th><th>J (kT)</th><th>P(same) preview</th><th></th></tr></thead><tbody>${inter}</tbody></table></div>
    <div class="fields"><button type="button" class="btn small" data-a="int-add">Add pair</button>
      ${S.sites.length ? `<button type="button" class="btn small" data-a="int-axes">Add ±a, ±b, ±c neighbours of ${esc(S.sites[0].label)}</button>` : ""}</div></div>
  <div class="block"><h3>Rigid-body disorder of molecules</h3><div class="fields">
    ${field("Translation σ (Å)", num("disorder.molecule_rigid.translation_sigma"))}${field("Libration σ (°)", num("disorder.molecule_rigid.libration_sigma_deg"))}</div>
    <p class="hint">Each molecule gets a Gaussian shift and a small rotation about a random axis on top of its preferred orientation. These are static (frozen) displacements, so they appear as diffuse scattering rather than as a Debye–Waller factor.</p></div>
  <div class="block"><h3>Atomic displacements</h3><div class="fields">
    ${sp.length ? sp.map(s => field(`σ ${esc(s)} (Å)`, `<input type="number" step="any" min="0" data-b="disorder.atomic_displacements.${esc(s)}" data-t="num" value="${esc(d.atomic_displacements[s] || 0)}">`)).join("") : '<span class="hint">No species yet.</span>'}</div>
    <p class="hint">Independent Gaussian displacements per atom (all three directions). Applied to molecule atoms too, which then distorts the molecules slightly.</p></div>
  <div class="block"><h3>Contact check</h3><div class="fields">${field("Warn below (Å)", num("disorder.clash_threshold", 'class="w-s"'))}</div>
    <p class="hint">The script reports atoms of different molecules or sites that come closer than this.</p></div>`;
};

// ---- 6 Monte Carlo (DISCUS mmc)
RENDER.mc = () => {
  const m = S.mmc, hasMol = S.sites.some(s => s.occupants.some(o => o.type === "molecule" && +o.p > 0));
  const siteOpts = S.sites.map((s, i) => [i + 1, `${i + 1} · ${s.label}`]);
  const sp = allSpecies().concat(S.sites.some(s => s.occupants.some(o => o.type === "vacancy")) ? ["VOID"] : []);
  const spOpts = sp.map(s => [s, s]);
  const nbs = m.neighborhoods.map((nb, k) => `<div class="site"><div class="site-head"><strong>Neighbourhood ${k + 1}</strong> ${chk(`mmc.neighborhoods.${k}.symmetric`, "add reversed vectors")}
      <span style="flex:1"></span><button type="button" class="btn ghost" data-a="nb-del" data-v="${k}">Delete</button></div>
    <div class="tbl-wrap"><table class="tbl"><thead><tr><th>From site</th><th>To site</th><th>Cell offset</th><th></th></tr></thead><tbody>
    ${nb.vectors.map((v, j) => `<tr><td>${sel(`mmc.neighborhoods.${k}.vectors.${j}.origin`, siteOpts)}</td><td>${sel(`mmc.neighborhoods.${k}.vectors.${j}.target`, siteOpts)}</td>
      <td>${trio(`mmc.neighborhoods.${k}.vectors.${j}.t`, 'step="1" class="w-s"')}</td><td><button type="button" class="btn ghost" data-a="vec-del" data-v="${k}.${j}" aria-label="Remove">✕</button></td></tr>`).join("")}</tbody></table></div>
    <div class="fields"><button type="button" class="btn small" data-a="vec-add" data-v="${k}">Add vector</button><button type="button" class="btn small" data-a="vec-axes" data-v="${k}">Add ±a, ±b, ±c of site 1</button></div></div>`).join("");
  const tg = m.targets.map((t, k) => {
    const p = `mmc.targets.${k}`;
    const extra = t.type === "displacement" ? field("Direction (Cartesian)", trio(`${p}.direction`))
      : t.type === "spring" ? field("Distance (Å)", num(`${p}.distance`, 'class="w-s"')) + field("k", num(`${p}.k`, 'class="w-s"'))
      : t.type === "lennard_jones" ? field("Distance (Å)", num(`${p}.distance`, 'class="w-s"')) + field("Depth", num(`${p}.depth`, 'class="w-s"')) : "";
    const valued = t.type === "occupational" || t.type === "displacement";
    return `<div class="fields">${field("Type", sel(`${p}.type`, [["occupational", "occupational correlation"], ["displacement", "displacement correlation"], ["spring", "spring (size effect)"], ["lennard_jones", "Lennard-Jones"]]))}
      ${field("Neighbourhood", sel(`${p}.neighborhood`, m.neighborhoods.map((x, i) => [i, String(i + 1)])))}${field("A", sel(`${p}.a`, spOpts))}${field("B", sel(`${p}.b`, spOpts))}
      ${valued ? field("Given as", sel(`${p}.mode`, [["correlation", "target correlation"], ["energy", "fixed energy (kT)"]])) + field("Value", num(`${p}.value`, 'class="w-s"')) : ""}
      ${extra}<button type="button" class="btn ghost" data-a="tg-del" data-v="${k}">Remove</button></div>`;
  }).join("");
  return head("DISCUS Monte Carlo", "Optional DISCUS mmc run on the crystal before the scattering calculation, through discus-python's MCConfig. It adds chemical short-range order, correlated displacements or size effects.") +
    (hasMol ? `<div class="msg warn">This model contains molecules. DISCUS mmc moves single atoms, so it is only offered for atom-only models. Use “Correlated occupants” under Disorder for molecule orientations.</div>` : "") + `
  <div class="block"><h3>${chk("mmc.enabled", "Run DISCUS mmc", hasMol ? "disabled" : "")}</h3>
    <p class="hint">Sites are numbered as in the Sites step (DISCUS site = position in the cell). Acceptance is heat-bath, not Metropolis. A vacancy is the atom type VOID.</p>
    <div class="fields">${field("Moves per atom", num("mmc.moves_per_atom", 'class="w-s"'))}${field("Feedback per atom", num("mmc.feedback_per_atom", 'class="w-s"'))}${field("kT", num("mmc.temperature", 'class="w-s"'))}</div></div>
  <div class="block"><h3>Neighbourhoods <span class="grow"></span><button type="button" class="btn small" data-a="nb-add">Add neighbourhood</button></h3><div class="cards">${nbs || '<p class="hint">None yet.</p>'}</div></div>
  <div class="block"><h3>Targets <span class="grow"></span><button type="button" class="btn small" data-a="tg-add">Add target</button></h3>${tg || '<p class="hint">None yet.</p>'}
    <p class="hint">At most one target may use a correlation (DISCUS feedback); give the others fixed energies. Each neighbourhood holds one occupational and one displacement target.</p></div>
  <div class="block"><h3>Moves</h3><div class="fields">${field("Swap probability", num("mmc.moves.swap", 'class="w-s"'))}${field("Pairing", sel("mmc.moves.pairing", [["all", "all"], ["local", "local"], ["slocal", "slocal"], ["site", "site"]]))}
    ${field("Shift probability", num("mmc.moves.shift", 'class="w-s"'))}${field("Swap + displacement", num("mmc.moves.swap_displacement", 'class="w-s"'))}</div>
    <div class="fields">${sp.filter(s => s !== "VOID").map(s => field(`Shift σ ${esc(s)} (Å)`, `<input type="number" step="any" min="0" data-b="mmc.moves.shift_sigma.${esc(s)}" data-t="num" value="${esc(m.moves.shift_sigma[s] || 0)}">`)).join("")}</div></div>`;
};

// ---- 7 calculations
RENDER.calc = () => {
  const c = S.calc, f = c.fourier, p = c.pdf, A = lattice(), n = S.supercell.map(Number);
  const shape = f.steps.map((s, k) => 2 * Math.round(+f.hmax[k] / +s) + 1);
  const vox = shape.reduce((a, b) => a * b, 1);
  const ext = f.steps.map(s => 1 / (2 * +s)), res = f.hmax.map(h => 1 / (2 * +h));
  return head("Calculations", "What the generated script asks DISCUS 6.23.01 to compute, with discus-python's verified macros. Grids are centred on hkl = 0 so the 3D-ΔPDF is an exact FFT.") + `
  <div class="block"><h3>Running DISCUS</h3><div class="fields">${field("Command that runs a macro", txt("calc.launcher", 'class="w-xl"'))}${field("Time limit (s)", num("calc.timeout_s", 'class="w-s"'))}</div>
    <div class="fields"><button type="button" class="btn small" data-a="launcher" data-v="linux">Linux / WSL shell</button><button type="button" class="btn small" data-a="launcher" data-v="wsl">Windows → WSL</button>
      ${chk("calc.debye_waller", "Debye–Waller factor from the B values (temp use)")}</div>
    <p class="hint">The script appends <code>run.mac</code> to this command and starts it in each run folder. You can also pass <code>--launcher</code> when you run the script.</p></div>
  <div class="block"><h3>${chk("calc.fourier.enabled", "3D diffuse scattering")}</h3>
    <div class="fields">${field("Step Δh Δk Δl (r.l.u.)", trio("calc.fourier.steps"))}<button type="button" class="btn small" data-a="step-commensurate">1 / supercell</button>
      ${field("Range ± h k l", trio("calc.fourier.hmax"))}</div>
    <div class="fields">${field("Bragg subtraction (set aver, %)", num("calc.fourier.aver_percent", 'min="0" max="100" class="w-s"'))}
      ${field("Friedel symmetry", sel("calc.fourier.symmetry", [["", "engine default"], ["friedel", "friedel"], ["ignore", "ignore"], ["apply", "apply"]]))}
      ${chk("calc.fourier.delta_pdf", "3D-ΔPDF from the diffuse volume")}${chk("calc.fourier.total", "Also total scattering + 3D-PDF (Patterson)")}</div>
    <dl class="facts"><dt>Grid</dt><dd>${shape.join(" × ")} = ${(vox / 1e6).toFixed(2)} M points</dd>
      <dt>3D-ΔPDF range</dt><dd>±${ext.map(x => fmt(x, 2)).join(", ± ")} cells${A ? ` (±${ext.map((x, k) => fmt(x * DB.norm(A[k]), 1)).join(", ±")} Å)` : ""}</dd>
      <dt>3D-ΔPDF voxel</dt><dd>${res.map(x => fmt(x, 3)).join(" × ")} cells${A ? ` (${res.map((x, k) => fmt(x * DB.norm(A[k]), 2)).join(" × ")} Å)` : ""}</dd>
      <dt>Radiation</dt><dd>X-ray, Mo Kα1 (the only setting verified for discus-python's 3D runs)</dd></dl>
    <details><summary>LOTS (average over sub-volumes)</summary><div class="fields" style="margin-top:8px">${chk("calc.fourier.lots.enabled", "use LOTS")}
      ${field("Shape", sel("calc.fourier.lots.shape", [["box", "box"], ["ellipsoid", "ellipsoid"]]))}${field("Size (cells)", trio("calc.fourier.lots.size", 'step="1" class="w-s"'))}
      ${field("Number", txt("calc.fourier.lots.count", 'class="w-s"'))}${chk("calc.fourier.lots.periodic", "periodic")}</div>
      <p class="hint">Number is a count or “all”. LOTS apply to the diffuse run only.</p></details></div>
  <div class="block"><h3>${chk("calc.pdf.enabled", "1D powder PDF")}</h3>
    <div class="fields">${field("Radiation", sel("calc.pdf.radiation", [["xray", "X-ray"], ["neutron", "neutron"]]))}${field("Qmin (Å⁻¹)", num("calc.pdf.qmin", 'class="w-s"'))}
      ${field("Qmax (Å⁻¹)", num("calc.pdf.qmax", 'class="w-s"'))}${field("ΔQ", num("calc.pdf.dq", 'class="w-s"'))}${field("Gaussian width (set delta)", num("calc.pdf.qsig", 'class="w-s"'))}
      ${field("λ (Å, empty = default)", txt("calc.pdf.wavelength", 'class="w-s"'))}</div>
    <div class="fields">${field("rmax (Å)", num("calc.pdf.rmax", 'class="w-s"'))}${field("Δr (Å)", num("calc.pdf.dr", 'class="w-s"'))}</div>
    <p class="hint">Computed by DISCUS's powder Debye sum over the disordered supercell (N² pairs), so long-range peaks are damped by the finite crystal size. DISCUS 6.23.01 caps Q at ≈ 12.3/λ.</p></div>`;
};

// ---- 8 script
RENDER.script = () => {
  const v = validation;
  let code = "", err = "";
  if (!v.errs.length) { try { code = DB.generateScript(withState(), TEMPLATE); } catch (e) { err = e.message; } }
  const lines = code ? code.split("\n").length : 0;
  return head("Python script", "One self-contained script: it rebuilds this model with numpy and discus-python, writes the DISCUS structure and runs the calculations. The model is embedded as JSON, so you can re-import the script here later.") + `
  ${msgs(v.errs, "err")}${msgs(v.warns, "warn")}${err ? msgs([err], "err") : ""}
  <div class="block"><h3>Export <span class="grow"></span></h3>
    <div class="fields"><button type="button" class="btn primary" data-a="dl-script" ${code ? "" : "disabled"}>Download ${esc(fileBase())}.py</button>
      <button type="button" class="btn" data-a="copy-script" ${code ? "" : "disabled"}>Copy script</button>
      <button type="button" class="btn" data-a="dl-json">Download model JSON</button>
      ${code ? `<span class="hint">${lines} lines</span>` : ""}</div>
    <div class="doc"><p>Run it where DISCUS_SUITE 6.23.01 and discus-python are installed:</p>
      <pre class="code">pip install "discuspdf[data] @ git+https://github.com/jmsweng/discus-python"
python ${esc(fileBase())}.py --dry-run        # build the crystal, stage the DISCUS inputs
python ${esc(fileBase())}.py --numpy-preview  # quick diffuse plane, no DISCUS needed
python ${esc(fileBase())}.py                  # full run</pre>
      <p>Results: <code>supercell.stru</code>, <code>supercell.xyz</code>, <code>diffuse.npz</code> and slice plots, <code>delta_pdf.npz</code> / <code>.vts</code>, <code>pdf_1d.txt</code>, and every DISCUS run folder with its macro and log.</p></div></div>
  ${code ? `<pre class="code" id="code">${esc(code)}</pre>` : ""}`;
};
const fileBase = () => (S.title || "model").toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "").slice(0, 40) || "model";
function withState() { return S; }

// ------------------------------------------------------------------------------------------ actions
const ACT = {
  "supercell": el => { const k = +el.dataset.v; S.supercell = [k, k, k]; commit(); },
  "example": el => { S = EX[el.dataset.v].make(); molSel = S.molecules[0] && S.molecules[0].id; grpSel = S.orientations[0] && S.orientations[0].id; previewSeed = 0; commit(); toast("Example loaded"); },
  "mol-add": () => { const id = uid("m"); S.molecules.push({ id, name: "New molecule", anchor: "centroid", atoms: [{ el: "C", xyz: [0, 0, 0] }] }); molSel = id; commit(); },
  "mol-select": el => { molSel = el.dataset.v; atomSel = -1; renderSection(); },
  "mol-dup": () => { const m = S.molecules.find(x => x.id === molSel); const c = JSON.parse(JSON.stringify(m)); c.id = uid("m"); c.name += " copy"; S.molecules.push(c); molSel = c.id; commit(); },
  "mol-del": () => {
    const used = S.orientations.some(g => g.molecule === molSel);
    if (used) { toast("Delete the orientation groups that use this molecule first"); return; }
    S.molecules = S.molecules.filter(x => x.id !== molSel); molSel = S.molecules[0] && S.molecules[0].id; commit();
  },
  "atom-add": () => { const m = S.molecules.find(x => x.id === molSel); m.atoms.push({ el: placeEl || "C", xyz: [0, 0, 0] }); atomSel = m.atoms.length - 1; commit(); },
  "atom-del": el => { const m = S.molecules.find(x => x.id === molSel); m.atoms.splice(+el.dataset.v, 1); if (typeof m.anchor === "number" && m.anchor >= m.atoms.length) m.anchor = "centroid"; atomSel = -1; commit(); },
  "atom-pick": el => { atomSel = +el.dataset.v; renderSection(); },
  "mol-centre": () => { const m = S.molecules.find(x => x.id === molSel); const loc = DB.moleculeLocal(m); m.atoms.forEach((a, i) => { a.xyz = loc.xyz[i].map(x => +x.toFixed(4)); }); commit(); },
  "xyz-import": () => {
    const m = S.molecules.find(x => x.id === molSel); const lines = $("#xyz-in").value.split("\n").map(l => l.trim()).filter(Boolean);
    const atoms = [];
    for (const l of lines) { const p = l.split(/[\s,]+/); if (p.length >= 4 && DB.elementOf(p[0]) && p.slice(1, 4).every(x => isFinite(+x))) atoms.push({ el: DB.elementOf(p[0]), xyz: p.slice(1, 4).map(Number) }); }
    if (!atoms.length) { toast("No “Element x y z” lines found"); return; }
    m.atoms = atoms; m.anchor = "centroid"; commit(); toast(`${atoms.length} atoms read`);
  },
  "place-toggle": () => { placing = !placing; const e = $("#place-el"); if (e) placeEl = e.value.trim() || "C"; renderSection(); },
  "mol-view": el => { const A = lattice(); if (!molViewer) return; if (el.dataset.v === "z") molViewer.viewAlong([0, 0, 1], [0, 1, 0]); else molViewer.viewAlong([1, 0, 0], [0, 0, 1]); window._molR = molViewer.R; void A; },
  "grp-add": () => {
    if (!S.molecules.length) { toast("Build a molecule first"); return; }
    const id = uid("g"); S.orientations.push({ id, name: "New group", molecule: (S.molecules.find(m => m.id === molSel) || S.molecules[0]).id, convention: "ZYZ", euler: [0, 0, 0], ops: "", frame: "lab", closure: true, weights: [] });
    grpSel = id; commit();
  },
  "grp-select": el => { grpSel = el.dataset.v; renderSection(); },
  "grp-dup": () => { const g = S.orientations.find(x => x.id === grpSel); const c = JSON.parse(JSON.stringify(g)); c.id = uid("g"); c.name += " copy"; S.orientations.push(c); grpSel = c.id; commit(); },
  "grp-del": () => {
    if (S.sites.some(s => s.occupants.some(o => o.type === "molecule" && o.group === grpSel))) { toast("A site still uses this group"); return; }
    S.orientations = S.orientations.filter(x => x.id !== grpSel); grpSel = S.orientations[0] && S.orientations[0].id; commit();
  },
  "op-add": el => { const g = S.orientations.find(x => x.id === grpSel); g.ops = (g.ops.trim() ? g.ops.trim() + "\n" : "") + el.dataset.v; commit(); },
  "mat-show": el => {
    const g = S.orientations.find(x => x.id === grpSel), mol = S.molecules.find(x => x.id === g.molecule);
    const R = DB.expandGroup(g, mol, lattice()).orientations[+el.dataset.v].R;
    const box = $("#matbox"); box.hidden = false;
    box.textContent = `Orientation #${+el.dataset.v + 1}, Cartesian rotation matrix (applied to the molecule's coordinates):\n` + R.map(r => r.map(x => (x >= 0 ? " " : "") + x.toFixed(6)).join("  ")).join("\n") + `\ndet = ${DB.det(R).toFixed(6)}`;
  },
  "site-add": () => { S.sites.push({ label: "S" + (S.sites.length + 1), frac: [0, 0, 0], b_iso: 1, draw: "random", occupants: [{ type: "atom", species: "C", p: 1, b: "" }] }); commit(); },
  "site-del": el => { const i = +el.dataset.v; S.sites.splice(i, 1);
    S.disorder.correlations.interactions = S.disorder.correlations.interactions.filter(it => +it.site_a !== i && +it.site_b !== i).map(it => ({ ...it, site_a: +it.site_a > i ? +it.site_a - 1 : +it.site_a, site_b: +it.site_b > i ? +it.site_b - 1 : +it.site_b }));
    commit(); },
  "occ-add": el => {
    const [si, type] = el.dataset.v.split("."); const s = S.sites[+si];
    const o = type === "atom" ? { type, species: "C", p: 0, b: "" } : type === "molecule" ? { type, group: (S.orientations[0] || {}).id || "", p: 0, b: "" } : { type, p: 0 };
    s.occupants.push(o); commit();
  },
  "occ-del": el => { const [si, oi] = el.dataset.v.split(".").map(Number); S.sites[si].occupants.splice(oi, 1); commit(); },
  "occ-norm": el => { const s = S.sites[+el.dataset.v]; const t = s.occupants.reduce((a, o) => a + (+o.p || 0), 0); if (t > 0) s.occupants.forEach(o => { o.p = +((+o.p || 0) / t).toFixed(6); }); commit(); },
  "int-add": () => { S.disorder.correlations.interactions.push({ site_a: 0, site_b: 0, t: [1, 0, 0], J: -0.5 }); S.disorder.correlations.enabled = true; commit(); },
  "int-axes": () => { [[1, 0, 0], [0, 1, 0], [0, 0, 1]].forEach(t => S.disorder.correlations.interactions.push({ site_a: 0, site_b: 0, t, J: -0.5 })); S.disorder.correlations.enabled = true; commit(); },
  "int-del": el => { S.disorder.correlations.interactions.splice(+el.dataset.v, 1); commit(); },
  "nb-add": () => { S.mmc.neighborhoods.push({ vectors: [{ origin: 1, target: 1, t: [1, 0, 0] }], symmetric: true }); commit(); },
  "nb-del": el => { S.mmc.neighborhoods.splice(+el.dataset.v, 1); S.mmc.targets = S.mmc.targets.filter(t => +t.neighborhood !== +el.dataset.v); commit(); },
  "vec-add": el => { S.mmc.neighborhoods[+el.dataset.v].vectors.push({ origin: 1, target: 1, t: [0, 1, 0] }); commit(); },
  "vec-axes": el => { const nb = S.mmc.neighborhoods[+el.dataset.v]; [[1, 0, 0], [0, 1, 0], [0, 0, 1]].forEach(t => nb.vectors.push({ origin: 1, target: 1, t })); commit(); },
  "vec-del": el => { const [k, j] = el.dataset.v.split(".").map(Number); S.mmc.neighborhoods[k].vectors.splice(j, 1); commit(); },
  "tg-add": () => { const sp = allSpecies(); S.mmc.targets.push({ type: "occupational", neighborhood: 0, a: sp[0] || "Cu", b: sp[1] || sp[0] || "Au", mode: "energy", value: 0.5, direction: [1, 0, 0], distance: 2.5, k: 1, depth: 1 }); commit(); },
  "tg-del": el => { S.mmc.targets.splice(+el.dataset.v, 1); commit(); },
  "launcher": el => { S.calc.launcher = el.dataset.v === "wsl" ? "wsl.exe -d Ubuntu-22.04 -- bash -lc 'discus_suite -noautorun -macro \"$0\"'" : "discus_suite -noautorun -macro"; commit(); },
  "step-commensurate": () => { S.calc.fourier.steps = S.supercell.map(n => +(1 / +n).toFixed(10)); commit(); },
  "dl-script": () => { try { download(`${fileBase()}.py`, DB.generateScript(withState(), TEMPLATE), "text/x-python"); } catch (e) { toast(e.message); } },
  "copy-script": () => copyText(DB.generateScript(withState(), TEMPLATE)),
  "dl-json": () => download(`${fileBase()}.json`, JSON.stringify({ schema: "discuspdf-disorder-builder-state", version: 1, state: S }, null, 1), "application/json"),
  "export": () => ACT["dl-json"](),
  "import": () => $("#file-in").click(),
  "reseed": () => { previewSeed++; updatePreview(); },
  "cells": el => { cellsShown = +el.dataset.v; drawCrystal(); $$('[data-a="cells"]').forEach(b => b.setAttribute("aria-pressed", String(+b.dataset.v === cellsShown))); },
  "color": el => { colorMode = el.dataset.v; drawCrystal(); $$('[data-a="color"]').forEach(b => b.setAttribute("aria-pressed", String(b.dataset.v === colorMode))); },
  "diffuse": () => runDiffuse(),
  "plane": el => { const v = el.dataset.v; const P = { hk0: [[1, 0, 0], [0, 1, 0]], h0l: [[1, 0, 0], [0, 0, 1]], "0kl": [[0, 1, 0], [0, 0, 1]], hhl: [[1, 1, 0], [0, 0, 1]] }[v];
    S.calc.preview.u = P[0]; S.calc.preview.v = P[1]; save(); renderBenchControls(); markDiffuseStale(); },
};
function download(name, text, type) {
  try {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = document.createElement("a"); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
    toast(`Downloading ${name}. If nothing arrives, use Copy instead.`);
  } catch (e) { toast("Download blocked here; use Copy instead"); }
}
function copyText(t) {
  const fallback = () => { const box = $("#code"); if (box) { const r = document.createRange(); r.selectNodeContents(box); const s = getSelection(); s.removeAllRanges(); s.addRange(r); } toast("Selected: press Ctrl/⌘+C to copy"); };
  try { navigator.clipboard.writeText(t).then(() => toast("Script copied"), fallback); } catch (e) { fallback(); }
}
let toastT;
function toast(m) { let el = $("#toast"); el.textContent = m; el.hidden = false; clearTimeout(toastT); toastT = setTimeout(() => { el.hidden = true; }, 2600); }

// ------------------------------------------------------------------------------------------ events
document.addEventListener("change", e => {
  const el = e.target;
  if (el.id === "file-in") { importFile(el.files[0]); el.value = ""; return; }
  if (el.id === "pv-range" || el.id === "pv-n" || el.id === "pv-d") { S.calc.preview.range = +$("#pv-range").value || 4; S.calc.preview.n = +$("#pv-n").value || 8; S.calc.preview.draws = +$("#pv-d").value || 4; save(); markDiffuseStale(); return; }
  const path = el.dataset.b; if (!path) return;
  let v = el.value;
  const t = el.dataset.t;
  if (t === "num") { v = el.value === "" ? 0 : +el.value; if (!isFinite(v)) v = 0; }
  else if (t === "numornull") v = el.value === "" ? "" : +el.value;
  else if (t === "bool") v = el.checked;
  if (path.endsWith(".type") && path.includes(".occupants.")) {
    const o = getPath(path.replace(/\.type$/, ""));
    if (v === "atom" && !o.species) o.species = "C";
    if (v === "molecule" && !o.group) o.group = (S.orientations[0] || {}).id || "";
  }
  if (path.endsWith(".group") && path.includes("disorder") === false) { /* group changes keep their weights */ }
  setPath(path, v);
  if (path === "orientations." + S.orientations.findIndex(g => g.id === grpSel) + ".ops" || /\.(closure|frame|euler\.\d|molecule|convention)$/.test(path)) {
    const g = S.orientations.find(x => x.id === grpSel); if (g && /orientations\.\d+\.(ops|closure|frame|molecule)$/.test(path)) g.weights = [];
  }
  commit();
});
document.addEventListener("click", e => {
  const s = e.target.closest("[data-sec]");
  if (s) { section = s.dataset.sec; placing = false; renderRail(); renderSection(); window.scrollTo(0, 0); return; }
  const a = e.target.closest("[data-a]");
  if (a && ACT[a.dataset.a] && !a.disabled) { e.preventDefault(); ACT[a.dataset.a](a); }
});
function importFile(file) {
  if (!file) return;
  const rd = new FileReader();
  rd.onload = () => {
    try {
      let text = String(rd.result), obj;
      const m = /CONFIG_JSON = r"""\s*([\s\S]*?)\s*"""/.exec(text);
      obj = JSON.parse(m ? m[1] : text);
      const st = obj.state || obj.builder_state || (obj.cell && obj.sites && obj.calc && obj.disorder && Array.isArray(obj.molecules) ? obj : null);
      if (!st) throw new Error("this file holds no builder model");
      S = upgrade(st); molSel = S.molecules[0] && S.molecules[0].id; grpSel = S.orientations[0] && S.orientations[0].id;
      commit(); toast("Model imported");
    } catch (err) { toast(`Could not import: ${err.message}`); }
  };
  rd.readAsText(file);
}

// ------------------------------------------------------------------------------------------ preview bench
let crystalViewer = null, previewTimer = null;
function schedulePreview() { clearTimeout(previewTimer); previewTimer = setTimeout(updatePreview, 250); }
function previewN() { const pn = +S.calc.preview.n || 8; return S.supercell.map(x => Math.max(1, Math.min(+x || 1, pn))); }
function updatePreview() {
  if (validation.errs.some(e => /cell|Supercell|probabilit|orientation group|not an element/i.test(e))) { $("#xstats").textContent = "Fix the problems listed under Script to see the preview."; return; }
  try { lastReal = DB.realize(S, { n: previewN(), seed: (+S.seed || 1) + 7919 * previewSeed, maxSweeps: 60 }); }
  catch (e) { $("#xstats").textContent = e.message; return; }
  drawCrystal();
  if (section === "disorder") { const y = window.scrollY; renderSection(); window.scrollTo(0, y); }
}
const ORI = ["#0d6a80", "#d9822b", "#5f9e3a", "#9b59b6", "#c0392b", "#2c7fb8", "#c9a227", "#16a085", "#e67e22", "#7f8c8d", "#8e44ad", "#27ae60"];
function drawCrystal() {
  if (!crystalViewer || !lastReal) return;
  const { atoms, A } = lastReal, k = cellsShown;
  const sub = atoms.filter(a => a.cell[0] < k && a.cell[1] < k && a.cell[2] < k);
  const groups = S.orientations.map(g => g.id);
  const color = a => {
    if (colorMode === "element") return DB.colorOf(a.el);
    if (colorMode === "site") return ORI[a.site % ORI.length];
    if (a.mol < 0) return "#9aa7ad";
    return ORI[(groups.indexOf(a.group) * 3 + a.orient) % ORI.length];
  };
  const mols = {};
  sub.forEach((a, i) => { if (a.mol >= 0) (mols[a.mol] = mols[a.mol] || []).push(i); });
  const bonds = [];
  Object.values(mols).forEach(ix => { for (let x = 0; x < ix.length; x++) for (let y = x + 1; y < ix.length; y++) {
    const p = sub[ix[x]], q = sub[ix[y]], d = DB.norm(DB.sub(p.p, q.p));
    if (d < DB.covOf(p.el) + DB.covOf(q.el) + 0.45) bonds.push([ix[x], ix[y]]);
  } });
  const lines = [], O = [0, 0, 0], a = DB.scale(A[0], k), b = DB.scale(A[1], k), c = DB.scale(A[2], k);
  const corners = [[O, a], [O, b], [O, c], [a, DB.add(a, b)], [a, DB.add(a, c)], [b, DB.add(a, b)], [b, DB.add(b, c)], [c, DB.add(a, c)], [c, DB.add(b, c)], [DB.add(a, b), DB.add(DB.add(a, b), c)], [DB.add(a, c), DB.add(DB.add(a, b), c)], [DB.add(b, c), DB.add(DB.add(a, b), c)]];
  corners.forEach(([p, q]) => lines.push([p, q, null, 0.5, 1]));
  crystalViewer.set({ atoms: sub.map(x => ({ p: x.p, color: color(x), r: 0.12 + 0.13 * DB.covOf(x.el) })), bonds, lines, axes: A,
    center: DB.scale(DB.add(DB.add(a, b), c), 0.5) });
  // legend + stats
  const comp = {};
  atoms.forEach(x => { comp[x.el] = (comp[x.el] || 0) + 1; });
  const nc = lastReal.n.reduce((p, q) => p * q, 1);
  let legend = "";
  if (colorMode === "element") legend = Object.keys(comp).map(e => `<span><span class="dot" style="background:${DB.colorOf(e)}"></span>${esc(e)}</span>`).join("");
  else if (colorMode === "site") legend = S.sites.map((s, i) => `<span><span class="dot" style="background:${ORI[i % ORI.length]}"></span>${esc(s.label)}</span>`).join("");
  else {
    const otab = DB.orientationTable(S, A);
    legend = Object.values(otab).flatMap((t, gi) => t.orientations.map((o, k2) => `<span><span class="dot" style="background:${ORI[(groups.indexOf(t.group.id) * 3 + k2) % ORI.length]}"></span>${esc(t.group.name)} #${k2 + 1}</span>`)).join("") + `<span><span class="dot" style="background:#9aa7ad"></span>atoms</span>`;
  }
  $("#xlegend").innerHTML = legend;
  $("#xstats").innerHTML = `${lastReal.n.join("×")} preview cells · ${(atoms.length / nc).toFixed(2)} atoms/cell · ${lastReal.nMolecules} molecules<br>` +
    Object.entries(comp).map(([e, c2]) => `${esc(e)} ${(c2 / nc).toFixed(3)}`).join(" · ") + " per cell";
}

// diffuse
let diffuseJob = 0;
function markDiffuseStale() { const s = $("#dstate"); if (s && !s.dataset.running) s.textContent = "Model changed: compute again to update."; }
function renderBenchControls() {
  const pv = S.calc.preview;
  const key = JSON.stringify([pv.u, pv.v]);
  const names = { hk0: [[1, 0, 0], [0, 1, 0]], h0l: [[1, 0, 0], [0, 0, 1]], "0kl": [[0, 1, 0], [0, 0, 1]], hhl: [[1, 1, 0], [0, 0, 1]] };
  $("#planes").innerHTML = Object.entries(names).map(([k, v]) => `<button type="button" class="btn small" data-a="plane" data-v="${k}" aria-pressed="${JSON.stringify(v) === key}">${k === "hhl" ? "hhl" : k}</button>`).join("");
  $("#pv-range").value = pv.range; $("#pv-n").value = pv.n || 8; $("#pv-d").value = pv.draws || 4;
}
function runDiffuse() {
  const job = ++diffuseJob;
  const pv = S.calc.preview, draws = Math.max(1, Math.min(16, +pv.draws || 4));
  let real, plane, P, sum, d = 0, a = 0, natoms = 0;
  const st = $("#dstate"), bar = $("#dprog span");
  const t0 = performance.now();
  const start = () => {
    real = DB.realize(S, { n: previewN(), seed: (+S.seed || 1) + 7919 * previewSeed + 104729 * d, maxSweeps: 60 });
    plane = DB.diffusePlaneRows(real, { u: pv.u, v: pv.v, range: +pv.range || 4 });
    if (!sum) { P = plane.P; sum = new Float64Array(P * P); }
    natoms = real.atoms.length; a = 0;
  };
  try { start(); } catch (e) { toast(e.message); return; }
  st.dataset.running = "1"; st.textContent = `Computing ${P}×${P} points, ${draws} random crystals of ${natoms} atoms…`;
  const step = () => {
    if (job !== diffuseJob) return;
    const tEnd = performance.now() + 30;
    while (a < P && performance.now() < tEnd) { const r = plane.row(a); for (let b = 0; b < P; b++) sum[a * P + b] += r[b]; a++; }
    bar.style.width = `${100 * (d * P + a) / (draws * P)}%`;
    if (a >= P) { d++; drawHeat(Float64Array.from(sum, x => x / d), P, plane.t); if (d < draws) { start(); } }
    if (d < draws) requestAnimationFrame(step);
    else { delete st.dataset.running; st.textContent = `Average of ${draws} random ${real.n.join("×")} crystals (${natoms} atoms each), ${((performance.now() - t0) / 1000).toFixed(1)} s. f = Z, no Debye–Waller: a shape check, not the DISCUS result.`; }
  };
  requestAnimationFrame(step);
}
function magma(t) {
  const st = [[0, 0, 4], [40, 11, 84], [101, 21, 110], [159, 42, 99], [212, 72, 66], [245, 125, 21], [250, 193, 39], [252, 253, 191]];
  t = Math.max(0, Math.min(1, t)) * (st.length - 1); const i = Math.min(st.length - 2, Math.floor(t)), f = t - i;
  return st[i].map((x, k) => x + (st[i + 1][k] - x) * f);
}
function drawHeat(img, P, t) {
  const cv = $("#dcv"), g = cv.getContext("2d");
  cv.width = P; cv.height = P;
  const vals = Array.from(img).filter(x => x > 0).sort((x, y) => x - y);
  const vmax = vals.length ? vals[Math.floor(0.995 * (vals.length - 1))] : 1;
  const id = g.createImageData(P, P);
  for (let a = 0; a < P; a++) for (let b = 0; b < P; b++) {
    const v = img[a * P + b], c = magma(Math.sqrt(Math.max(0, v) / (vmax || 1)));
    const o = 4 * ((P - 1 - b) * P + a);
    id.data[o] = c[0]; id.data[o + 1] = c[1]; id.data[o + 2] = c[2]; id.data[o + 3] = 255;
  }
  g.putImageData(id, 0, 0);
  const pv = S.calc.preview, lab = v => `[${v.join(" ")}]`;
  $("#daxes").innerHTML = `<span>${fmt(t[0], 2)}</span><span>→ ${lab(pv.u)} (r.l.u.) · ↑ ${lab(pv.v)}</span><span>${fmt(t[t.length - 1], 2)}</span>`;
}

// ------------------------------------------------------------------------------------------ boot
function boot() {
  renderRail();
  try { validation = DB.validate(S); } catch (e) { validation = { errs: [e.message], warns: [] }; }
  renderStatus();
  renderSection();
  crystalViewer = new Viewer($("#xcv"), { aspect: 0.85 });
  $$('[data-a="cells"]').forEach(b => b.setAttribute("aria-pressed", String(+b.dataset.v === cellsShown)));
  $$('[data-a="color"]').forEach(b => b.setAttribute("aria-pressed", String(b.dataset.v === colorMode)));
  renderBenchControls();
  updatePreview();
  setTimeout(runDiffuse, 400);
  if (window.matchMedia) window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { drawCrystal(); if (molViewer) drawMolecule(); });
}
boot();
})();
