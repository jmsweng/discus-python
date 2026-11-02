/* Disorder Builder core: lattice maths, orientations, random realisation, diffuse preview and
   script generation. No DOM access, so it runs in Node for tests as well as in the page. */
(function (root) {
"use strict";

// ---------------------------------------------------------------------------------- elements
const SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr " +
  "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W " +
  "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm").split(" ");
const Z = {};
SYMBOLS.forEach((s, i) => { Z[s] = i + 1; });
// covalent radii (Å, Cordero et al. 2008) for the common elements; others fall back to 1.5 Å
const COV = { H: .31, He: .28, Li: 1.28, Be: .96, B: .84, C: .76, N: .71, O: .66, F: .57, Ne: .58, Na: 1.66, Mg: 1.41,
  Al: 1.21, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02, Ar: 1.06, K: 2.03, Ca: 1.76, Sc: 1.70, Ti: 1.60, V: 1.53, Cr: 1.39,
  Mn: 1.39, Fe: 1.32, Co: 1.26, Ni: 1.24, Cu: 1.32, Zn: 1.22, Ga: 1.22, Ge: 1.20, As: 1.19, Se: 1.20, Br: 1.20,
  Kr: 1.16, Rb: 2.20, Sr: 1.95, Y: 1.90, Zr: 1.75, Nb: 1.64, Mo: 1.54, Ag: 1.45, Cd: 1.44, In: 1.42, Sn: 1.39,
  Sb: 1.39, Te: 1.38, I: 1.39, Xe: 1.40, Cs: 2.44, Ba: 2.15, La: 2.07, Pt: 1.36, Au: 1.36, Hg: 1.32, Pb: 1.46, Bi: 1.48,
  U: 1.96 };
const COLORS = { H: "#e8ecef", C: "#5b6670", N: "#3d6fd6", O: "#e0412f", F: "#7fd34e", Cl: "#2fbf5a", Br: "#a3392b",
  I: "#7b3fa8", S: "#e6c229", P: "#f08a24", Na: "#8a5cf0", K: "#9f6ae6", Cs: "#5a2d9c", Li: "#c086ff", Mg: "#4fa36b",
  Ca: "#57b884", Si: "#d9b38c", Al: "#b9a6a6", Fe: "#d2722e", Cu: "#c8803c", Ag: "#b9c2cc", Au: "#e3b43c",
  Zn: "#7d80b0", Ni: "#4fbf6a", Co: "#e07ba0", Ti: "#9aa3ab", Zr: "#7fbfbf", Pb: "#57595f", Sn: "#7b8c99",
  B: "#e9a3a3", Ba: "#2fae8e", Sr: "#3fbf9f", Mn: "#a67fd6", Cr: "#8a9bd6", V: "#a6a6ab", La: "#5fb3d6", U: "#2f8fd6" };
function elementOf(sp) { const m = /^([A-Z][a-z]?)/.exec(sp || ""); return m && Z[m[1]] ? m[1] : null; }
function zOf(sp) { const e = elementOf(sp); return e ? Z[e] : 0; }
function colorOf(sp) {
  const e = elementOf(sp);
  if (e && COLORS[e]) return COLORS[e];
  const z = e ? Z[e] : 7; return `hsl(${(z * 47) % 360} 45% 55%)`;
}
function covOf(sp) { const e = elementOf(sp); return (e && COV[e]) || 1.5; }

// ---------------------------------------------------------------------------------- 3x3 maths
const I3 = () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
function mul(A, B) {
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
  return C;
}
function tr(A) { return [[A[0][0], A[1][0], A[2][0]], [A[0][1], A[1][1], A[2][1]], [A[0][2], A[1][2], A[2][2]]]; }
function mv(A, v) { return [A[0][0] * v[0] + A[0][1] * v[1] + A[0][2] * v[2], A[1][0] * v[0] + A[1][1] * v[1] + A[1][2] * v[2], A[2][0] * v[0] + A[2][1] * v[1] + A[2][2] * v[2]]; }
function det(A) {
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}
function inv(A) {
  const d = det(A);
  if (Math.abs(d) < 1e-14) throw new Error("singular matrix");
  const m = [[A[1][1] * A[2][2] - A[1][2] * A[2][1], A[0][2] * A[2][1] - A[0][1] * A[2][2], A[0][1] * A[1][2] - A[0][2] * A[1][1]],
    [A[1][2] * A[2][0] - A[1][0] * A[2][2], A[0][0] * A[2][2] - A[0][2] * A[2][0], A[0][2] * A[1][0] - A[0][0] * A[1][2]],
    [A[1][0] * A[2][1] - A[1][1] * A[2][0], A[0][1] * A[2][0] - A[0][0] * A[2][1], A[0][0] * A[1][1] - A[0][1] * A[1][0]]];
  return m.map(r => r.map(x => x / d));
}
function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
function norm(a) { return Math.hypot(a[0], a[1], a[2]); }
function maxDiff(A, B) { let m = 0; for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) m = Math.max(m, Math.abs(A[i][j] - B[i][j])); return m; }
function isOrthogonal(R, tol = 1e-6) { return maxDiff(mul(R, tr(R)), I3()) < tol; }
function clean(R) { return R.map(r => r.map(x => (Math.abs(x) < 1e-12 ? 0 : x))); }

// rows a, b, c in the discuspdf canonical orientation (a ∥ x, b in xy, right-handed)
function latticeMatrix(cell) {
  const d = Math.PI / 180, { a, b, c } = cell;
  const ca = Math.cos(cell.alpha * d), cb = Math.cos(cell.beta * d), cg = Math.cos(cell.gamma * d), sg = Math.sin(cell.gamma * d);
  const cx = c * cb, cy = c * (ca - cb * cg) / sg, cz2 = c * c - cx * cx - cy * cy;
  if (!(cz2 > 0) || !(a > 0 && b > 0 && c > 0)) throw new Error("invalid unit cell");
  return [[a, 0, 0], [b * cg, b * sg, 0], [cx, cy, Math.sqrt(cz2)]];
}
const frac2cart = (A, f) => [f[0] * A[0][0] + f[1] * A[1][0] + f[2] * A[2][0], f[0] * A[0][1] + f[1] * A[1][1] + f[2] * A[2][1], f[0] * A[0][2] + f[1] * A[1][2] + f[2] * A[2][2]];

function axisAngle(axis, deg) {
  const n = norm(axis);
  if (!(n > 0)) throw new Error("zero rotation axis");
  const [x, y, z] = scale(axis, 1 / n), t = deg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t), C = 1 - c;
  return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s], [y * x * C + z * s, c + y * y * C, y * z * C - x * s], [z * x * C - y * s, z * y * C + x * s, c + z * z * C]];
}
const Rx = d => axisAngle([1, 0, 0], d), Ry = d => axisAngle([0, 1, 0], d), Rz = d => axisAngle([0, 0, 1], d);
// active rotations; intrinsic convention (R = R1(α) R2(β) R3(γ))
function eulerMatrix(conv, a, b, g) {
  if (conv === "ZXZ") return mul(mul(Rz(a), Rx(b)), Rz(g));
  if (conv === "XYZ") return mul(mul(Rx(a), Ry(b)), Rz(g));
  return mul(mul(Rz(a), Ry(b)), Rz(g)); // ZYZ
}

// ---------------------------------------------------------------------------------- symmetry operations
// Syntax (one per line):
//   -y,x,z            Jones symbol in crystal fractions (translations ignored)
//   C4 [0 0 1]        n-fold rotation about a lattice direction [u v w] (or Cartesian <x y z>)
//   S4 [0 0 1]        rotoreflection;  m [1 -1 0]  mirror with that normal;  -1  inversion;  E  identity
//   R 30 [1 1 1]      rotation by an angle in degrees
//   M a b c d e f g h i   explicit Cartesian 3x3 matrix (rows)
function parseDirection(txt, A) {
  const m = /([\[<])\s*([-+\d.eE\s,]+)\s*[\]>]/.exec(txt);
  if (!m) throw new Error(`no direction in "${txt}" (write [u v w] or <x y z>)`);
  const v = m[2].trim().split(/[\s,]+/).map(Number);
  if (v.length !== 3 || v.some(x => !isFinite(x))) throw new Error(`direction needs three numbers in "${txt}"`);
  return m[1] === "[" ? frac2cart(A, v) : v;
}
function parseJones(txt, A) {
  const parts = txt.split(",").map(s => s.trim().toLowerCase());
  if (parts.length !== 3) throw new Error(`"${txt}" needs three comma-separated parts`);
  const M = parts.map(p => {
    const row = [0, 0, 0];
    const re = /([+-]?)\s*(\d*\.?\d*(?:\/\d+)?)\s*\*?\s*([xyz])?/g;
    let m, any = false;
    const s = p.replace(/\s+/g, "");
    while ((m = re.exec(s)) !== null) {
      if (m[0] === "") { re.lastIndex++; continue; }
      if (!m[3]) continue; // translation part, ignored
      any = true;
      let coef = 1;
      if (m[2]) coef = m[2].includes("/") ? Number(m[2].split("/")[0]) / Number(m[2].split("/")[1]) : Number(m[2]);
      if (m[1] === "-") coef = -coef;
      row["xyz".indexOf(m[3])] += coef;
    }
    if (!any) throw new Error(`"${p}" has no x, y or z`);
    return row;
  });
  const At = tr(A);
  return mul(mul(At, M), inv(At));
}
function parseOp(line, A) {
  const txt = line.trim();
  if (!txt || txt.startsWith("#")) return null;
  let R, kind;
  let m;
  if (/^(E|1|id|identity)$/i.test(txt)) { R = I3(); kind = "E"; }
  else if (/^(-1|i|inv|inversion)$/i.test(txt)) { R = scale3(I3(), -1); kind = "-1"; }
  else if ((m = /^C(\d+)\s*(.*)$/i.exec(txt))) { R = axisAngle(parseDirection(m[2], A), 360 / Number(m[1])); kind = `C${m[1]}`; }
  else if ((m = /^S(\d+)\s*(.*)$/i.exec(txt))) {
    const d = parseDirection(m[2], A);
    R = mul(mirror(d), axisAngle(d, 360 / Number(m[1]))); kind = `S${m[1]}`;
  }
  else if ((m = /^m\s*(.*)$/i.exec(txt)) && /[\[<]/.test(m[1])) { R = mirror(parseDirection(m[1], A)); kind = "m"; }
  else if ((m = /^R\s+([-+\d.eE]+)\s*(.*)$/i.exec(txt))) { R = axisAngle(parseDirection(m[2], A), Number(m[1])); kind = "R"; }
  else if ((m = /^M\s+(.*)$/i.exec(txt))) {
    const v = m[1].trim().split(/[\s,]+/).map(Number);
    if (v.length !== 9 || v.some(x => !isFinite(x))) throw new Error("M needs nine numbers");
    R = [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)]; kind = "M";
  }
  else if (/[xyz]/i.test(txt) && txt.includes(",")) { R = parseJones(txt, A); kind = "jones"; }
  else throw new Error(`cannot read "${txt}"`);
  if (!isOrthogonal(R, 1e-5)) {
    throw new Error(kind === "jones"
      ? `"${txt}" is not a symmetry of this lattice (the Cartesian matrix is not orthogonal)`
      : `"${txt}" is not orthogonal`);
  }
  return { R: clean(R), text: txt, improper: det(R) < 0 };
}
function scale3(M, s) { return M.map(r => r.map(x => x * s)); }
function mirror(n) {
  const l = norm(n), u = scale(n, 1 / l);
  return [[1 - 2 * u[0] * u[0], -2 * u[0] * u[1], -2 * u[0] * u[2]], [-2 * u[1] * u[0], 1 - 2 * u[1] * u[1], -2 * u[1] * u[2]], [-2 * u[2] * u[0], -2 * u[2] * u[1], 1 - 2 * u[2] * u[2]]];
}
function closure(gens, limit = 240) {
  const G = [I3()];
  const has = M => G.some(X => maxDiff(X, M) < 1e-6);
  let grew = true;
  while (grew) {
    grew = false;
    for (const X of G.slice()) for (const g of gens) {
      const Y = clean(mul(g, X));
      if (!has(Y)) { G.push(Y); grew = true; if (G.length > limit) throw new Error(`group closure exceeds ${limit} operations; check the generators`); }
    }
  }
  return G;
}

// ---------------------------------------------------------------------------------- molecules and orientations
function moleculeLocal(mol) {
  const xyz = mol.atoms.map(a => a.xyz.map(Number));
  let o = [0, 0, 0];
  if (mol.anchor === "com") {
    let mt = 0;
    mol.atoms.forEach((a, i) => { const m = massOf(a.el); mt += m; o = add(o, scale(xyz[i], m)); });
    o = scale(o, 1 / mt);
  } else if (typeof mol.anchor === "number") o = xyz[mol.anchor] || [0, 0, 0];
  else { xyz.forEach(p => { o = add(o, p); }); o = scale(o, 1 / Math.max(1, xyz.length)); }
  return { species: mol.atoms.map(a => a.el), xyz: xyz.map(p => sub(p, o)) };
}
function massOf(sp) { const z = zOf(sp); return z === 1 ? 1.008 : 2 * z + (z > 20 ? 0.5 : 0); }
function bondsOf(species, xyz, tol = 0.45) {
  const out = [];
  for (let i = 0; i < xyz.length; i++) for (let j = i + 1; j < xyz.length; j++) {
    const d = norm(sub(xyz[i], xyz[j]));
    if (d > 0.4 && d < covOf(species[i]) + covOf(species[j]) + tol) out.push([i, j]);
  }
  return out;
}
// same molecule placement? (rotated atom sets coincide, species respected)
function samePlacement(species, xyz, R1, R2, tol = 0.02) {
  const P1 = xyz.map(p => mv(R1, p)), P2 = xyz.map(p => mv(R2, p));
  const used = new Array(P2.length).fill(false);
  for (let i = 0; i < P1.length; i++) {
    let found = -1;
    for (let j = 0; j < P2.length; j++) if (!used[j] && species[j] === species[i] && norm(sub(P1[i], P2[j])) < tol) { found = j; break; }
    if (found < 0) return false;
    used[found] = true;
  }
  return true;
}
// expand an orientation group: base Euler orientation, then symmetry operations
function expandGroup(group, mol, A) {
  const local = moleculeLocal(mol);
  const R0 = eulerMatrix(group.convention || "ZYZ", +group.euler[0] || 0, +group.euler[1] || 0, +group.euler[2] || 0);
  const errors = [], ops = [];
  (group.ops || "").split("\n").forEach((line, k) => {
    try { const op = parseOp(line, A); if (op) ops.push(op); } catch (e) { errors.push(`line ${k + 1}: ${e.message}`); }
  });
  let G;
  if (group.closure) {
    try { G = closure(ops.map(o => o.R)); } catch (e) { errors.push(e.message); G = [I3()]; }
  } else G = [I3(), ...ops.map(o => o.R)];
  const out = [];
  for (const g of G) {
    const R = clean(group.frame === "molecule" ? mul(R0, g) : mul(g, R0));
    const dup = out.findIndex(o => samePlacement(local.species, local.xyz, o.R, R));
    if (dup >= 0) { out[dup].multiplicity++; continue; }
    out.push({ R, improper: det(R) < 0, multiplicity: 1 });
  }
  return { orientations: out, errors, nOps: G.length };
}

// ---------------------------------------------------------------------------------- random numbers
function rng32(seed) {
  let a = (seed >>> 0) || 1;
  const next = () => { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  next.normal = () => { let u = 0; while (u === 0) u = next(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * next()); };
  next.int = n => Math.floor(next() * n);
  return next;
}

// ---------------------------------------------------------------------------------- model → classes
function orientationTable(state, A) {
  const tab = {};
  for (const g of state.orientations) {
    const mol = state.molecules.find(m => m.id === g.molecule);
    if (!mol) continue;
    const ex = expandGroup(g, mol, A);
    const w = ex.orientations.map((o, k) => {
      const v = g.weights && g.weights[k] != null ? +g.weights[k] : 1;
      return isFinite(v) && v >= 0 ? v : 0;
    });
    tab[g.id] = { group: g, mol, orientations: ex.orientations, weights: w, errors: ex.errors, nOps: ex.nOps };
  }
  return tab;
}
function siteClasses(site, otab) {
  const cls = [];
  for (const o of site.occupants) {
    const p = +o.p;
    if (!(p > 0)) continue;
    if (o.type === "atom") cls.push({ type: "atom", species: o.species, p, b: o.b, key: "atom:" + o.species });
    else if (o.type === "vacancy") cls.push({ type: "vacancy", p, key: "vacancy" });
    else if (o.type === "molecule") {
      const t = otab[o.group];
      if (!t) continue;
      const ws = t.weights.reduce((x, y) => x + y, 0) || 1;
      t.orientations.forEach((or, k) => {
        if (t.weights[k] > 0) cls.push({ type: "molecule", group: o.group, index: k, R: or.R, mol: t.mol, p: p * t.weights[k] / ws, b: o.b, key: `mol:${o.group}:${k}` });
      });
    }
  }
  return cls;
}

// ---------------------------------------------------------------------------------- realisation (preview)
function realize(state, opts = {}) {
  const n = opts.n || state.supercell.map(Number);
  const A = latticeMatrix(state.cell);
  const otab = orientationTable(state, A);
  const rnd = rng32(opts.seed != null ? opts.seed : +state.seed || 1);
  const classes = state.sites.map(s => siteClasses(s, otab));
  const ncell = n[0] * n[1] * n[2];
  const idx = (x, y, z) => x + n[0] * (y + n[1] * z);
  const choice = classes.map((cl, s) => {
    const arr = new Int32Array(ncell);
    const p = cl.map(c => c.p), tot = p.reduce((a, b) => a + b, 0);
    if (!cl.length) return arr;
    if (state.sites[s].draw === "exact") {
      const counts = p.map(x => Math.floor(x / tot * ncell));
      let rest = ncell - counts.reduce((a, b) => a + b, 0);
      const order = p.map((x, i) => [x / tot * ncell - counts[i] + rnd() * 1e-9, i]).sort((a, b) => b[0] - a[0]);
      for (let k = 0; k < rest; k++) counts[order[k][1]]++;
      let pos = 0; counts.forEach((c, i) => { for (let k = 0; k < c; k++) arr[pos++] = i; });
      for (let i = ncell - 1; i > 0; i--) { const j = rnd.int(i + 1); const t = arr[i]; arr[i] = arr[j]; arr[j] = t; }
    } else {
      for (let c = 0; c < ncell; c++) { let u = rnd() * tot, k = 0; while (k < p.length - 1 && u >= p[k]) { u -= p[k]; k++; } arr[c] = k; }
    }
    return arr;
  });
  // correlated occupancies (same model as the Python script; separate random stream)
  const corr = state.disorder.correlations;
  const report = [];
  if (corr && corr.enabled) {
    const inter = corr.interactions.filter(it => +it.J !== 0);
    const byGroup = corr.match === "group";
    const ids = {};
    const key = classes.map(cl => Int32Array.from(cl.map(c => { const k = byGroup && c.group ? c.group : c.key; if (!(k in ids)) ids[k] = Object.keys(ids).length; return ids[k]; })));
    const nb = state.sites.map(() => []);
    inter.forEach(it => { nb[+it.site_a].push([+it.site_b, it.t.map(Number), +it.J]); nb[+it.site_b].push([+it.site_a, it.t.map(x => -x), +it.J]); });
    const wrap = (v, m) => ((v % m) + m) % m;
    const cellOf = c => [c % n[0], Math.floor(c / n[0]) % n[1], Math.floor(c / (n[0] * n[1]))];
    const E = (s, c, cls) => {
      const [x, y, z] = cellOf(c); let e = 0; const k = key[s][cls];
      for (const [o, t, J] of nb[s]) if (key[o][choice[o][idx(wrap(x + t[0], n[0]), wrap(y + t[1], n[1]), wrap(z + t[2], n[2]))]] === k) e += J;
      return e;
    };
    const measure = () => inter.map(it => {
      const a = +it.site_a, b = +it.site_b, t = it.t.map(Number); let same = 0;
      for (let c = 0; c < ncell; c++) { const [x, y, z] = cellOf(c); if (key[a][choice[a][c]] === key[b][choice[b][idx(wrap(x + t[0], n[0]), wrap(y + t[1], n[1]), wrap(z + t[2], n[2]))]]) same++; }
      return same / ncell;
    });
    const before = measure(), kT = +corr.kT || 1;
    const sweeps = Math.min(+corr.sweeps || 0, opts.maxSweeps || 1e9);
    for (let sw = 0; sw < sweeps; sw++) for (let s = 0; s < classes.length; s++) {
      if (!nb[s].length || classes[s].length < 2) continue;
      for (let it = 0; it < ncell; it++) {
        const c1 = rnd.int(ncell), c2 = rnd.int(ncell), x1 = choice[s][c1], x2 = choice[s][c2];
        if (key[s][x1] === key[s][x2]) continue;
        const e0 = E(s, c1, x1) + E(s, c2, x2);
        choice[s][c1] = x2; choice[s][c2] = x1;
        const dE = E(s, c1, x2) + E(s, c2, x1) - e0;
        if (!(dE <= 0 || rnd() < Math.exp(-dE / kT))) { choice[s][c1] = x1; choice[s][c2] = x2; }
      }
    }
    const after = measure();
    inter.forEach((it, i) => report.push({ it, before: before[i], after: after[i] }));
  }
  // place atoms
  const rig = state.disorder.molecule_rigid || {}, tS = +rig.translation_sigma || 0, lS = +rig.libration_sigma_deg || 0;
  const disp = state.disorder.atomic_displacements || {};
  const atoms = [];
  let molCount = 0;
  for (let z = 0; z < n[2]; z++) for (let y = 0; y < n[1]; y++) for (let x = 0; x < n[0]; x++) {
    const c = idx(x, y, z);
    state.sites.forEach((site, s) => {
      const cl = classes[s][choice[s][c]];
      if (!cl) return;
      const centre = frac2cart(A, add([x, y, z], site.frac.map(Number)));
      if (cl.type === "atom") atoms.push({ el: cl.species, p: centre, cell: [x, y, z], site: s, cls: choice[s][c], mol: -1 });
      else if (cl.type === "molecule") {
        const loc = moleculeLocal(cl.mol);
        let R = cl.R;
        if (lS > 0) { const ax = [rnd.normal(), rnd.normal(), rnd.normal()]; R = mul(R, axisAngle(ax, rnd.normal() * lS)); }
        const sh = tS > 0 ? [rnd.normal() * tS, rnd.normal() * tS, rnd.normal() * tS] : [0, 0, 0];
        loc.xyz.forEach((p, j) => atoms.push({ el: loc.species[j], p: add(add(centre, sh), mv(R, p)), cell: [x, y, z], site: s, cls: choice[s][c], mol: molCount, orient: cl.index, group: cl.group }));
        molCount++;
      }
    });
  }
  for (const a of atoms) { const s = +disp[a.el] || 0; if (s > 0) a.p = add(a.p, [rnd.normal() * s, rnd.normal() * s, rnd.normal() * s]); }
  return { atoms, n, A, report, classes, nMolecules: molCount };
}

// ---------------------------------------------------------------------------------- diffuse preview (f = Z)
// I_diff(h) = (|Σ f e^{2πi h·x}|² − |Σ_cells e^{2πi h·R} ⟨F_cell⟩|²) / N_cells, x in unit-cell fractions.
function diffusePlaneRows(real, plane) {
  const { atoms, n, A } = real;
  const Ainv = inv(A), ncell = n[0] * n[1] * n[2];
  const loc = new Float64Array(atoms.length * 3), f = new Float64Array(atoms.length), cid = new Int32Array(atoms.length);
  atoms.forEach((a, i) => {
    const fr = [a.p[0] * Ainv[0][0] + a.p[1] * Ainv[1][0] + a.p[2] * Ainv[2][0], a.p[0] * Ainv[0][1] + a.p[1] * Ainv[1][1] + a.p[2] * Ainv[2][1], a.p[0] * Ainv[0][2] + a.p[1] * Ainv[1][2] + a.p[2] * Ainv[2][2]];
    loc[3 * i] = fr[0] - a.cell[0]; loc[3 * i + 1] = fr[1] - a.cell[1]; loc[3 * i + 2] = fr[2] - a.cell[2];
    f[i] = zOf(a.el); cid[i] = a.cell[0] + n[0] * (a.cell[1] + n[1] * a.cell[2]);
  });
  const step = 1 / Math.max(...n), m = Math.round(plane.range / step), t = [];
  for (let k = -m; k <= m; k++) t.push(k * step);
  const P = t.length, u = plane.u, v = plane.v;
  const Fr = new Float64Array(ncell * P), Fi = new Float64Array(ncell * P);
  const rows = [];
  function row(a) {
    Fr.fill(0); Fi.fill(0);
    const ta = t[a];
    for (let i = 0; i < atoms.length; i++) {
      const base = cid[i] * P, x = loc[3 * i], y = loc[3 * i + 1], z = loc[3 * i + 2];
      const p0 = ta * (u[0] * x + u[1] * y + u[2] * z), dp = v[0] * x + v[1] * y + v[2] * z;
      for (let b = 0; b < P; b++) { const ph = 2 * Math.PI * (p0 + t[b] * dp); Fr[base + b] += f[i] * Math.cos(ph); Fi[base + b] += f[i] * Math.sin(ph); }
    }
    const out = new Float64Array(P);
    for (let b = 0; b < P; b++) {
      const H = [ta * u[0] + t[b] * v[0], ta * u[1] + t[b] * v[1], ta * u[2] + t[b] * v[2]];
      let sr = 0, si = 0, mr = 0, mi = 0, lr = 0, li = 0;
      for (let c = 0; c < ncell; c++) {
        const cx = c % n[0], cy = Math.floor(c / n[0]) % n[1], cz = Math.floor(c / (n[0] * n[1]));
        const ph = 2 * Math.PI * (H[0] * cx + H[1] * cy + H[2] * cz), co = Math.cos(ph), sn = Math.sin(ph);
        const fr = Fr[c * P + b], fi = Fi[c * P + b];
        sr += fr * co - fi * sn; si += fr * sn + fi * co;
        mr += fr; mi += fi; lr += co; li += sn;
      }
      mr /= ncell; mi /= ncell;
      const ar = mr * lr - mi * li, ai = mr * li + mi * lr;
      out[b] = (sr * sr + si * si - ar * ar - ai * ai) / ncell;
    }
    return out;
  }
  return { t, P, row };
}

// ---------------------------------------------------------------------------------- validation
function validate(state) {
  const errs = [], warns = [];
  let A;
  try { A = latticeMatrix(state.cell); } catch (e) { errs.push(e.message); return { errs, warns }; }
  const n = state.supercell.map(Number);
  if (n.some(x => !(Number.isInteger(x) && x >= 1))) errs.push("Supercell needs three whole numbers ≥ 1.");
  if (!state.sites.length) errs.push("Add at least one site.");
  const otab = orientationTable(state, A);
  Object.values(otab).forEach(t => t.errors.forEach(e => errs.push(`Orientation group “${t.group.name}”: ${e}`)));
  state.molecules.forEach(m => {
    if (!m.atoms.length) errs.push(`Molecule “${m.name}” has no atoms.`);
    m.atoms.forEach((a, i) => { if (!elementOf(a.el)) errs.push(`Molecule “${m.name}”, atom ${i + 1}: “${a.el}” is not an element symbol.`); });
  });
  const labels = new Set();
  state.sites.forEach(s => {
    if (!/^[A-Za-z0-9_]{1,10}$/.test(s.label)) errs.push(`Site label “${s.label}”: use 1–10 letters, digits or _.`);
    if (labels.has(s.label)) errs.push(`Site label “${s.label}” is used twice.`);
    labels.add(s.label);
    const tot = s.occupants.reduce((a, o) => a + (+o.p || 0), 0);
    if (Math.abs(tot - 1) > 1e-6) errs.push(`Site ${s.label}: probabilities add up to ${tot.toFixed(4)}, not 1.`);
    s.occupants.forEach(o => {
      if (o.type === "atom") {
        if (!elementOf(o.species)) errs.push(`Site ${s.label}: “${o.species}” is not an element symbol.`);
        else if (o.species.length > 4) errs.push(`Site ${s.label}: species names longer than 4 characters break DISCUS.`);
      }
      if (o.type === "molecule") {
        const t = otab[o.group];
        if (!t) errs.push(`Site ${s.label}: pick an orientation group for the molecule.`);
        else if (!t.weights.some(w => w > 0)) errs.push(`Site ${s.label}: group “${t.group.name}” has no orientation with weight > 0.`);
        else if (t.orientations.some((x, k) => x.improper && t.weights[k] > 0)) warns.push(`Group “${t.group.name}” uses improper operations: those orientations are mirror images of the molecule.`);
      }
      if (+o.p < 0) errs.push(`Site ${s.label}: negative probability.`);
    });
  });
  const hasMol = state.sites.some(s => s.occupants.some(o => o.type === "molecule" && +o.p > 0));
  if (state.mmc.enabled && hasMol) errs.push("DISCUS Monte Carlo works on atom-only models; switch it off or remove molecules.");
  if (state.mmc.enabled) {
    const nfb = state.mmc.targets.filter(t => t.mode === "correlation" && (t.type === "occupational" || t.type === "displacement")).length;
    if (nfb > 1) errs.push("Monte Carlo: only one target may use a correlation (feedback); give the others energies.");
    if (n.some(x => x < 2)) errs.push("Monte Carlo needs at least 2 cells along every axis.");
    if (!state.mmc.targets.length) errs.push("Monte Carlo: add a target.");
  }
  const f3 = state.calc.fourier;
  if (f3.enabled) {
    f3.steps.forEach((s, k) => {
      const per = 1 / +s;
      if (!(+s > 0)) errs.push("3D grid steps must be positive.");
      else if (Math.abs(per - Math.round(per)) > 1e-6 || n[k] % Math.round(per) !== 0)
        warns.push(`3D grid step along ${"hkl"[k]} (${(+s).toPrecision(4)}) is not 1/m with m dividing the supercell (${n[k]}); the diffuse pattern will show finite-size ripples. Use 1/${n[k]} or 1/${n[k] / 2 || 1}.`);
    });
    const shape = f3.steps.map((s, k) => 2 * Math.round(+f3.hmax[k] / +s) + 1);
    const vox = shape.reduce((a, b) => a * b, 1);
    if (vox > 3e7) warns.push(`3D grid has ${(vox / 1e6).toFixed(0)} million points; the HDF5 file and the FFT will be large.`);
    if (f3.lots.enabled && f3.lots.size.some((s, k) => +s > n[k])) errs.push("LOTS larger than the crystal (DISCUS error −101).");
  }
  const p = state.calc.pdf;
  if (p.enabled) {
    const lam = +p.wavelength || (p.radiation === "neutron" ? 1.54 : 1.5406);
    if (+p.qmax * lam > 12) warns.push(`1D PDF: DISCUS 6.23.01 cannot use Q above ≈ ${(12.3 / lam).toFixed(1)} Å⁻¹ at λ = ${lam} Å; Qmax ${p.qmax} behaves like ${(12.3 / lam).toFixed(1)}.`);
    const r = +p.rmax / +p.dr;
    if (Math.abs(r - Math.round(r)) > 1e-7) errs.push("1D PDF: rmax must be a whole number of dr steps.");
    const L = Math.min(...n.map((x, k) => x * norm(A[k])));
    if (+p.rmax > L / 2) warns.push(`1D PDF: rmax ${p.rmax} Å exceeds half the smallest crystal dimension (${(L / 2).toFixed(1)} Å); DISCUS's Debye sum treats the supercell as a finite particle.`);
  }
  if (state.calc.debye_waller) {
    const bad = state.sites.some(s => (+s.b_iso || 0) < 0.001) || state.sites.some(s => s.occupants.some(o => o.b != null && o.b !== "" && +o.b < 0.001));
    if (bad) errs.push("Debye–Waller on: every B must be ≥ 0.001 Å² (DISCUS replaces smaller values by 1).");
  }
  return { errs, warns };
}

// ---------------------------------------------------------------------------------- export
function buildConfig(state) {
  const A = latticeMatrix(state.cell);
  const otab = orientationTable(state, A);
  const molecules = {}, orientations = {};
  state.molecules.forEach(m => { molecules[m.id] = { name: m.name, anchor: m.anchor, atoms: m.atoms.map(a => ({ el: a.el, xyz: a.xyz.map(Number) })) }; });
  Object.entries(otab).forEach(([gid, t]) => {
    orientations[gid] = {
      molecule: t.mol.id, name: t.group.name,
      matrices: t.orientations.map(o => o.R.map(r => r.map(x => +x.toFixed(12)))),
      weights: t.weights,
      definition: { convention: t.group.convention, euler_deg: t.group.euler.map(Number), ops: t.group.ops, frame: t.group.frame, closure: !!t.group.closure },
    };
  });
  const num = v => (v === "" || v == null ? null : +v);
  return {
    schema: "discuspdf-disorder-builder", version: 1, title: state.title, seed: +state.seed,
    output_dir: state.output_dir || "discus_model_out",
    cell: Object.fromEntries(["a", "b", "c", "alpha", "beta", "gamma"].map(k => [k, +state.cell[k]])),
    supercell: state.supercell.map(Number), molecules, orientations,
    sites: state.sites.map(s => ({
      label: s.label, frac: s.frac.map(Number), b_iso: +s.b_iso || 1, draw: s.draw,
      occupants: s.occupants.map(o => o.type === "atom" ? { type: "atom", species: o.species, p: +o.p, b: num(o.b) }
        : o.type === "molecule" ? { type: "molecule", group: o.group, p: +o.p, b: num(o.b) } : { type: "vacancy", p: +o.p }),
    })),
    disorder: {
      atomic_displacements: Object.fromEntries(Object.entries(state.disorder.atomic_displacements).filter(([, v]) => +v > 0).map(([k, v]) => [k, +v])),
      molecule_rigid: { translation_sigma: +state.disorder.molecule_rigid.translation_sigma || 0, libration_sigma_deg: +state.disorder.molecule_rigid.libration_sigma_deg || 0 },
      clash_threshold: +state.disorder.clash_threshold || 0.8,
      correlations: { ...state.disorder.correlations, sweeps: +state.disorder.correlations.sweeps, kT: +state.disorder.correlations.kT,
        interactions: state.disorder.correlations.interactions.map(i => ({ site_a: +i.site_a, site_b: +i.site_b, t: i.t.map(Number), J: +i.J })) },
    },
    mmc: state.mmc.enabled ? JSON.parse(JSON.stringify(state.mmc)) : { enabled: false },
    calc: JSON.parse(JSON.stringify(state.calc)),
    builder_state: JSON.parse(JSON.stringify(state)),
  };
}
function generateScript(state, template) {
  const cfg = buildConfig(state);
  const json = JSON.stringify(cfg, null, 1);
  if (!template.includes('CONFIG_JSON = r"""__CONFIG__"""')) throw new Error("script template is damaged");
  return template.replace('CONFIG_JSON = r"""__CONFIG__"""', () => 'CONFIG_JSON = r"""\n' + json + '\n"""');
}

const api = { SYMBOLS, Z, elementOf, zOf, colorOf, covOf, I3, mul, tr, mv, det, inv, add, sub, scale, norm, latticeMatrix,
  frac2cart, axisAngle, eulerMatrix, parseOp, closure, moleculeLocal, bondsOf, samePlacement, expandGroup, rng32,
  orientationTable, siteClasses, realize, diffusePlaneRows, validate, buildConfig, generateScript, isOrthogonal };
if (typeof module !== "undefined" && module.exports) module.exports = api; else root.DB = api;
})(typeof window !== "undefined" ? window : globalThis);
