/* Small orthographic ball-and-stick renderer on a 2D canvas (no WebGL needed). */
(function (root) {
"use strict";
const { mul, tr, mv, axisAngle, sub, add, scale } = root.DB;

function cssVar(name, fallback) {
  try { return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback; } catch (e) { return fallback; }
}
function shade(hex, f) {
  const c = toRgb(hex);
  if (!c) return hex;
  const m = f >= 0 ? c.map(x => Math.round(x + (255 - x) * f)) : c.map(x => Math.round(x * (1 + f)));
  return `rgb(${m[0]},${m[1]},${m[2]})`;
}
const rgbCache = {};
function toRgb(col) {
  if (rgbCache[col]) return rgbCache[col];
  let out = null;
  if (/^#([0-9a-f]{6})$/i.test(col)) out = [0, 2, 4].map(i => parseInt(col.slice(1 + i, 3 + i), 16));
  else {
    try {
      const cv = document.createElement("canvas"); cv.width = cv.height = 1;
      const g = cv.getContext("2d"); g.fillStyle = col; g.fillRect(0, 0, 1, 1);
      const d = g.getImageData(0, 0, 1, 1).data; out = [d[0], d[1], d[2]];
    } catch (e) { out = [128, 128, 128]; }
  }
  rgbCache[col] = out;
  return out;
}

class Viewer {
  constructor(canvas, opts = {}) {
    this.cv = canvas; this.opts = opts;
    this.R = mul(axisAngle([1, 0, 0], -62), axisAngle([0, 0, 1], -28));
    this.zoom = 1; this.scene = { atoms: [], bonds: [], lines: [] };
    this.pan = [0, 0];
    if (!opts.static) this._bind();
    this.resize();
    if (typeof ResizeObserver !== "undefined" && !opts.static) new ResizeObserver(() => { this.resize(); this.draw(); }).observe(canvas);
  }
  resize() {
    const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
    const w = this.cv.clientWidth || this.cv.width, h = this.opts.aspect ? w * this.opts.aspect : (this.cv.clientHeight || this.cv.height);
    this.cv.width = Math.max(10, Math.round(w * dpr)); this.cv.height = Math.max(10, Math.round(h * dpr));
    if (this.opts.aspect) this.cv.style.height = h + "px";
    this.dpr = dpr;
  }
  set(scene, { keepView = true } = {}) {
    this.scene = scene;
    const P = scene.atoms.map(a => a.p).concat(scene.lines.flatMap(l => [l[0], l[1]]));
    if (!P.length) { this.center = [0, 0, 0]; this.extent = 2; }
    else {
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      P.forEach(p => { for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); } });
      this.center = scene.center || scale(add(lo, hi), 0.5);
      this.extent = Math.max(1.2, ...P.map(p => Math.hypot(...sub(p, this.center)))) + Math.max(0, ...scene.atoms.map(a => a.r || 0.3)) + 0.15;
    }
    if (!keepView) { this.zoom = 1; this.pan = [0, 0]; }
    this.draw();
  }
  scaleFactor() { return Math.min(this.cv.width, this.cv.height) * 0.46 / this.extent * this.zoom; }
  project(p) {
    const q = mv(this.R, sub(p, this.center)), s = this.scaleFactor();
    return [this.cv.width / 2 + this.pan[0] + s * q[0], this.cv.height / 2 + this.pan[1] - s * q[1], q[2]];
  }
  unproject(x, y) { // screen point -> world point in the view plane through the centre
    const s = this.scaleFactor();
    const q = [(x * this.dpr - this.cv.width / 2 - this.pan[0]) / s, -(y * this.dpr - this.cv.height / 2 - this.pan[1]) / s, 0];
    return add(mv(tr(this.R), q), this.center);
  }
  draw() {
    const g = this.cv.getContext("2d"), W = this.cv.width, H = this.cv.height;
    g.fillStyle = cssVar("--view", "#f7f9f9"); g.fillRect(0, 0, W, H);
    const s = this.scaleFactor(), ink = cssVar("--view-ink", "#223");
    const lineCol = cssVar("--muted", "#789");
    g.lineCap = "round";
    for (const l of this.scene.lines) {
      const a = this.project(l[0]), b = this.project(l[1]);
      g.strokeStyle = l[2] || lineCol; g.globalAlpha = l[3] || 0.55; g.lineWidth = Math.max(1, this.dpr * (l[4] || 1));
      g.beginPath(); g.moveTo(a[0], a[1]); g.lineTo(b[0], b[1]); g.stroke();
    }
    g.globalAlpha = 1;
    const P = this.scene.atoms.map(a => this.project(a.p));
    const order = P.map((p, i) => i).sort((i, j) => P[i][2] - P[j][2]);
    const bondsBy = {};
    for (const [i, j] of this.scene.bonds) { (bondsBy[i] = bondsBy[i] || []).push(j); (bondsBy[j] = bondsBy[j] || []).push(i); }
    const zmin = Math.min(...P.map(p => p[2]), 0), zmax = Math.max(...P.map(p => p[2]), 0), zr = (zmax - zmin) || 1;
    const fog = this.scene.atoms.length > 300;
    for (const i of order) {
      const a = this.scene.atoms[i], p = P[i];
      for (const j of bondsBy[i] || []) {
        if (P[j][2] > p[2] || (P[j][2] === p[2] && j > i)) continue;
        const m = [(p[0] + P[j][0]) / 2, (p[1] + P[j][1]) / 2];
        g.strokeStyle = shade(this.scene.atoms[i].color, -0.15); g.lineWidth = Math.max(1.5, s * 0.16);
        g.beginPath(); g.moveTo(p[0], p[1]); g.lineTo(m[0], m[1]); g.stroke();
      }
      const r = Math.max(1.2 * this.dpr, s * (a.r || 0.3));
      g.globalAlpha = fog ? 0.45 + 0.55 * (p[2] - zmin) / zr : 1;
      const grd = g.createRadialGradient(p[0] - r * .35, p[1] - r * .35, r * .1, p[0], p[1], r);
      grd.addColorStop(0, shade(a.color, 0.55)); grd.addColorStop(0.55, a.color); grd.addColorStop(1, shade(a.color, -0.35));
      g.fillStyle = grd; g.beginPath(); g.arc(p[0], p[1], r, 0, 2 * Math.PI); g.fill();
      if (a.ring) { g.strokeStyle = a.ring; g.lineWidth = 2 * this.dpr; g.stroke(); }
      for (const j of bondsBy[i] || []) {
        if (!(P[j][2] > p[2] || (P[j][2] === p[2] && j > i))) continue;
        const m = [(p[0] + P[j][0]) / 2, (p[1] + P[j][1]) / 2];
        g.strokeStyle = shade(a.color, -0.15); g.lineWidth = Math.max(1.5, s * 0.16);
        g.beginPath(); g.moveTo(p[0], p[1]); g.lineTo(m[0], m[1]); g.stroke();
      }
      if (a.label) { g.globalAlpha = 1; g.fillStyle = ink; g.font = `${11 * this.dpr}px ui-monospace, monospace`; g.fillText(a.label, p[0] + r + 2, p[1] - r); }
    }
    g.globalAlpha = 1;
    if (this.scene.axes) this._axes(g);
  }
  _axes(g) {
    const L = 22 * this.dpr, o = [26 * this.dpr, this.cv.height - 24 * this.dpr];
    const cols = ["#d0533a", "#3f9a5a", "#3d6fd6"];
    this.scene.axes.forEach((v, k) => {
      const n = Math.hypot(...v); const q = mv(this.R, scale(v, 1 / n));
      g.strokeStyle = cols[k]; g.lineWidth = 2 * this.dpr; g.beginPath(); g.moveTo(o[0], o[1]); g.lineTo(o[0] + L * q[0], o[1] - L * q[1]); g.stroke();
      g.fillStyle = cols[k]; g.font = `${11 * this.dpr}px ui-monospace, monospace`; g.fillText("abc"[k], o[0] + (L + 6) * q[0] - 3, o[1] - (L + 6) * q[1] + 4);
    });
  }
  _bind() {
    let drag = null;
    this.cv.addEventListener("pointerdown", e => {
      if (this.opts.onClick && this.opts.placing && this.opts.placing()) {
        const r = this.cv.getBoundingClientRect();
        this.opts.onClick(this.unproject(e.clientX - r.left, e.clientY - r.top), e);
        return;
      }
      drag = { x: e.clientX, y: e.clientY, shift: e.shiftKey }; this.cv.setPointerCapture(e.pointerId);
    });
    this.cv.addEventListener("pointermove", e => {
      if (!drag) return;
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY;
      if (drag.shift) { this.pan = [this.pan[0] + dx * this.dpr, this.pan[1] + dy * this.dpr]; }
      else this.R = mul(mul(axisAngle([0, 1, 0], dx * 0.6), axisAngle([1, 0, 0], dy * 0.6)), this.R);
      this.draw();
    });
    const end = () => { drag = null; };
    this.cv.addEventListener("pointerup", end); this.cv.addEventListener("pointercancel", end);
    this.cv.addEventListener("wheel", e => { e.preventDefault(); this.zoom *= Math.exp(-e.deltaY * 0.0015); this.zoom = Math.min(20, Math.max(0.2, this.zoom)); this.draw(); }, { passive: false });
  }
  viewAlong(v, up) { // look down v (Cartesian), with up roughly 'up'
    const z = scale(v, 1 / Math.hypot(...v));
    let x = [up[1] * z[2] - up[2] * z[1], up[2] * z[0] - up[0] * z[2], up[0] * z[1] - up[1] * z[0]];
    if (Math.hypot(...x) < 1e-6) x = [1, 0, 0];
    x = scale(x, 1 / Math.hypot(...x));
    const y = [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]];
    this.R = [x, y, z]; this.pan = [0, 0]; this.draw();
  }
}
root.Viewer = Viewer;
root.toRgb = toRgb;
})(typeof window !== "undefined" ? window : globalThis);
