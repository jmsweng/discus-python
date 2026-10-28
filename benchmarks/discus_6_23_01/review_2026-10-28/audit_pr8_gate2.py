"""Read-only audit of the PR #8 Gate 2 evidence (reproduces every number in ../REVIEW_2026-10-28.md).

Usage: python audit_pr8_gate2.py > audit_output.txt
Reads only committed files; writes nothing else. Needs numpy and h5py.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

B6 = Path(__file__).resolve().parents[1]
WK = json.loads((B6 / "gate2_prereg" / "wk_coefficients.json").read_text())["elements"]
IT92 = {"GA": ((15.2354, 6.7006, 4.3591, 2.9623), (3.0669, 0.2412, 10.7805, 61.4135), 1.7189),
        "AS": ((16.6723, 6.0701, 3.4313, 4.2779), (2.6345, 0.2647, 12.9479, 47.7972), 2.531)}
RUNS = {"corr": ("correlated_run/corr_model.stru", "correlated_run/corr_total.h5",
                 "correlated_run/corr_diffuse.h5", "correlated_run/corr_total_3dpdf.h5"),
        "sr": ("short_range_run/short_range_model.stru", "short_range_run/sr_total.h5", None,
               "short_range_run/sr_total_3dpdf.h5")}


def f_of(el, s, table):
    if table == "it92":
        a, b, c = IT92[el]
    else:
        w = WK[el]["discus" if table == "wk-discus" else "dabax"]; a, b, c = w["a"], w["b"], w["c"]
    return c + sum(ai * np.exp(-bi * s * s) for ai, bi in zip(a, b))


def stru(p):
    sp, x, B, cell = [], [], [], None
    for line in (B6 / p).read_text().splitlines():
        t = line.replace(",", " ").split()
        if t and t[0].lower() == "cell":
            cell = [float(v) for v in t[1:7]]
        if t and t[0].upper() in ("GA", "AS"):
            sp.append(t[0].upper()); x.append([float(v) for v in t[1:4]]); B.append(float(t[4]))
    return np.array(sp), np.array(x), np.array(B), cell


def h5(p):
    with h5py.File(B6 / p) as f:
        d = f["entry/data"]
        return {k: d[k][()] for k in d if d[k].ndim <= 3}


def txt(v):
    return v.decode() if isinstance(v, bytes) else v


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


print("## 1. Inventory of committed engine output")
for tag, (_, tot, dif, pdf) in RUNS.items():
    for p in (tot, dif, pdf):
        if p is None:
            print(f"{tag}: no diffuse (set aver,100) run exists"); continue
        d = h5(p); v = d["data_values"]
        print(f"{p}: sha256 {hashlib.sha256((B6 / p).read_bytes()).hexdigest()} shape {v.shape} "
              f"min {v.min():.4g} max {v.max():.4g} bragg={txt(d['data_type_with_bragg'])} "
              f"axes={txt(d['data_type_axes'])} corner={d['data_corner'].tolist()} "
              f"step={np.diag(d['data_increment_vector']).tolist()}")

print("\n## 2. Independent forward model of TOTAL intensity (finite atom list, B on F)")
for tag, (sp_path, tot, _, _) in RUNS.items():
    sp, x, B, cell = stru(sp_path); a = cell[0]
    d = h5(tot); v = d["data_values"]
    ijk = np.indices(v.shape).reshape(3, -1).T.astype(float)
    hkl = d["data_corner"] + ijk @ d["data_increment_vector"].T
    s_exact = np.linalg.norm(hkl / a, axis=1) / 2
    for table, quant in (("it92", False), ("wk", False), ("wk", True), ("wk-discus", True)):
        s = np.floor(s_exact / 0.001 + 0.5) * 0.001 if quant else s_exact
        F = np.zeros(len(hkl), complex)
        for el in ("GA", "AS"):
            m = sp == el
            F += f_of(el, s, table) * np.exp(-B[m][0] * s * s) * np.exp(2j * np.pi * hkl @ x[m].T).sum(1)
        print(f"{tag}: {table:9s} s_rounded_to_0.001={quant!s:5s} relL2 = {rel((np.abs(F)**2).reshape(v.shape), v):.3e}")
    print(f"{tag}: d*_max = {np.max(np.linalg.norm(hkl, axis=1)) / a:.4f} 1/A (d_min {a / np.max(np.linalg.norm(hkl, axis=1)):.2f} A)")

print("\n## 3. The committed sign check (compare_with_engine.py) is a tautology on the wrong line")
preds = {"corr": [-1, -1, 1, 1, -1, -1], "sr": [-1, -1, 1, 1, 1, 1]}
for tag, (_, tot, _, _) in RUNS.items():
    d = h5(tot); v = d["data_values"]; n = v.shape[0]; a = float(d["unit_cell_lengths"][0]); dh = 0.125
    pdf = np.fft.ifftn(v).real; c = n // 2; du = a / (n * dh)
    got = []
    for ax in range(3):
        line = [pdf[:, c, c], pdf[c, :, c], pdf[c, c, :]][ax]   # exactly the committed script's lines
        got += [np.sign(line[np.argmax(line)]), np.sign(line[np.argmin(line)])]
    match = sum(int(g == p) for g, p in zip(got, preds[tag]))
    print(f"{tag}: '+dir' := argmax (always >= 0), '-dir' := argmin (always <= 0) -> match {match}/6; "
          f"any +-symmetric prediction scores exactly 3/6")
    print(f"{tag}: raw ifftn index {c} on the other two axes is u = {c * du:.2f} A, not 0; "
          f"real-space spacing {du:.3f} A; nearest grid points to 2.5 A: 0 and {du:.2f} A")
    h = -2 + dh * np.arange(n)
    for r in (0.0, 2.5, 5.0, 10.0):
        print(f"{tag}: 1D window response at {r:4.1f} A / origin = {np.mean(np.cos(2 * np.pi * h * r / a)):+.3f}")
    for tau in ((2.5, 0, 0), (0, 2.5, 0), (0, 0, 2.5), (0, 0, 0)):
        u = np.array(tau) / a
        E = [np.exp(-2j * np.pi * u[k] * h) for k in range(3)]
        rho = np.einsum("i,j,k,ijk->", *E, v).real
        print(f"{tag}: exact windowed transform at tau={tau} A: {rho:.6e}")
    un = np.fft.ifftshift(v); opp = un[np.ix_(*[(-np.arange(m)) % m for m in v.shape])]
    print(f"{tag}: Friedel error {np.linalg.norm(un - opp) / np.linalg.norm(un):.1e} "
          f"(log: 'Friedel pairs merged' -> imposed by the engine, uninformative)")

print("\n## 4. Pair correlations at the predicted displacements: open vs periodic boundaries")
for tag, (sp_path, *_rest) in RUNS.items():
    sp, x, _, cell = stru(sp_path); a = cell[0]; L = int(round(a / 2.5))
    s = np.where(sp == "GA", 1, -1)
    key = {tuple(np.round(xi * a / 2.5).astype(int)): i for i, xi in enumerate(x)}
    for tau in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (4, 0, 0)):
        op, pe = [], []
        for g, i in key.items():
            g2 = tuple(np.add(g, tau)); g2p = tuple(np.mod(g2, L))
            if g2p in key:
                pe.append(s[i] * s[key[g2p]])
                if g2 == g2p:
                    op.append(s[i] * s[key[g2]])
        print(f"{tag}: tau={tuple(2.5 * np.array(tau))} A open n={len(op)} <ss>={np.mean(op):+.4f} | "
              f"periodic n={len(pe)} <ss>={np.mean(pe):+.4f} (SE ~ {1 / np.sqrt(len(pe)):.3f})")

print("\n## 5. Why set aver,100 returned exactly zero for corr_model.stru")
sp, x, _, _ = stru(RUNS["corr"][0])
ext = x.max(0) - x.min(0)
print(f"no 'ncell' line; DISCUS sets cr_icc = INT(extent + 1) = {np.floor(ext + 1).astype(int).tolist()} "
      f"and cr_ncatoms = {len(sp)} (structur.f90 do_readstru_disk) -> one unit cell, <F> = F, diffuse = 0")

print("\n## 6. DISCUS native 3DPDF vs direct transform of the total volume (Gate 3 note)")
for tag, (_, tot, _, pdf) in RUNS.items():
    I = h5(tot)["data_values"]; dP = h5(pdf); P = dP["data_values"]
    h = -2 + 0.125 * np.arange(33); u = dP["data_corner"][0] + dP["data_increment_vector"][0, 0] * np.arange(33)
    E = np.exp(-2j * np.pi * np.outer(u, h))
    D = np.einsum("ai,bj,ck,ijk->abc", E, E, E, I, optimize=True).real
    (sc, off), *_ = np.linalg.lstsq(np.c_[D.ravel(), np.ones(D.size)], P.ravel(), rcond=None)
    print(f"{tag}: native = {sc:.6e} * direct + {off:.4e}; residual {np.linalg.norm(sc * D + off - P) / np.linalg.norm(P - P.mean()):.2e} of variation")
