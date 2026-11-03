"""Monte Carlo on the occupant choices of a disorder-builder model.

Two energies, both optional: a Potts-like neighbour energy J along given cell vectors (correlated
occupants), and penalties for atom pairs of neighbouring occupants closer than a rule's minimum distance
(forbidden contacts, evaluated on the ideal atom positions). Moves swap the occupants of two cells of one
site (composition kept) and, with contacts, redraw an orientation within its molecule group. Moved from the
disorder builder's ``engine.py`` (builder @ 92fe94e) unchanged, so the same model and seed give the same
configuration. Works for molecules and vacancies, which DISCUS ``mmc`` (:mod:`discuspdf.mc`) cannot handle.
"""
from __future__ import annotations

import math

import numpy as np

from ..atom_names import element_of

__all__ = ["contact_model", "ordering_mc"]


def _rule_hits(r: str, sp: str) -> bool:
    """A contact rule names an atom name or an element: 'Ni' covers Ni and Ni62, 'H' covers H and D."""
    return r == "any" or str(r).upper() == str(sp).upper() or r == element_of(sp)


def contact_model(cfg: dict, classes, A: np.ndarray, mols: dict) -> dict:
    """Penalty and contact-count tables for every pair of neighbouring site positions (ideal atom positions,
    no translation / libration / displacement noise). Same construction as the web page."""
    con = cfg.get("disorder", {}).get("contacts") or {}
    rules = [(r["a"], r["b"], float(r["dmin"]) ** 2, float(r["penalty"])) for r in con.get("rules", [])
             if float(r.get("dmin", 0)) > 0 and float(r.get("penalty", 0)) > 0]
    sites = cfg["sites"]
    entries = [[] for _ in sites]
    pairs = []
    if not con.get("enabled") or not rules:
        return {"entries": entries, "pairs": pairs, "rules": rules}

    def shape(c):
        if c["type"] == "atom":
            return [c["species"]], np.zeros((1, 3))
        if c["type"] == "molecule":
            sp, loc = mols[c["molecule"]]
            return sp, loc @ np.asarray(c["R"]).T
        return [], np.zeros((0, 3))

    shapes = [[shape(c) for c in cl] for cl in classes]
    R = [max([0.0] + [float(np.linalg.norm(p, axis=1).max()) for _, p in sh if len(p)]) for sh in shapes]
    dmax = math.sqrt(max(r[2] for r in rules))
    recip = np.linalg.norm(np.linalg.inv(A), axis=0)
    match = lambda r, e1, e2: (_rule_hits(r[0], e1) and _rule_hits(r[1], e2)) or (_rule_hits(r[0], e2) and _rule_hits(r[1], e1))
    for i, s1 in enumerate(sites):
        for j, s2 in enumerate(sites):
            cut = R[i] + R[j] + dmax
            m = [int(math.ceil(cut * r)) + 1 for r in recip]
            for tx in range(-m[0], m[0] + 1):
                for ty in range(-m[1], m[1] + 1):
                    for tz in range(-m[2], m[2] + 1):
                        if i == j and not (tx or ty or tz):
                            continue
                        off = (np.array([tx, ty, tz], float) + np.asarray(s2["frac"], float) - np.asarray(s1["frac"], float)) @ A
                        dist = float(np.linalg.norm(off))
                        if dist > cut + 1e-9:
                            continue
                        K, L = len(classes[i]), len(classes[j])
                        table = np.zeros((K, L))
                        count = np.zeros((K, L), dtype=np.int64)
                        for k in range(K):
                            ea, pa = shapes[i][k]
                            for l in range(L):
                                eb, pb = shapes[j][l]
                                if not len(pa) or not len(pb):
                                    continue
                                d2 = (((off + pb)[None, :, :] - pa[:, None, :]) ** 2).sum(axis=2)
                                e = 0.0
                                hits = np.zeros(d2.shape, dtype=bool)
                                for r in rules:
                                    mk = np.array([[match(r, ea[u], eb[v]) for v in range(len(eb))] for u in range(len(ea))])
                                    h = mk & (d2 < r[2])
                                    e += r[3] * int(h.sum())
                                    hits |= h
                                table[k, l] = e
                                count[k, l] = int(hits.sum())
                        if table.any():
                            entries[i].append({"o": j, "t": (tx, ty, tz), "table": table, "count": count})
                            if i <= j:
                                pairs.append({"i": i, "j": j, "t": (tx, ty, tz), "dist": dist,
                                              "n_forbidden_pairs": int((table > 0).sum()), "n_pairs": int(table.size)})
    return {"entries": entries, "pairs": pairs, "rules": rules}


def ordering_mc(cfg: dict, classes, choice: np.ndarray, A: np.ndarray, mols: dict, rng, log) -> dict:
    """Monte Carlo on the occupant choices (same model as the web page):
    - neighbour energies J added when two listed neighbours hold the same occupant (Correlated occupants), and
    - contact penalties for atom pairs of neighbouring occupants closer than a rule's minimum distance
      (Forbidden contacts).
    Moves: swaps of occupants between cells of one site (composition kept) and, with contacts, a new orientation
    drawn with the weights from the same molecule group (change = "orientation") or from every molecule of the
    site (change = "group"); change = "none" keeps swaps only. Metropolis acceptance; optional annealing from
    10 kT down to kT over the first 60 % of the sweeps."""
    dis = cfg.get("disorder", {})
    corr = dis.get("correlations") or {}
    con = dis.get("contacts") or {}
    potts = [i for i in corr.get("interactions", []) if float(i.get("J", 0.0)) != 0.0] if corr.get("enabled") else []
    cm = contact_model(cfg, classes, A, mols)
    use_contacts = bool(cm["pairs"])
    if not potts and not use_contacts:
        return {}
    n = choice.shape[1:]
    nx, ny, nz = (int(x) for x in n)
    ncell = nx * ny * nz
    flat = [list(choice[s].transpose(2, 1, 0).ravel()) for s in range(len(classes))]   # cell = x + nx (y + ny z)
    by_group = corr.get("match") == "group"
    ids: dict = {}
    key = [[ids.setdefault(c.get("group", c["key"]) if by_group else c["key"], len(ids)) for c in cl] for cl in classes]
    change = con.get("change") or ("none" if con.get("reorient") is False else "orientation")
    # occupants a run may turn into each other: orientations of one group, or every molecule of the site
    fam = [[(("mol" if change == "group" else "mol:" + c["group"]) if c["type"] == "molecule" else c["key"]) for c in cl]
           for cl in classes]
    maps: dict = {}

    def map_of(t):
        t = tuple(int(x) for x in t)
        if t not in maps:
            x = np.arange(nx)[None, None, :]
            y = np.arange(ny)[None, :, None]
            z = np.arange(nz)[:, None, None]
            idx = ((x + t[0]) % nx) + nx * (((y + t[1]) % ny) + ny * ((z + t[2]) % nz))
            maps[t] = list(np.broadcast_to(idx, (nz, ny, nx)).ravel())
        return maps[t]

    nb = [[] for _ in classes]
    for it in potts:
        a, b, t, J = int(it["site_a"]), int(it["site_b"]), [int(x) for x in it["t"]], float(it["J"])
        nb[a].append((b, map_of(t), J, None))
        nb[b].append((a, map_of([-x for x in t]), J, None))
    for i, lst in enumerate(cm["entries"]):
        for e in lst:
            nb[i].append((e["o"], map_of(e["t"]), 0.0, e["table"].tolist()))

    def term(s, e, k, k2):
        o, _, J, table = e
        if table is not None:
            return table[k][k2]
        return J if key[s][k] == key[o][k2] else 0.0

    def local(s, c, k, excl):
        E = 0.0
        for e in nb[s]:
            c2 = e[1][c]
            if e[0] == s and c2 == excl:
                continue
            E += term(s, e, k, flat[e[0]][c2])
        return E

    def bond(s, c1, k1, c2, k2):
        return sum(term(s, e, k1, k2) for e in nb[s] if e[0] == s and e[1][c1] == c2)

    def measure():
        out = []
        for it in potts:
            a, b, mp = int(it["site_a"]), int(it["site_b"]), map_of(it["t"])
            out.append(sum(key[a][flat[a][c]] == key[b][flat[b][mp[c]]] for c in range(ncell)) / ncell)
        return out

    def n_contacts():
        tot = 0
        for i, lst in enumerate(cm["entries"]):
            for e in lst:
                mp, cnt = map_of(e["t"]), e["count"]
                tot += int(sum(cnt[flat[i][c], flat[e["o"]][mp[c]]] for c in range(ncell)))
        return tot // 2

    before, c_before = measure(), (n_contacts() if use_contacts else 0)
    kT = float(corr.get("kT", 1.0))
    sweeps = int(corr.get("sweeps", 50))
    reorient = use_contacts and change != "none"
    anneal = use_contacts and con.get("anneal", True)
    fam_idx = []
    for s, cl in enumerate(classes):
        g: dict = {}
        for k in range(len(cl)):
            g.setdefault(fam[s][k], []).append(k)
        fam_idx.append([g[fam[s][k]] for k in range(len(cl))])
    rnd = rng.random
    for sw in range(sweeps):
        T = kT * 10 ** max(0.0, 1 - sw / max(1.0, 0.6 * sweeps)) if anneal else kT
        for s, cls in enumerate(classes):
            if not nb[s] or len(cls) < 2:
                continue
            fl = flat[s]
            for _ in range(ncell):
                if reorient:
                    c = int(rnd() * ncell)
                    x = fl[c]
                    opts = fam_idx[s][x]
                    if len(opts) > 1:
                        u = rnd() * sum(cls[k]["p"] for k in opts)
                        y = opts[0]
                        for k in opts:
                            u -= cls[k]["p"]
                            if u <= 0:
                                y = k
                                break
                        if y != x:
                            dE = local(s, c, y, -1) - local(s, c, x, -1)
                            if dE <= 0 or rnd() < math.exp(-dE / T):
                                fl[c] = y
                c1, c2 = int(rnd() * ncell), int(rnd() * ncell)
                x1, x2 = fl[c1], fl[c2]
                if c1 == c2 or (fam[s][x1] == fam[s][x2] and (reorient or key[s][x1] == key[s][x2])):
                    continue
                if not use_contacts and key[s][x1] == key[s][x2]:
                    continue
                e0 = local(s, c1, x1, c2) + local(s, c2, x2, c1) + bond(s, c1, x1, c2, x2)
                fl[c1], fl[c2] = x2, x1
                dE = local(s, c1, x2, c2) + local(s, c2, x1, c1) + bond(s, c1, x2, c2, x1) - e0
                if not (dE <= 0 or rnd() < math.exp(-dE / T)):
                    fl[c1], fl[c2] = x1, x2
    for s in range(len(classes)):
        choice[s] = np.array(flat[s]).reshape(nz, ny, nx).transpose(2, 1, 0)
    after = measure()
    report = {"sweeps": sweeps, "kT": kT, "interactions": potts,
              "fraction_same_before": before, "fraction_same_after": after}
    for it, b0, b1 in zip(potts, before, after):
        log(f"  correlation site {it['site_a']}→{it['site_b']} t={tuple(it['t'])} J={it['J']}: "
            f"P(same) {b0:.3f} → {b1:.3f}")
    if use_contacts:
        c_after = n_contacts()
        report["contacts"] = {"rules": [{"a": r[0], "b": r[1], "dmin": math.sqrt(r[2]), "penalty": r[3]} for r in cm["rules"]],
                              "neighbour_pairs": cm["pairs"], "before": c_before, "after": c_after,
                              "change": change, "anneal": anneal}
        log(f"  forbidden contacts (ideal positions): {c_before} → {c_after} over {len(cm['pairs'])} neighbour vectors")
        for s, cl in enumerate(classes):
            h = np.bincount(np.asarray(flat[s]), minlength=len(cl)) / ncell
            groups: dict = {}
            for k, c in enumerate(cl):
                groups[c.get("group", c["key"])] = groups.get(c.get("group", c["key"]), 0.0) + h[k]
            log(f"  site {cfg['sites'][s]['label']}: " + ", ".join(f"{g} {v:.3f}" for g, v in groups.items()))
    return report
