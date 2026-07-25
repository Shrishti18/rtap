#!/usr/bin/env python3
"""Scan the local ICSD CIF cache for MX4 tetrahedra and classify them."""
import csv, os, re, sys, warnings
warnings.filterwarnings("ignore")

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xtal

CIFDIR = "/home/user/rtap/data/priority3_cifs/cif"
META = "/home/user/rtap/data/priority2_flatbands/flatband_icsd_metadata.csv"
EXTRA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extra_cif")

PRIORITY = [("Co", {"As", "P"}, 1), ("Co", {"S", "Se"}, 1),
            ("Ni", {"As", "P", "C", "N"}, 2), ("Cu", {"As", "P", "C", "N"}, 2),
            ("Fe", {"As", "P", "S", "Se"}, 3),
            ("Mn", {"P", "As"}, 4)]
DIMNAME = {0: "0D (isolated)", 1: "1D chain/ladder", 2: "2D layer", 3: "3D network"}


def analyse(path, icsd, formula, spg, spgsym, src):
    try:
        cell, ops, sites = xtal.read_cif(path)
        A = xtal.cart_matrix(cell)
        atoms = xtal.expand(cell, ops, sites)
    except Exception as e:
        return [("ERR", "%s: %s" % (type(e).__name__, e))]
    els = {a["el"] for a in atoms}
    out = []
    for M, XS, prio in PRIORITY:
        if M not in els:
            continue
        X = XS & els
        if not X:
            continue
        midx = [i for i, a in enumerate(atoms) if a["el"] == M]
        rej = {}
        bonded = {(M, x) for x in X} | {(x, M) for x in X}
        adj = None
        seenkey = set()
        for i in midx:
            tet, shape = xtal.tetrahedron(atoms, A, i, X)
            if tet is None:
                rej.setdefault((M, "/".join(sorted(X))), shape)
                continue
            ang = xtal.angles(tet)
            gr = xtal.group_angles(ang)
            d = sorted(round(r, 3) for r, _, _, _ in tet)
            key = (tuple(round(a, 1) for a in ang), tuple(d))
            if key in seenkey:
                continue
            seenkey.add(key)
            if adj is None:
                adj = xtal.build_adjacency(atoms, A, bonded)
                dim, nsite = xtal.dimensionality(atoms, A, midx, X, adj=adj)
            xel = sorted({atoms[j]["el"] for _, j, _, _ in tet})
            out.append(dict(
                icsd=icsd, source=src, formula=formula, spg=spg, spg_symbol=spgsym,
                M=M, X="/".join(xel), priority=prio,
                dim=dim, motif=DIMNAME.get(dim, "?"),
                n_M_sites=len(midx),
                angles_raw=" ".join("%.2f" % a for a in ang),
                angles_grouped="; ".join("%.2f x%d" % (v, n) for v, n in gr),
                angle_min=round(min(ang), 2), angle_max=round(max(ang), 2),
                spread=round(max(ang) - min(ang), 2),
                d_MX=" ".join("%.3f" % x for x in d),
                d_MX_mean=round(float(np.mean(d)), 3),
                regime=xtal.regime(gr)))
        if not any(o["M"] == M for o in out) and rej:
            out.append(dict(icsd=icsd, source=src, formula=formula, spg=spg,
                            spg_symbol=spgsym, M=M, X="/".join(sorted(X)),
                            priority=prio, dim=-1, motif="EXCLUDED",
                            n_M_sites=len(midx), angles_raw="", angles_grouped="",
                            angle_min="", angle_max="", spread="", d_MX="",
                            d_MX_mean="", regime="not-tetrahedral: " + list(rej.values())[0]))
    return out


def main():
    meta = {r["icsd"]: r for r in csv.DictReader(open(META))}
    jobs = []
    for icsd, m in meta.items():
        e = set(re.findall(r"([A-Z][a-z]?)", m["formula"] or ""))
        if any(M in e and (XS & e) for M, XS, _ in PRIORITY):
            jobs.append((os.path.join(CIFDIR, "%s.cif" % icsd), icsd,
                         m["formula"], m["space_group_number"],
                         m["space_group_symbol"], "ICSD (local cache)"))
    if os.path.isdir(EXTRA):
        for fn in sorted(os.listdir(EXTRA)):
            if fn.endswith(".cif"):
                jobs.append((os.path.join(EXTRA, fn), fn[:-4], "", "", "", "COD"))
    print("structures to analyse: %d" % len(jobs), flush=True)
    rows, errs = [], []
    for k, j in enumerate(jobs, 1):
        for r in analyse(*j):
            if isinstance(r, tuple):
                errs.append((j[1], r[1]))
            else:
                rows.append(r)
        if k % 100 == 0:
            print("  %d/%d  tetrahedra so far: %d" % (k, len(jobs), len(rows)), flush=True)
    order = {"SPLIT-A": 0, "SPLIT-B": 1, "REGULAR": 2, "OTHER": 3}
    rows.sort(key=lambda r: (r["priority"], order.get(r["regime"], 9),
                             r["dim"], r["formula"]))
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mx4_tetrahedra.csv")
    with open(p, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("-> %s (%d tetrahedra, %d structures, %d errors)"
          % (p, len(rows), len({r["icsd"] for r in rows}), len(errs)))
    print("MX4 SCAN DONE", flush=True)


if __name__ == "__main__":
    main()
