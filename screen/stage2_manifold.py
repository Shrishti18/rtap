#!/usr/bin/env python3
"""Stage 2 (v5) — two-band manifolds, not single isolated bands.

auto_bands() required one band gapped from both neighbours, which an n_phi=2
doublet can never satisfy: a 2D site irrep forces the pair degenerate at Gamma.
This driver uses harvest.manifold_bands / descriptors_manifold instead.

Run in two phases, because the expensive steps must not be spent on every
manifold:

  PHASE A (this script, default): every material, every nb-band manifold,
    descriptors + tc_max on the NAIVE metric. No Powell, no fine grid.
    Since M_min <= M_naive, Tc_max(naive) is an upper bound on Tc_max(M_min)
    whenever the band is stiffness-bound, and exact when it is amplitude-bound.
    Ranking on it therefore cannot hide a winner.

  PHASE B (--refine): the globally top-ranked manifolds get the
    gauge-invariant minimal metric and the fine-grid d_iso re-check.
"""
import csv, os, sys, time, traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harvest, fastops

HERE = os.path.dirname(os.path.abspath(__file__))
HR = os.path.join(HERE, "hr")
SMALL = "/home/user/rtap/data/priority1_wannier_jarvis/extracted_small"

NB = 2
WMAX, ISO_MIN = 1.0, 0.05
R_VALUES = [1.00, 0.85, 0.55]
NK_WANT = 40
MEM = 8.0e9


def poscar_lattice(jid):
    L = open(os.path.join(SMALL, jid, "POSCAR")).read().splitlines()
    s = float(L[1].split()[0])
    return np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s


def frac_centres(jid, dim):
    p = os.path.join(HR, "%s.wout" % jid)
    if not os.path.exists(p):
        return None
    c = harvest.read_wout_centres(p)
    if c is None:
        return None
    return (c @ np.linalg.inv(poscar_lattice(jid)))[:, :dim]


def pick_nk(norb, dim, want=NK_WANT):
    nk = want
    while nk > 12 and 2.5 * (nk ** dim) * norb * norb * 16 > MEM:
        nk -= 4
    return nk


def manifolds_from_E(E, nb=NB, wmax=WMAX, iso_min=ISO_MIN):
    """harvest.manifold_bands, reusing an eigendecomposition already done."""
    nbnd = E.shape[-1]
    out = []
    for b in range(nbnd - nb + 1):
        sub = E[..., b:b + nb]
        lo, hi = sub.min(), sub.max()
        gl = lo - E[..., b - 1].max() if b > 0 else np.inf
        gu = E[..., b + nb].min() - hi if b + nb < nbnd else np.inf
        iso = min(gl, gu)
        if (hi - lo) <= wmax and iso >= iso_min:
            out.append((tuple(range(b, b + nb)), float(hi - lo), float(iso)))
    return out


def phase_a():
    recs = list(csv.DictReader(open(os.path.join(HERE, "stage2_matched.csv"))))
    byj = {}
    for r in recs:
        byj.setdefault(r["jid"], r)
    jids = sorted(byj, key=lambda s: int(s.split("-")[1]))
    rows, errs = [], []
    for i, jid in enumerate(jids, 1):
        rec = byj[jid]
        dim = 2 if rec["dimensionality"] == "2D" else 3
        try:
            t0 = time.time()
            R, H, deg = harvest.read_hr(os.path.join(HR, "%s_hr.dat" % jid))
            norb = H.shape[1]
            nk = pick_nk(norb, dim)
            Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
            E, V = np.linalg.eigh(Hk)
            del Hk
            mans = manifolds_from_E(E)
            for bands, W, iso in mans:
                d = fastops.descriptors_manifold_fast(E, V, dks, bands)
                sp = harvest.spec(d, Tc_K=300, dim=dim, R=1.0)
                row = dict(jid=jid, formula=rec["jarvis_formula"],
                           spg=rec["jarvis_spg"], dim=dim, norb=norb, nk=nk,
                           bands="-".join(map(str, bands)), b0=bands[0],
                           nb=len(bands), M_trg=d["M_trg"], lam=d["lam"],
                           nphi=d["nphi"], W=d["W"], d_iso=d["d_iso"],
                           unif=d["unif"], Emid=d["Emid"],
                           U_req=sp["U_req"], iso_req=sp["iso_req"],
                           pass_iso=sp["pass_iso"])
                for Rv in R_VALUES:
                    tc = harvest.tc_max(d, dim=dim, R=Rv)
                    row["Tc_max_K_R%.2f" % Rv] = tc["Tc_max_K"]
                    row["binds_R%.2f" % Rv] = tc["binds"]
                row["U_used"] = harvest.tc_max(d, dim=dim, R=1.0)["U_used"]
                rows.append(row)
            del E, V
            best = max([r["Tc_max_K_R1.00"] for r in rows if r["jid"] == jid],
                       default=0.0)
            print("[%2d/%d] %-13s %-10s dim=%d norb=%-3d nk=%-3d manifolds=%-3d "
                  "bestTc=%7.1f K  %.0fs"
                  % (i, len(jids), jid, rec["jarvis_formula"], dim, norb, nk,
                     len(mans), best, time.time() - t0), flush=True)
        except Exception as e:
            errs.append((jid, "%s: %s" % (type(e).__name__, e)))
            print("[%2d/%d] %-13s ERROR %s" % (i, len(jids), jid, e), flush=True)
            traceback.print_exc()
    if rows:
        p = os.path.join(HERE, "stage2_manifolds.csv")
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("-> %s (%d manifold rows)" % (p, len(rows)))
    if errs:
        with open(os.path.join(HERE, "stage2_manifold_errors.csv"), "w", newline="") as f:
            csv.writer(f).writerows([["jid", "error"]] + errs)
    print("PHASE A DONE", flush=True)


if __name__ == "__main__":
    phase_a()
