#!/usr/bin/env python3
"""Phase B on the Phase-A top bands.

Two things Phase A could not settle:
  * M_trg at nk=24 is a LOWER bound (~2% low on a smooth band, much worse near
    a touching). Re-run on a fine grid.
  * d_iso at nk=24 is an UPPER bound -- a coarse grid can miss the touching
    point entirely, which is the unsafe direction. Re-run on the same fine grid.
  * M_naive is only an upper bound on the gauge-invariant metric. Run the
    minimal metric with the .wout centres seeded (both signs).

The fine pass never forms the whole (nk^dim, norb, norb) array: it walks k3
slice by slice and keeps only the eigenvector of the bands of interest, so
memory is nk^3 * norb * 16 bytes rather than nk^3 * norb^2 * 16.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "4"

import csv, sys, time, zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harvest, fastops

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
SMALL = "/home/user/rtap/data/priority1_wannier_jarvis/extracted_small"
RS = [1.00, 0.85, 0.55]
NK_FINE_MAX, NK_COARSE = 64, 24


def fine_nk(norb, dim):
    nk = NK_FINE_MAX
    while nk > 24 and (nk ** dim) * norb * 16 > 1.2e9:
        nk -= 8
    return nk


def fine_pass(R, H, deg, dim, nk, bands):
    """Per-band eigenvector on a fine grid plus the full band extrema,
    computed slice by slice so the projector is never materialised."""
    n = H.shape[1]
    Hf = (H / deg[:, None, None]).reshape(H.shape[0], n * n)
    ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    lo = np.full(n, np.inf)
    hi = np.full(n, -np.inf)
    if dim == 2:
        m = np.meshgrid(ks, ks, indexing="ij")
        kf = np.stack([x.ravel() for x in m], 1)
        Hk = (np.exp(1j * (kf @ R[:, :2].T.astype(float))) @ Hf).reshape(-1, n, n)
        E, V = np.linalg.eigh(Hk)
        lo, hi = E.min(0), E.max(0)
        U = {b: V[:, :, b].reshape(nk, nk, n) for b in bands}
        return lo, hi, U, [2 * np.pi / nk] * 2
    m1, m2 = np.meshgrid(ks, ks, indexing="ij")
    base = np.stack([m1.ravel(), m2.ravel()], 1)
    U = {b: np.empty((nk, nk, nk, n), complex) for b in bands}
    for i3, k3 in enumerate(ks):
        kf = np.concatenate([base, np.full((base.shape[0], 1), k3)], 1)
        Hk = (np.exp(1j * (kf @ R[:, :3].T.astype(float))) @ Hf).reshape(-1, n, n)
        E, V = np.linalg.eigh(Hk)
        lo = np.minimum(lo, E.min(0))
        hi = np.maximum(hi, E.max(0))
        for b in bands:
            U[b][:, :, i3, :] = V[:, :, b].reshape(nk, nk, n)
    return lo, hi, U, [2 * np.pi / nk] * 3


def frac_centres(jid, dim, norb):
    p = os.path.join(ZIPS, "%s.zip" % jid)
    with zipfile.ZipFile(p) as z:
        if "wannier90.wout" not in z.namelist():
            return None
        txt = z.read("wannier90.wout").decode("utf8", "replace")
    import re
    blocks = re.findall(r"Final State(.*?)(?:Sum of centres|\Z)", txt, re.S)
    if not blocks:
        return None
    rows = re.findall(r"WF centre and spread\s+\d+\s*\(\s*([-\d.]+)\s*,\s*"
                      r"([-\d.]+)\s*,\s*([-\d.]+)\s*\)", blocks[-1])
    if not rows or len(rows) != norb:
        return None
    c = np.array([[float(a) for a in r] for r in rows])
    L = open(os.path.join(SMALL, jid, "POSCAR")).read().splitlines()
    s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    return (c @ np.linalg.inv(A))[:, :dim]


def main():
    work = list(csv.DictReader(open(os.path.join(HERE, "phaseB_worklist.csv"))))
    bygroup = {}
    for r in work:
        bygroup.setdefault((r["jid"], r["group"]), []).append(r)
    print("Phase B: %d bands in %d groups on %d materials"
          % (len(work), len(bygroup), len({k[0] for k in bygroup})), flush=True)
    out = []
    for i, ((jid, group), rs) in enumerate(sorted(bygroup.items()), 1):
        t0 = time.time()
        dim = int(rs[0]["dim"])
        bands = sorted(int(r["band"]) for r in rs)
        gb = [int(x) for x in group.split("-")]
        try:
            with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
                nm = next(n for n in z.namelist() if n.endswith("_hr.dat"))
                R, H, deg = fastops.read_hr_bytes(z.read(nm))
            norb = H.shape[1]
            nk = fine_nk(norb, dim)
            lo, hi, U, dks = fine_pass(R, H, deg, dim, nk, bands)
            g0, g1 = gb[0], gb[-1]
            gl = lo[g0] - hi[g0 - 1] if g0 > 0 else np.inf
            gu = lo[g1 + 1] - hi[g1] if g1 + 1 < len(lo) else np.inf
            d_iso_fine = float(min(gl, gu))
            # coarse grid for the (expensive) minimal-metric optimisation
            Hc, dksc = fastops.hk_grid_fast(R, H, deg, NK_COARSE, dim=dim)
            Ec, Vc = np.linalg.eigh(Hc)
            del Hc
            centres = frac_centres(jid, dim, norb)
            for r in rs:
                b = int(r["band"])
                m_fine = fastops.trg_from_u(U[b], dks)
                w_fine = float(hi[b] - lo[b])
                mm = fastops.minimal_metric_fast(None, dksc, b, dim,
                                                 centres=centres, restarts=2,
                                                 u=Vc[..., :, b], maxfev=30000)
                rec = dict(r)
                rec.update(nk_fine=nk, M_trg_fine=m_fine,
                           M_trg_coarse=float(r["M_trg"]),
                           M_ratio_fine_coarse=m_fine / max(float(r["M_trg"]), 1e-12),
                           W_band_fine=w_fine, d_iso_fine=d_iso_fine,
                           d_iso_coarse=float(r["d_iso"]),
                           M_naive_c=mm["M_naive"], M_min_c=mm["M_min"],
                           gauge_ratio=(mm["M_naive"] / mm["M_min"]) if mm["M_min"] > 0 else np.inf)
                # best estimate: fine metric scaled by the gauge correction,
                # with the fine (smaller, safer) isolation
                gauge = (mm["M_min"] / mm["M_naive"]) if mm["M_naive"] > 0 else 1.0
                d = dict(M_trg=m_fine * gauge, M=0.0, lam=float(r["lam"]),
                         d_iso=d_iso_fine)
                for Rv in RS:
                    tc = harvest.tc_max(d, dim=dim, R=Rv)
                    rec["TcB_R%.2f" % Rv] = tc["Tc_max_K"]
                    rec["bindsB_R%.2f" % Rv] = tc["binds"]
                    rec["capB_R%.2f" % Rv] = tc["U_capped_by"]
                rec["M_trg_best"] = m_fine * gauge
                out.append(rec)
            del Ec, Vc, U
            print("[%2d/%d] %-13s %-10s group %-7s nk_fine=%-3d  %.0fs"
                  % (i, len(bygroup), jid, rs[0]["formula"], group, nk,
                     time.time() - t0), flush=True)
        except Exception as e:
            print("[%2d/%d] %-13s ERROR %s" % (i, len(bygroup), jid, e), flush=True)
    if out:
        p = os.path.join(HERE, "phaseB_results.csv")
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print("-> %s (%d bands)" % (p, len(out)))
    print("PHASE B DONE", flush=True)


if __name__ == "__main__":
    main()
