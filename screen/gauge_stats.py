#!/usr/bin/env python3
"""M_naive vs M_min over the bands that still matter after the E_F gate.

The original Phase B was re-checking the Phase-A top 20 on an nk=64 grid, but
the E_F gate eliminated all of those, so that work no longer bears on any
conclusion. What is still wanted is the gauge statistic -- how optimistic
M_naive is in general -- and that needs only the coarse grid, since it is a
ratio of two quantities computed on the SAME grid.

Sampled here: every band whose enclosing group lies within 0.5 eV of E_F (the
only bands with any physical claim), plus the single E_F-gate survivor.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import csv, sys, time, zipfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
SMALL = "/home/user/rtap/data/priority1_wannier_jarvis/extracted_small"
NEAR_EV = 0.5
MAX_NORB = 60      # Powell cost scales with (norb-1)*dim parameters
MAXFEV = 3000      # bounded refinement from the .wout start, not a full search


def frac_centres(jid, dim, norb, z):
    import re
    if "wannier90.wout" not in z.namelist():
        return None
    txt = z.read("wannier90.wout").decode("utf8", "replace")
    blocks = re.findall(r"Final State(.*?)(?:Sum of centres|\Z)", txt, re.S)
    if not blocks:
        return None
    rows = re.findall(r"WF centre and spread\s+\d+\s*\(\s*([-\d.]+)\s*,\s*"
                      r"([-\d.]+)\s*,\s*([-\d.]+)\s*\)", blocks[-1])
    if not rows or len(rows) != norb:
        return None
    c = np.array([[float(a) for a in r] for r in rows])
    p = os.path.join(SMALL, jid, "POSCAR")
    if not os.path.exists(p):
        return None
    L = open(p).read().splitlines()
    s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    return (c @ np.linalg.inv(A))[:, :dim]


def one(item):
    import fastops
    jid, dim, nk, bands, meta = item
    try:
        with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
            nm = next(n for n in z.namelist() if n.endswith("_hr.dat"))
            R, H, deg = fastops.read_hr_bytes(z.read(nm))
            norb = H.shape[1]
            centres = frac_centres(jid, dim, norb, z)
        Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
        E, V = np.linalg.eigh(Hk)
        del Hk
        out = []
        for b in bands:
            mm = fastops.minimal_metric_fast(None, dks, b, dim, centres=centres,
                                             restarts=0, u=V[..., :, b],
                                             maxfev=MAXFEV)
            out.append(dict(meta[b], M_naive=mm["M_naive"], M_min=mm["M_min"],
                            gauge_ratio=(mm["M_naive"] / mm["M_min"])
                            if mm["M_min"] > 1e-12 else np.inf,
                            had_centres=bool(centres is not None)))
        del E, V
        return jid, out, None
    except Exception as e:
        return jid, [], "%s: %s" % (type(e).__name__, e)


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(HERE, "full_bands_ef.csv")))
            if r["efermi"] != ""]
    sel = []
    for r in rows:
        ef, lo, hi = float(r["efermi"]), float(r["group_Emin"]), float(r["group_Emax"])
        dist = 0.0 if lo < ef < hi else (lo - ef if ef < lo else ef - hi)
        if (dist <= NEAR_EV or r["ef_status"] == "PASS") and int(r["norb"]) <= MAX_NORB:
            r["ef_distance"] = dist
            sel.append(r)
    print("bands near E_F (<= %.1f eV, norb <= %d): %d on %d materials"
          % (NEAR_EV, MAX_NORB, len(sel), len({r["jid"] for r in sel})), flush=True)
    jobs = {}
    for r in sel:
        k = r["jid"]
        jobs.setdefault(k, [int(r["dim"]), int(r["nk"]), [], {}])
        jobs[k][2].append(int(r["band"]))
        jobs[k][3][int(r["band"])] = r
    items = [(j, d, nk, sorted(set(bs)), meta) for j, (d, nk, bs, meta) in jobs.items()]
    out, errs = [], []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for i, (jid, rs, err) in enumerate(ex.map(one, items, chunksize=1), 1):
            if err:
                errs.append((jid, err))
            out.extend(rs)
            if i % 4 == 0:
                print("  %d/%d materials" % (i, len(items)), flush=True)
    if out:
        p = os.path.join(HERE, "gauge_stats.csv")
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()), extrasaction="ignore")
            w.writeheader()
            w.writerows(out)
        print("-> %s (%d bands)" % (p, len(out)))
    if errs:
        print("errors: %d %s" % (len(errs), errs[:3]))
    print("GAUGE STATS DONE", flush=True)


if __name__ == "__main__":
    main()
