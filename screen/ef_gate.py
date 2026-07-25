#!/usr/bin/env python3
"""Fermi-level gate over the Phase-A candidates.

A band isolated by several eV in a wide-gap insulator is completely full or
completely empty, and an unfillable band cannot superconduct at any Tc. The
Phase-A screen never checked that E_F lies inside the isolated group.

For each material carrying a group:
  * E_F is read from the OUTCAR shipped in the material's own archive, so it is
    on the same absolute scale as the _hr.dat eigenvalues.
  * a group passes only if min(E_group) < E_F < max(E_group), strictly.
  * per band, the filling f = fraction of the BZ with E_band(k) < E_F is
    recorded; f outside [0.05, 0.95] is discarded, since the theory needs
    f(1-f) finite and it peaks at f = 1/2.

Only eigenVALUES are needed, so this is far cheaper than the Phase-A sweep.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import csv, re, sys, zipfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
EF_RE = re.compile(r"E-fermi\s*:\s*([-\d.]+)")


def one(item):
    import fastops
    jid, dim, nk, groups = item
    try:
        with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
            nm = next(n for n in z.namelist() if n.endswith("_hr.dat"))
            R, H, deg = fastops.read_hr_bytes(z.read(nm))
            ef = None
            if "OUTCAR" in z.namelist():
                m = EF_RE.findall(z.read("OUTCAR").decode("utf8", "replace"))
                if m:
                    ef = float(m[-1])
        if ef is None:
            return jid, None, "no E-fermi in OUTCAR"
        Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
        E = np.linalg.eigvalsh(Hk)
        del Hk
        Ef = E.reshape(-1, E.shape[-1])
        out = {}
        for g in groups:
            bs = [int(x) for x in g.split("-")]
            sub = Ef[:, bs]
            gmin, gmax = float(sub.min()), float(sub.max())
            rec = dict(efermi=ef, group_Emin=gmin, group_Emax=gmax,
                       ef_inside=bool(gmin < ef < gmax))
            for b in bs:
                rec["fill_%d" % b] = float((Ef[:, b] < ef).mean())
                rec["Emin_%d" % b] = float(Ef[:, b].min())
                rec["Emax_%d" % b] = float(Ef[:, b].max())
            out[g] = rec
        return jid, out, None
    except Exception as e:
        return jid, None, "%s: %s" % (type(e).__name__, e)


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "full_bands.csv"))))
    jobs = {}
    for r in rows:
        jobs.setdefault(r["jid"], (r["jid"], int(r["dim"]), int(r["nk"]), set()))[3].add(r["group"])
    items = [(j, d, nk, sorted(gs)) for (j, d, nk, gs) in jobs.values()]
    print("materials to gate: %d" % len(items), flush=True)
    res, errs = {}, []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for i, (jid, out, err) in enumerate(ex.map(one, items, chunksize=1), 1):
            if err:
                errs.append((jid, err))
            else:
                res[jid] = out
            if i % 25 == 0:
                print("  %d/%d" % (i, len(items)), flush=True)

    keep = []
    for r in rows:
        g = res.get(r["jid"], {}).get(r["group"])
        if g is None:
            r["efermi"] = ""
            r["ef_inside"] = ""
            r["filling"] = ""
            r["ef_status"] = "no_efermi"
            keep.append(r)
            continue
        b = int(r["band"])
        f = g["fill_%d" % b]
        r["efermi"] = g["efermi"]
        r["group_Emin"] = g["group_Emin"]
        r["group_Emax"] = g["group_Emax"]
        r["ef_inside"] = g["ef_inside"]
        r["filling"] = f
        r["band_Emin"] = g["Emin_%d" % b]
        r["band_Emax"] = g["Emax_%d" % b]
        ok_ef = g["ef_inside"]
        ok_f = 0.05 <= f <= 0.95
        r["ef_status"] = ("PASS" if (ok_ef and ok_f) else
                          "fail_ef_outside" if not ok_ef else "fail_filling")
        r["lamM"] = float(r["lam"]) * float(r["M_trg"])
        keep.append(r)

    p = os.path.join(HERE, "full_bands_ef.csv")
    cols = list(keep[0].keys())
    for r in keep:
        for c in cols:
            r.setdefault(c, "")
    with open(p, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(keep)
    print("-> %s" % p)
    if errs:
        print("errors: %d  %s" % (len(errs), errs[:5]))
    print("EF GATE DONE", flush=True)


if __name__ == "__main__":
    main()
