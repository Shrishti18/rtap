#!/usr/bin/env python3
"""Phase A over ALL 1,772 JARVIS Wannier materials — the flat-band catalogue is
ignored for this pass, so the selection criterion is ours.

The 69.5 GB problem does not need batching: each wannier90_hr.dat is read
straight out of its zip into memory, used, and dropped. Nothing is written to
disk, so peak disk cost is zero and peak memory is one material per worker.

manifold_bands is used ONLY TO LOCATE isolated groups. The metric is then
taken PER BAND, because when a group's bands jointly span a fixed orbital
subspace the projector is a subspace identity and M_manifold vanishes
identically even though each band has a healthy metric (verified: 2.3e-14 for
the group against 1.62 per band). M_manifold is kept as a diagnostic column
only.

nk is deliberately coarse (24 in 3D, ~2% low for a smooth band). M is
underestimated on a coarse grid but the ORDERING is stable, and Phase B
re-runs the survivors at nk=64.
"""
import os
# BLAS reads these at import time, so they must be set before numpy loads;
# each worker gets one thread and parallelism comes from the process pool.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import csv, sys, time, traceback, zipfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
MANIFEST = "/home/user/rtap/data/priority1_wannier_jarvis/jarvis_wtb_manifest.csv"

NB = 2
WMAX, ISO_MIN = 0.5, 0.05      # wmax now gates the WIDEST INDIVIDUAL BAND
R_VALUES = [1.00, 0.85, 0.55]
NK3D, NK2D = 24, 40
MEM_PER_WORKER = 1.6e9


def pick_nk(norb, dim):
    nk = NK2D if dim == 2 else NK3D
    while nk > 8 and 2.5 * (nk ** dim) * norb * norb * 16 > MEM_PER_WORKER:
        nk -= 4
    return nk


def manifolds_from_E(E, nb=NB, wmax=WMAX, iso_min=ISO_MIN):
    """harvest.manifold_bands (v7) reusing an eigendecomposition already done.
    Gates the widest INDIVIDUAL band, not the manifold span."""
    nbnd = E.shape[-1]
    widths = [float(E[..., j].max() - E[..., j].min()) for j in range(nbnd)]
    out = []
    for b in range(nbnd - nb + 1):
        sub = E[..., b:b + nb]
        lo, hi = sub.min(), sub.max()
        gl = lo - E[..., b - 1].max() if b > 0 else np.inf
        gu = E[..., b + nb].min() - hi if b + nb < nbnd else np.inf
        iso = min(gl, gu)
        wb = max(widths[b:b + nb])
        if wb <= wmax and iso >= iso_min:
            out.append((tuple(range(b, b + nb)), wb, float(iso), float(hi - lo)))
    return out


def one(rec):
    import harvest, fastops
    jid = rec["jid"]
    dim = 2 if rec["dimensionality"] == "2D" else 3
    norb = int(rec["num_wann"] or 0)
    t0 = time.time()
    try:
        with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
            name = next(n for n in z.namelist() if n.endswith("_hr.dat"))
            R, H, deg = fastops.read_hr_bytes(z.read(name))
        norb = H.shape[1]
        nk = pick_nk(norb, dim)
        Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
        # Locate first with eigenVALUES only; most materials have no isolated
        # group at all, and eigvalsh is far cheaper than a full eigh.
        E = np.linalg.eigvalsh(Hk)
        mans = manifolds_from_E(E)
        if not mans:
            del Hk, E
            return jid, [], None, time.time() - t0, norb, nk
        E, V = np.linalg.eigh(Hk)
        del Hk
        rows = []
        for bands, wb, iso, wspan in mans:
            dm = fastops.descriptors_manifold_fast(E, V, dks, bands)
            for b in bands:
                d = fastops.descriptors_fast(E, V, dks, b)
                # isolation of the ENCLOSING group, not of the single band
                d_iso_group = iso
                dd = dict(d, d_iso=d_iso_group)
                row = dict(jid=jid, formula=rec["formula"], spg=rec["spg_number"],
                           dim=dim, norb=norb, nk=nk,
                           group="-".join(map(str, bands)), band=b,
                           rank=len(bands),
                           W_band=d["W"], W_group_max=wb, W_span=wspan,
                           d_iso=d_iso_group,
                           M_trg=d["M_trg"], lam=d["lam"], nphi=d["nphi"],
                           unif=d["unif"], Emid=d["Emid"],
                           M_manifold_diag=dm["M_trg"],
                           lam_manifold_diag=dm["lam"])
                for Rv in R_VALUES:
                    tc = harvest.tc_max(dd, dim=dim, R=Rv)
                    row["Tc_R%.2f" % Rv] = tc["Tc_max_K"]
                    row["binds_R%.2f" % Rv] = tc["binds"]
                    row["cap_R%.2f" % Rv] = tc["U_capped_by"]
                    row["U_R%.2f" % Rv] = tc["U_used"]
                rows.append(row)
        del E, V
        return jid, rows, None, time.time() - t0, norb, nk
    except Exception as e:
        return jid, [], "%s: %s" % (type(e).__name__, e), time.time() - t0, norb, 0


def main():
    recs = list(csv.DictReader(open(MANIFEST)))
    todo = sys.argv[1:]
    if todo:
        recs = [r for r in recs if r["jid"] in set(todo)]
    print("materials to scan: %d" % len(recs), flush=True)
    rows, errs, done, withman = [], [], 0, 0
    nw = max(1, min(4, (os.cpu_count() or 2)))
    out_p = os.path.join(HERE, "full_bands.csv")
    fh = None
    wr = None
    with ProcessPoolExecutor(max_workers=nw) as ex:
        for jid, rs, err, dt, norb, nk in ex.map(one, recs, chunksize=1):
            done += 1
            if err:
                errs.append((jid, err))
            if rs:
                withman += 1
                rows.extend(rs)
                if fh is None:      # checkpoint as we go
                    fh = open(out_p, "w", newline="", encoding="utf8")
                    wr = csv.DictWriter(fh, fieldnames=list(rs[0].keys()))
                    wr.writeheader()
                wr.writerows(rs)
                fh.flush()
            if done % 25 == 0 or err:
                best = max([r["Tc_R1.00"] for r in rows], default=0.0)
                print("  %4d/%d  %-13s norb=%-3d nk=%-3d man=%-2d  "
                      "materials_with_manifold=%d  bestTc=%.1f K  %s"
                      % (done, len(recs), jid, norb, nk, len(rs), withman, best,
                         err or ""), flush=True)
    if fh:
        fh.close()
        print("-> %s (%d band rows, %d materials)"
              % (out_p, len(rows), len({r["jid"] for r in rows})))
    with open(os.path.join(HERE, "full_scan_errors.csv"), "w", newline="") as f:
        csv.writer(f).writerows([["jid", "error"]] + errs)
    print("materials scanned: %d | with >=1 manifold: %d | errors: %d"
          % (done, withman, len(errs)), flush=True)
    print("FULL SCAN DONE", flush=True)


if __name__ == "__main__":
    main()
