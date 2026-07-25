"""
Redo the abundance claim under the self-consistency rule.

THE CLAIM UNDER TEST: "geometry is abundant -- median M_trg = 0.95, 545/728
bands above the 0.339 threshold." That is SINGLE-BAND M throughout, and it has
been cited in every strategic decision since, including the geometry-first
reframe.

THE RULE (now standing, not a per-case repair). A quantity must be computed on
a projector whose rank matches the degeneracy structure of what it projects
onto, and which rank that is depends on Tc, not on convention:

  * rank-r MANIFOLD treatment is valid iff the intra-group splitting delta is
    SMALL compared with Tc -- the bands then pair coherently as one manifold.
      U <= d_iso/2  (group isolation),  M = M_manifold,  lam = lam_manifold
  * SINGLE-BAND treatment is valid iff delta is LARGE compared with Tc -- the
    partner is then spectrally inert.
      U <= delta_nearest/2 (the band's OWN nearest gap),  M = M_trg,  lam = lam

Using one treatment's U with the other's M is what inflated LiPbAu2 from 9 K to
186 K, and it is the form the D1 error took in the primary data product.

delta_nearest and delta_intra are recomputed here from the hr.dat files as
MINIMA over the BZ; M, lam, M_manifold, lam_manifold and d_iso are taken from
the stored full-scan.
"""
import os, sys, csv, time, zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastops

ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
HERE = os.path.dirname(os.path.abspath(__file__))
KB = 0.086173e-3
C3 = 0.55 * 0.67


def pick_nk(norb, dim):
    nk = 40 if dim == 2 else 24
    while nk > 8 and 2.5 * (nk ** dim) * norb * norb * 16 > 1.6e9:
        nk -= 4
    return nk


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "full_bands.csv"))))
    seen = {}
    for r in rows:
        seen[(r["jid"], r["band"])] = r
    rows = list(seen.values())
    M0 = np.array([float(r["M_trg"]) for r in rows])
    print(f"STORED (single-band) over {len(rows)} unique bands:")
    print(f"   median M_trg = {np.median(M0):.4f}   "
          f"above 0.339: {(M0 > 0.339).sum()}/{len(M0)} "
          f"({100*(M0>0.339).mean():.1f}%)\n")

    byj = {}
    for r in rows:
        byj.setdefault((r["jid"], int(r["dim"])), []).append(r)
    groups = {}
    for r in rows:
        groups.setdefault((r["jid"], r["group"]), []).append(int(r["band"]))

    out, t0, nfail = [], time.time(), 0
    for i, ((jid, dim), rs) in enumerate(sorted(byj.items())):
        try:
            with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
                nm = next(n for n in z.namelist() if n.endswith("_hr.dat"))
                R, H, deg = fastops.read_hr_bytes(z.read(nm))
            nk = pick_nk(H.shape[1], dim)
            Hk, _ = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
            E = np.linalg.eigvalsh(Hk).reshape(-1, H.shape[1])
            del Hk
        except Exception as e:
            nfail += 1
            print(f"  FAIL {jid}: {type(e).__name__}: {e}", flush=True)
            continue
        nb = E.shape[1]
        for r in rs:
            b = int(r["band"])
            g = sorted(groups[(jid, r["group"])])
            # the band's OWN nearest gap, over the BZ
            cand = []
            if b > 0:
                cand.append(E[:, b] - E[:, b - 1])
            if b < nb - 1:
                cand.append(E[:, b + 1] - E[:, b])
            dnear = float(np.min(np.stack(cand), axis=0).min()) if cand else np.inf
            # intra-group splitting: largest adjacent gap INSIDE the group
            dintra = 0.0
            for x, y in zip(g[:-1], g[1:]):
                dintra = max(dintra, float((E[:, y] - E[:, x]).min()))
            out.append((r, dnear, dintra, len(g)))
        if (i + 1) % 40 == 0:
            print(f"  {i+1}/{len(byj)}, {time.time()-t0:.0f}s", flush=True)
    print(f"\n{len(out)} bands, {nfail} materials failed, {time.time()-t0:.0f}s\n")

    res = []
    for r, dnear, dintra, r_grp in out:
        Us = max(dnear, 0.0) / 2
        Um = float(r["d_iso"]) / 2
        Ms, lms = float(r["M_trg"]), float(r["lam"])
        Mm = float(r["M_manifold_diag"])
        lmm = float(r["lam_manifold_diag"])
        tcs = min(Us * lms / 4, C3 * Us * Ms) / KB
        tcm = min(Um * lmm / 4, C3 * Um * Mm) / KB
        # self-consistency: manifold valid if the intra-group splitting is below
        # the Tc it predicts; single valid if it is well above.
        man_ok = (r_grp > 1) and (dintra <= KB * tcm)
        sgl_ok = (r_grp == 1) or (dintra >= KB * tcs)
        if man_ok and not sgl_ok:
            which, M, tc = "manifold", Mm, tcm
        elif sgl_ok and not man_ok:
            which, M, tc = "single", Ms, tcs
        elif man_ok and sgl_ok:
            which, M, tc = "either", min(Ms, Mm), min(tcs, tcm)
        else:
            which, M, tc = "neither", min(Ms, Mm), min(tcs, tcm)
        res.append(dict(jid=r["jid"], formula=r["formula"], band=r["band"],
                        rank=r_grp, d_near=dnear, d_intra=dintra,
                        M_single=Ms, M_manifold=Mm, Tc_single=tcs,
                        Tc_manifold=tcm, treatment=which, M_applied=M, Tc=tc))

    with open(os.path.join(HERE, "selfconsist.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
        w.writeheader(); w.writerows(res)

    Ma = np.array([x["M_applied"] for x in res])
    Ms = np.array([x["M_single"] for x in res])
    tr = np.array([x["treatment"] for x in res])
    print("Which treatment is self-consistent:")
    for k in ["single", "manifold", "either", "neither"]:
        print(f"   {k:<10}{(tr==k).sum():>5} / {len(tr)}")
    print(f"\nABUNDANCE CLAIM, before and after:")
    print(f"   stored  (single-band M): median {np.median(Ms):.4f}   "
          f"above 0.339: {(Ms>0.339).sum()}/{len(Ms)} ({100*(Ms>0.339).mean():.1f}%)")
    print(f"   applied (self-consistent): median {np.median(Ma):.4f}   "
          f"above 0.339: {(Ma>0.339).sum()}/{len(Ma)} ({100*(Ma>0.339).mean():.1f}%)")
    for q in [10, 25, 50, 75, 90, 99]:
        print(f"      {q:>3}th pct:  stored {np.percentile(Ms,q):8.4f}   "
              f"applied {np.percentile(Ma,q):8.4f}")
    tc = np.array([x["Tc"] for x in res])
    print(f"\n   Tc under the applied treatment: median {np.median(tc):.1f} K   "
          f"max {tc.max():.0f} K   above 100 K: {(tc>100).sum()}/{len(tc)}")


if __name__ == "__main__":
    main()
