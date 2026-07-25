"""
Does any REAL material already have the winding structure?

The winding result: if a two-band manifold is H(k) = d0(k)*I + d(k).sigma with
|d(k)| CONSTANT but dhat(k) winding, then the gap is 2|d| = const everywhere,
both bands are exactly flat, and M_trg = |n|^2/4 independent of the gap. That
is the structure that turns the 8.1 K terminus into 334 K.

The objection is that it needs H(R) nonzero ONLY at R = +-n -- an n-th
neighbour hop with all shorter-range hops absent -- which nobody has shown a
mechanism for.

But the signature is directly measurable in band data: a flat band whose gap to
its partner is NEARLY CONSTANT across the BZ. This screens all 728 candidate
bands (179 materials) from the full JARVIS-WTB sweep for it. Eigenvalues only,
so it is cheap.

    unif_gap = min_k Delta(k) / max_k Delta(k)      on the partner side

    winding model      -> 1.000  (exactly)
    atomic-SOC isolated -> much less than 1
"""
import os, sys, csv, time, zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastops

ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
HERE = os.path.dirname(os.path.abspath(__file__))
NK3D, NK2D = 24, 40
MEM = 1.6e9


def pick_nk(norb, dim):
    nk = NK2D if dim == 2 else NK3D
    while nk > 8 and 2.5 * (nk ** dim) * norb * norb * 16 > MEM:
        nk -= 4
    return nk


def gaps_for(jid, dim, bands):
    """Return {band: (unif_lo, unif_hi, mean_lo, mean_hi, W)}."""
    with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
        name = next(n for n in z.namelist() if n.endswith("_hr.dat"))
        R, H, deg = fastops.read_hr_bytes(z.read(name))
    norb = H.shape[1]
    nk = pick_nk(norb, dim)
    Hk, _ = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
    E = np.linalg.eigvalsh(Hk)
    del Hk
    E = E.reshape(-1, E.shape[-1])
    out = {}
    for b in bands:
        if b <= 0 or b >= E.shape[1] - 1:
            continue
        glo = E[:, b] - E[:, b - 1]
        ghi = E[:, b + 1] - E[:, b]
        f = lambda g: float(g.min() / g.max()) if g.max() > 1e-12 else 0.0
        out[b] = (f(glo), f(ghi), float(glo.mean()), float(ghi.mean()),
                  float(np.ptp(E[:, b])))
    return out


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "full_bands.csv"))))
    byj = {}
    for r in rows:
        byj.setdefault((r["jid"], int(r["dim"])), []).append(r)
    print(f"{len(rows)} bands in {len(byj)} materials; eigenvalues only\n")
    res, t0, nfail = [], time.time(), 0
    for i, ((jid, dim), rs) in enumerate(sorted(byj.items())):
        try:
            g = gaps_for(jid, dim, sorted({int(r["band"]) for r in rs}))
        except Exception as e:
            nfail += 1
            print(f"  FAIL {jid}: {type(e).__name__}: {e}", flush=True)
            continue
        for r in rs:
            b = int(r["band"])
            if b not in g:
                continue
            ulo, uhi, mlo, mhi, W = g[b]
            # the partner is the side with the smaller mean gap
            u, mg = (ulo, mlo) if mlo <= mhi else (uhi, mhi)
            res.append(dict(jid=jid, formula=r["formula"], spg=r["spg"],
                            dim=dim, band=b, W_band=float(r["W_band"]),
                            M_trg=float(r["M_trg"]), d_iso=float(r["d_iso"]),
                            nphi=float(r["nphi"]), gap_mean=mg, unif_gap=u))
        if (i + 1) % 40 == 0:
            print(f"  {i+1}/{len(byj)} materials, {time.time()-t0:.0f}s",
                  flush=True)
    print(f"\n{len(res)} bands measured, {nfail} materials failed, "
          f"{time.time()-t0:.0f}s")

    with open(os.path.join(HERE, "constd.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
        w.writeheader()
        w.writerows(res)

    u = np.array([r["unif_gap"] for r in res])
    M = np.array([r["M_trg"] for r in res])
    print("\nDistribution of unif_gap = min_k Delta / max_k Delta:")
    for q in [50, 75, 90, 95, 99, 100]:
        print(f"   {q:3d}th percentile: {np.percentile(u, q):.4f}")
    print(f"\n   unif_gap > 0.5 : {int((u > 0.5).sum()):4d} / {len(u)}")
    print(f"   unif_gap > 0.7 : {int((u > 0.7).sum()):4d} / {len(u)}")
    print(f"   unif_gap > 0.9 : {int((u > 0.9).sum()):4d} / {len(u)}")
    sel = (u > 0.5) & (M > 0.339)
    print(f"\n   BOTH unif_gap > 0.5 AND M_trg > 0.339 : {int(sel.sum())} / {len(u)}")
    print(f"   correlation(unif_gap, M_trg) = {np.corrcoef(u, M)[0,1]:+.4f}")

    print("\nTop 20 by unif_gap:")
    print(f"   {'formula':<16}{'spg':>5}{'jid':>12}{'b':>4}{'W_band':>9}"
          f"{'gap':>8}{'unif':>8}{'M_trg':>9}{'nphi':>7}")
    for r in sorted(res, key=lambda r: -r["unif_gap"])[:20]:
        print(f"   {r['formula'][:15]:<16}{r['spg']:>5}{r['jid']:>12}"
              f"{r['band']:>4}{r['W_band']:>9.4f}{r['gap_mean']:>8.3f}"
              f"{r['unif_gap']:>8.4f}{r['M_trg']:>9.4f}{r['nphi']:>7.3f}")


if __name__ == "__main__":
    main()
