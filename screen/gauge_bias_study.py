#!/usr/bin/env python3
"""Item 4: how much was the v3 phase-then-roll boundary bug biasing the metric?

At shifts = 0 the two forms are identical by construction (V = I, so there is
no phase to roll and no commutator), so M_naive(new)/M_naive(old) is exactly 1
and carries no information. The bias lives entirely in the SHIFTED objective,
i.e. in M_min and in the value at the physical Wannier centres. Both are
reported here.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harvest, fastops
import stage2_harvest as S

HERE = os.path.dirname(os.path.abspath(__file__))


def trg_old(P, dks, kgrids, shifts):
    """The v3 phase-then-roll form, kept here for the comparison only."""
    ph = sum(kgrids[i][..., None] * shifts[:, i] for i in range(len(dks)))
    V = np.exp(1j * ph)
    Ps = P * V[..., :, None] * np.conj(V)[..., None, :]
    trg = np.zeros(Ps.shape[:-2])
    for i in range(len(dks)):
        dP = (np.roll(Ps, -1, i) - np.roll(Ps, 1, i)) / (2 * dks[i])
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    return float(trg.mean())


def study(jid, dim, nk, bands, rng):
    R, H, deg = harvest.read_hr(os.path.join(S.HR, "%s_hr.dat" % jid))
    Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
    E, V = np.linalg.eigh(Hk)
    ks = [2 * np.pi * (np.arange(nk) + 0.5) / nk for _ in range(dim)]
    kg = np.meshgrid(*ks, indexing="ij")
    centres = S.frac_centres(jid, dim)
    out = []
    for b in bands:
        u = V[..., :, b]
        P = np.einsum('...i,...j->...ij', u, u.conj())
        norb = u.shape[-1]
        z = np.zeros((norb, dim))
        new0 = fastops.trg_from_u(u, dks)
        old0 = trg_old(P, dks, kg, z)
        rec = dict(jid=jid, band=b, M_naive_new=new0, M_naive_old=old0,
                   ratio_naive=new0 / old0 if old0 else np.nan)
        if centres is not None:
            c = np.asarray(centres, float)
            s = np.zeros((norb, dim))
            s[1:] = -(c[1:] - c[0])          # v4 sign convention
            rec["at_centres_new"] = fastops.trg_from_u(u, dks, s)
            rec["at_centres_old"] = trg_old(P, dks, kg, s)
            rec["ratio_at_centres"] = rec["at_centres_new"] / rec["at_centres_old"]
        vals = []
        for _ in range(5):
            s = rng.uniform(-1, 1, (norb, dim))
            s[0] = 0
            vals.append((fastops.trg_from_u(u, dks, s), trg_old(P, dks, kg, s)))
        vn = np.array([v[0] for v in vals])
        vo = np.array([v[1] for v in vals])
        rec["ratio_random_mean"] = float((vn / vo).mean())
        rec["ratio_random_min"] = float((vn / vo).min())
        rec["ratio_random_max"] = float((vn / vo).max())
        out.append(rec)
    return out


if __name__ == "__main__":
    import csv
    rng = np.random.default_rng(0)
    rows = list(csv.DictReader(open(os.path.join(HERE, "stage2_bands.csv"))))
    # a handful of real bands, smallest orbital counts first to keep it cheap
    sel = sorted(rows, key=lambda r: int(r["norb"]))[:6]
    res = []
    for r in sel:
        dim = int(r["dim"])
        nk = 16 if dim == 3 else 40
        res += study(r["jid"], dim, nk, [int(r["band"])], rng)
        print("done %s band %s" % (r["jid"], r["band"]), flush=True)
    print("\n%-13s %5s | %12s %12s %7s | %12s %12s %7s | %s"
          % ("jid", "band", "naive_new", "naive_old", "ratio",
             "centres_new", "centres_old", "ratio", "random-shift ratio (min/mean/max)"))
    for r in res:
        print("%-13s %5d | %12.8f %12.8f %7.4f | %12s %12s %7s | %.3f / %.3f / %.3f"
              % (r["jid"], r["band"], r["M_naive_new"], r["M_naive_old"],
                 r["ratio_naive"],
                 ("%.8f" % r["at_centres_new"]) if "at_centres_new" in r else "-",
                 ("%.8f" % r["at_centres_old"]) if "at_centres_old" in r else "-",
                 ("%.4f" % r["ratio_at_centres"]) if "ratio_at_centres" in r else "-",
                 r["ratio_random_min"], r["ratio_random_mean"], r["ratio_random_max"]))
    with open(os.path.join(HERE, "gauge_bias.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader()
        w.writerows(res)
    print("-> gauge_bias.csv")
