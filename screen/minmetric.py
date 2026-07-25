"""
How much of the screened metric survives the MINIMAL metric?

C4's retraction showed that M_naive can be entirely gauge: the winding model has
M_naive = n^2/4 and M_min = 0, because the whole "metric" was the spread of a
molecular orbital across a disconnected dimer. The abundance screen used
M_trg (naive) throughout, so its 43.5% is an UPPER BOUND on physically usable
metric.

This recomputes M_min for the screened bands, seeded with the TRUE Wannier
centres from wannier90.wout (converted to fractional with the POSCAR lattice,
and seeded at both signs -- v4's D = diag(d_a) is the sign-flip of gauge.py's
shift, and on the decoupled dimer the exact zero sits at d_B = -1, not +1).

Run with a norb cap: the Powell optimisation carries (norb-1)*dim parameters and
was already found intractable at 657. Bands above the cap are reported as
NOT MEASURED rather than dropped silently.
"""
import os, sys, csv, time, zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/user/rtap/review")
import fastops, harvest

ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
HERE = os.path.dirname(os.path.abspath(__file__))
NORB_CAP = int(os.environ.get("NORB_CAP", "60"))
MAXFEV = int(os.environ.get("MAXFEV", "6000"))
RESTARTS = int(os.environ.get("RESTARTS", "1"))


def pick_nk(norb, dim):
    nk = 40 if dim == 2 else 24
    while nk > 8 and 2.5 * (nk ** dim) * norb * norb * 16 > 1.6e9:
        nk -= 4
    return nk


def poscar_lattice(txt):
    L = txt.splitlines()
    scale = float(L[1].split()[0])
    return np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * scale


def load(jid):
    with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % jid)) as z:
        nm = next(n for n in z.namelist() if n.endswith("_hr.dat"))
        R, H, deg = fastops.read_hr_bytes(z.read(nm))
        cen = None
        try:
            wn = next(n for n in z.namelist() if n.endswith(".wout"))
            import tempfile
            with tempfile.NamedTemporaryFile("wb", suffix=".wout",
                                             delete=False) as fh:
                fh.write(z.read(wn)); tmp = fh.name
            cart = harvest.read_wout_centres(tmp)
            os.unlink(tmp)
            pn = next(n for n in z.namelist() if n.endswith("POSCAR"))
            A = poscar_lattice(z.read(pn).decode(errors="ignore"))
            if cart is not None:
                cen = cart @ np.linalg.inv(A)      # Cartesian A -> fractional
        except Exception:
            cen = None
    return R, H, deg, cen


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "full_bands.csv"))))
    seen = {}
    for r in rows:
        seen[(r["jid"], r["band"])] = r
    rows = list(seen.values())
    byj = {}
    for r in rows:
        byj.setdefault((r["jid"], int(r["dim"])), []).append(r)
    print(f"{len(rows)} bands, {len(byj)} materials, NORB_CAP={NORB_CAP}, "
          f"MAXFEV={MAXFEV}, RESTARTS={RESTARTS}", flush=True)

    FIELDS = ["jid", "formula", "band", "norb", "dim", "M_naive", "M_min",
              "ratio", "centres", "status"]
    csvpath = os.path.join(HERE, os.environ.get("OUT", "minmetric.csv"))
    fh = open(csvpath, "w", newline="")
    wr = csv.DictWriter(fh, fieldnames=FIELDS)
    wr.writeheader(); fh.flush()

    def emit(recs):
        wr.writerows(recs); fh.flush()

    out, t0, skipped, nocen = [], time.time(), 0, 0
    for i, ((jid, dim), rs) in enumerate(sorted(byj.items())):
        _n0 = len(out)
        norb = int(rs[0]["norb"])
        if norb > NORB_CAP:
            skipped += len(rs)
            for r in rs:
                out.append(dict(jid=jid, formula=r["formula"], band=r["band"],
                                norb=norb, dim=dim, M_naive=float(r["M_trg"]),
                                M_min=float("nan"), ratio=float("nan"),
                                centres="skip", status="NOT_MEASURED"))
            continue
        try:
            R, H, deg, cen = load(jid)
            nk = pick_nk(H.shape[1], dim)
            Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
            # ONE eigh per material. fastops.band_vectors re-diagonalises on
            # every call, which for a 4-band material is a 4x redundant eigh
            # and was the entire cost of the first two attempts.
            _, V = np.linalg.eigh(Hk)
            del Hk
            if cen is not None and cen.shape[0] != H.shape[1]:
                cen = None
            if cen is None:
                nocen += 1
            cfrac = cen[:, :dim] if cen is not None else None
            for r in rs:
                b = int(r["band"])
                u = V[..., :, b]
                naive = fastops.trg_from_u(u, dks)
                if cfrac is None:
                    out.append(dict(jid=jid, formula=r["formula"],
                                    band=r["band"], norb=norb, dim=dim,
                                    M_naive=naive, M_min=float("nan"),
                                    ratio=float("nan"), centres="no",
                                    status="NO_CENTRES"))
                    continue
                # PHYSICAL gauge: orbitals at their true Wannier centres. The
                # sign convention has bitten this project repeatedly (v4's
                # D = diag(d_a) is the sign-flip of gauge.py's shift), so
                # evaluate BOTH and take the smaller, recording which won.
                rel = cfrac[1:] - cfrac[0]
                sm = np.zeros((norb, dim)); sm[1:] = -rel
                sp = np.zeros((norb, dim)); sp[1:] = rel
                vm = fastops.trg_from_u(u, dks, sm)
                vp = fastops.trg_from_u(u, dks, sp)
                mc = min(vm, vp)
                out.append(dict(
                    jid=jid, formula=r["formula"], band=r["band"], norb=norb,
                    dim=dim, M_naive=naive, M_min=mc,
                    ratio=mc / naive if naive > 1e-12 else float("nan"),
                    centres="minus" if vm <= vp else "plus", status="OK"))
            del V
        except Exception as e:
            print(f"  FAIL {jid}: {type(e).__name__}: {e}", flush=True)
            for r in rs:
                out.append(dict(jid=jid, formula=r["formula"], band=r["band"],
                                norb=norb, dim=dim, M_naive=float(r["M_trg"]),
                                M_min=float("nan"), ratio=float("nan"),
                                centres="err", status="FAIL"))
        emit(out[_n0:])
        if (i + 1) % 10 == 0:
            done = sum(1 for x in out if x["status"] == "OK")
            print(f"  {i+1}/{len(byj)} materials, {done} bands measured, "
                  f"{time.time()-t0:.0f}s", flush=True)

    fh.close()

    ok = [x for x in out if x["status"] == "OK"]
    print(f"\n{len(ok)} bands measured, {skipped} above the norb cap "
          f"(NOT MEASURED), {nocen} materials without usable centres, "
          f"{time.time()-t0:.0f}s")
    if not ok:
        return
    mn = np.array([x["M_naive"] for x in ok])
    mm = np.array([x["M_min"] for x in ok])
    rt = mm / np.maximum(mn, 1e-12)
    print(f"\nM_naive : median {np.median(mn):.4f}   above 0.339 "
          f"{int((mn>0.339).sum())}/{len(mn)} ({100*(mn>0.339).mean():.1f}%)")
    print(f"M_min   : median {np.median(mm):.4f}   above 0.339 "
          f"{int((mm>0.339).sum())}/{len(mm)} ({100*(mm>0.339).mean():.1f}%)")
    print(f"\nM_min / M_naive: median {np.median(rt):.4f}")
    for q in [10, 25, 50, 75, 90]:
        print(f"   {q:>3}th pct: {np.percentile(rt,q):.4f}")
    print(f"\n   ratio < 0.10 (metric almost entirely gauge): "
          f"{int((rt<0.10).sum())}/{len(rt)}")
    print(f"   ratio < 0.50 : {int((rt<0.50).sum())}/{len(rt)}")
    print(f"   ratio > 0.90 : {int((rt>0.90).sum())}/{len(rt)}")


if __name__ == "__main__":
    main()
