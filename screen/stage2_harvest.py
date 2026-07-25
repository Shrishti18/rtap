#!/usr/bin/env python3
"""Stage 2b — harvest the PASS materials that have a Wannier Hamiltonian.

Numerics go through fastops, which is validated against the shipped harvest
routines to ~1e-14 (see validate_fastops.py). The shipped hk_grid and
_trg_shifted are correct but allocate a full (nk^dim, norb, norb) temporary per
R-vector / per objective evaluation, which makes Powell intractable here.

The three convergence traps are handled as instructed:
  1. M UNDERestimated on coarse grids -> raw M / M_min / M_naive always
     reported, never only the boolean.
  2. d_iso OVERestimated on coarse grids -> every band with pass_iso=True is
     re-checked on a much finer grid, computed slice-by-slice.
  3. pair_criterion `ratio` diverging when eps0 == 0 -> flagged as
     "eps0 ~ 0, ideal" and excluded from the finite-ratio statistics.
"""
import csv, os, sys, time, traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harvest, fastops

HERE = os.path.dirname(os.path.abspath(__file__))
HR = os.path.join(HERE, "hr")
SMALL = "/home/user/rtap/data/priority1_wannier_jarvis/extracted_small"

M_THRESH = {2: 0.84, 3: 0.34}      # user-supplied M_min thresholds
WMAX, ISO_MIN = 0.5, 0.2
NK_WANT, NK_FINE_3D, NK_FINE_2D = 40, 96, 400
MEM = 8.0e9                        # peak bytes for Hk + V


def poscar_lattice(jid):
    L = open(os.path.join(SMALL, jid, "POSCAR")).read().splitlines()
    s = float(L[1].split()[0])
    return np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s


def frac_centres(jid, dim):
    """wannier90.wout centres are Cartesian Angstrom; the shift phase is
    k_i * d_i with k in units of 2*pi per lattice vector, so they must be
    converted to fractional coordinates first."""
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


def auto_bands_from_E(E, wmax=WMAX, iso_min=ISO_MIN):
    """harvest.auto_bands, but reusing an eigendecomposition already done."""
    nb = E.shape[-1]
    out = []
    for b in range(nb):
        Eb = E[..., b]
        W = Eb.max() - Eb.min()
        gl = Eb.min() - E[..., b - 1].max() if b > 0 else np.inf
        gu = E[..., b + 1].min() - Eb.max() if b < nb - 1 else np.inf
        iso = min(gl, gu)
        if W <= wmax and iso >= iso_min:
            out.append((b, float(W), float(iso)))
    return out


def band_extrema_fine(R, H, deg, nk, dim):
    """Per-band min/max on a fine grid, slice by slice (eigenvalues only)."""
    n = H.shape[1]
    Hf = (H / deg[:, None, None]).reshape(H.shape[0], n * n)
    ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    lo = np.full(n, np.inf)
    hi = np.full(n, -np.inf)
    if dim == 2:
        m = np.meshgrid(ks, ks, indexing="ij")
        kf = np.stack([x.ravel() for x in m], 1)
        E = np.linalg.eigvalsh((np.exp(1j * (kf @ R[:, :2].T.astype(float))) @ Hf)
                               .reshape(-1, n, n))
        return E.min(0), E.max(0)
    m1, m2 = np.meshgrid(ks, ks, indexing="ij")
    base = np.stack([m1.ravel(), m2.ravel()], 1)
    for k3 in ks:
        kf = np.concatenate([base, np.full((base.shape[0], 1), k3)], 1)
        E = np.linalg.eigvalsh((np.exp(1j * (kf @ R[:, :3].T.astype(float))) @ Hf)
                               .reshape(-1, n, n))
        lo = np.minimum(lo, E.min(0))
        hi = np.maximum(hi, E.max(0))
    return lo, hi


def iso_from_extrema(lo, hi, b):
    n = len(lo)
    gl = lo[b] - hi[b - 1] if b > 0 else np.inf
    gu = lo[b + 1] - hi[b] if b < n - 1 else np.inf
    return float(min(gl, gu)), float(hi[b] - lo[b])


def run_one(rec, min_budget=2):
    jid = rec["jid"]
    dim = 2 if rec["dimensionality"] == "2D" else 3
    R, H, deg = harvest.read_hr(os.path.join(HR, "%s_hr.dat" % jid))
    norb = H.shape[1]
    nk = pick_nk(norb, dim)
    t0 = time.time()
    Hk, dks = fastops.hk_grid_fast(R, H, deg, nk, dim=dim)
    E, V = np.linalg.eigh(Hk)
    del Hk
    cands = auto_bands_from_E(E)
    centres = frac_centres(jid, dim)
    nk_fine = NK_FINE_2D if dim == 2 else NK_FINE_3D
    lo = hi = None
    rows = []
    for b, W, iso in cands:
        d = fastops.descriptors_fast(E, V, dks, b)
        sp = harvest.spec(d, Tc_K=300, dim=dim)
        row = dict(jid=jid, formula=rec["jarvis_formula"], spg=rec["jarvis_spg"],
                   dim=dim, norb=norb, nk=nk, band=b,
                   M_grid=d["M"], lam=d["lam"], W=d["W"], d_iso=d["d_iso"],
                   unif=d["unif"], nphi=d["nphi"], Emid=d["Emid"],
                   U_req=sp["U_req"], iso_req=sp["iso_req"], M_req_spec=sp["M_req"],
                   pass_iso=sp["pass_iso"], pass_M_spec=sp["pass_M"],
                   pass_unif=sp["pass_unif"], calibrated=sp["calibrated"],
                   M_min="", M_naive="", M_ratio="", pass_M_min="",
                   d_iso_fine="", W_fine="", nk_fine="", pass_iso_fine="")
        if sp["pass_iso"]:                       # trap 2 — the unsafe direction
            if lo is None:
                lo, hi = band_extrema_fine(R, H, deg, nk_fine, dim)
            iso_f, W_f = iso_from_extrema(lo, hi, b)
            row.update(d_iso_fine=iso_f, W_fine=W_f, nk_fine=nk_fine,
                       pass_iso_fine=bool(iso_f >= sp["iso_req"]))
        rows.append((row, d))

    # gauge-invariant metric on the most promising bands only
    order = sorted(range(len(rows)),
                   key=lambda i: (not rows[i][0]["pass_iso"],
                                  rows[i][0]["unif"]))
    for i in order[:min_budget]:
        row, d = rows[i]
        mm = fastops.minimal_metric_fast(None, dks, row["band"], dim,
                                         centres=centres, restarts=2,
                                         u=V[..., :, row["band"]], maxfev=40000)
        row.update(M_min=mm["M_min"] * (2 * np.pi) ** (dim - 1),
                   M_naive=mm["M_naive"] * (2 * np.pi) ** (dim - 1))
        row["M_ratio"] = (row["M_naive"] / row["M_min"]) if row["M_min"] > 0 else np.inf
        row["pass_M_min"] = bool(row["M_min"] >= M_THRESH[dim])
    del V, E
    return [r for r, _ in rows], [d for _, d in rows], cands, nk, dim, time.time() - t0


def main():
    recs = list(csv.DictReader(open(os.path.join(HERE, "stage2_matched.csv"))))
    byj = {}
    for r in recs:
        byj.setdefault(r["jid"], r)
    jids = sorted(byj, key=lambda s: int(s.split("-")[1]))
    only = sys.argv[1:] or jids
    allrows, prows, errs = [], [], []
    for i, jid in enumerate(only, 1):
        rec = byj[jid]
        try:
            rows, dsc, cands, nk, dim, dt = run_one(rec)
            allrows.extend(rows)
            npi = sum(1 for r in rows if r["pass_iso"])
            print("[%2d/%d] %-13s %-10s dim=%d norb=%-3d nk=%-3d bands=%-2d "
                  "pass_iso=%d  %.0fs"
                  % (i, len(only), jid, rec["jarvis_formula"], dim, rows[0]["norb"]
                     if rows else int(rec["num_wann"]), nk, len(rows), npi, dt),
                  flush=True)
            if cands:
                R, H, deg = harvest.read_hr(os.path.join(HR, "%s_hr.dat" % jid))
                idx = [c[0] for c in cands]
                for b in idx:
                    if b + 1 in idx:
                        pc = harvest.pair_criterion(R, H, deg, [b, b + 1],
                                                    nk=min(nk, 24), dim=dim)
                        pc.update(jid=jid, i1=b, i2=b + 1, band_ref=b,
                                  basis="band-index pair (as instructed)")
                        prows.append(pc)
                for (row, d) in zip(rows, dsc):
                    o1, o2 = (int(x) for x in np.argsort(-d["w"])[:2])
                    pc = harvest.pair_criterion(R, H, deg, [o1, o2],
                                                nk=min(nk, 24), dim=dim)
                    pc.update(jid=jid, i1=o1, i2=o2, band_ref=row["band"],
                              basis="dominant Wannier orbitals of the band")
                    prows.append(pc)
        except Exception as e:
            errs.append((jid, "%s: %s" % (type(e).__name__, e)))
            print("[%2d/%d] %-13s ERROR %s" % (i, len(only), jid, e), flush=True)
            traceback.print_exc()

    if allrows:
        p = os.path.join(HERE, "stage2_bands.csv")
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
            w.writeheader()
            w.writerows(allrows)
        print("-> %s (%d band rows)" % (p, len(allrows)))
    if prows:
        p = os.path.join(HERE, "stage2_pairs.csv")
        keys = sorted({k for r in prows for k in r})
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(prows)
        print("-> %s (%d pair rows)" % (p, len(prows)))
    if errs:
        with open(os.path.join(HERE, "stage2_errors.csv"), "w", newline="") as f:
            csv.writer(f).writerows([["jid", "error"]] + errs)
        print("errors: %d" % len(errs))
    print("STAGE2 DONE", flush=True)


if __name__ == "__main__":
    main()
