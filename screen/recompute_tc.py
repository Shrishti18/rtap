#!/usr/bin/env python3
"""Recompute tc_max from the stored Phase-A descriptors under harvest v6.

Phase A does not need rerunning: lam and U enter only tc_max, and the stored
`lam` column is lam_RAW (the max eigenvalue of A before the r^2 normalisation),
so lam_eff = lam_raw / rank^2 recovers v6's convention exactly. M_trg, d_iso,
nphi and W are unchanged by v6.

Scans R in {1.00, 0.85, 0.55} x U_cap in {1.0, 1.5, 2.5} eV, and reports which
of the three caps on U is binding for each candidate.
"""
import csv, os, sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harvest

HERE = os.path.dirname(os.path.abspath(__file__))
RS = [1.00, 0.85, 0.55]
UCAPS = [1.0, 1.5, 2.5]

rows = list(csv.DictReader(open(os.path.join(HERE, "stage2_manifolds.csv"))))
print("Phase-A manifolds loaded: %d  (materials: %d)"
      % (len(rows), len({r["jid"] for r in rows})))

out = []
for r in rows:
    rank = int(r["nb"])
    lam_raw = float(r["lam"])
    d = dict(M_trg=float(r["M_trg"]), M=float(r["M_trg"]) * (2 * np.pi) ** (int(r["dim"]) - 1),
             d_iso=float(r["d_iso"]), lam=lam_raw / rank ** 2)
    rec = dict(jid=r["jid"], formula=r["formula"], spg=r["spg"], dim=int(r["dim"]),
               norb=int(r["norb"]), nk=int(r["nk"]), bands=r["bands"], rank=rank,
               M_trg=d["M_trg"], M_trg_per_band=d["M_trg"] / rank,
               lam_raw=lam_raw, lam_eff=d["lam"], nphi=float(r["nphi"]),
               W=float(r["W"]), d_iso=d["d_iso"], Emid=float(r["Emid"]),
               Tc_old_R1=float(r.get("Tc_max_K_R1.00", "nan") or "nan"))
    for R in RS:
        for U in UCAPS:
            t = harvest.tc_max(d, dim=rec["dim"], R=R, U_cap=U)
            tag = "R%.2f_U%.1f" % (R, U)
            rec["Tc_" + tag] = t["Tc_max_K"]
            rec["U_" + tag] = t["U_used"]
            rec["cap_" + tag] = t["U_capped_by"]
            rec["binds_" + tag] = t["binds"]
    out.append(rec)

p = os.path.join(HERE, "stage2_tc_v6.csv")
with open(p, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
print("-> %s" % p)

print("\n" + "#" * 118)
print("HEADLINE  (v6: lam_eff = lam_raw/r^2 ; U <= min(d_iso/2, U_cap, 2*w0))")
print("#" * 118)
print("  highest Tc_max, at U_cap = 1.5 eV (default):")
for R in RS:
    k = "Tc_R%.2f_U1.5" % R
    best = max(out, key=lambda x: x[k])
    print("     R=%.2f : %9.1f K   %-13s %-9s bands %-7s  U=%.2f (%s)  binds=%s"
          % (R, best[k], best["jid"], best["formula"], best["bands"],
             best["U_R%.2f_U1.5" % R], best["cap_R%.2f_U1.5" % R],
             best["binds_R%.2f_U1.5" % R]))

print("\n  sensitivity to U_cap (the parameter now doing the most work):")
print("     %-10s %12s %12s %12s" % ("R", "U_cap=1.0", "U_cap=1.5", "U_cap=2.5"))
for R in RS:
    vals = [max(x["Tc_R%.2f_U%.1f" % (R, U)] for x in out) for U in UCAPS]
    print("     R=%.2f     %12.1f %12.1f %12.1f" % (R, *vals))

print("\n  distribution of U_capped_by (all %d manifolds):" % len(out))
for R in RS:
    for U in UCAPS:
        c = Counter(x["cap_R%.2f_U%.1f" % (R, U)] for x in out)
        b = Counter(x["binds_R%.2f_U%.1f" % (R, U)] for x in out)
        print("     R=%.2f U_cap=%.1f : caps %-42s | binds %s"
              % (R, U, dict(c), dict(b)))

print("\n  effect of the two v6 fixes on the old headline:")
old = max(out, key=lambda x: x["Tc_old_R1"])
print("     old max (lam unnormalised, iso cap only) : %9.1f K  (%s %s)"
      % (old["Tc_old_R1"], old["jid"], old["formula"]))
print("     v6  max at R=1.00, U_cap=1.5             : %9.1f K"
      % max(x["Tc_R1.00_U1.5"] for x in out))
print("     v6  max at R=0.55, U_cap=1.5             : %9.1f K"
      % max(x["Tc_R0.55_U1.5"] for x in out))

print("\n" + "-" * 118)
print("TOP 20 MANIFOLDS by Tc_max (R=0.55, U_cap=1.5 -- the most conservative)")
print("-" * 118)
print("%-13s %-9s %3s %5s %-7s | %9s %8s %8s %6s %7s | %8s %8s %8s | %5s %6s"
      % ("jid", "formula", "dim", "norb", "bands", "M_trg", "lam_eff", "lam_raw",
         "nphi", "d_iso", "Tc R1.00", "Tc R0.85", "Tc R0.55", "binds", "cap"))
for x in sorted(out, key=lambda x: -x["Tc_R0.55_U1.5"])[:20]:
    print("%-13s %-9s %3d %5d %-7s | %9.5f %8.4f %8.4f %6.2f %7.3f | %8.1f %8.1f %8.1f | %5s %6s"
          % (x["jid"], x["formula"], x["dim"], x["norb"], x["bands"],
             x["M_trg"], x["lam_eff"], x["lam_raw"], x["nphi"], x["d_iso"],
             x["Tc_R1.00_U1.5"], x["Tc_R0.85_U1.5"], x["Tc_R0.55_U1.5"],
             x["binds_R0.55_U1.5"], x["cap_R0.55_U1.5"]))

print("\n" + "-" * 118)
print("lam_eff DISTRIBUTION (target: lam_eff ~ 0.5 for a genuine n_phi=2 doublet)")
print("-" * 118)
le = np.array([x["lam_eff"] for x in out])
np_ = np.array([x["nphi"] for x in out])
print("  lam_eff : min %.4f  median %.4f  max %.4f" % (le.min(), np.median(le), le.max()))
print("  nphi    : min %.3f  median %.3f  max %.3f" % (np_.min(), np.median(np_), np_.max()))
edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
for lo, hi in zip(edges[:-1], edges[1:]):
    n = int(((le >= lo) & (le < hi)).sum())
    tag = "  <- n_phi ~ 2 target" if lo in (0.4, 0.5) else ""
    print("    [%.1f, %.1f) : %4d %-50s%s" % (lo, hi, n, "#" * int(round(50.0 * n / max(len(le), 1))), tag))
SUB = [x for x in out if 0.40 <= x["lam_eff"] <= 0.60]
print("\n  lam_eff in [0.40, 0.60]: %d of %d manifolds on %d materials"
      % (len(SUB), len(out), len({x["jid"] for x in SUB})))
if SUB:
    print("    max M_trg in subset      : %.6f" % max(x["M_trg"] for x in SUB))
    for R in RS:
        b = max(SUB, key=lambda x: x["Tc_R%.2f_U1.5" % R])
        print("    max Tc (R=%.2f, U_cap=1.5): %8.1f K  %s %s bands %s"
              % (R, b["Tc_R%.2f_U1.5" % R], b["jid"], b["formula"], b["bands"]))
