#!/usr/bin/env python3
"""Report the full 1,772-material Phase-A sweep. Raw distributions first."""
import csv, os, sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RS = [1.00, 0.85, 0.55]
MANIFEST = "/home/user/rtap/data/priority1_wannier_jarvis/jarvis_wtb_manifest.csv"

rows = list(csv.DictReader(open(os.path.join(HERE, "full_bands.csv"))))
allmat = list(csv.DictReader(open(MANIFEST)))
errs = []
ep = os.path.join(HERE, "full_scan_errors.csv")
if os.path.exists(ep):
    errs = [r for r in csv.DictReader(open(ep))]

F = lambda x: float(x)
mats = {r["jid"] for r in rows}
groups = {(r["jid"], r["group"]) for r in rows}

print("=" * 122)
print("FULL SWEEP — Phase A over every JARVIS Wannier material")
print("=" * 122)
print("  materials scanned                     : %d" % len(allmat))
print("  materials with >=1 isolated 2-group   : %d  (%.1f%%)"
      % (len(mats), 100.0 * len(mats) / len(allmat)))
print("  isolated groups found                 : %d" % len(groups))
print("  candidate BANDS (metric taken per band): %d" % len(rows))
print("  scan errors                           : %d" % len(errs))


def hist(name, vals, edges, note=""):
    v = np.array(vals, float)
    print("\n  %s : min %.5f  median %.5f  mean %.5f  max %.5f %s"
          % (name, v.min(), np.median(v), v.mean(), v.max(), note))
    for lo, hi in zip(edges[:-1], edges[1:]):
        n = int(((v >= lo) & (v < hi)).sum())
        print("      [%8.3f, %8.3f) : %5d %s"
              % (lo, hi, n, "#" * int(round(55.0 * n / max(len(v), 1)))))


print("\n" + "#" * 122)
print("DISTRIBUTIONS over the %d candidate bands" % len(rows))
print("#" * 122)
hist("W_band (eV)", [F(r["W_band"]) for r in rows],
     [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.51])
hist("d_iso  (eV, of enclosing group)", [F(r["d_iso"]) for r in rows],
     [0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 100.0])
hist("lam (per band; 1/n_phi, target ~0.5)", [F(r["lam"]) for r in rows],
     [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.001])
hist("M_trg (per band, nk coarse -> lower bound)", [F(r["M_trg"]) for r in rows],
     [0, 0.05, 0.1, 0.2, 0.339, 0.5, 1.0, 2.0, 1e9])

print("\n" + "#" * 122)
print("Tc_max  (v6: lam_eff, U <= min(d_iso/2, U_cap, 2*w0))")
print("#" * 122)
for R in RS:
    v = np.array([F(r["Tc_R%.2f" % R]) for r in rows])
    best = max(rows, key=lambda r: F(r["Tc_R%.2f" % R]))
    print("  R=%.2f : max %9.1f K   median %7.2f   >300K: %4d   >77K: %4d"
          % (R, v.max(), np.median(v), int((v > 300).sum()), int((v > 77).sum())))
    print("           top = %-13s %-10s group %-7s band %-4s  "
          "M_trg %.4f lam %.4f d_iso %.3f  U %.2f (%s)  binds %s"
          % (best["jid"], best["formula"], best["group"], best["band"],
             F(best["M_trg"]), F(best["lam"]), F(best["d_iso"]),
             F(best["U_R%.2f" % R]), best["cap_R%.2f" % R], best["binds_R%.2f" % R]))

print("\n  U_capped_by distribution:")
for R in RS:
    print("     R=%.2f : caps %-46s binds %s"
          % (R, dict(Counter(r["cap_R%.2f" % R] for r in rows)),
             dict(Counter(r["binds_R%.2f" % R] for r in rows))))

print("\n" + "-" * 122)
print("TOP 20 BANDS by Tc_max at R=0.55 (most conservative mediator)")
print("-" * 122)
print("%-13s %-10s %3s %5s %-7s %4s | %8s %7s %6s %7s %7s | %9s %9s %9s | %5s %5s"
      % ("jid", "formula", "dim", "norb", "group", "band", "M_trg", "lam", "nphi",
         "W_band", "d_iso", "Tc R1.00", "Tc R0.85", "Tc R0.55", "binds", "cap"))
top = sorted(rows, key=lambda r: -F(r["Tc_R0.55"]))[:20]
for r in top:
    print("%-13s %-10s %3s %5s %-7s %4s | %8.4f %7.4f %6.2f %7.4f %7.3f | "
          "%9.1f %9.1f %9.1f | %5s %5s"
          % (r["jid"], r["formula"], r["dim"], r["norb"], r["group"], r["band"],
             F(r["M_trg"]), F(r["lam"]), F(r["nphi"]), F(r["W_band"]), F(r["d_iso"]),
             F(r["Tc_R1.00"]), F(r["Tc_R0.85"]), F(r["Tc_R0.55"]),
             r["binds_R0.55"], r["cap_R0.55"]))

print("\n" + "-" * 122)
print("TOP 20 BANDS by Tc_max at R=1.00 (instantaneous/electronic pairing)")
print("-" * 122)
top1 = sorted(rows, key=lambda r: -F(r["Tc_R1.00"]))[:20]
for r in top1:
    print("%-13s %-10s %3s %5s %-7s %4s | %8.4f %7.4f %6.2f %7.4f %7.3f | "
          "%9.1f %9.1f %9.1f | %5s %5s"
          % (r["jid"], r["formula"], r["dim"], r["norb"], r["group"], r["band"],
             F(r["M_trg"]), F(r["lam"]), F(r["nphi"]), F(r["W_band"]), F(r["d_iso"]),
             F(r["Tc_R1.00"]), F(r["Tc_R0.85"]), F(r["Tc_R0.55"]),
             r["binds_R1.00"], r["cap_R1.00"]))

SUB = [r for r in rows if 0.40 <= F(r["lam"]) <= 0.60]
print("\n" + "-" * 122)
print("GENUINE n_phi ~ 2 SUBSET: lam in [0.40, 0.60]")
print("-" * 122)
print("  %d of %d bands (%.1f%%) on %d materials"
      % (len(SUB), len(rows), 100.0 * len(SUB) / max(len(rows), 1),
         len({r["jid"] for r in SUB})))
if SUB:
    print("  max M_trg in subset : %.5f" % max(F(r["M_trg"]) for r in SUB))
    for R in RS:
        b = max(SUB, key=lambda r: F(r["Tc_R%.2f" % R]))
        print("  max Tc (R=%.2f)      : %9.1f K  %-13s %-10s group %s band %s"
              % (R, F(b["Tc_R%.2f" % R]), b["jid"], b["formula"], b["group"], b["band"]))
    print("\n  top 10 in subset by Tc at R=0.55:")
    for r in sorted(SUB, key=lambda r: -F(r["Tc_R0.55"]))[:10]:
        print("    %-13s %-10s group %-7s band %-4s M_trg %8.4f lam %.4f "
              "d_iso %6.3f  Tc %8.1f / %8.1f / %8.1f"
              % (r["jid"], r["formula"], r["group"], r["band"], F(r["M_trg"]),
                 F(r["lam"]), F(r["d_iso"]), F(r["Tc_R1.00"]), F(r["Tc_R0.85"]),
                 F(r["Tc_R0.55"])))

# write the Phase-B worklist: union of the two top-20s
wl = {(r["jid"], r["group"], r["band"]): r for r in top + top1}
p = os.path.join(HERE, "phaseB_worklist.csv")
with open(p, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(wl.values())
print("\n-> %s (%d bands for Phase B)" % (p, len(wl)))
