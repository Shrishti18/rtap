#!/usr/bin/env python3
"""Stage-2 report. Raw numbers first; Tc_max in place of pass/fail."""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RS = [1.00, 0.85, 0.55]


def f(x, d=None):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return d


bands = list(csv.DictReader(open(os.path.join(HERE, "stage2_bands.csv"))))
pp = os.path.join(HERE, "stage2_pairs.csv")
pairs = list(csv.DictReader(open(pp))) if os.path.exists(pp) else []
matched = list(csv.DictReader(open(os.path.join(HERE, "stage2_matched.csv"))))
jids = sorted({m["jid"] for m in matched})

print("=" * 124)
print("STAGE 2 — Stage-1 PASS materials that have a Wannier Hamiltonian")
print("=" * 124)
print("  JARVIS materials screened                            : %d" % len(jids))
print("  materials with >=1 candidate band (W<=0.5, iso>=0.2)  : %d"
      % len({b["jid"] for b in bands}))
print("  candidate bands total                                : %d" % len(bands))

if not bands:
    print("\n  No candidate band found in any material.")
    raise SystemExit

key = lambda R: "Tc_max_K_Mmin_R%.2f" % R
keyn = lambda R: "Tc_max_K_R%.2f" % R

print("\n" + "#" * 124)
print("HEADLINE")
print("#" * 124)
print("  1. highest Tc_max in the set (on the gauge-invariant M_min):")
for R in RS:
    v = [(f(b[key(R)], -1), b) for b in bands if b.get(key(R))]
    if not v:
        continue
    top, b = max(v, key=lambda t: t[0])
    tn = f(b[keyn(R)], float("nan"))
    print("       R=%.2f : %8.2f K   %s band %s (%s, SG %s, %sD)   [naive M: %.2f K]  binds=%s"
          % (R, top, b["jid"], b["band"], b["formula"], b["spg"], b["dim"], tn,
             b.get("binds_Mmin_R%.2f" % R, "?")))
print("     same, on the naive metric (upper bound):")
for R in RS:
    v = [(f(b[keyn(R)], -1), b) for b in bands if b.get(keyn(R))]
    top, b = max(v, key=lambda t: t[0])
    print("       R=%.2f : %8.2f K   %s band %s (%s)"
          % (R, top, b["jid"], b["band"], b["formula"]))

print("\n  2. distribution of M_trg and M_min_trg:")
mt = np.array([f(b["M_trg"], np.nan) for b in bands])
mm = np.array([f(b["M_min_trg"], np.nan) for b in bands if f(b["M_min_trg"]) is not None])
print("       M_trg     : max %.6f  median %.6f  min %.6f" % (np.nanmax(mt), np.nanmedian(mt), np.nanmin(mt)))
if len(mm):
    print("       M_min_trg : max %.6f  median %.6f  min %.6f" % (mm.max(), np.median(mm), mm.min()))

print("\n  3. how often the gauge fix matters (M_naive / M_min):")
rat = np.array([f(b["M_ratio"], np.nan) for b in bands])
rat = rat[np.isfinite(rat)]
if len(rat):
    print("       min %.3f  median %.3f  mean %.3f  max %.3f" % (rat.min(), np.median(rat), rat.mean(), rat.max()))
    for thr in (1.1, 1.5, 2.0, 5.0):
        print("       ratio > %.1f : %d / %d bands (%.1f%%)"
              % (thr, int((rat > thr).sum()), len(rat), 100.0 * (rat > thr).mean()))

print("\n" + "-" * 124)
print("TOP 20 BANDS BY Tc_max (R=1.00, on M_min)")
print("-" * 124)
print("%-13s %-9s %3s %5s %4s | %9s %9s %6s | %6s %8s %6s %5s | %8s %8s %8s | %5s"
      % ("jid", "formula", "dim", "norb", "band", "M_min_trg", "M_trg", "ratio",
         "lam", "d_iso", "W", "nphi", "Tc R1.00", "Tc R0.85", "Tc R0.55", "binds"))
rank = sorted(bands, key=lambda b: -(f(b.get(key(1.0)), -1) or -1))
for b in rank[:20]:
    print("%-13s %-9s %3s %5s %4s | %9.6f %9.6f %6s | %6.4f %8.4f %6.3f %5.2f | "
          "%8.2f %8.2f %8.2f | %5s"
          % (b["jid"], b["formula"], b["dim"], b["norb"], b["band"],
             f(b["M_min_trg"], np.nan), f(b["M_trg"], np.nan),
             ("%.2f" % f(b["M_ratio"])) if f(b["M_ratio"]) is not None
             and np.isfinite(f(b["M_ratio"])) else "-",
             f(b["lam"], np.nan), f(b["d_iso"], np.nan), f(b["W"], np.nan),
             f(b["nphi"], np.nan),
             f(b.get(key(1.0)), np.nan), f(b.get(key(0.85)), np.nan),
             f(b.get(key(0.55)), np.nan), b.get("binds_Mmin_R1.00", "?")))

npi = [b for b in bands if b["pass_iso"] == "True"]
print("\n  d_iso >= 2*U_req on the coarse grid : %d / %d bands" % (len(npi), len(bands)))
conf = [b for b in npi if b["pass_iso_fine"] == "True"]
print("  ...confirmed on the fine grid (trap 2): %d  (false positives: %d)"
      % (len(conf), len(npi) - len(conf)))
for b in npi[:25]:
    print("     %-13s band %-3s  d_iso %.4f (nk=%s) -> %.4f (nk=%s)   %s"
          % (b["jid"], b["band"], f(b["d_iso"], np.nan), b["nk"],
             f(b["d_iso_fine"], np.nan), b["nk_fine"],
             "CONFIRMED" if b["pass_iso_fine"] == "True" else "FALSE POSITIVE"))

binds = {}
for b in bands:
    binds[b.get("binds_Mmin_R1.00", "?")] = binds.get(b.get("binds_Mmin_R1.00", "?"), 0) + 1
print("\n  what limits Tc (R=1.00, on M_min): %s" % binds)

print("\n" + "=" * 124)
print("PAIR CRITERION — `ratio` = inter-orbital / intra-orbital hopping")
print("=" * 124)
if not pairs:
    print("  no pairs evaluated")
else:
    fin, ideal = [], []
    for p in pairs:
        r = f(p["ratio"])
        (ideal if (r is not None and r > 1e6) else fin).append((r, p))
    print("  pair evaluations: %d   (finite %d, eps0 ~ 0 'ideal' %d)"
          % (len(pairs), len(fin), len(ideal)))
    vals = [r for r, _ in fin if r is not None]
    if vals:
        v = np.array(vals)
        print("  finite ratio: min %.4f  median %.4f  mean %.4f  MAX %.4f"
              % (v.min(), np.median(v), v.mean(), v.max()))
        for q in (50, 75, 90, 95, 99):
            print("      p%-3d %.4f" % (q, np.percentile(v, q)))
        print("  ratio >= 4 (the target): %d of %d finite (%.2f%%)"
              % (int((v >= 4).sum()), len(v), 100.0 * (v >= 4).mean()))
        print("\n  highest finite ratios:")
        print("    %-13s %-40s %9s %9s %8s %8s %8s"
              % ("jid", "basis", "ratio", "gap_min", "W_eps0", "d_var", "dzhat"))
        for r, p in sorted(fin, key=lambda t: -(t[0] or -1))[:12]:
            print("    %-13s %-40s %9.4f %9.4f %8.4f %8.3f %8.4f"
                  % (p["jid"], p["basis"][:40], r, f(p["gap_min"], np.nan),
                     f(p["W_eps0"], np.nan), f(p["d_var"], np.nan),
                     f(p["dzhat"], np.nan)))
    if ideal:
        print("\n  eps0 ~ 0 (perfectly flat identity channel, the GOOD case): %d" % len(ideal))
        for r, p in ideal[:10]:
            print("    %-13s %-40s ratio %.3e  gap_min %.4f"
                  % (p["jid"], p["basis"][:40], r, f(p["gap_min"], np.nan)))
