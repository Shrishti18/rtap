#!/usr/bin/env python3
"""Stage-2 report: raw numbers first, booleans alongside."""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def f(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


bands = list(csv.DictReader(open(os.path.join(HERE, "stage2_bands.csv"))))
pairs = []
pp = os.path.join(HERE, "stage2_pairs.csv")
if os.path.exists(pp):
    pairs = list(csv.DictReader(open(pp)))

matched = list(csv.DictReader(open(os.path.join(HERE, "stage2_matched.csv"))))
jids = sorted({m["jid"] for m in matched})

print("=" * 118)
print("STAGE 2 — PASS materials with a Wannier Hamiltonian")
print("=" * 118)
print("  JARVIS materials screened      : %d" % len(jids))
print("  materials yielding >=1 candidate band (W<=0.5 eV, iso>=0.2 eV): %d"
      % len({b["jid"] for b in bands}))
print("  candidate bands total          : %d" % len(bands))

if not bands:
    print("\n  No candidate band in any of the %d materials." % len(jids))
else:
    print("\n" + "-" * 118)
    print("PER BAND (M is the naive grid metric; M_min is gauge-invariant)")
    print("-" * 118)
    h = ("%-13s %-9s %4s %5s %4s | %8s %6s %6s %7s %6s %5s | %4s %4s %4s | %8s %8s %6s %4s"
         % ("jid", "formula", "dim", "norb", "band", "M_grid", "lam", "W", "d_iso",
            "unif", "nphi", "iso", "M", "unif", "M_min", "M_naive", "ratio", "Mmin"))
    print(h)
    for b in sorted(bands, key=lambda r: (r["jid"], int(r["band"]))):
        mm, mn = f(b["M_min"]), f(b["M_naive"])
        rt = f(b["M_ratio"])
        print("%-13s %-9s %4s %5s %4s | %8.4f %6.3f %6.3f %7.4f %6.3f %5.2f | "
              "%4s %4s %4s | %8s %8s %6s %4s"
              % (b["jid"], b["formula"], b["dim"], b["norb"], b["band"],
                 f(b["M_grid"]), f(b["lam"]), f(b["W"]), f(b["d_iso"]),
                 f(b["unif"]), f(b["nphi"]),
                 "OK" if b["pass_iso"] == "True" else "no",
                 "OK" if b["pass_M_spec"] == "True" else "no",
                 "OK" if b["pass_unif"] == "True" else "no",
                 ("%.4f" % mm) if mm is not None else "-",
                 ("%.4f" % mn) if mn is not None else "-",
                 ("%.2f" % rt) if rt is not None and np.isfinite(rt) else "-",
                 b["pass_M_min"] or "-"))

    npi = [b for b in bands if b["pass_iso"] == "True"]
    print("\n  pass_iso (coarse grid) : %d / %d bands" % (len(npi), len(bands)))
    conf = [b for b in npi if b["pass_iso_fine"] == "True"]
    print("  ...still passing on the fine grid (trap 2): %d" % len(conf))
    for b in npi:
        print("     %-13s band %-3s d_iso %.4f (nk=%s) -> %.4f (nk=%s)  %s"
              % (b["jid"], b["band"], f(b["d_iso"]), b["nk"],
                 f(b["d_iso_fine"], float("nan")), b["nk_fine"],
                 "CONFIRMED" if b["pass_iso_fine"] == "True" else "FALSE POSITIVE"))

    have = [b for b in bands if f(b["M_min"]) is not None]
    if have:
        rs = [f(b["M_ratio"]) for b in have if f(b["M_ratio"]) is not None
              and np.isfinite(f(b["M_ratio"]))]
        print("\n  M_min computed on %d bands" % len(have))
        if rs:
            print("  M_naive / M_min : min %.3f  median %.3f  max %.3f"
                  % (min(rs), float(np.median(rs)), max(rs)))
        fp = [b for b in have if b["pass_M_spec"] == "True" and b["pass_M_min"] == "False"]
        print("  naive metric clears but M_min does not (false positives): %d / %d"
              % (len(fp), len(have)))
        print("  M_min passing threshold (0.34 3D / 0.84 2D): %d"
              % sum(1 for b in have if b["pass_M_min"] == "True"))

    allp = [b for b in bands if b["pass_iso"] == "True" and b["pass_unif"] == "True"]
    print("\n  bands passing BOTH d_iso and unif: %d" % len(allp))

print("\n" + "=" * 118)
print("PAIR CRITERION — `ratio` = inter-orbital / intra-orbital hopping")
print("=" * 118)
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
        print("  finite ratio: min %.4f  median %.4f  mean %.4f  max %.4f"
              % (v.min(), np.median(v), v.mean(), v.max()))
        for q in (50, 75, 90, 95, 99):
            print("      p%-3d %.4f" % (q, np.percentile(v, q)))
        print("  ratio >= 4 (the target): %d of %d finite (%.2f%%)"
              % (int((v >= 4).sum()), len(v), 100.0 * (v >= 4).mean()))
        top = sorted(fin, key=lambda t: -(t[0] or -1))[:12]
        print("\n  highest finite ratios:")
        print("    %-13s %-38s %8s %9s %8s %8s %8s"
              % ("jid", "basis", "ratio", "gap_min", "W_eps0", "d_var", "dzhat"))
        for r, p in top:
            print("    %-13s %-38s %8.3f %9.4f %8.4f %8.3f %8.4f"
                  % (p["jid"], p["basis"][:38], r, f(p["gap_min"]), f(p["W_eps0"]),
                     f(p["d_var"]), f(p["dzhat"])))
    if ideal:
        print("\n  eps0 ~ 0 (perfectly flat identity channel, the GOOD case): %d" % len(ideal))
        for r, p in ideal[:10]:
            print("    %-13s %-38s ratio %.3e  gap_min %.4f"
                  % (p["jid"], p["basis"][:38], r, f(p["gap_min"])))
