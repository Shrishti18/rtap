#!/usr/bin/env python3
"""Stage-1 analysis: headline PASS fraction plus the requested breakdowns."""
import csv, collections, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
META = "/home/user/rtap/data/priority2_flatbands/flatband_icsd_metadata.csv"

sites = list(csv.DictReader(open(os.path.join(HERE, "stage1_sites.csv"))))
meta = {r["icsd"]: r for r in csv.DictReader(open(META))}
fails = list(csv.DictReader(open(os.path.join(HERE, "stage1_parse_failures.csv"))))

by_icsd = collections.defaultdict(list)
for s in sites:
    by_icsd[s["icsd"]].append(s)

ALL = set(meta)
CUR = {i for i, m in meta.items() if m["curated_flatband"] == "1"}
BEST = {i for i, m in meta.items() if m["best_flatband"] == "1"}
ATOM = {i for i, m in meta.items() if m["atomic_flatband"] == "1"}

passing = {i for i, ss in by_icsd.items() if any(s["verdict"] == "PASS" for s in ss)}
dilute = {i for i, ss in by_icsd.items()
          if i not in passing and any(s["verdict"] == "dilute" for s in ss)}


def mat(icsds):
    return {meta[i]["material_id"] for i in icsds if i in meta}


def line(label, subset):
    n = len(subset)
    p = len(subset & passing)
    d = len(subset & dilute)
    ms, mp = mat(subset), mat(subset & passing)
    print("  %-26s %5d ICSD | PASS %5d (%5.1f%%) | dilute %5d (%4.1f%%) "
          "| materials %5d -> %5d (%5.1f%%)"
          % (label, n, p, 100.0 * p / n if n else 0, d, 100.0 * d / n if n else 0,
             len(ms), len(mp), 100.0 * len(mp) / len(ms) if ms else 0))


print("=" * 108)
print("STAGE 1 HEADLINE — at least one PASS site (n_phi == 2)")
print("=" * 108)
line("ALL harvested CIFs", ALL)
line("curated (2,379 materials)", CUR)
line("top candidates (345 mat.)", BEST)
line("atomic (1,102 materials)", ATOM)
print("\n  parse failures: %d" % len(fails))
print("  site rows: %d over %d CIFs" % (len(sites), len(by_icsd)))

psites = [s for s in sites if s["verdict"] == "PASS"]
print("\n  PASS sites: %d, on %d distinct CIFs" % (len(psites), len(passing)))

print("\n" + "=" * 108)
print("(a) SHELL of the doublet-bearing PASS site")
print("=" * 108)
c = collections.Counter(s["shell"] for s in psites)
shell_mat = collections.defaultdict(set)
for s in psites:
    if s["icsd"] in meta:
        shell_mat[s["shell"]].add(meta[s["icsd"]]["material_id"])
for k in ["3d", "4d", "5d", "-"]:
    flag = "  <-- FLAG (4d/5d)" if k in ("4d", "5d") else ""
    print("  %-6s %6d PASS sites | %5d materials%s"
          % (k if k != "-" else "other", c.get(k, 0), len(shell_mat.get(k, ())), flag))

print("\n" + "=" * 108)
print("(b) DOUBLET TYPE at PASS sites (a site may carry more than one)")
print("=" * 108)
kc = collections.Counter()
km = collections.defaultdict(set)
for s in psites:
    for k in s["doublet_kinds"].split("|"):
        if not k:
            continue
        kc[k] += 1
        if s["icsd"] in meta:
            km[k].add(meta[s["icsd"]]["material_id"])
for k, v in kc.most_common():
    print("  %-22s %6d PASS sites | %5d materials" % (k, v, len(km[k])))

print("\n" + "=" * 108)
print("(c) CHEMISTRY FLAG — 6s^2 lone pair (Bi,Pb,Tl,Sb,Sn,Te) or JT ion")
print("      (Cu2+,Mn3+,Ni3+,Cr2+), intersected with PASS")
print("=" * 108)
chem = {i for i, ss in by_icsd.items() if ss and ss[0]["chem_flag"] == "1"}
lone = {i for i, ss in by_icsd.items() if ss and ss[0]["lone_pair"]}
jt = {i for i, ss in by_icsd.items() if ss and ss[0]["jt_ions"]}
for label, s in [("chem flag (either)", chem), ("  lone pair only", lone), ("  JT ion only", jt)]:
    inter = s & passing
    print("  %-22s %5d ICSD | PASS&flag %5d | materials %5d"
          % (label, len(s), len(inter), len(mat(inter))))

print("\n  --- PASS & chem flag, by subset ---")
for label, subset in [("ALL", ALL), ("curated", CUR), ("top-345", BEST), ("atomic", ATOM)]:
    inter = subset & passing & chem
    print("    %-10s %5d ICSD | %5d materials" % (label, len(inter), len(mat(inter))))

# the survivors, written out for stage 2
surv = sorted(passing, key=int)
with open(os.path.join(HERE, "stage1_pass_icsds.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["icsd", "material_id", "formula", "space_group_number", "curated",
                "best", "atomic", "pass_sites", "shells", "doublet_kinds",
                "chem_flag", "lone_pair", "jt_ions"])
    for i in surv:
        ss = [s for s in by_icsd[i] if s["verdict"] == "PASS"]
        m = meta.get(i, {})
        w.writerow([i, m.get("material_id", ""), m.get("formula", ""),
                    m.get("space_group_number", ""), m.get("curated_flatband", ""),
                    m.get("best_flatband", ""), m.get("atomic_flatband", ""),
                    "|".join(s["site_label"] for s in ss),
                    "|".join(sorted({s["shell"] for s in ss})),
                    "|".join(sorted({k for s in ss for k in s["doublet_kinds"].split("|") if k})),
                    ss[0]["chem_flag"], ss[0]["lone_pair"], ss[0]["jt_ions"]])
print("\n-> stage1_pass_icsds.csv (%d rows, %d materials)" % (len(surv), len(mat(passing))))
