#!/usr/bin/env python3
"""Stage 1 — symmetry pre-filter over every harvested CIF.

Runs cifscan.scan() on all 8,168 ICSD CIFs and emits one CSV row per
(ICSD, site). Parse failures are recorded with their exception, never skipped.

  n_phi == 2  -> PASS    (one doublet-bearing site per primitive cell)
  n_phi  > 2  -> dilute  (Tc falls as 1/n_phi)
  no doublet  -> FAIL    (no 3-, 4- or 6-fold axis at the site)
"""
import csv, glob, os, re, sys, traceback
from concurrent.futures import ProcessPoolExecutor

import cifscan

CIFDIR = "/home/user/rtap/data/priority3_cifs/cif"
OUT = os.path.dirname(os.path.abspath(__file__))

# 6s^2 lone-pair cations (the classic stereochemically-active-lone-pair set)
LONEPAIR = {"Bi", "Pb", "Tl", "Sb", "Sn", "Te"}
# Jahn-Teller-active ions, matched on the CIF's own oxidation-state labels
JT_IONS = {"Cu2+", "Mn3+", "Ni3+", "Cr2+"}

DOUBLET_KINDS = {
    frozenset(("z2", "x2-y2")): "e_g{z2,x2-y2}",
    frozenset(("xz", "yz")): "{xz,yz}",
    frozenset(("xy", "x2-y2")): "{xy,x2-y2}",
    frozenset(("xy", "yz")): "{xy,yz}",
    frozenset(("xy", "xz")): "{xy,xz}",
    frozenset(("x", "y")): "p{x,y}",
}


def doublet_kind(label):
    return DOUBLET_KINDS.get(frozenset(label.split("+")), "other:" + label)


def atom_types(path):
    """Element set and oxidation-state-labelled ion set from the CIF."""
    txt = open(path, errors="replace").read()
    ions = set(re.findall(r"^\s*([A-Z][a-z]?\d*[+-])\s+-?[\d.]+\s*$", txt, re.M))
    elems = {re.match(r"([A-Z][a-z]?)", i).group(1) for i in ions}
    for m in re.finditer(r"_chemical_formula_(?:sum|structural)\s+'([^']*)'", txt):
        elems |= set(re.findall(r"([A-Z][a-z]?)", m.group(1)))
    return elems, ions


def one(path):
    icsd = os.path.splitext(os.path.basename(path))[0]
    try:
        rows = cifscan.scan(path, verbose=False)
        elems, ions = atom_types(path)
        lone = sorted(elems & LONEPAIR)
        jt = sorted(ions & JT_IONS)
        out = []
        for r in rows:
            d2 = r["d2"] or []
            p2 = r["p2"] or []
            kinds = sorted({doublet_kind(x) for x in d2} | {doublet_kind(x) for x in p2})
            nphi = r["nphi"]
            verdict = ("PASS" if nphi == 2 else "dilute" if nphi else "FAIL")
            out.append(dict(
                icsd=icsd, site_label=r["label"], element=cifscan._elem(r["label"]),
                shell=r["shell"], site_sym_order=r["order"], max_rot_order=r["maxrot"],
                orbit_conv=r["mult"], orbit_prim=int(round(r["prim"])),
                d_doublets="|".join(d2), p_doublets="|".join(p2),
                doublet_kinds="|".join(kinds), n_phi=("" if nphi is None else nphi),
                verdict=verdict, d_multiplets=r["d"],
                lone_pair="|".join(lone), jt_ions="|".join(jt),
                chem_flag=int(bool(lone or jt)),
            ))
        return icsd, out, None
    except Exception as e:
        return icsd, [], "%s: %s" % (type(e).__name__, e)


def main():
    files = sorted(glob.glob(os.path.join(CIFDIR, "*.cif")),
                   key=lambda p: int(os.path.splitext(os.path.basename(p))[0]))
    print("CIFs to scan: %d" % len(files), flush=True)
    rows, fails = [], []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        for i, (icsd, out, err) in enumerate(ex.map(one, files, chunksize=16), 1):
            if err:
                fails.append((icsd, err))
            rows.extend(out)
            if i % 1000 == 0:
                print("  scanned %d/%d" % (i, len(files)), flush=True)

    cols = ["icsd", "site_label", "element", "shell", "site_sym_order",
            "max_rot_order", "orbit_conv", "orbit_prim", "d_doublets",
            "p_doublets", "doublet_kinds", "n_phi", "verdict", "d_multiplets",
            "lone_pair", "jt_ions", "chem_flag"]
    p = os.path.join(OUT, "stage1_sites.csv")
    with open(p, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("-> %s  (%d site rows)" % (p, len(rows)))

    pf = os.path.join(OUT, "stage1_parse_failures.csv")
    with open(pf, "w", newline="", encoding="utf8") as f:
        w = csv.writer(f)
        w.writerow(["icsd", "exception"])
        w.writerows(fails)
    print("-> %s  (%d failures)" % (pf, len(fails)))
    for icsd, err in fails[:20]:
        print("   FAIL %s  %s" % (icsd, err))


if __name__ == "__main__":
    main()
