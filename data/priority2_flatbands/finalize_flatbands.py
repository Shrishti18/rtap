#!/usr/bin/env python3
"""Turn the fetched per-ICSD records into the deliverable candidate lists, and
validate the CIFs harvested alongside them.

Materials are grouped from the `otherICSDs` field (the site's own grouping of
ICSD entries that describe the same material), so the material-level lists do
not depend on the PDF table parse. The SI Table XI parse is then repaired
against this grouping.
"""
import ast, csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, "records")
CIF = os.path.abspath(os.path.join(HERE, "..", "priority3_cifs", "cif"))
SYMOP = ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz")


def strip_tags(s):
    return re.sub(r"<[^>]+>", "", s or "").strip()


def composition(f):
    """Element -> count, so that 'NaAlSi3O8' and 'Al1 Na1 O8 Si3' compare equal.
    Parenthesised groups are expanded, e.g. 'Sr3(PbNiO6)' and 'K2(W2O7)'."""
    f = (f or "").replace(" ", "")
    if not f:
        return {}

    def parse(s, i=0):
        counts, n = defaultdict(float), len(s)
        while i < n:
            if s[i] == "(":
                sub, i = parse(s, i + 1)
                m = re.match(r"\d+(?:\.\d+)?", s[i:])
                mult = float(m.group(0)) if m else 1.0
                i += len(m.group(0)) if m else 0
                for k, v in sub.items():
                    counts[k] += v * mult
            elif s[i] == ")":
                return counts, i + 1
            else:
                m = re.match(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)?", s[i:])
                if not m:
                    i += 1
                    continue
                counts[m.group(1)] += float(m.group(2)) if m.group(2) else 1.0
                i += len(m.group(0))
        return counts, i

    c = parse(f)[0]
    return {k: v for k, v in c.items() if v}


def sig(formula, sgnum):
    """Canonical (reduced composition, space group) key for matching a table row
    against a site material, independent of how the formula was written."""
    c = composition(formula)
    if not c or not str(sgnum).strip():
        return None
    ints = all(abs(v - round(v)) < 1e-6 for v in c.values())
    if ints:
        vals = [int(round(v)) for v in c.values()]
        g = 0
        for v in vals:
            g = math.gcd(g, v)
        g = g or 1
        items = tuple(sorted((k, int(round(v)) // g) for k, v in c.items()))
    else:
        m = min(c.values())
        items = tuple(sorted((k, round(v / m, 3)) for k, v in c.items()))
    return (items, str(sgnum).strip())


def same_composition(a, b):
    ca, cb = composition(a), composition(b)
    if not ca or not cb or set(ca) != set(cb):
        return False
    # compare up to an overall multiplicative factor (Z differences)
    k0 = sorted(ca)[0]
    r = cb[k0] / ca[k0]
    return all(abs(cb[k] - ca[k] * r) < 1e-6 for k in ca)


def load():
    recs = {}
    for p in glob.glob(os.path.join(REC, "*.json")):
        try:
            d = json.load(open(p, encoding="utf8"))
        except Exception:
            continue
        recs[str(d["ICSD"])] = d
    return recs


def group(recs):
    """Union ICSD entries that belong to the same material."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for icsd, d in recs.items():
        find(icsd)
        others = d.get("otherICSDs") or "[]"
        if isinstance(others, str):
            others = re.findall(r"\d+", others)
        for o in others:
            o = str(o)
            if o in recs:
                union(icsd, o)
    g = defaultdict(list)
    for icsd in recs:
        g[find(icsd)].append(icsd)
    return {min(v, key=int): sorted(v, key=int) for v in g.values()}


def main():
    recs = load()
    print("records loaded: %d" % len(recs))
    groups = group(recs)
    print("materials (grouped by otherICSDs): %d" % len(groups))

    # ---------------- per-ICSD metadata ---------------------------------
    rows = []
    for icsd, d in sorted(recs.items(), key=lambda kv: int(kv[0])):
        cifp = os.path.join(CIF, "%s.cif" % icsd)
        has_cif = os.path.exists(cifp)
        symops = 0
        if has_cif:
            txt = open(cifp, encoding="utf8", errors="replace").read()
            if any(k in txt for k in SYMOP):
                # count operation lines inside the symmetry loop
                m = re.search(r"(%s)(.*?)(?=\nloop_|\n_[a-zA-Z])" % "|".join(SYMOP),
                              txt, re.S)
                if m:
                    symops = len(re.findall(r"['\"][^'\"]*[xyz][^'\"]*['\"]", m.group(2)))
        rows.append({
            "icsd": icsd,
            "formula": d.get("CompoundChem", ""),
            "space_group_number": d.get("SG", ""),
            "space_group_symbol": strip_tags(d.get("SGname", "")),
            "point_group": d.get("PGname", ""),
            "material_id": next(k for k, v in groups.items() if icsd in v),
            "curated_flatband": d.get("curatedflatband", 0),
            "best_flatband": d.get("bestflatband", 0),
            "atomic_flatband": d.get("atomicflatband", 0),
            "high_quality": d.get("highquality", 0),
            "topology_soc": d.get("toposubclass", ""),
            "topology_nosoc": d.get("nosoctoposubclass", ""),
            "superconductor": d.get("superconductor", 0),
            "nonzero_magnetisation": d.get("nonzeromag", 0),
            "theory_entry": d.get("theory", 0),
            "n_electrons": d.get("nbrelectrons", ""),
            "structure_type": d.get("structuretype", ""),
            "a": d.get("celllengtha", ""), "b": d.get("celllengthb", ""),
            "c": d.get("celllengthc", ""),
            "alpha": d.get("cellanglealpha", ""), "beta": d.get("cellanglebeta", ""),
            "gamma": d.get("cellanglegamma", ""),
            "n_kagome": d.get("nbrrigorouskagomes", 0),
            "n_pyrochlore": d.get("nbrpyrochlores", 0),
            "n_lieb": d.get("nbrrigorousliebs", 0),
            "n_bipartite": d.get("nbrbipartites", 0),
            "n_split": d.get("nbrsplits", 0),
            "materials_project": d.get("linkmaterialsproject", ""),
            "cif_file": "cif/%s.cif" % icsd if has_cif else "",
            "cif_symmetry_ops": symops,
        })
    out = os.path.join(HERE, "flatband_icsd_metadata.csv")
    with open(out, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("-> %s (%d rows)" % (out, len(rows)))

    by_icsd = {r["icsd"]: r for r in rows}

    # ---------------- material-level lists ------------------------------
    def material_rows(flag):
        out = []
        for mid, members in sorted(groups.items(), key=lambda kv: int(kv[0])):
            hit = [m for m in members if int(by_icsd[m][flag] or 0) == 1]
            if not hit:
                continue
            rep = by_icsd[mid] if mid in by_icsd else by_icsd[hit[0]]
            withcif = [m for m in members if by_icsd[m]["cif_file"]]
            out.append({
                "material_id": mid,
                "formula": rep["formula"],
                "space_group_number": rep["space_group_number"],
                "space_group_symbol": rep["space_group_symbol"],
                "topology_soc": rep["topology_soc"],
                "high_quality": rep["high_quality"],
                "superconductor": rep["superconductor"],
                "nonzero_magnetisation": rep["nonzero_magnetisation"],
                "is_best_flatband": int(any(int(by_icsd[m]["best_flatband"] or 0)
                                            for m in members)),
                "is_curated_flatband": int(any(int(by_icsd[m]["curated_flatband"] or 0)
                                               for m in members)),
                "is_atomic_flatband": int(any(int(by_icsd[m]["atomic_flatband"] or 0)
                                              for m in members)),
                "n_icsd": len(members),
                "icsd_ids": " ".join(members),
                "n_icsd_with_cif": len(withcif),
                "materials_project": rep["materials_project"],
            })
        return out

    stats = {}
    for flag, name in (("curated_flatband", "materials_curated_flatband"),
                       ("best_flatband", "materials_best_flatband"),
                       ("atomic_flatband", "materials_atomic_flatband")):
        mr = material_rows(flag)
        p = os.path.join(HERE, "site_%s.csv" % name)
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(mr[0].keys()))
            w.writeheader()
            w.writerows(mr)
        nic = sum(1 for r in rows if int(r[flag] or 0) == 1)
        stats[flag] = (len(mr), nic)
        print("-> %s : %d materials, %d ICSD entries" % (p, len(mr), nic))

    # ---------------- repair the SI Table XI ICSD assignment -------------
    xi = os.path.join(HERE, "tableXI_curated_flatband_materials.csv")
    if os.path.exists(xi):
        trows = list(csv.DictReader(open(xi)))
        # Index the site's material groups by (reduced composition, space group).
        index = defaultdict(list)
        for mid in groups:
            rep = by_icsd[mid]
            k = sig(rep["formula"], rep["space_group_number"])
            if k:
                index[k].append(mid)
        claimed = set()

        # Pass 1: unambiguous composition+space-group match, preferring the
        # candidate that the row's own parsed ICSD ids point at.
        for r in trows:
            r["material_id"] = ""
            r["verification"] = "unmatched"
            ids = [i for i in r["icsd_ids"].split() if i in by_icsd]
            cands = [m for m in index.get(sig(r["formula"], r["space_group_number"]) or (), [])
                     if m not in claimed]
            if not cands:
                continue
            pick = None
            if len(cands) > 1 and ids:
                want = {by_icsd[i]["material_id"] for i in ids}
                pref = [m for m in cands if m in want]
                pick = pref[0] if len(pref) == 1 else None
            if pick is None:
                pick = cands[0] if len(cands) == 1 else None
            if pick is None:
                continue
            claimed.add(pick)
            r["material_id"] = pick
            r["verification"] = "matched"

        # Pass 2: rows still unmatched fall back to their parsed ICSD ids, but
        # only when those ids agree with the row's own formula.
        for r in trows:
            if r["material_id"]:
                continue
            ids = [i for i in r["icsd_ids"].split() if i in by_icsd]
            votes = defaultdict(int)
            for i in ids:
                votes[by_icsd[i]["material_id"]] += 1
            for mid in sorted(votes, key=lambda m: -votes[m]):
                if mid in claimed:
                    continue
                if same_composition(r["formula"], by_icsd[mid]["formula"]):
                    claimed.add(mid)
                    r["material_id"] = mid
                    r["verification"] = "matched_by_icsd"
                    break

        for r in trows:
            mid = r["material_id"]
            if mid:
                rep = by_icsd[mid]
                r["icsd_ids_verified"] = " ".join(groups[mid])
                r["n_icsd_verified"] = len(groups[mid])
                r["site_formula"] = rep["formula"]
                r["site_space_group_number"] = rep["space_group_number"]
            else:
                r["icsd_ids_verified"] = ""
                r["n_icsd_verified"] = 0
                r["site_formula"] = ""
                r["site_space_group_number"] = ""
        cols = list(trows[0].keys())
        p = os.path.join(HERE, "tableXI_curated_flatband_materials_verified.csv")
        with open(p, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(trows)
        import collections as _c
        vc = _c.Counter(r["verification"] for r in trows)
        ok = vc["matched"]
        tot = sum(int(r["n_icsd_verified"]) for r in trows)
        print("-> %s : %d/%d rows cleanly matched, %d ICSD entries"
              % (p, ok, len(trows), tot))
        for k, v in sorted(vc.items()):
            print("     %-30s %d" % (k, v))

    # ---------------- CIF validation -------------------------------------
    cifs = glob.glob(os.path.join(CIF, "*.cif"))
    withops = sum(1 for r in rows if r["cif_symmetry_ops"] > 0)
    print("\nCIFs on disk            : %d" % len(cifs))
    print("CIFs with explicit symops: %d" % withops)
    print("CIFs WITHOUT symops      : %d"
          % sum(1 for r in rows if r["cif_file"] and r["cif_symmetry_ops"] == 0))
    print("ICSD entries with no CIF : %d" % sum(1 for r in rows if not r["cif_file"]))


if __name__ == "__main__":
    main()
