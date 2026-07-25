#!/usr/bin/env python3
"""Stage 2a — cross-match the Stage-1 PASS list against the 1,772 JARVIS
Wannier materials, by ICSD id where available, else reduced composition +
space group."""
import csv, math, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
MAN = "/home/user/rtap/data/priority1_wannier_jarvis/jarvis_wtb_manifest.csv"


def composition(f):
    f = (f or "").replace(" ", "")
    if not f:
        return {}

    def parse(s, i=0):
        c, n = collections.defaultdict(float), len(s)
        while i < n:
            if s[i] == "(":
                sub, i = parse(s, i + 1)
                m = re.match(r"\d+(?:\.\d+)?", s[i:])
                mult = float(m.group(0)) if m else 1.0
                i += len(m.group(0)) if m else 0
                for k, v in sub.items():
                    c[k] += v * mult
            elif s[i] == ")":
                return c, i + 1
            else:
                m = re.match(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)?", s[i:])
                if not m:
                    i += 1
                    continue
                c[m.group(1)] += float(m.group(2)) if m.group(2) else 1.0
                i += len(m.group(0))
        return c, i
    return {k: v for k, v in parse(f)[0].items() if v}


def redsig(f):
    c = composition(f)
    if not c:
        return None
    if all(abs(v - round(v)) < 1e-6 for v in c.values()):
        g = 0
        for v in c.values():
            g = math.gcd(g, int(round(v)))
        g = g or 1
        return tuple(sorted((k, int(round(v)) // g) for k, v in c.items()))
    m = min(c.values())
    return tuple(sorted((k, round(v / m, 3)) for k, v in c.items()))


jar = list(csv.DictReader(open(MAN)))
by_icsd, by_sig = {}, collections.defaultdict(list)
for j in jar:
    for i in (j["icsd"] or "").replace(";", ",").split(","):
        i = i.strip()
        if i.isdigit():
            by_icsd.setdefault(i, j)
    s = redsig(j["formula"])
    if s:
        by_sig[(s, str(j["spg_number"]).strip())].append(j)

rows = list(csv.DictReader(open(os.path.join(HERE, "stage1_pass_icsds.csv"))))
print("Stage-1 PASS: %d ICSD entries, %d materials"
      % (len(rows), len({r["material_id"] for r in rows})))

out, hit_icsd, hit_sig = [], 0, 0
for r in rows:
    j = by_icsd.get(r["icsd"])
    how = "icsd" if j else ""
    if not j:
        cands = by_sig.get((redsig(r["formula"]), str(r["space_group_number"]).strip()), [])
        if len(cands) >= 1:
            j = cands[0]
            how = "composition+spg"
    if j:
        hit_icsd += how == "icsd"
        hit_sig += how == "composition+spg"
        out.append(dict(r, jid=j["jid"], match=how, jarvis_formula=j["formula"],
                        jarvis_spg=j["spg_number"], dimensionality=j["dimensionality"],
                        num_wann=j["num_wann"], hr_bytes=j["hr_dat_bytes"],
                        has_wout=j["has_wout"]))

print("matched: %d ICSD rows (%d by ICSD id, %d by composition+spg)"
      % (len(out), hit_icsd, hit_sig))
jids = sorted({o["jid"] for o in out}, key=lambda s: int(s.split("-")[1]))
mats = {o["material_id"] for o in out}
print("distinct JARVIS materials with a Wannier Hamiltonian: %d" % len(jids))
print("distinct flat-band materials covered: %d" % len(mats))
print("of which top-345: %d | curated: %d | atomic: %d"
      % (len({o["material_id"] for o in out if o["best"] == "1"}),
         len({o["material_id"] for o in out if o["curated"] == "1"}),
         len({o["material_id"] for o in out if o["atomic"] == "1"})))
tot = sum(int(o["hr_bytes"]) for o in {o["jid"]: o for o in out}.values())
print("total wannier90_hr.dat to extract: %.2f GB" % (tot / 1e9))

cols = list(out[0].keys())
p = os.path.join(HERE, "stage2_matched.csv")
with open(p, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(out)
print("-> %s" % p)
with open(os.path.join(HERE, "stage2_jids.txt"), "w") as f:
    f.write("\n".join(jids))
print("-> stage2_jids.txt (%d ids)" % len(jids))
