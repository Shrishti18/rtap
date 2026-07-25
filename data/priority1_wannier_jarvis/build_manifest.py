#!/usr/bin/env python3
"""Verify every JARVIS-WTB zip, extract the small per-material files, and join
the JARVIS-DFT metadata to produce the JVASP id -> formula/space-group map."""
import csv, json, os, re, zipfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = os.path.join(HERE, "zips")
SMALL = os.path.join(HERE, "extracted_small")
os.makedirs(SMALL, exist_ok=True)

# small files worth having unpacked next to the archives
KEEP = ("POSCAR", "wannier90.win", "INCAR", "KPOINTS", "OSZICAR")


def scan(fn):
    jid = fn[:-4]
    p = os.path.join(ZIPS, fn)
    try:
        z = zipfile.ZipFile(p)
        names = z.namelist()
    except Exception as e:
        return {"jid": jid, "zip": fn, "error": str(e)}
    info = {i.filename: i for i in z.infolist()}
    hr = next((n for n in names if n.endswith("_hr.dat")), "")
    d = os.path.join(SMALL, jid)
    os.makedirs(d, exist_ok=True)
    for k in KEEP:
        if k in info:
            with open(os.path.join(d, k), "wb") as f:
                f.write(z.read(k))
    rec = {
        "jid": jid,
        "zip": fn,
        "zip_bytes": os.path.getsize(p),
        "n_files_in_zip": len(names),
        "has_hr_dat": bool(hr),
        "hr_dat_name": hr,
        "hr_dat_bytes": info[hr].file_size if hr else 0,
        "has_win": "wannier90.win" in info,
        "has_wout": "wannier90.wout" in info,
        "has_eig": "wannier90.eig" in info,
        "has_poscar": "POSCAR" in info,
        "files": ";".join(sorted(names)),
        "error": "",
    }
    # wannier90.win records the projections / number of Wannier functions
    if "wannier90.win" in info:
        win = z.read("wannier90.win").decode("utf8", "replace")
        m = re.search(r"num_wann\s*[=:]\s*(\d+)", win, re.I)
        rec["num_wann"] = m.group(1) if m else ""
        m = re.search(r"num_bands\s*[=:]\s*(\d+)", win, re.I)
        rec["num_bands"] = m.group(1) if m else ""
    else:
        rec["num_wann"] = rec["num_bands"] = ""
    # the JVASP-*.json 'path' field carries the Materials Project id
    jf = "%s.json" % jid
    rec["mp_id"] = ""
    if jf in info:
        try:
            path = json.loads(z.read(jf)).get("path", "")
            m = re.search(r"(mp-\d+)", path)
            rec["mp_id"] = m.group(1) if m else ""
        except Exception:
            pass
    return rec


fns = sorted(os.listdir(ZIPS))
with ThreadPoolExecutor(max_workers=8) as ex:
    recs = list(ex.map(scan, fns))

# ---- join JARVIS-DFT metadata -------------------------------------------
meta = {}
for name, dim in (("jarvis_dft_3d.json", "3D"), ("jarvis_dft_2d.json", "2D")):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        continue
    for e in json.load(open(p)):
        meta.setdefault(e["jid"], {
            "formula": e.get("formula", ""),
            "spg_number": e.get("spg_number", ""),
            "spg_symbol": e.get("spg_symbol", ""),
            "icsd": e.get("icsd", ""),
            "jarvis_dimensionality": dim,
            "optb88vdw_bandgap": e.get("optb88vdw_bandgap", ""),
        })

article = {}
for f in json.load(open(os.path.join(HERE, "figshare_filelist.json"))):
    article[f["name"][:-4]] = f["article_title"]

for r in recs:
    m = meta.get(r["jid"], {})
    r["formula"] = m.get("formula", "")
    r["spg_number"] = m.get("spg_number", "")
    r["spg_symbol"] = m.get("spg_symbol", "")
    r["icsd"] = m.get("icsd", "")
    r["optb88vdw_bandgap"] = m.get("optb88vdw_bandgap", "")
    r["figshare_article"] = article.get(r["jid"], "")
    r["dimensionality"] = "2D" if "2D materials" in r["figshare_article"] else "3D"

cols = ["jid", "formula", "spg_number", "spg_symbol", "dimensionality", "icsd",
        "mp_id", "num_wann", "num_bands", "optb88vdw_bandgap", "zip",
        "zip_bytes", "n_files_in_zip", "has_hr_dat", "hr_dat_name",
        "hr_dat_bytes", "has_win", "has_wout", "has_eig", "has_poscar",
        "figshare_article", "error"]

out = os.path.join(HERE, "jarvis_wtb_manifest.csv")
with open(out, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(sorted(recs, key=lambda r: int(r["jid"].split("-")[1])))

# compact id -> formula/spacegroup map, which is what most users want
out2 = os.path.join(HERE, "jvasp_id_formula_spacegroup.csv")
with open(out2, "w", newline="", encoding="utf8") as f:
    w = csv.writer(f)
    w.writerow(["jid", "formula", "spg_number", "spg_symbol", "dimensionality",
                "icsd", "mp_id", "num_wann"])
    for r in sorted(recs, key=lambda r: int(r["jid"].split("-")[1])):
        w.writerow([r["jid"], r["formula"], r["spg_number"], r["spg_symbol"],
                    r["dimensionality"], r["icsd"], r["mp_id"], r["num_wann"]])

n = len(recs)
print("zips scanned      : %d" % n)
print("with _hr.dat      : %d" % sum(1 for r in recs if r["has_hr_dat"]))
print("with .win/.wout   : %d / %d" % (sum(1 for r in recs if r["has_win"]),
                                       sum(1 for r in recs if r["has_wout"])))
print("errors            : %d" % sum(1 for r in recs if r["error"]))
print("matched metadata  : %d (formula non-empty)"
      % sum(1 for r in recs if r["formula"]))
print("3D / 2D           : %d / %d"
      % (sum(1 for r in recs if r["dimensionality"] == "3D"),
         sum(1 for r in recs if r["dimensionality"] == "2D")))
print("hr.dat uncompressed total: %.1f GB"
      % (sum(r["hr_dat_bytes"] for r in recs) / 1e9))
print("->", out, "\n->", out2)
