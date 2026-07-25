#!/usr/bin/env python3
"""Extract Wannier90 files from the JARVIS-WTB archives.

The 1,772 archives hold 69.5 GB of uncompressed `wannier90_hr.dat` (12.6 GB
zipped), which is why they are shipped here as-is. Use this to unpack all of
them, or a subset, once you have the disk for it.

  python3 extract_hr_dat.py --out /big/disk/wtb                # every material
  python3 extract_hr_dat.py --out ./sub --jid JVASP-5 JVASP-8  # named ones
  python3 extract_hr_dat.py --out ./sub --formula FeSe Bi2Te3  # by formula
  python3 extract_hr_dat.py --out ./sub --jid JVASP-5 --all-files

By default only `wannier90_hr.dat` is written, as `<out>/<JVASP-id>_hr.dat`.
"""
import argparse, csv, os, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ZIPS = os.path.join(HERE, "zips")
MAP = os.path.join(HERE, "jvasp_id_formula_spacegroup.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--jid", nargs="*", default=[])
    ap.add_argument("--formula", nargs="*", default=[])
    ap.add_argument("--all-files", action="store_true",
                    help="extract the whole archive, not just wannier90_hr.dat")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(MAP)))
    sel = set(a.jid)
    if a.formula:
        want = {f.lower() for f in a.formula}
        sel |= {r["jid"] for r in rows if r["formula"].lower() in want}
    if not sel:
        sel = {r["jid"] for r in rows}

    os.makedirs(a.out, exist_ok=True)
    n = 0
    for jid in sorted(sel, key=lambda j: int(j.split("-")[1])):
        p = os.path.join(ZIPS, "%s.zip" % jid)
        if not os.path.exists(p):
            print("missing archive: %s" % jid, file=sys.stderr)
            continue
        z = zipfile.ZipFile(p)
        if a.all_files:
            z.extractall(os.path.join(a.out, jid))
        else:
            name = next((x for x in z.namelist() if x.endswith("_hr.dat")), None)
            if not name:
                print("no _hr.dat in %s" % jid, file=sys.stderr)
                continue
            with open(os.path.join(a.out, "%s_hr.dat" % jid), "wb") as f:
                f.write(z.read(name))
        n += 1
        if n % 100 == 0:
            print("extracted %d/%d" % (n, len(sel)), flush=True)
    print("done: %d materials -> %s" % (n, a.out))


if __name__ == "__main__":
    main()
