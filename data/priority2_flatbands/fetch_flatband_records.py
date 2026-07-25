#!/usr/bin/env python3
"""Fetch per-ICSD records from the Materials Flatband Database
(https://www.topologicalquantumchemistry.fr/flatbands/showmaterial.php?ICSD=...).

Each response is a JSON blob containing the material's identifiers, the
curated / best / atomic flat-band flags, the sibling ICSD entries of the same
material, the original ICSD CIF (with explicit symmetry operations) and a
POSCAR.

Seeds come from the SI table parse plus the site's own Curated/Atomic search
listings; the ICSD set is then closed under the `otherICSDs` field so that no
entry of a curated material is missed.

Writes:
  records/<ICSD>.json      trimmed record (band-structure plot data stripped)
  ../priority3_cifs/cif/<ICSD>.cif
  poscar/<ICSD>.POSCAR
  flatband_icsd_metadata.csv
"""
import csv, json, os, re, sys, threading, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, "records")
POS = os.path.join(HERE, "poscar")
CIF = os.path.abspath(os.path.join(HERE, "..", "priority3_cifs", "cif"))
for d in (REC, POS, CIF):
    os.makedirs(d, exist_ok=True)

URL = "https://www.topologicalquantumchemistry.fr/flatbands/showmaterial.php?ICSD=%s"
# large plot payloads we do not need to keep on disk
DROP = ("bsdata", "dosdata", "nosocbsdata", "nosocdosdata",
        "socbandsets", "nosocbandsets", "cif", "POSCAR")

lock = threading.Lock()
seen = set()
queue = []


def fetch_one(icsd):
    p = os.path.join(REC, "%s.json" % icsd)
    if os.path.exists(p):
        try:
            return json.load(open(p, encoding="utf8"))
        except Exception:
            pass
    last = None
    for attempt in range(6):
        try:
            req = urllib.request.Request(URL % icsd, headers={"User-Agent": "curl/8"})
            # short timeout: the server stalls rather than refuses when busy,
            # so it is better to give up quickly and retry than to block.
            raw = urllib.request.urlopen(req, timeout=45).read().decode("utf8", "replace")
            d = json.loads(raw)
            cif, poscar = d.get("cif"), d.get("POSCAR")
            if cif and "_atom_site" in cif:
                with open(os.path.join(CIF, "%s.cif" % icsd), "w", encoding="utf8") as f:
                    f.write(cif)
            if poscar:
                with open(os.path.join(POS, "%s.POSCAR" % icsd), "w", encoding="utf8") as f:
                    f.write(poscar)
            trim = {k: v for k, v in d.items() if k not in DROP}
            trim["_has_cif"] = bool(cif and "_atom_site" in cif)
            trim["_has_poscar"] = bool(poscar)
            with open(p, "w", encoding="utf8") as f:
                json.dump(trim, f)
            return trim
        except Exception as e:
            last = e
            time.sleep(min(2 ** attempt, 20))
    print("FAIL %s %s" % (icsd, last), flush=True)
    return None


def worker(icsd):
    d = fetch_one(icsd)
    if not d:
        return
    others = d.get("otherICSDs") or "[]"
    if isinstance(others, str):
        others = re.findall(r"\d+", others)
    new = []
    with lock:
        for o in others:
            o = str(o)
            if o not in seen:
                seen.add(o)
                new.append(o)
    return new


def main(seeds):
    global queue
    with lock:
        for s in seeds:
            if s not in seen:
                seen.add(s)
    queue = sorted(seen)
    done = 0
    with ThreadPoolExecutor(max_workers=10) as ex:
        while queue:
            batch, queue = queue, []
            for new in ex.map(worker, batch):
                done += 1
                if new:
                    queue.extend(new)
                if done % 100 == 0:
                    print("fetched %d, pending %d" % (done, len(queue)), flush=True)
            queue = sorted(set(queue))
    print("FETCH COMPLETE: %d ICSD records" % done, flush=True)


if __name__ == "__main__":
    seeds = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    print("seeds: %d" % len(seeds), flush=True)
    main(seeds)
