#!/usr/bin/env python3
"""Download all JARVIS-WTB zips from figshare with md5 verification + resume."""
import json, os, sys, hashlib, urllib.request, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.join(HERE, "zips")
os.makedirs(DEST, exist_ok=True)
files = json.load(open(os.path.join(HERE, "figshare_filelist.json")))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def get(rec):
    p = os.path.join(DEST, rec["name"])
    if os.path.exists(p) and os.path.getsize(p) == rec["size"]:
        if not rec["md5"] or md5(p) == rec["md5"]:
            return ("cached", rec["name"])
    for attempt in range(5):
        try:
            req = urllib.request.Request(rec["url"], headers={"User-Agent": "curl/8"})
            with urllib.request.urlopen(req, timeout=600) as r, open(p, "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(p) != rec["size"]:
                raise IOError("size mismatch %d != %d" % (os.path.getsize(p), rec["size"]))
            if rec["md5"] and md5(p) != rec["md5"]:
                raise IOError("md5 mismatch")
            return ("ok", rec["name"])
        except Exception as e:
            if attempt == 4:
                return ("FAIL:%s" % e, rec["name"])
            time.sleep(2 ** attempt)


done = 0
with ThreadPoolExecutor(max_workers=8) as ex:
    for status, name in ex.map(get, files):
        done += 1
        if status.startswith("FAIL"):
            print("FAIL %s %s" % (name, status), flush=True)
        if done % 50 == 0:
            print("progress %d/%d" % (done, len(files)), flush=True)
print("DOWNLOAD COMPLETE %d files" % done, flush=True)
