#!/usr/bin/env python3
"""Parse Supplementary Tables XI and XII of Regnault et al., Nature 603, 824 (2022)
from the Nature SI PDF (41586_2022_4519_MOESM1_ESM.pdf) using pdftotext word boxes.

Table XI  (pp. 48-94):  2,379 curated flat-band materials      (10 columns)
Table XII (pp. 95-114): 1,102 curated flat ATOMIC band mats.   ( 9 columns, no "best")

Two things make -layout output unusable and are handled here:
  * column x-positions are re-typeset on every page, so cells are found by
    intra-row x-gap clustering rather than fixed x ranges;
  * the ICSD cell of a row is a multi-line block vertically *centred* on the
    row, so ICSD text lines are clustered by vertical gap and each cluster is
    assigned to its nearest row anchor.
"""
import csv, os, re, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "SI_41586_2022_4519_MOESM1_ESM.pdf")

WORD_RE = re.compile(
    r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>')
ENT = {"&amp;": "&", "&lt;": "<", "&gt;": ">", "&quot;": '"', "&apos;": "'"}

SG_NUM_RE = re.compile(r"^\(\d{1,3}\)$")
ICSD_RE = re.compile(r"^(\d{2,7})(\([^)]*\))?$")
DASHES = {"—", "–", "-", "—", "–"}

CELL_GAP = 6.0        # x-gap (pt) that separates two table cells within a row
ICSD_Y_OFFSET = 1.4   # ICSD text lines sit ~1.4pt below their row anchor
ICSD_CLUSTER_GAP = 6.6  # intra-cell line pitch ~6.0pt, inter-cell gap >7pt
ICSD_MIN_X = 300.0    # ICSD column is always well right of the formula/SG


def unescape(s):
    for k, v in ENT.items():
        s = s.replace(k, v)
    return s


def page_words(page):
    xml = subprocess.run(["pdftotext", "-f", str(page), "-l", str(page),
                          "-bbox-layout", PDF, "-"],
                         capture_output=True, text=True).stdout
    out = []
    for x0, y0, x1, y1, t in WORD_RE.findall(xml):
        t = unescape(t).strip()
        if t:
            out.append((float(x0), float(y0), float(x1), float(y1), t))
    return out


def split_cells(words):
    """Group x-sorted words of one row into cells using the x-gap."""
    cells, cur = [], []
    for w in words:
        if cur and w[0] - max(c[2] for c in cur) > CELL_GAP:
            cells.append(cur)
            cur = []
        cur.append(w)
    if cur:
        cells.append(cur)
    return cells


def cell_text(cell):
    return "".join(w[4] for w in sorted(cell, key=lambda w: w[0]))


def parse_page(page, expect_next, ncat):
    """ncat = number of categorical cells after the space-group cell."""
    words = page_words(page)
    anchors, anchor_words = [], []
    for w in sorted((w for w in words if w[0] < 95 and w[1] > 40),
                    key=lambda w: w[1]):
        if re.fullmatch(r"\d{1,4}", w[4]) and int(w[4]) == expect_next + len(anchors):
            anchors.append((int(w[4]), w[1]))
            anchor_words.append(w)
    if not anchors:
        return [], expect_next
    anchor_set = set(anchor_words)

    # ---- ICSD blocks: cluster by y, assign each cluster to nearest anchor ----
    by_line = {}
    for w in words:
        by_line.setdefault(round(w[1], 2), []).append(w)
    # Restrict to the vertical extent of the table body, so that the printed
    # page number (top-right) and any body text above the table are excluded.
    ytop = anchors[0][1] - 5.0
    ybot = anchors[-1][1] + 15.0
    icsd_words = []
    for y, lw in by_line.items():
        lw.sort(key=lambda w: w[0])
        for i, w in enumerate(lw):
            if w[0] < ICSD_MIN_X or not ICSD_RE.match(w[4]):
                continue
            if not (ytop <= w[1] <= ybot):
                continue
            # "Fig. 43" -> 43 is a figure number, not an ICSD id
            if i and lw[i - 1][4].rstrip(".").lower() == "fig":
                continue
            icsd_words.append(w)

    ylines = sorted({round(w[1], 2) for w in icsd_words})
    clusters, cur = [], []
    for y in ylines:
        if cur and y - cur[-1] > ICSD_CLUSTER_GAP:
            clusters.append(cur)
            cur = []
        cur.append(y)
    if cur:
        clusters.append(cur)

    icsd_for = {n: [] for n, _ in anchors}
    for cl in clusters:
        centre = (cl[0] + cl[-1]) / 2.0 - ICSD_Y_OFFSET
        n = min(anchors, key=lambda a: abs(a[1] - centre))[0]
        ys = set(cl)
        for w in sorted((w for w in icsd_words if round(w[1], 2) in ys),
                        key=lambda w: (w[1], w[0])):
            icsd_for[n].append(w[4])

    # ---- row cells ---------------------------------------------------------
    rows, problems = [], []
    icsd_ids = {round(w[1], 2) for w in icsd_words}
    for n, y in anchors:
        # Long formulas overflow left of the Num column, so exclude the row
        # anchors themselves rather than everything left of a fixed x.
        band = [w for w in words
                if y - 2.0 <= w[1] <= y + 6.5
                and w not in anchor_set
                and not (w[0] >= ICSD_MIN_X and ICSD_RE.match(w[4])
                         and w in icsd_words)]
        band.sort(key=lambda w: w[0])
        cells = split_cells(band)
        sg_i = next((i for i, c in enumerate(cells)
                     if any(SG_NUM_RE.match(w[4]) for w in c)), None)
        if sg_i is None:
            problems.append(n)
            continue
        formula = "".join(cell_text(c) for c in cells[:sg_i])
        sg = cell_text(cells[sg_i])
        cat = [cell_text(c) for c in cells[sg_i + 1:]][:ncat]
        cat += [""] * (ncat - len(cat))
        rows.append((n, formula, sg, cat, icsd_for[n]))
    return (rows, problems), anchors[-1][0] + 1


def clean(n, formula, sg, cat, icsd, has_best):
    m = re.search(r"\((\d{1,3})\)$", sg)
    sg_num = int(m.group(1)) if m else ""
    sg_sym = re.sub(r"\(\d{1,3}\)$", "", sg).strip()
    d = lambda s: "" if s.strip() in DASHES else s.strip()
    if has_best:
        topo, best, sub, mag, sc, hq = cat
    else:
        topo, sub, mag, sc, hq = cat
        best = ""
    ids, labels = [], []
    for tok in icsd:
        mm = ICSD_RE.match(tok)
        ids.append(mm.group(1))
        if mm.group(2):
            labels.append(mm.group(1) + mm.group(2))
    best = d(best)
    rec = {
        "num": n,
        "formula": formula,
        "space_group_symbol": sg_sym,
        "space_group_number": sg_num,
        "topology_at_Ef": d(topo),
        "sublattices": d(sub).strip("()"),
        "magnetic": d(mag),
        "superconductor": d(sc),
        "high_quality": d(hq),
        "n_icsd": len(ids),
        "icsd_ids": " ".join(ids),
        "icsd_extra_sublattice_labels": " ".join(labels),
    }
    if has_best:
        rec["is_top_candidate"] = "yes" if best.lower().startswith("fig") else "no"
        rec["best_figure"] = best
    return rec


def run(first, last, expect_rows, ncat, has_best, outfile, label):
    allrows, allprob, nxt = [], [], 1
    for p in range(first, last + 1):
        (rows, prob), nxt = parse_page(p, nxt, ncat)
        allrows.extend(rows)
        allprob.extend(prob)
    recs = [clean(*r, has_best=has_best) for r in allrows]
    nums = [r["num"] for r in recs]
    assert nums == list(range(1, len(nums) + 1)), "row numbering broke: %s" % allprob
    order = ["num", "formula", "space_group_symbol", "space_group_number",
             "topology_at_Ef"]
    if has_best:
        order += ["is_top_candidate", "best_figure"]
    order += ["sublattices", "magnetic", "superconductor", "high_quality",
              "n_icsd", "icsd_ids", "icsd_extra_sublattice_labels"]
    with open(outfile, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=order)
        w.writeheader()
        w.writerows(recs)
    tot = sum(r["n_icsd"] for r in recs)
    uniq = len({i for r in recs for i in r["icsd_ids"].split()})
    print("%s\n   rows=%d (expect %d)  ICSD=%d (%d unique)  unparsed=%d"
          % (label, len(recs), expect_rows, tot, uniq, len(allprob)))
    if has_best:
        print("   top candidates flagged: %d"
              % sum(1 for r in recs if r["is_top_candidate"] == "yes"))
    return recs


if __name__ == "__main__":
    xi = run(48, 94, 2379, 6, True,
             os.path.join(HERE, "tableXI_curated_flatband_materials.csv"),
             "Table XI  - curated flat-band materials")
    xii = run(95, 114, 1102, 5, False,
              os.path.join(HERE, "tableXII_curated_flat_atomic_band_materials.csv"),
              "Table XII - curated flat ATOMIC band materials")
    top = [r for r in xi if r["is_top_candidate"] == "yes"]
    out = os.path.join(HERE, "top_candidates_best_flatbands.csv")
    with open(out, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(xi[0].keys()))
        w.writeheader()
        w.writerows(top)
    print("Top-candidate subset: %d rows -> %s" % (len(top), out))
