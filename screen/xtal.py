"""CIF -> MX4 tetrahedra, X-M-X angles, and connectivity dimensionality.

Everything is computed from the CIF itself: the asymmetric unit is expanded by
the symmetry operations in the file, neighbours are found over periodic images,
and the dimensionality of the M-X network is obtained from the rank of the
lattice translations spanned by one connected component.
"""
import itertools, re

import numpy as np

# Cordero et al. (2008) covalent radii, Angstrom; high-spin values for the
# 3d metals so the neighbour cutoff is not accidentally tight.
RCOV = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.61, "Fe": 1.52, "Co": 1.50,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38,
    "I": 1.39, "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03,
    "Nd": 2.01, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Th": 2.06,
    "U": 1.96,
}
SYMLOOP = re.compile(
    r"(?:_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz)\s*\n"
    r"(.*?)(?=\nloop_|\n\s*_[a-zA-Z]|\Z)", re.S)
OPCHARS = re.compile(r"^[-+0-9xyzXYZ/,.\s]+$")


def _num(s):
    return float(re.sub(r"\(.*?\)", "", s).strip())


def _parse_op(s):
    W = np.zeros((3, 3))
    w = np.zeros(3)
    for i, comp in enumerate(s.replace(" ", "").lower().split(",")):
        for m in re.finditer(r"([+-]?)(\d+/\d+|\d*\.?\d*)\*?([xyz]?)", comp):
            sgn, num, var = m.group(1), m.group(2), m.group(3)
            if not num and not var:
                continue
            if num == "":
                v = 1.0
            elif "/" in num:
                p, q = num.split("/")
                v = float(p) / float(q)
            else:
                v = float(num)
            if sgn == "-":
                v = -v
            if var:
                W[i, "xyz".index(var)] += v
            else:
                w[i] += v
    return W, w


def read_cif(path):
    txt = open(path, errors="replace").read()
    cell = {}
    for k in ("a", "b", "c"):
        m = re.search(r"_cell_length_%s\s+(\S+)" % k, txt)
        cell[k] = _num(m.group(1))
    for k, nm in (("al", "alpha"), ("be", "beta"), ("ga", "gamma")):
        m = re.search(r"_cell_angle_%s\s+(\S+)" % nm, txt)
        cell[k] = _num(m.group(1))
    ops = []
    for m in SYMLOOP.finditer(txt):
        for line in m.group(1).splitlines():
            line = line.strip()
            if not line:
                continue
            q = re.findall(r"'([^']*)'|\"([^\"]*)\"", line)
            cand = [a or b for a, b in q] or [re.sub(r"^\d+\s+", "", line)]
            for s in cand:
                if s.count(",") == 2 and OPCHARS.match(s):
                    ops.append(_parse_op(s))
    if not ops:
        ops = [(np.eye(3), np.zeros(3))]
    # atom site loop
    sites = []
    lines = txt.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            hdr, j = [], i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                hdr.append(lines[j].strip())
                j += 1
            if any("_atom_site_fract_x" in h for h in hdr):
                idx = {h: k for k, h in enumerate(hdr)}
                cx, cy, cz = (idx["_atom_site_fract_x"], idx["_atom_site_fract_y"],
                              idx["_atom_site_fract_z"])
                clab = idx.get("_atom_site_label", 0)
                ctyp = idx.get("_atom_site_type_symbol", clab)
                cocc = idx.get("_atom_site_occupancy")
                while (j < len(lines) and lines[j].strip()
                       and not lines[j].strip().startswith(("_", "loop_", "#"))):
                    t = lines[j].split()
                    if len(t) > max(cx, cy, cz):
                        el = re.match(r"([A-Z][a-z]?)", t[ctyp])
                        occ = _num(t[cocc]) if (cocc is not None and len(t) > cocc) else 1.0
                        sites.append(dict(label=t[clab],
                                          el=el.group(1) if el else t[ctyp],
                                          f=np.array([_num(t[cx]), _num(t[cy]), _num(t[cz])]),
                                          occ=occ))
                    j += 1
            i = j
        else:
            i += 1
    return cell, ops, sites


def cart_matrix(c):
    a, b, cc = c["a"], c["b"], c["c"]
    al, be, ga = np.radians([c["al"], c["be"], c["ga"]])
    v1 = np.array([a, 0, 0])
    v2 = np.array([b * np.cos(ga), b * np.sin(ga), 0])
    cx = cc * np.cos(be)
    cy = cc * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    v3 = np.array([cx, cy, np.sqrt(max(cc ** 2 - cx ** 2 - cy ** 2, 0.0))])
    return np.array([v1, v2, v3])          # ROWS = lattice vectors


def expand(cell, ops, sites, tol=1e-4):
    """Asymmetric unit -> all atoms in one cell (deduplicated)."""
    out = []
    for s in sites:
        for W, w in ops:
            y = (W @ s["f"] + w) % 1.0
            if not any(o["el"] == s["el"] and
                       np.allclose(np.minimum(np.abs(y - o["f"]),
                                              1 - np.abs(y - o["f"])), 0, atol=tol)
                       for o in out):
                out.append(dict(s, f=y))
    return out


def neighbours(atoms, A, i, want_els, rmax=3.8):
    """Neighbours of atom i drawn from `want_els`, over periodic images."""
    res = []
    fi = atoms[i]["f"]
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    for j, at in enumerate(atoms):
        if at["el"] not in want_els:
            continue
        d = at["f"] + shifts - fi
        cart = d @ A
        r = np.linalg.norm(cart, axis=1)
        for k, rr in enumerate(r):
            if 0.3 < rr <= rmax:
                res.append((rr, j, tuple(shifts[k]), cart[k]))
    res.sort(key=lambda t: t[0])
    return res


def tetrahedron(atoms, A, i, want_els, gap=1.18, slack=1.15):
    """The 4 nearest `want_els` neighbours, if they form a genuine MX4
    tetrahedron. Returns (neighbours, shape) or (None, reason).

    Two rejections matter in practice. A distance cutoff of 1.15*(r_M + r_X)
    keeps out shells that are merely the 4 closest atoms of a species that is
    not bonded at all (Co-Se at 3.3 A in an oxide, say). And any shell with an
    angle above 150 deg is square-planar or trans-coordinated rather than
    tetrahedral -- a real tetrahedron cannot have a trans pair.
    """
    nb = neighbours(atoms, A, i, want_els)
    if len(nb) < 4:
        return None, "fewer than 4 %s neighbours" % "/".join(sorted(want_els))
    rM = RCOV.get(atoms[i]["el"], 1.4)
    d4 = nb[3][0]
    lim = slack * (rM + max(RCOV.get(e, 1.2) for e in want_els))
    if d4 > lim:
        return None, "4th neighbour at %.2f A > cutoff %.2f A (not bonded)" % (d4, lim)
    if len(nb) > 4 and nb[4][0] < gap * d4:
        return None, "5th neighbour at %.2f A (coordination > 4)" % nb[4][0]
    tet = nb[:4]
    amax = max(angles(tet))
    if amax > 150.0:
        return None, "square-planar / trans (max angle %.1f deg)" % amax
    return tet, "tetrahedral"


def angles(tet):
    """All six X-M-X angles, in degrees, sorted ascending."""
    out = []
    for a, b in itertools.combinations(range(4), 2):
        u, v = tet[a][3], tet[b][3]
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return sorted(out)


def group_angles(vals, tol=1.5):
    """Collapse near-equal angles into (value, multiplicity)."""
    out = []
    for v in vals:
        if out and abs(v - out[-1][0]) <= tol:
            n = out[-1][1] + 1
            out[-1] = ((out[-1][0] * out[-1][1] + v) / n, n)
        else:
            out.append((v, 1))
    return out


def build_adjacency(atoms, A, bonded, rmax=3.8):
    """i -> [(j, shift)] for every pair whose elements are in `bonded`.

    Built once per structure: the dimensionality walk would otherwise redo the
    O(natoms x 27) image scan at every step.
    """
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    F = np.array([a["f"] for a in atoms])
    els = [a["el"] for a in atoms]
    adj = {i: [] for i in range(len(atoms))}
    cart_sh = shifts @ A
    for i in range(len(atoms)):
        d = (F - F[i]) @ A
        for k in range(len(shifts)):
            r = np.linalg.norm(d + cart_sh[k], axis=1)
            for j in np.nonzero((r > 0.3) & (r <= rmax))[0]:
                if (els[i], els[j]) in bonded:
                    lim = 1.35 * (RCOV.get(els[i], 1.4) + RCOV.get(els[j], 1.4))
                    if r[j] <= lim:
                        adj[i].append((int(j), tuple(shifts[k])))
    return adj


def dimensionality(atoms, A, m_idx, want_els, adj=None, rmax=3.8):
    """Rank of the lattice translations spanned by one connected M-X component.

    Nodes are (atom index, lattice translation). Whenever the same atom is
    reached at two different translations, their difference is a lattice vector
    of the component; the rank of those vectors is the dimensionality.
    """
    if adj is None:
        mel = {atoms[k]["el"] for k in m_idx}
        bonded = {(a, b) for a in mel for b in want_els}
        bonded |= {(b, a) for a, b in bonded}
        adj = build_adjacency(atoms, A, bonded, rmax)
    # The component is infinite, so the walk is bounded to a small translation
    # box; +-2 cells is ample to expose the rank of the spanned lattice.
    BOX = 2
    seen = {(m_idx[0], (0, 0, 0))}
    first = {m_idx[0]: (0, 0, 0)}
    stack = [(m_idx[0], (0, 0, 0))]
    vecs = []
    while stack:
        i, ti = stack.pop()
        for j, sh in adj[i]:
            tj = (ti[0] + sh[0], ti[1] + sh[1], ti[2] + sh[2])
            if max(abs(v) for v in tj) > BOX:
                continue
            if (j, tj) in seen:
                continue
            seen.add((j, tj))
            if j in first:
                d = tuple(a - b for a, b in zip(tj, first[j]))
                if any(d):
                    vecs.append(d)
            else:
                first[j] = tj
            stack.append((j, tj))
    if not vecs:
        return 0, len(first)
    return int(np.linalg.matrix_rank(np.array(vecs), tol=1e-8)), len(first)


def regime(gr):
    """Tag the angle pattern. gr = [(angle, multiplicity), ...] ascending."""
    vals = [v for v, n in gr]
    if all(abs(v - 109.47) <= 3.0 for v in vals):
        return "REGULAR"
    lo = [(v, n) for v, n in gr if v < 109.47]
    hi = [(v, n) for v, n in gr if v >= 109.47]
    nlo = sum(n for _, n in lo)
    nhi = sum(n for _, n in hi)
    if nhi == 2 and nlo == 4:
        op = max(v for v, _ in hi)
        cl = [v for v, _ in lo]
        if 118 <= op <= 132 and all(98 <= v <= 106 for v in cl):
            return "SPLIT-A"
    if nlo == 2 and nhi == 4:
        cl = min(v for v, _ in lo)
        op = [v for v, _ in hi]
        if 88 <= cl <= 102 and all(112 <= v <= 122 for v in op):
            return "SPLIT-B"
    return "OTHER"
