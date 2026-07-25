"""
Structure-file pre-filter for the derived spec.

PASS requires, from symmetry alone (no DFT):
  * a site whose site-symmetry group carries a 2-fold orbital multiplet
    (equivalently: a 3-, 4-, or 6-fold axis; all triclinic/monoclinic/
    orthorhombic sites are excluded outright)
  * exactly ONE such site per PRIMITIVE cell, so n_phi = 2
  * cubic sites qualify only via e_g {d_z2, d_x2-y2}; t2g and p are 3-fold
    multiplets -> n_phi = 3, diluted

Orbital multiplets are obtained by decomposing the real p (l=1) and d (l=2)
representations over the actual Cartesian site-symmetry group. Real matrices
on real orbitals => the decomposition returned is already TRS-glued.
"""
import re
import numpy as np

DNAMES = ['xy', 'yz', 'xz', 'x2-y2', 'z2']
PNAMES = ['x', 'y', 'z']


# ---------------- CIF parsing ----------------

def _num(s):
    return float(re.sub(r'\(.*?\)', '', s).strip())


def _parse_xyz(s):
    """'-y+1/2, x, z' -> (W 3x3, w 3)"""
    W = np.zeros((3, 3))
    w = np.zeros(3)
    for i, comp in enumerate(s.replace(' ', '').lower().split(',')):
        for m in re.finditer(r'([+-]?)(\d+/\d+|\d*\.?\d*)\*?([xyz]?)', comp):
            sgn, num, var = m.group(1), m.group(2), m.group(3)
            if not num and not var:
                continue
            if num == '':
                v = 1.0
            elif '/' in num:
                p, q = num.split('/')
                v = float(p) / float(q)
            else:
                v = float(num)
            if sgn == '-':
                v = -v
            if var:
                W[i, 'xyz'.index(var)] += v
            else:
                w[i] += v
    return W, w


SYMLOOP_RE = re.compile(
    r'(?:_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz)\s*\n'
    r'(.*?)(?=\nloop_|\n\s*_[a-zA-Z]|\Z)', re.S)
_OPCHARS = re.compile(r'^[-+0-9xyzXYZ/,.\s]+$')


def symmetry_ops(txt):
    """Symmetry operation strings, taken ONLY from the symmetry loop.

    The original `re.findall(r"'([^']*[xyz][^']*)'", txt)` picked up any quoted
    string containing an x, y or z -- publication titles and chemical names
    among them -- which both corrupts the stabiliser and raises on strings like
    'Nb2 Se ... data'. Restricting to the loop body and requiring the token to
    look like an operation fixes it.
    """
    out = []
    for m in SYMLOOP_RE.finditer(txt):
        for line in m.group(1).splitlines():
            line = line.strip()
            if not line:
                continue
            q = re.findall(r"'([^']*)'|\"([^\"]*)\"", line)
            cand = [a or b for a, b in q]
            if not cand:                       # unquoted: strip a leading index
                cand = [re.sub(r'^\d+\s+', '', line)]
            for s in cand:
                if s.count(',') == 2 and _OPCHARS.match(s):
                    out.append(s)
    return out


def parse_cif(path):
    txt = open(path, errors='replace').read()
    cell = {}
    for key in ['a', 'b', 'c']:
        m = re.search(r'_cell_length_%s\s+(\S+)' % key, txt)
        cell[key] = _num(m.group(1))
    for key, nm in [('al', 'alpha'), ('be', 'beta'), ('ga', 'gamma')]:
        m = re.search(r'_cell_angle_%s\s+(\S+)' % nm, txt)
        cell[key] = _num(m.group(1))
    ops = [_parse_xyz(s) for s in symmetry_ops(txt)]
    if not ops:
        raise ValueError('no symmetry operations found')
    sites = []
    lines = txt.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == 'loop_':
            hdr, j = [], i + 1
            while j < len(lines) and lines[j].strip().startswith('_'):
                hdr.append(lines[j].strip())
                j += 1
            if any('_atom_site_fract_x' in h for h in hdr):
                cx = [k for k, h in enumerate(hdr) if h == '_atom_site_fract_x'][0]
                cy = [k for k, h in enumerate(hdr) if h == '_atom_site_fract_y'][0]
                cz = [k for k, h in enumerate(hdr) if h == '_atom_site_fract_z'][0]
                cl = [k for k, h in enumerate(hdr) if 'label' in h or 'type_symbol' in h]
                cl = cl[0] if cl else 0
                while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith(('_', 'loop_', '#')):
                    t = lines[j].split()
                    if len(t) > max(cx, cy, cz):
                        sites.append((t[cl], np.array([_num(t[cx]), _num(t[cy]), _num(t[cz])])))
                    j += 1
            i = j
        else:
            i += 1
    return cell, ops, sites


def cart_matrix(c):
    a, b, cc = c['a'], c['b'], c['c']
    al, be, ga = np.radians([c['al'], c['be'], c['ga']])
    v1 = np.array([a, 0, 0])
    v2 = np.array([b * np.cos(ga), b * np.sin(ga), 0])
    cx = cc * np.cos(be)
    cy = cc * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    v3 = np.array([cx, cy, np.sqrt(max(cc**2 - cx**2 - cy**2, 0.0))])
    return np.array([v1, v2, v3]).T          # columns = lattice vectors


# ---------------- symmetry analysis ----------------

def stabilizer(x0, ops, tol=1e-4):
    """Ops fixing x0 modulo lattice translation -> rotational parts."""
    out = []
    for W, w in ops:
        d = W @ x0 + w - x0
        if np.allclose(d - np.round(d), 0, atol=tol):
            if not any(np.allclose(W, X, atol=tol) for X in out):
                out.append(W)
    return out


def orbit(x0, ops, tol=1e-4):
    pts = []
    for W, w in ops:
        y = (W @ x0 + w) % 1.0
        if not any(np.allclose(np.abs(y - z) % 1.0, 0, atol=tol) or
                   np.allclose(np.abs(y - z) % 1.0, 1, atol=tol) for z in pts):
            pts.append(y)
    return len(pts)


def n_centering(ops, tol=1e-6):
    return sum(1 for W, w in ops if np.allclose(W, np.eye(3), atol=tol))


def to_cart(Ws, A):
    Ai = np.linalg.inv(A)
    return [A @ W @ Ai for W in Ws]


# ---------------- orbital multiplets ----------------

def _dbasis():
    s2, s6 = np.sqrt(2), np.sqrt(6)
    return [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / s2,
            np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / s2,
            np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / s2,
            np.diag([1.0, -1.0, 0.0]) / s2,
            np.diag([-1.0, -1.0, 2.0]) / s6]


def _drep(R):
    B = _dbasis()
    M = np.zeros((5, 5))
    for j, Q in enumerate(B):
        Qp = R @ Q @ R.T
        for i, P in enumerate(B):
            M[i, j] = np.sum(P * Qp)
    return M


def multiplets(G, kind='d'):
    dim, rep = (5, _drep) if kind == 'd' else (3, lambda R: R)
    names = DNAMES if kind == 'd' else PNAMES
    rng = np.random.default_rng(0)
    A = rng.normal(size=(dim, dim))
    A = A + A.T
    M = sum(rep(g).T @ A @ rep(g) for g in G) / len(G)
    wv, V = np.linalg.eigh(M)
    out, cur = [], [0]
    for i in range(1, dim):
        if abs(wv[i] - wv[i - 1]) < 1e-7:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    res = []
    for gr in out:
        wgt = (V[:, gr] ** 2).sum(1)
        res.append((len(gr), '+'.join(names[i] for i in np.argsort(-wgt)[:len(gr)])))
    return res


def max_rot_order(G, tol=1e-6):
    best = 1
    for R in G:
        if np.linalg.det(R) < 0:
            continue
        c = (np.trace(R) - 1) / 2
        for n in (6, 4, 3, 2):
            if abs(c - np.cos(2 * np.pi / n)) < 1e-5:
                best = max(best, n)
                break
    return best


# ---------------- driver ----------------

TM3D = set('Sc Ti V Cr Mn Fe Co Ni Cu Zn'.split())
TM4D = set('Y Zr Nb Mo Tc Ru Rh Pd Ag Cd'.split())
TM5D = set('Lu Hf Ta W Re Os Ir Pt Au Hg'.split())


def _elem(lab):
    m = re.match(r'([A-Z][a-z]?)', lab)
    return m.group(1) if m else lab


def _shell(lab):
    e = _elem(lab)
    return '3d' if e in TM3D else '4d' if e in TM4D else '5d' if e in TM5D else '-'


def scan(path, verbose=True):
    cell, ops, sites = parse_cif(path)
    A = cart_matrix(cell)
    ncent = max(n_centering(ops), 1)
    rows = []
    for lab, x0 in sites:
        G = to_cart(stabilizer(x0, ops), A)
        norb = orbit(x0, ops)
        nprim = norb / ncent
        dm = multiplets(G, 'd')
        pm = multiplets(G, 'p')
        d2 = [l for n, l in dm if n == 2]
        p2 = [l for n, l in pm if n == 2]
        nphi = 2 * nprim if (d2 or p2) else None
        rows.append(dict(label=lab, order=len(G), maxrot=max_rot_order(G),
                         mult=norb, prim=nprim, d2=d2, p2=p2, nphi=nphi,
                         shell=_shell(lab),
                         d=' '.join(f"{n}({l})" for n, l in dm)))
    if verbose:
        print(f"  cell {cell['a']:.3f} {cell['b']:.3f} {cell['c']:.3f}  "
              f"centering x{ncent}  ops {len(ops)}")
        print(f"  {'site':>6} {'sh':>3} {'|G|':>4} {'rot':>4} {'prim':>5} "
              f"{'d-doublet':<14} {'p-doublet':<10} {'n_phi':>6} {'':>9}")
        for r in rows:
            ok = 'FAIL'
            if r['nphi'] == 2:
                ok = 'PASS' + ('*' if r['shell'] in ('4d', '5d') and r['d2'] else '')
            elif r['nphi']:
                ok = 'dilute'
            print(f"  {r['label']:>6} {r['shell']:>3} {r['order']:>4} {r['maxrot']:>4} "
                  f"{r['prim']:>5.0f} {(','.join(r['d2']) or '-'):<14} "
                  f"{(','.join(r['p2']) or '-'):<10} {str(r['nphi'] or '-'):>6} {ok:>9}")
    return rows
