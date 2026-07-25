"""
Flat-band bipolaron ED: Holstein (site/density) vs bond-Peierls (SSH) coupling.

Same lattice as ed.py: H(k)=D[cos k sx + sin k sz], two exactly flat bands at +-D,
<rho>=(1/2,1/2), tr g = 1/4 uniform (kappa=1), M_min = 0.25.

Holstein : H_ep = g sum_s  n_s X_s              (2N modes, one per site)
SSH      : H_ep = g sum_i  X_i T_i / D          (N modes, one per link)
           T_i = link-i hopping operator; /D makes g a fractional modulation.

The two conventions induce different effective attractions, so comparing at
equal g is meaningless. Instead we compare at EQUAL PAIR BINDING ENERGY:
  E_b = E(2 electrons) - 2 E(1 electron)
find U_inst reproducing the same E_b in the instantaneous (phononless)
attractive-Hubbard model, and report

  R = stiffness(phonon model) / stiffness(instantaneous model at same E_b)

R ~ 1        -> no mass penalty; the pair carries its cloud for free
R ~ exp(-l)  -> Lang-Firsov suppression survives
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from itertools import product
from scipy.optimize import brentq

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


W_DISP = [0.0]   # intra-orbital hopping; 0 = exactly flat, >0 = dispersive


def HR1(D):
    return (D / 2) * (SX - 1j * SZ) + (W_DISP[0] / 2) * np.eye(2)


def hop_matrix(N, D, theta):
    M = np.zeros((2 * N, 2 * N), complex)
    h = HR1(D)
    for i in range(N):
        j = (i + 1) % N
        ph = np.exp(1j * theta / N)
        for a in range(2):
            for b in range(2):
                M[2 * i + a, 2 * j + b] += h[a, b] * ph
                M[2 * j + b, 2 * i + a] += np.conj(h[a, b] * ph)
    return M


def link_op(N, D, i, theta):
    """Hopping operator across link i<->i+1 only (Hermitian), normalised by D."""
    M = np.zeros((2 * N, 2 * N), complex)
    h = HR1(D)
    j = (i + 1) % N
    ph = np.exp(1j * theta / N)
    for a in range(2):
        for b in range(2):
            M[2 * i + a, 2 * j + b] += h[a, b] * ph
            M[2 * j + b, 2 * i + a] += np.conj(h[a, b] * ph)
    return M / D


def ph_basis(nmode, nmax):
    cfgs = [c for c in product(range(nmax + 1), repeat=nmode) if sum(c) <= nmax]
    return cfgs, {c: i for i, c in enumerate(cfgs)}


def ph_ops(cfgs, index, nmode):
    X = []
    for s in range(nmode):
        r, c, v = [], [], []
        for i, cfg in enumerate(cfgs):
            n = cfg[s]
            up = list(cfg); up[s] = n + 1
            t = tuple(up)
            if t in index:
                r.append(index[t]); c.append(i); v.append(np.sqrt(n + 1))
            if n > 0:
                dn = list(cfg); dn[s] = n - 1
                t = tuple(dn)
                r.append(index[t]); c.append(i); v.append(np.sqrt(n))
        X.append(sp.coo_matrix((v, (r, c)), shape=(len(cfgs),) * 2).tocsr())
    return X, sp.diags([float(sum(c)) for c in cfgs])


def lift(O, ns, nel):
    """Single-particle matrix -> many-body matrix (nel = 1 or 2 distinguishable)."""
    if nel == 1:
        return sp.csr_matrix(O)
    idx = {(i, j): i * ns + j for i in range(ns) for j in range(ns)}
    r, c, v = [], [], []
    for (i, j), e in idx.items():
        for ip in range(ns):
            if abs(O[ip, i]) > 1e-14:
                r.append(idx[(ip, j)]); c.append(e); v.append(O[ip, i])
        for jp in range(ns):
            if abs(O[jp, j]) > 1e-14:
                r.append(idx[(i, jp)]); c.append(e); v.append(O[jp, j])
    return sp.coo_matrix((v, (r, c)), shape=(ns * ns,) * 2).tocsr()


def E0(N, D, theta, g, w0, nmax, kind, nel):
    ns = 2 * N
    ne = ns if nel == 1 else ns * ns
    Hel = lift(hop_matrix(N, D, theta), ns, nel)
    if kind == 'inst':
        if nel == 1:
            return float(np.linalg.eigvalsh(hop_matrix(N, D, theta))[0])
        dia = np.array([-g if i == j else 0.0
                        for i in range(ns) for j in range(ns)])   # g carries U here
        H = Hel + sp.diags(dia)
        return float(spla.eigsh(H, k=1, which='SA', maxiter=8000)[0][0])
    nmode = ns if kind == 'holstein' else N
    cfgs, ix = ph_basis(nmode, nmax)
    X, ntot = ph_ops(cfgs, ix, nmode)
    Iel = sp.identity(ne, format='csr')
    Iph = sp.identity(len(cfgs), format='csr')
    H = sp.kron(Hel, Iph) + w0 * sp.kron(Iel, ntot)
    for s in range(nmode):
        if kind == 'holstein':
            if nel == 1:
                dens = np.array([1.0 if i == s else 0.0 for i in range(ns)])
            else:
                dens = np.array([float((i == s) + (j == s))
                                 for i in range(ns) for j in range(ns)])
            Oel = sp.diags(dens)
        else:
            Oel = lift(link_op(N, D, s, theta), ns, nel)
        H = H + g * sp.kron(Oel, X[s])
    return float(spla.eigsh(H.tocsr(), k=1, which='SA', maxiter=12000)[0][0])


def stiff(N, D, g, w0, nmax, kind, ths=(0.0, 0.10, 0.20)):
    es = [E0(N, D, t, g, w0, nmax, kind, 2) for t in ths]
    return 2 * np.polyfit(np.array(ths), np.array(es), 2)[0]


def binding(N, D, g, w0, nmax, kind):
    return E0(N, D, 0.0, g, w0, nmax, kind, 2) - 2 * E0(N, D, 0.0, g, w0, nmax, kind, 1)


def match_U(N, D, Eb, lo=1e-4, hi=12.0):
    f = lambda U: binding(N, D, U, 0, 0, 'inst') - Eb
    return brentq(f, lo, hi, xtol=1e-8)
