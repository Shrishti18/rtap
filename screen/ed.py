"""
Flat-band bipolaron ED.

Lattice: N cells x 2 orbitals, H(k) = D[cos k sx + sin k sz]
         -> two EXACTLY flat bands at +-D, gap 2D, <rho>=(1/2,1/2), tr g = 1/4 uniform.
Two electrons (one up, one down), Holstein phonons (one mode per site).
Flux theta threaded through the ring: H(R) -> H(R) exp(i theta R / N).

Stiffness = d^2 E / d theta^2.  Pair moves ONLY via the geometric channel
(the band is flat), so this directly tests whether phonon dressing suppresses
metric-derived pair motion.

Reference: instantaneous limit w0 -> inf at fixed U = 2 g^2 / w0, which is
plain attractive Hubbard -U on the same lattice with no phonons.
If Lang-Firsov applied we would see stiffness ratio ~ exp(-g^2/w0^2).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from itertools import product

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def hop_matrix(N, D, theta):
    """Single-particle hopping on 2N sites, site = 2*cell + orbital."""
    HR = {1: (D / 2) * (SX - 1j * SZ), -1: (D / 2) * (SX + 1j * SZ)}
    M = np.zeros((2 * N, 2 * N), complex)
    for i in range(N):
        for R, h in HR.items():
            j = (i + R) % N
            ph = np.exp(1j * theta * R / N)
            for a in range(2):
                for b in range(2):
                    M[2 * i + a, 2 * j + b] += h[a, b] * ph
    return M


def phonon_basis(nsite, nmax):
    cfgs = [c for c in product(range(nmax + 1), repeat=nsite) if sum(c) <= nmax]
    return cfgs, {c: i for i, c in enumerate(cfgs)}


def phonon_ops(cfgs, index, nsite):
    """X_s = b_s + b_s^dag, and the total-number diagonal."""
    X = []
    for s in range(nsite):
        r, c, v = [], [], []
        for i, cfg in enumerate(cfgs):
            n = cfg[s]
            up = list(cfg); up[s] = n + 1; up = tuple(up)
            if up in index:
                r.append(index[up]); c.append(i); v.append(np.sqrt(n + 1))
            if n > 0:
                dn = list(cfg); dn[s] = n - 1; dn = tuple(dn)
                r.append(index[dn]); c.append(i); v.append(np.sqrt(n))
        X.append(sp.coo_matrix((v, (r, c)), shape=(len(cfgs),) * 2).tocsr())
    ntot = sp.diags([float(sum(c)) for c in cfgs])
    return X, ntot


def elec_space(nsite):
    """Two distinguishable electrons (up, down): index = i*nsite + j."""
    return [(i, j) for i in range(nsite) for j in range(nsite)]


def ground_energy(N, D, theta, g, w0, nmax, U_inst=None):
    ns = 2 * N
    ES = elec_space(ns)
    ne = len(ES)
    hop = hop_matrix(N, D, theta)
    r, c, v = [], [], []
    for e, (i, j) in enumerate(ES):
        for ip in range(ns):
            if abs(hop[ip, i]) > 1e-14:
                r.append(ES.index((ip, j))); c.append(e); v.append(hop[ip, i])
        for jp in range(ns):
            if abs(hop[jp, j]) > 1e-14:
                r.append(ES.index((i, jp))); c.append(e); v.append(hop[jp, j])
    Hel = sp.coo_matrix((v, (r, c)), shape=(ne, ne)).tocsr()

    if U_inst is not None:                      # instantaneous reference, no phonons
        dia = np.array([-U_inst if i == j else 0.0 for i, j in ES])
        H = Hel + sp.diags(dia)
        return float(spla.eigsh(H, k=1, which='SA', maxiter=5000)[0][0])

    cfgs, idx = phonon_basis(ns, nmax)
    X, ntot = phonon_ops(cfgs, idx, ns)
    npn = len(cfgs)
    Iel = sp.identity(ne, format='csr')
    Iph = sp.identity(npn, format='csr')
    H = sp.kron(Hel, Iph) + w0 * sp.kron(Iel, ntot)
    for s in range(ns):
        dens = np.array([float((i == s) + (j == s)) for i, j in ES])
        H = H + g * sp.kron(sp.diags(dens), X[s])
    H = H.tocsr()
    return float(spla.eigsh(H, k=1, which='SA', maxiter=10000)[0][0])


def stiffness(N, D, g, w0, nmax, U_inst=None, ths=(0.0, 0.10, 0.20)):
    es = [ground_energy(N, D, t, g, w0, nmax, U_inst) for t in ths]
    return 2 * np.polyfit(np.array(ths), np.array(es), 2)[0], es[0]
