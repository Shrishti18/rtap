"""
Resolve the ED-vs-BdG discrepancy on a model that is NOT pathological.

Sawtooth chain: A on a chain (hopping t), B bridging A_i and A_{i+1} (hopping t').
    H_AA(k) = 2t cos k ,  H_AB(k) = t'(1 + e^{-ik})
t' = sqrt(2) t gives an EXACTLY FLAT band at E = -2t.

Why this model settles it:
  * hoppings are REAL  -> TRS, so singlet BdG at (k,-k) is the right structure
  * H(R=1) = [[t,0],[0,0]] is NOT nilpotent -> a pair CAN circulate the ring,
    unlike the two earlier test lattices where H(R=1)^2 = 0 forced zero
    transport and made the ED/BdG comparison meaningless
  * the flat band is isolated: gap runs 2t (at k=pi) to 6t (at k=0)

Normalisation. BdG returns Draw = d^2 e/dq^2 with e = energy PER CELL and pair
momentum 2q. Flux theta on an N-cell ring shifts each single-particle k by
theta/N, so the pair momentum shifts by 2 theta/N, giving q = theta/N and
    S_ED = d^2 E_total/dtheta^2 = N * Draw / N^2 = Draw / N.
(An earlier note of mine said Draw/(4N); that used q = theta/(2N), which is
wrong. Both are printed below so the data decides.)
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize_scalar

T = 1.0
TP = np.sqrt(2.0) * T


def hk(k):
    k = np.atleast_1d(k)
    H = np.zeros(k.shape + (2, 2), complex)
    H[..., 0, 0] = 2 * T * np.cos(k)
    H[..., 0, 1] = TP * (1 + np.exp(-1j * k))
    H[..., 1, 0] = np.conj(H[..., 0, 1])
    return H


def hop_real(N, theta):
    """2N x 2N single-particle matrix. site = 2*i + (0=A, 1=B)."""
    M = np.zeros((2 * N, 2 * N), complex)
    ph = np.exp(1j * theta / N)
    for i in range(N):
        j = (i + 1) % N
        A, B, Aj = 2 * i, 2 * i + 1, 2 * j
        M[A, Aj] += T * ph;        M[Aj, A] += T * np.conj(ph)   # A_i - A_{i+1}
        M[A, B] += TP;             M[B, A] += TP                 # A_i - B_i
        M[Aj, B] += TP * np.conj(ph)                             # A_{i+1} - B_i
        M[B, Aj] += TP * ph
    return M


def metric(nk=4000):
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    E, V = np.linalg.eigh(hk(k))
    b = int(np.argmin([E[:, j].max() - E[:, j].min() for j in range(2)]))
    u = V[:, :, b]
    P = np.einsum('ki,kj->kij', u, u.conj())
    dP = (np.roll(P, -1, 0) - np.roll(P, 1, 0)) / (2 * (2 * np.pi / nk))
    trg = 0.5 * np.real(np.einsum('kij,kji->k', dP, dP))
    rho = np.abs(u) ** 2
    A = (rho.T @ rho) / nk
    gap = float(np.min(np.abs(E[:, 1 - b] - E[:, b])))
    return dict(band=b, Eflat=float(E[:, b].mean()),
                W=float(E[:, b].max() - E[:, b].min()),
                M_trg=float(trg.mean()), lam=float(np.linalg.eigvalsh(A)[-1]),
                rho=rho.mean(0), gap=gap)


def ed_stiff(N, U, mu, ths=(0.0, 0.12, 0.24)):
    ns = 2 * N
    idx = {(i, j): i * ns + j for i in range(ns) for j in range(ns)}

    def E0(theta):
        h = hop_real(N, theta) - mu * np.eye(ns)
        r, c, v = [], [], []
        for (i, j), e in idx.items():
            for ip in range(ns):
                if abs(h[ip, i]) > 1e-14:
                    r.append(idx[(ip, j)]); c.append(e); v.append(h[ip, i])
            for jp in range(ns):
                if abs(h[jp, j]) > 1e-14:
                    r.append(idx[(i, jp)]); c.append(e); v.append(h[jp, j])
        H = sp.coo_matrix((v, (r, c)), shape=(ns * ns,) * 2).tocsr()
        H = H + sp.diags([-U if i == j else 0.0 for i in range(ns) for j in range(ns)])
        return float(spla.eigsh(H, k=1, which='SA', maxiter=20000)[0][0])

    es = np.array([E0(t) for t in ths])
    return 2 * np.polyfit(np.array(ths), es, 2)[0]


def bdg_draw(U, mu, nk=2000, qs=(0.0, 0.015, 0.030)):
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk

    def E0(dd, q):
        M = np.zeros((nk, 4, 4), complex)
        M[:, :2, :2] = hk(k + q) - mu * np.eye(2)
        M[:, 2:, 2:] = -(hk(k - q) - mu * np.eye(2))
        M[:, :2, 2:] = dd * np.eye(2)
        M[:, 2:, :2] = dd * np.eye(2)
        return -0.5 * np.abs(np.linalg.eigvalsh(M)).sum(-1).mean() + 2 * dd ** 2 / U

    es = [minimize_scalar(E0, bracket=(1e-4, 0.4), args=(q,)).fun for q in qs]
    return 2 * np.polyfit(np.array(qs), np.array(es), 2)[0]
