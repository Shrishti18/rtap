"""
Cs2CoS2 chain, built from the MEASURED crystallography.

Geometry check first. Cs2CoS2 reports X-Co-X = 96.12 x2, 114.92 x2, 118.17 x2.
For an SiS2-type edge-sharing tetrahedral chain with Co at (0,0,0),(0,0,h) and
shared-edge S at (+-u, 0, h/2), the two shared-edge ligands subtend 2*theta
with tan(theta)=u/(h/2), and every cross angle is arccos(-cos^2 theta).
    2*theta = 96.12  ->  theta = 48.06 deg
    arccos(-cos^2 48.06) = 116.5 deg  vs measured mean of 114.92 and 118.17
                                       = 116.5 deg.   EXACT.
So the measured tetrahedron IS the ideal edge-sharing chain with theta=48.06,
split only by a small orthorhombic distortion. The model below uses that.

Electronic model:
  * 2 Co per cell (edge orientation alternates x, y -> period 2h)
  * 5 d orbitals each -> 10x10 H(k)
  * on-site: AOM crystal field from the 4 ligands (e_sig, e_pi)
  * inter-site: superexchange-like d-d hopping through each BRIDGING S,
        H_eff = -J * w(n_A->S) w(n_B->S)^T ,  J = t_pd^2 / Delta_inv
    with w(n) the sigma-bonding d combination for ligand direction n, and
    t_pd from Harrison, Delta_inv = eps_p(S) - eps_d(Co) = +4.54 eV.
"""
import numpy as np

NAMES = ['z2', 'x2-y2', 'xy', 'xz', 'yz']
HBM = 7.6199
RD_CO = 1.045
EPS_INV = 4.54            # S 3p above Co 3d, eV
S2, S6 = np.sqrt(2), np.sqrt(6)
BASIS = [np.diag([-1, -1, 2]) / S6,
         np.diag([1, -1, 0]) / S2,
         np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / S2,
         np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / S2,
         np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / S2]


def dmat(n):
    """Rotation of the real-d basis taking z to direction n."""
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    th = np.arccos(np.clip(n[2], -1, 1))
    ph = np.arctan2(n[1], n[0])
    Ry = np.array([[np.cos(th), 0, np.sin(th)], [0, 1, 0], [-np.sin(th), 0, np.cos(th)]])
    Rz = np.array([[np.cos(ph), -np.sin(ph), 0], [np.sin(ph), np.cos(ph), 0], [0, 0, 1]])
    R = Rz @ Ry
    D = np.zeros((5, 5))
    for j, Q in enumerate(BASIS):
        Qr = R @ Q @ R.T
        for i, P in enumerate(BASIS):
            D[i, j] = np.sum(P * Qr)
    return D


def wsig(n):
    """sigma-bonding d combination for a ligand along n (first column)."""
    return dmat(n)[:, 0]


def wpi(n):
    """the two pi combinations."""
    D = dmat(n)
    return D[:, 3], D[:, 4]


def geometry(dCoS=2.365, theta_deg=48.06):
    th = np.radians(theta_deg)
    u = dCoS * np.sin(th)
    h = 2 * dCoS * np.cos(th)          # Co-Co spacing along the chain
    A = np.array([0.0, 0.0, 0.0])
    B = np.array([0.0, 0.0, h])
    SAB = [np.array([u, 0.0, h / 2]), np.array([-u, 0.0, h / 2])]      # shared A-B
    SBA = [np.array([0.0, u, 3 * h / 2]), np.array([0.0, -u, 3 * h / 2])]  # shared B-A(+c)
    return A, B, SAB, SBA, h, u


def hk_chain(kz, dCoS=2.365, theta_deg=48.06, epi_ratio=0.25, tscale=1.0):
    A, B, SAB, SBA, h, u = geometry(dCoS, theta_deg)
    c = 2 * h
    t_pd = tscale * 2.95 * HBM * RD_CO ** 1.5 / dCoS ** 3.5
    J = t_pd ** 2 / EPS_INV
    esig, epi = 1.0, epi_ratio
    # on-site crystal field (same for both Co by symmetry, up to axis swap)
    def cf(site, ligs):
        H = np.zeros((5, 5))
        for L in ligs:
            D = dmat(L - site)
            H += D @ np.diag([esig, 0, 0, epi, epi]) @ D.T
        return H
    ligA = SAB + [np.array([0.0, u, -h / 2]), np.array([0.0, -u, -h / 2])]
    ligB = SBA + SAB
    # scale CF so e_sig is a physical ~1.0 eV for a Co-S sigma bond
    CFA, CFB = cf(A, ligA), cf(B, ligB)
    kz = np.atleast_1d(kz)
    H = np.zeros(kz.shape + (10, 10), complex)
    H[..., :5, :5] = CFA
    H[..., 5:, 5:] = CFB
    # A -> B through the two shared S (intra-cell)
    tAB = np.zeros((5, 5))
    for Sp in SAB:
        tAB += -J * np.outer(wsig(Sp - A), wsig(Sp - B))
    # B -> A(+c) through the other two shared S (inter-cell)
    tBA = np.zeros((5, 5))
    for Sp in SBA:
        tBA += -J * np.outer(wsig(Sp - B), wsig(Sp - (A + np.array([0, 0, c]))))
    ph = np.exp(1j * kz * c)
    H[..., :5, 5:] += tAB
    H[..., 5:, :5] += tAB.T
    H[..., 5:, :5] += tBA.T * ph[..., None, None]
    H[..., :5, 5:] += tBA * np.conj(ph)[..., None, None]
    return H, c, t_pd, J
