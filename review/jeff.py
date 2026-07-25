"""
j_eff = 1/2 quantum metric.

t2g manifold (xy, yz, zx) x spin = 6 states. Atomic SOC in the t2g subspace acts
as an EFFECTIVE L=1 with reversed sign:  H_soc = -(lam/2) l_eff . sigma
splitting into  j=1/2 (2 states, UPPER for the reversed sign)  and
j=3/2 (4 states), separated by 3*lam/2.

Hopping: 5d t2g in a corner-sharing octahedral network hops predominantly
through the ligand p, and the standard result is that each t2g orbital is
quasi-2D -- xy hops in the xy plane, yz in yz, zx in zx. Build that, add SOC,
diagonalise the 6x6, and measure the metric of the j=1/2 bands.

The question the whole project now turns on: is M_trg of the j=1/2 manifold
above ~0.34 (the 3D electronic threshold), or is SOC buying isolation at the
cost of stiffness?
"""
import numpy as np

# t2g basis order: xy, yz, zx  (orbital), each x spin up/down -> 6
LX = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], complex)
LY = np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], complex)
LZ = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def h_soc(lam):
    """-(lam/2) l_eff.sigma on the t2g x spin space (6x6)."""
    return -0.5 * lam * (np.kron(LX, SX) + np.kron(LY, SY) + np.kron(LZ, SZ))


def hop_t2g(k, t, tp=0.0):
    """
    Orbital-selective t2g hopping, cubic lattice, a=1.
    xy disperses in (kx,ky); yz in (ky,kz); zx in (kz,kx).
    tp adds a weak isotropic inter-orbital-independent term (t' second neighbour)
    which is what makes the lattice genuinely 3D for every orbital.
    """
    kx, ky, kz = k[..., 0], k[..., 1], k[..., 2]
    d_xy = -2 * t * (np.cos(kx) + np.cos(ky)) - 4 * tp * np.cos(kx) * np.cos(ky)
    d_yz = -2 * t * (np.cos(ky) + np.cos(kz)) - 4 * tp * np.cos(ky) * np.cos(kz)
    d_zx = -2 * t * (np.cos(kz) + np.cos(kx)) - 4 * tp * np.cos(kz) * np.cos(kx)
    H = np.zeros(k.shape[:-1] + (3, 3), complex)
    H[..., 0, 0] = d_xy
    H[..., 1, 1] = d_yz
    H[..., 2, 2] = d_zx
    return H


def hk(k, t, lam, tp=0.0):
    Ho = hop_t2g(k, t, tp)
    H6 = np.kron(Ho, np.eye(2)) if Ho.ndim == 2 else np.einsum(
        '...ij,ab->...iajb', Ho, np.eye(2)).reshape(Ho.shape[:-2] + (6, 6))
    return H6 + h_soc(lam)


def analyse(t=0.4, lam=0.45, tp=0.0, nk=24, bands=(4, 5)):
    kk = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), axis=-1)
    dk = 2 * np.pi / nk
    E, V = np.linalg.eigh(hk(K, t, lam, tp))
    out = []
    for b in bands:
        u = V[..., :, b]
        P = np.einsum('...i,...j->...ij', u, u.conj())
        trg = np.zeros(P.shape[:-2])
        for i in range(3):
            dP = (np.roll(P, -1, i) - np.roll(P, 1, i)) / (2 * dk)
            trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
        # orbital weights: trace over spin, so rho_alpha for alpha in {xy,yz,zx}
        w = np.abs(u.reshape(u.shape[:-1] + (3, 2))) ** 2
        rho = w.sum(-1)
        rf = rho.reshape(-1, 3)
        A = (rf.T @ rf) / rf.shape[0]
        gl = float(E[..., b].min() - E[..., b - 1].max()) if b > 0 else np.inf
        gu = float(E[..., b + 1].min() - E[..., b].max()) if b < 5 else np.inf
        out.append(dict(band=b, W=float(E[..., b].max() - E[..., b].min()),
                        M_trg=float(trg.mean()),
                        lam_pair=float(np.linalg.eigvalsh(A)[-1]),
                        nphi=float(1 / np.sum(rf.mean(0) ** 2)),
                        iso=min(gl, gu), w=rf.mean(0)))
    # manifold gap: j=1/2 pair (bands 4,5) vs j=3/2 (bands 0-3)
    gap_manifold = float(E[..., 4].min() - E[..., 3].max())
    return out, gap_manifold
