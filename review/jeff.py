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
    """Metric of the j=1/2 doublet.

    Bands 4 and 5 are a KRAMERS PAIR, degenerate by time reversal to machine
    precision (max|E5-E4| ~ 3e-15 across the BZ). A rank-1 projector built from
    one member is therefore not a property of the Hamiltonian: the eigenvector
    LAPACK returns inside the degenerate subspace is arbitrary, it re-randomises
    between adjacent k-points, and the resulting <tr g> DIVERGES as nk^2
    (2.50 at nk=16 -> 37.2 at nk=64) while also shifting under a random U(2)
    rotation inside the pair. The earlier "M_trg = 5.3" was that artefact.

    The rank-2 projector onto the whole doublet is the gauge-invariant object
    and it converges: 0.4965 / 0.5312 / 0.5447 / 0.5585 at nk = 16 / 24 / 32 / 64,
    Richardson (q=2) limit M_trg = 0.563.

    lam_pair and n_phi are unaffected -- the j=1/2 doublet is orbitally isotropic
    (exactly 1/3 on each t2g orbital), so any mixing inside the pair preserves
    them. They come out at exactly 1/3 and 3 either way.
    """
    kk = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), axis=-1)
    dk = 2 * np.pi / nk
    E, V = np.linalg.eigh(hk(K, t, lam, tp))
    bl = list(bands)

    # rank-len(bands) projector: gauge invariant under mixing inside the manifold
    Um = V[..., :, bl]
    P = np.einsum('...ia,...ja->...ij', Um, Um.conj())
    trg = np.zeros(P.shape[:-2])
    for i in range(3):
        dP = (np.roll(P, -1, i) - np.roll(P, 1, i)) / (2 * dk)
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))

    r = len(bl)
    rho = np.real(np.einsum('...ii->...i', P)).reshape(-1, 6)
    orb = rho.reshape(-1, 3, 2).sum(-1)            # trace over spin -> xy,yz,zx
    Am = (orb.T @ orb) / orb.shape[0]
    w = orb.mean(0)
    b0, b1 = bl[0], bl[-1]
    gl = float(E[..., bl].min() - E[..., b0 - 1].max()) if b0 > 0 else np.inf
    gu = float(E[..., b1 + 1].min() - E[..., bl].max()) if b1 < 5 else np.inf

    out = [dict(bands=tuple(bl), rank=r,
                W=float(max(E[..., j].max() - E[..., j].min() for j in bl)),
                M_trg=float(trg.mean()),
                lam_pair=float(np.linalg.eigvalsh(Am)[-1]) / r ** 2,
                nphi=float(r ** 2 / np.sum(w ** 2)),
                iso=float(min(gl, gu)), w=w)]
    gap_manifold = float(E[..., 4].min() - E[..., 3].max())
    return out, gap_manifold
