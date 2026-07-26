"""
CALIBRATION: does the framework reproduce the ONE measured flat-band
superconductor?

Magic-angle twisted bilayer graphene is the only system where flat-band,
quantum-geometric superconductivity is believed to operate AND Tc is measured:
**Tc ~ 1.7 K**. Everything this project has predicted (348 K, 330 K, 8.1 K, ...)
comes from a formula that has never once reproduced a measured transition
temperature. This is the test.

    Tc = min( U*lam/4 , 0.369*U*M_min )

Built from scratch (Bistritzer-MacDonald) rather than taking a literature metric,
because a published "quantum metric" may be in the M = (2pi)^(d-1) * M_trg
normalisation, a 39.5x error in 3D and 6.28x in 2D, and that exact confusion has
already cost this project real work.

MODEL. One valley, plane-wave expansion on the moire reciprocal lattice:
    h_l(k) = -hbar v_F (k - K_l) . sigma_theta_l        (Dirac, rotated)
    T_1 = [[wAA, wAB],[wAB, wAA]]
    T_2 = [[wAA, wAB/w],[wAB*w, wAA]],  T_3 = [[wAA, wAB*w],[wAB/w, wAA]]
    w = exp(2 pi i / 3)
Units: moire lattice constant = 1, so <tr g> is directly comparable with every
other number in this project (the 0.84 2D threshold).
"""
import numpy as np
from scipy.optimize import minimize

KB = 8.617333e-5           # eV/K
S0 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
W3 = np.exp(2j * np.pi / 3)

# --- physical parameters -----------------------------------------------------
A_GR = 2.46                 # graphene lattice constant, Angstrom
HVF = 5.944                 # hbar v_F in eV*Angstrom  (v_F = 0.9e6 m/s)
WAA = 0.0797                # interlayer AA tunnelling, eV (relaxed value)
WAB = 0.0975                # interlayer AB tunnelling, eV
THETA_DEG = 1.05            # magic angle


def build(theta_deg=THETA_DEG, wAA=WAA, wAB=WAB, N=4):
    """Return (kfun, dim, kD, ktheta). N = shell cutoff in moire G vectors."""
    th = np.deg2rad(theta_deg)
    kD = 4 * np.pi / (3 * A_GR)
    kt = 2 * kD * np.sin(th / 2)                    # moire wavevector
    q1 = kt * np.array([0.0, -1.0])
    q2 = kt * np.array([np.sqrt(3) / 2, 0.5])
    q3 = kt * np.array([-np.sqrt(3) / 2, 0.5])
    G1, G2 = q2 - q1, q3 - q1
    T = [np.array([[wAA, wAB], [wAB, wAA]], complex),
         np.array([[wAA, wAB / W3], [wAB * W3, wAA]], complex),
         np.array([[wAA, wAB * W3], [wAB / W3, wAA]], complex)]
    qs = [q1, q2, q3]

    Gs = [m * G1 + n * G2 for m in range(-N, N + 1) for n in range(-N, N + 1)
          if abs(m + n) <= N]
    Gs = np.array(Gs)
    nG = len(Gs)
    idx = {(round(g[0], 6), round(g[1], 6)): i for i, g in enumerate(Gs)}

    def dirac(kx, ky, sgn):
        """-hvF (k) . sigma rotated by +-theta/2."""
        c, s = np.cos(sgn * th / 2), np.sin(sgn * th / 2)
        kxr, kyr = c * kx + s * ky, -s * kx + c * ky
        return -HVF * (kxr * SX + kyr * SY)

    def H(kx, ky):
        dim = 4 * nG                                # layer x sublattice x G
        M = np.zeros((dim, dim), complex)
        for i, g in enumerate(Gs):
            # layer 1 block (rows 0:2 of each G), layer 2 block (rows 2:4)
            M[4*i:4*i+2, 4*i:4*i+2] = dirac(kx + g[0], ky + g[1], +1)
            M[4*i+2:4*i+4, 4*i+2:4*i+4] = dirac(kx + g[0] - q1[0],
                                                ky + g[1] - q1[1], -1)
            for j in range(3):
                dg = q1 - qs[j]
                key = (round(g[0] + dg[0], 6), round(g[1] + dg[1], 6))
                if key in idx:
                    k2 = idx[key]
                    M[4*i:4*i+2, 4*k2+2:4*k2+4] += T[j]
                    M[4*k2+2:4*k2+4, 4*i:4*i+2] += T[j].conj().T
        return M

    return H, 4 * nG, kt, (G1, G2)


def bands_and_metric(theta_deg=THETA_DEG, nk=24, N=4, **kw):
    H, dim, kt, (G1, G2) = build(theta_deg, N=N, **kw)
    # moire BZ sampled on the G1,G2 reciprocal basis
    f = (np.arange(nk) + 0.5) / nk
    F1, F2 = np.meshgrid(f, f, indexing='ij')
    KX = F1 * G1[0] + F2 * G2[0]
    KY = F1 * G1[1] + F2 * G2[1]
    nflat = dim // 2                                # two flat bands: nflat-1, nflat
    E = np.zeros((nk, nk, dim))
    V = np.zeros((nk, nk, dim, dim), complex)
    for i in range(nk):
        for j in range(nk):
            e, v = np.linalg.eigh(H(KX[i, j], KY[i, j]))
            E[i, j], V[i, j] = e, v
    lo, hi = nflat - 1, nflat
    Wf = float(max(np.ptp(E[..., lo]), np.ptp(E[..., hi])))
    gap_up = float(E[..., hi + 1].min() - E[..., hi].max())
    gap_dn = float(E[..., lo].min() - E[..., lo - 1].max())

    # rank-2 projector onto the two flat bands (they touch, so rank-1 is
    # meaningless -- the D1 lesson)
    Um = V[..., :, [lo, hi]]
    P = np.einsum('...ia,...ja->...ij', Um, Um.conj())
    # d/dk along the two moire basis directions, in units where a_moire = 1
    # (fractional coords f1,f2 run 0..1, so d/df = nk * finite difference)
    trg = np.zeros((nk, nk))
    for ax in (0, 1):
        dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * (2 * np.pi / nk))
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    return dict(W=Wf, gap_up=gap_up, gap_dn=gap_dn, M_trg=float(trg.mean()),
                E=E, lo=lo, kt=kt, dim=dim)


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)
    print("1. BAND STRUCTURE at the magic angle (convergence in cutoff N)")
    print(f"   {'N':>3}{'dim':>6}{'W_flat(meV)':>13}{'gap up(meV)':>13}"
          f"{'gap dn(meV)':>13}{'M_trg':>9}")
    for N in [2, 3, 4]:
        r = bands_and_metric(nk=18, N=N)
        print(f"   {N:>3}{r['dim']:>6}{r['W']*1000:>13.2f}{r['gap_up']*1000:>13.2f}"
              f"{r['gap_dn']*1000:>13.2f}{r['M_trg']:>9.4f}")

    print("\n2. GRID CONVERGENCE of the metric (N = 4)")
    for nk in [12, 18, 24, 30]:
        r = bands_and_metric(nk=nk, N=4)
        print(f"   nk = {nk:>3}: W = {r['W']*1000:>6.2f} meV   M_trg = {r['M_trg']:.4f}")

    r = bands_and_metric(nk=30, N=4)
    M = r['M_trg']
    iso = min(r['gap_up'], r['gap_dn'])
    print(f"\n3. THE CALIBRATION.   measured Tc(MATBG) ~ 1.7 K")
    print(f"   flat-band width      W       = {r['W']*1000:.2f} meV")
    print(f"   isolation gap        Delta   = {iso*1000:.2f} meV")
    print(f"   quantum metric       M_trg   = {M:.4f}   (2D threshold 0.84)")
    print(f"   pairing eigenvalue   lam     = 0.5 (two Wannier orbitals/valley/spin)")
    print()
    lam = 0.5
    print(f"   {'U (meV)':>9}{'source':<34}{'amp(K)':>9}{'geo(K)':>9}"
          f"{'Tc pred':>9}{'vs 1.7 K':>10}")
    for U, src in [(iso / 2, "the projection cap Delta/2"),
                   (0.025, "Coulomb e^2/(eps L_M), eps=10"),
                   (0.010, "screened estimate"),
                   (0.005, "strongly screened")]:
        amp = U * lam / 4 / KB
        geo = 0.369 * U * M / KB
        tc = min(amp, geo)
        print(f"   {U*1000:>9.1f}{src:<34}{amp:>9.1f}{geo:>9.1f}{tc:>9.1f}"
              f"{tc/1.7:>9.1f}x")


# ===========================================================================
# RESULT, AND A THIRD BZ-BOUNDARY BUG CORRECTED
# ===========================================================================
#
# THE CALIBRATION (theta = 1.05 deg, N = 4):
#     W_flat     = 11.63 meV        (literature MATBG: ~5-15 meV)  OK
#     Delta_iso  = 22.45 meV        (literature: ~20-30 meV)       OK
#     lam        = 0.5              (two moire Wannier orbitals/valley/spin)
#     U          = min(Coulomb, Delta_iso/2) = 11.22 meV
#
#     Tc(amplitude cap) = U*lam/4        = 16.3 K   <- BINDS
#     Tc(geometric cap) = 0.369*U*M_trg  = 71.9 K
#     Tc predicted      = 16.3 K
#     Tc MEASURED       =  1.7 K
#     ------------------------------------------------
#     OVERESTIMATE      =  9.6x
#
# The binding cap is the AMPLITUDE cap, which uses only lam and U -- NOT the
# metric. So this calibration is independent of the metric bug below.
#
# SECOND CALIBRATION POINT, already in ceiling.py: BKBO predicted 235 K,
# measured 30 K -> 7.7x. Two unrelated systems (moire flat bands; s-orbital
# valence skipping), two overestimates, 7.7x and 9.6x.
#
# METRIC BUG (third BZ-boundary occurrence in this project). The BM Hamiltonian
# in the plane-wave basis satisfies H(k+G) = V(G) H(k) V(G)^dag where V(G)
# PERMUTES the G-vector labels -- it is NOT periodic. np.roll across the BZ
# boundary therefore differences two mismatched bases. Symptom: M_trg climbing
# 0.836 -> 1.166 -> 1.451 at nk = 12/18/24 despite a POSITIVE 22 meV gap.
# Recomputed with finite offsets that never cross the boundary, using
# gauge-invariant projectors:
#     theta = 0.95:  M_trg = 0.603 - 0.619   (converged in nk AND in offset)
#     theta = 1.05:  M_trg = 0.575 - 0.583
# So the true MATBG flat-band metric is ~0.58-0.62, not 1.5. It is BELOW the
# 2D threshold of 0.84 -- consistent with MATBG being a 1.7 K superconductor
# rather than a room-temperature one.
#
# CONSEQUENCE FOR EVERY NUMBER IN THIS PROJECT: divide predicted Tc by ~8-10.
#   reaching 300 K requires the formula to predict 2400-3000 K
#   the Kramers-Creutz proposal's 330 K  ->  ~35 K
#   the required U rises from 158 meV    ->  1.3-1.6 eV, i.e. bismuthate scale
# The correction does not kill the programme; it selects hard for the largest
# available U, which is the s2 valence-skipping mechanism.
