"""
Valence-skipping s orbital on a two-site 3D DIAMOND net.

THE IDEA, and why it defeats my own exclusion. I ruled out every s^2 valence
skipper (Bi, Pb, Tl, Sn, Sb, In, Ge) because a single s orbital gives n_phi = 1,
a k-independent projector, and tr g == 0 identically. That theorem is about
ORBITAL degeneracy on one atom. Put the same s orbital on TWO SYMMETRY-RELATED
SITES and n_phi = 2 comes from the SUBLATTICE index instead. The metric is then
not forced to zero, and the ~1.9 eV bismuthate-scale negative U is retained.

WHY THE PRIZE IS LARGE. The 0.34 metric threshold was derived assuming
U ~ 0.2 eV. The requirements scale as 1/U:
    Tc = min(U*lam/4, 0.369*U*M_min) >= 25.9 meV
    at U = 1.90 eV  ->  lam >= 0.055  AND  M_min >= 0.037
    at U = 0.24 eV  ->  lam >= 0.43   AND  M_min >= 0.292
So a bismuthate-scale U makes the geometric conditions nearly free. The s^2
skippers were excluded on a technicality that a two-site net removes.

THE OBSTRUCTION TESTED HERE. Diamond is bipartite with the two sites related by
inversion about the bond centre. Writing H = d0 + d.sigma in the sublattice
space, a d_z term is ODD under that inversion and therefore FORBIDDEN while the
sites are equivalent. So dhat lies on the equator and the gap is 2|f(k)| with
    f(k) = t (1 + e^{-i q1} + e^{-i q2} + e^{-i q3}).
f = 0 is two real equations in three variables -> a NODAL LINE. Equivalent
sublattices force a gapless semimetal. This is the C1-vs-C3 tension in its
sharpest form: <rho> = 1/2 requires equivalence, equivalence forbids the gap.

THE ESCAPE, which our own C1' already licenses. C1' does NOT require balanced
weights -- it says maximise min(lam/4, 0.369 M_min), and unbalanced weights give
lam > 1/2, which is BETTER. So allow a breathing distortion (a sublattice energy
difference D, exactly what BaBiO3 does) to open the gap, and pay in rho balance.
Scanned below.
"""
import numpy as np
from scipy.optimize import minimize

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
KB = 0.086173e-3


def dvec(Q1, Q2, Q3, D=0.0, t=1.0):
    """s orbitals on diamond; D = sublattice (breathing) energy difference."""
    f = t * (1 + np.exp(-1j * Q1) + np.exp(-1j * Q2) + np.exp(-1j * Q3))
    return np.stack([f.real, -f.imag, np.full_like(f.real, D / 2)], -1)


def grid(nk):
    q = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    return np.meshgrid(q, q, q, indexing='ij'), 2 * np.pi / nk


def analyse(D, nk=32, t=1.0):
    (Q1, Q2, Q3), dk = grid(nk)
    d = dvec(Q1, Q2, Q3, D, t)
    n = np.linalg.norm(d, axis=-1)
    dh = d / n[..., None]
    H = (d[..., 0, None, None] * SX + d[..., 1, None, None] * SY
         + d[..., 2, None, None] * SZ)
    E, V = np.linalg.eigh(H)
    u = V[..., :, 0]

    def trg(delta):
        ang = delta[0] * Q1 + delta[1] * Q2 + delta[2] * Q3
        us = u * np.stack([np.ones_like(ang), np.exp(1j * ang)], -1)
        P = np.einsum('...i,...j->...ij', us, us.conj())
        tot = 0.0
        for ax in (0, 1, 2):
            dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * dk)
            tot = tot + 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
        return float(tot.mean())

    best = trg(np.zeros(3))
    naive = best
    for x0 in [np.zeros(3), [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0],
               [1, 1, 1], [-1, -1, -1], [.5, .5, .5]]:
        r = minimize(trg, np.asarray(x0, float), method="Powell",
                     options=dict(xtol=1e-8, ftol=1e-10))
        best = min(best, float(r.fun))
    rho = (np.abs(u) ** 2).reshape(-1, 2)
    lam = float(np.linalg.eigvalsh((rho.T @ rho) / rho.shape[0])[-1])
    return dict(gap=float((E[..., 1] - E[..., 0]).min()), naive=naive,
                M_min=best, lam=lam, w=rho.mean(0),
                W=float(np.ptp(E[..., 0])))


print(__doc__)
print("=" * 78)
print("1. THRESHOLDS SCALE AS 1/U -- what a bismuthate-scale U buys")
print(f"   {'U (eV)':>8}{'lam needed':>13}{'M_min needed':>15}   source")
for U, src in [(1.90, "BaBiO3 valence skipping (measured scale)"),
               (0.95, "half of it"),
               (0.24, "the Au2+ d9 estimate"),
               (0.158, "the e_g route requirement")]:
    print(f"   {U:>8.2f}{4*0.0259/U:>13.3f}{0.0259/(0.369*U):>15.3f}   {src}")

print("\n2. DIAMOND NET, s ORBITALS, vs sublattice breathing D")
print("   D = 0 is the symmetric (equivalent-site) case.")
print(f"   {'D/t':>6}{'gap/t':>8}{'W/t':>8}{'rho_A/rho_B':>14}{'lam':>8}"
      f"{'M_naive':>10}{'M_min':>9}{'M_min/naive':>13}")
rows = []
for D in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]:
    r = analyse(D)
    rows.append((D, r))
    print(f"   {D:>6.2f}{r['gap']:>8.4f}{r['W']:>8.3f}"
          f"   {r['w'][0]:.3f}/{r['w'][1]:.3f}{r['lam']:>8.4f}"
          f"{r['naive']:>10.4f}{r['M_min']:>9.4f}{r['M_min']/r['naive']:>13.4f}")

print("\n3. Tc AT BISMUTHATE-SCALE U = 1.90 eV")
print("   (t is the free scale: the model is dimensionless, so any t works;")
print("    what matters is lam and M_min, both t-independent)")
U = 1.90
print(f"   {'D/t':>6}{'lam':>8}{'M_min':>9}{'amp (K)':>10}{'geo (K)':>10}"
      f"{'Tc (K)':>9}{'300K?':>7}")
for D, r in rows:
    amp = U * r['lam'] / 4 / KB
    geo = 0.369 * U * r['M_min'] / KB
    tc = min(amp, geo)
    print(f"   {D:>6.2f}{r['lam']:>8.4f}{r['M_min']:>9.4f}{amp:>10.0f}"
          f"{geo:>10.0f}{tc:>9.0f}{'PASS' if tc >= 300 else 'fail':>7}")


# ---------------------------------------------------------------------------
# CORRECTED ANALYSIS. The section above applies U = 1.90 eV WITHOUT enforcing
# the projection condition U <= Delta_iso/2, and gets 746-2856 K. That is
# wrong. Enforcing it, and reporting the flatness ratio:
#
#   D/t  gap(eV)  W(eV)  U<=iso/2  lam    M_min   amp(K) geo(K) Tc(K)  W/Tc
#   0.00   0.196  3.898     0.098  0.500  0.5221     142    219   142  317.8
#   0.25   0.318  3.839     0.159  0.505  0.4652     233    317   233  191.3
#   0.50   0.537  3.736     0.269  0.518  0.3882     404    446   404  107.3
#   1.00   1.019  3.518     0.510  0.557  0.2846     824    621   621   65.7
#                                                    (t = 1.0 eV)
#
# W/Tc IS A LATTICE INVARIANT. Tc = min(lam/4, 0.369 M_min) * gap/2, and gap
# and W both scale linearly with t, so W/Tc is INDEPENDENT of t. Measured:
#   D/t = 0.00  gap/W = 0.050  Tc/t = 0.0123 eV  W/Tc = 317.8
#   D/t = 0.50  gap/W = 0.144  Tc/t = 0.0348 eV  W/Tc = 107.3
# No choice of hopping puts the diamond s-net into the flat-band regime. The
# formula Tc = U*lam/4 requires W <~ Tc; here the band is 65-318x wider.
#
# DESIGN CRITERION that falls out, and it is general:
#     flat-band validity needs  gap/W  >~  2 / min(lam/4, 0.369 M_min)  ~  20
#   diamond s-net           gap/W = 0.05-1.6  -> W/Tc = 66-318   FAILS
#   Kramers-Creutz proposal gap/W = 5.7       -> W/Tc = 4.2      marginal
#   normalised |d| = const  gap/W = infinite  -> W = 0 exactly   PASSES
#
# EMPIRICAL ANCHOR, already in ceiling.py and decisive here. BKBO is this exact
# physics -- s orbital, breathing distortion, ~1.9 eV valence-skipping negative
# U -- and the framework's own ceiling for it is 235 K against a MEASURED 30 K:
# 13% of ceiling, a 7.7x overestimate. BKBO gets no geometric help (n_phi = 1,
# tr g == 0). A diamond net WOULD get geometric help, which is the real
# improvement here, but the starting point for extrapolation is 30 K, not 300.
