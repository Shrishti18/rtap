"""
Can an e_g doublet on a triangular lattice reach the reference model?

MOTIVATION. The reference model needs |C| >= 2, i.e. an e^{2i theta} winding in
a two-orbital space. The e_g doublet {d_z2, d_x2-y2} carries exactly that: under
threefold rotation it transforms with l_eff = 2, so the two-centre Slater-Koster
hopping along a bond at angle theta is the orbital-compass form

    T(theta) = t0 * [ a*I + b*( cos(2 theta) tau_z + sin(2 theta) tau_x ) ]

That is the chemistry the (sin kx + i sin ky)^C construction was asking for, and
it comes with the valence-skipping negative U of a d9 ion (Au2+: d9 -> d10 + d8)
on the SAME two orbitals -- which is what BaBiO3 could not do, its 1.9 eV
negative U sitting on 6s where n_phi = 1 and tr g == 0 identically.

THE OBSTRUCTION TESTED HERE. The compass hopping is REAL: it generates only
tau_z and tau_x, never tau_y. So dhat is confined to a great circle of S^2, the
map T^2 -> S^2 has degree 0, and C = 0 no matter how large the winding. A mass
term is required, and within e_g on-site SOC CANNOT supply it -- e_g has quenched
orbital angular momentum, so <e_g| L |e_g> = 0 and lambda_SOC does not act at
first order. tau_y is TRS-odd, so it must come from broken time reversal:
Haldane-like complex second-neighbour hopping (orbital loop currents), or
magnetic order.

Tested: (1) compass alone -> C = 0? (2) with a Haldane tau_y term -> |C| = 2?
(3) M_min and lambda of the result, against the reference model.
"""
import numpy as np
from scipy.optimize import minimize

TX = np.array([[0, 1], [1, 0]], complex)
TY = np.array([[0, -1j], [1j, 0]], complex)
TZ = np.array([[1, 0], [0, -1]], complex)

# triangular lattice: NN at 0, 120, 240 deg ; NNN at 60, 180, 300 deg
NN = [np.array([np.cos(t), np.sin(t)]) for t in np.deg2rad([0, 120, 240])]
NNN = [np.array([np.cos(t), np.sin(t)]) * np.sqrt(3)
       for t in np.deg2rad([60, 180, 300])]


def hk(KX, KY, a=0.0, b=1.0, t2=0.0, m=0.0):
    """e_g compass on triangular + Haldane-like tau_y + on-site tau_z mass."""
    K = np.stack([KX, KY], -1)
    d = np.zeros(KX.shape + (3,))
    for dl in NN:                       # real compass hopping, l_eff = 2
        th = np.arctan2(dl[1], dl[0])
        c = 2 * np.cos(K @ dl)
        d[..., 0] += c * b * np.sin(2 * th)      # tau_x
        d[..., 2] += c * b * np.cos(2 * th)      # tau_z
    for dl in NNN:                      # TRS-breaking tau_y
        d[..., 1] += t2 * 2 * np.sin(K @ dl)
    d[..., 2] += m
    return d


def analyse(nk=120, **kw):
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KX, KY = np.meshgrid(k, k, indexing='ij')
    dk = 2 * np.pi / nk
    d = hk(KX, KY, **kw)
    n = np.linalg.norm(d, axis=-1)
    H = (d[..., 0, None, None] * TX + d[..., 1, None, None] * TY
         + d[..., 2, None, None] * TZ)
    E, V = np.linalg.eigh(H)
    u = V[..., :, 0]

    def trg_of(sh):
        us = u * np.exp(1j * (KX[..., None] * sh[:, 0] + KY[..., None] * sh[:, 1]))
        P = np.einsum('...i,...j->...ij', us, us.conj())
        t = np.zeros(P.shape[:-2])
        for ax in (0, 1):
            dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * dk)
            t += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
        return float(t.mean())

    naive = trg_of(np.zeros((2, 2)))
    best = naive
    for x0 in [np.zeros(2), [.5, .5], [-.5, -.5], [1., 0.], [0., 1.]]:
        r = minimize(lambda x: trg_of(np.array([[0., 0.], [x[0], x[1]]])),
                     np.asarray(x0, float), method="Powell",
                     options=dict(xtol=1e-8, ftol=1e-10))
        best = min(best, float(r.fun))

    def link(a_, ax):
        o = np.einsum('...i,...i->...', a_.conj(), np.roll(a_, -1, ax))
        return o / np.abs(o)
    U1, U2 = link(u, 0), link(u, 1)
    C = np.angle(U1 * np.roll(U2, -1, 0) / np.roll(U1, -1, 1) / U2).sum() / (2 * np.pi)
    rho = (np.abs(u) ** 2).reshape(-1, 2)
    lam = float(np.linalg.eigvalsh((rho.T @ rho) / rho.shape[0])[-1])
    return dict(C=C, gap=float((E[..., 1] - E[..., 0]).min()),
                W=float(np.ptp(E[..., 0])), naive=naive, M_min=best,
                lam=lam, w=rho.mean(0))


print(__doc__)
print("=" * 78)
print("1. COMPASS ALONE (t2 = 0): is the map degree-0 as predicted?")
r = analyse(t2=0.0, m=0.0)
print(f"   C = {r['C']:+.4f}   gap = {r['gap']:.4f}   M_naive = {r['naive']:.4f}"
      f"   M_min = {r['M_min']:.4f}")
print("   -> real hopping, dhat confined to the tau_x-tau_z great circle.")

print("\n2. WITH TRS-BREAKING tau_y (Haldane-like NNN, amplitude t2)")
print(f"   {'t2':>6}{'m':>6}{'C':>9}{'gap':>8}{'W_band':>9}{'M_naive':>10}"
      f"{'M_min':>9}{'lam':>8}{'w1/w2':>12}")
for t2 in [0.15, 0.35, 0.6]:
    for m in [0.0, 1.0]:
        r = analyse(t2=t2, m=m)
        print(f"   {t2:>6.2f}{m:>6.2f}{r['C']:>9.4f}{r['gap']:>8.4f}"
              f"{r['W']:>9.3f}{r['naive']:>10.4f}{r['M_min']:>9.4f}"
              f"{r['lam']:>8.4f}   {r['w'][0]:.2f}/{r['w'][1]:.2f}")

print("""
READING -- and the first line is a trap I set for myself and fell into.

  * "C = +2.0000" at t2 = 0 comes with gap = 0.0012, i.e. the bands TOUCH. The
    Chern number of a gapless band is not defined and the number is numerical
    noise at the touching points, where dhat is singular. It is not a C = 2
    state.
  * Reading the table by GAP instead: every gapped row (m = 1, gap 0.30-0.56)
    has C = 0.0000, and every row with C != 0 has gap ~ 0.002. GAPPED => C = 0.
    That is the predicted result: the compass hopping is REAL, generating only
    tau_z and tau_x, so dhat is confined to a great circle of S^2 and any
    gapped map has degree zero. l_eff = 2 supplies the WINDING but not the
    DEGREE.
  * On-site SOC cannot rescue it: e_g has QUENCHED orbital angular momentum,
    <e_g|L|e_g> = 0, so lambda_SOC does not act at first order within the
    doublet. tau_y is TRS-odd and must come from broken time reversal. The
    Haldane-like NNN term tried here does not gap it either -- it vanishes at
    the touching points.
  * THE BIND. |C| != 0 on an e_g doublet requires MAGNETISM or orbital loop
    currents. But the valence-skipping negative U (d9 + d9 -> d10 + d8) moves a
    SPIN-SINGLET pair and needs both spin species. A ferromagnetic solution is
    self-defeating.
  * THE RESOLUTION, and the thing to look for: a QSH-like state, spin-up with
    C = +2 and spin-down with C = -2, total C = 0 so time reversal is intact
    and singlet pairing is allowed, with each spin sector carrying
    M_min >= 2/(2pi) = 0.318. Note the gapped rows already show M_min =
    0.66-0.87, ABOVE the 0.34 threshold, and lam = 0.56-0.58 -- so the metric
    and pairing eigenvalue are not the problem here. Only the degree is.
""")
