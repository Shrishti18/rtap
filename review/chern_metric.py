"""
Does a CHERN winding survive what killed C4?

C4's winding was a labelling convention: d(k) = D[cos(nk)sx + sin(nk)sy] gives
H_12(R) = D delta_{R,n}, one bond, a disconnected dimer lattice, and
    M_naive = n^2/4    M_min = 0    D_s = 0.
The whole "metric" was the Wannier spread across a dimer, removable by putting
the orbitals at their true centres.

The proposed alternative is winding carried by an INTERNAL index (two orbitals /
sublattice / spin), which is what Hermiticity demands and what dipolar m = +-1
exchange, Slater-Koster complex orbitals, and circular Floquet drive all supply.
The question is whether that evades the C4 failure or merely dresses it.

It evades it, and for a reason that is a theorem rather than a model detail:

    tr g >= |Omega|  pointwise, in ANY gauge      (Omega = Berry curvature)
    => <tr g> >= (1/(2pi)^2) |int Omega| = |C| / (2pi)

The Chern number is a topological invariant: it cannot be removed by orbital
shifts, so this LOWER bound survives the minimal-metric operation that took C4
to zero. A Chern band has a floor on its gauge-invariant metric.

    M_trg >= |C| / (2pi) = 0.1592 |C|

Thresholds: 0.339 (3D) needs |C| >= 2.13; 0.839 (2D) needs |C| >= 5.27 -- from
the FLOOR alone. Actual metric can exceed it, so these are not requirements,
but they say a C = 1 band is not automatically enough.

Tested on Qi-Wu-Zhang, the minimal two-band Chern model:
    H(k) = sin kx sx + sin ky sy + (m + cos kx + cos ky) sz
    C = -1 for 0 < m < 2, C = +1 for -2 < m < 0, C = 0 for |m| > 2.
"""
import numpy as np
from scipy.optimize import minimize

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def hk(kx, ky, m):
    return (np.sin(kx)[..., None, None] * SX
            + np.sin(ky)[..., None, None] * SY
            + (m + np.cos(kx) + np.cos(ky))[..., None, None] * SZ)


def analyse(m, nk=96):
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KX, KY = np.meshgrid(k, k, indexing='ij')
    dk = 2 * np.pi / nk
    E, V = np.linalg.eigh(hk(KX, KY, m))
    u = V[..., :, 0]
    P = np.einsum('...i,...j->...ij', u, u.conj())

    def trg_of(sh):
        """<tr g> with orbital shifts sh (shape (2,2): orbital, dim)."""
        ph = np.exp(1j * (KX[..., None] * sh[:, 0] + KY[..., None] * sh[:, 1]))
        us = u * ph
        Q = np.einsum('...i,...j->...ij', us, us.conj())
        t = np.zeros(Q.shape[:-2])
        for ax in (0, 1):
            dQ = (np.roll(Q, -1, ax) - np.roll(Q, 1, ax)) / (2 * dk)
            t += 0.5 * np.real(np.einsum('...ij,...ji->...', dQ, dQ))
        return float(t.mean())

    naive = trg_of(np.zeros((2, 2)))
    obj = lambda x: trg_of(np.array([[0.0, 0.0], [x[0], x[1]]]))
    best = naive
    for x0 in [np.zeros(2), np.array([0.5, 0.5]), np.array([-0.5, -0.5]),
               np.array([1.0, 0.0]), np.array([0.0, 1.0])]:
        r = minimize(obj, x0, method="Powell",
                     options=dict(xtol=1e-8, ftol=1e-10))
        best = min(best, float(r.fun))

    # Berry curvature by the gauge-invariant plaquette (Fukui-Hatsugai-Suzuki)
    def link(a, ax):
        b = np.roll(a, -1, ax)
        o = np.einsum('...i,...i->...', a.conj(), b)
        return o / np.abs(o)
    U1, U2 = link(u, 0), link(u, 1)
    F = np.angle(U1 * np.roll(U2, -1, 0) / np.roll(U1, -1, 1) / U2)
    C = F.sum() / (2 * np.pi)
    gap = float((E[..., 1] - E[..., 0]).min())
    return dict(naive=naive, M_min=best, C=C, gap=gap,
                floor=abs(C) / (2 * np.pi))


print(__doc__)
print("=" * 78)
print(f"{'m':>7}{'C':>7}{'gap':>8}{'M_naive':>10}{'M_min':>10}"
      f"{'M_min/naive':>13}{'floor |C|/2pi':>15}{'ok':>5}")
for m in [-3.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 3.0]:
    r = analyse(m)
    ok = "yes" if r["M_min"] >= r["floor"] - 1e-6 else "NO"
    print(f"{m:>7.1f}{r['C']:>7.2f}{r['gap']:>8.3f}{r['naive']:>10.4f}"
          f"{r['M_min']:>10.4f}{r['M_min']/r['naive']:>13.4f}"
          f"{r['floor']:>15.4f}{ok:>5}")

print("""
READING.
  * The Chern floor is real and gauge-invariant: M_min >= |C|/2pi is respected
    at every m, and for C = +-1 the metric CANNOT be taken to zero. That is
    precisely what C4 lacked -- its winding carried no flux, so the whole
    metric was removable.
  * CAVEAT on this fixture: QWZ's two orbitals are CO-LOCATED (both at the site
    origin), so there is no orbital-shift freedom to exploit and M_min =
    M_naive identically, including at C = 0. This model therefore demonstrates
    the floor but does NOT demonstrate removability in the trivial case. C4's
    collapse needed orbitals at 0 and n. A displaced-orbital Chern model would
    be needed to test that half.
  * THE NUMBERS ARE THE PROBLEM. C = +-1 gives M_trg = 0.19-0.29, which is
    BELOW the 3D threshold 0.339 and 3x below the 2D threshold 0.839. The
    topological floor protects the metric from being gauge; it does not supply
    enough of it. M_min/floor = 1.79 at best here, so clearing 3D needs C >= 2
    and 2D needs C >= 3 if that ratio persists -- i.e. HIGH-Chern flat bands,
    not the C = 1 models that dominate the literature.
""")
