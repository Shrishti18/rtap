"""
How fine-tuned is the Kramers-Creutz construction?

The proposed H0 has A-B hopping ONLY at |R| = 2 (cos 2qx, sin 2qx), and EXACTLY
zero at |R| = 0 and |R| = 1. Inversion symmetry does not force that -- for two
orbitals exchanged by inversion, P H(q) P = H(-q) constrains H_AB(q) to be
real-even plus imaginary-odd, which permits nonzero hopping at every shell. So
the vanishing short-range A-B hopping is a fine-tuning assumption, and this
measures how much contamination the construction tolerates.

Two symmetry-allowed perturbations are added to d_x:
    e     * cos(qx)   -- nearest-cell A-B hopping, |R| = 1
    intra * 1         -- intra-cell A-B hopping,   |R| = 0

Reported: spread of |d| (0 = the ideal constant-length texture), the resulting
band width, the isolation gap, and the MINIMISED metric.
"""
import numpy as np
from scipy.optimize import minimize
from kramers_creutz import grid, T, T0, U

nk = 48
(QX, QY, QZ), dk = grid(nk)


def dfull(e=0.0, intra=0.0):
    return np.stack([np.cos(2 * QX) * np.cos(QZ) + e * np.cos(QX) + intra,
                     np.sin(2 * QX),
                     np.cos(2 * QX) * np.sin(QZ)], -1)


def metrics(e=0.0, intra=0.0):
    d = dfull(e, intra)
    n = np.linalg.norm(d, axis=-1)
    dh = d / n[..., None]

    def trg(delta):
        ang = delta[0] * QX + delta[1] * QY + delta[2] * QZ
        c, s = np.cos(ang), np.sin(ang)
        dr = np.stack([dh[..., 0] * c - dh[..., 1] * s,
                       dh[..., 0] * s + dh[..., 1] * c, dh[..., 2]], -1)
        tot = 0.0
        for ax in (0, 1, 2):
            g = (np.roll(dr, -1, ax) - np.roll(dr, 1, ax)) / (2 * dk)
            tot = tot + 0.25 * np.sum(g ** 2, -1)
        return float(tot.mean())

    best = trg(np.zeros(3))
    for x0 in [np.zeros(3), [1, 0, 0], [-1, 0, 0], [0, 0, 1], [2, 0, 0], [-2, 0, 0]]:
        r = minimize(trg, np.asarray(x0, float), method="Powell",
                     options=dict(xtol=1e-7, ftol=1e-9))
        best = min(best, float(r.fun))
    eps0 = 2 * T0 * (np.cos(QX) + np.cos(QY) + np.cos(QZ))
    Elo, Ehi = eps0 - 2 * T * n, eps0 + 2 * T * n
    W = float(np.ptp(Elo))
    return dict(M=best, W=W, iso=float((Ehi - Elo).min()) - W, dn=float(np.ptp(n)))


if __name__ == "__main__":
    print(__doc__)
    print("=" * 74)
    for lbl, key, vals in [("nearest-cell A-B hopping, |R| = 1", "e",
                            [0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50]),
                           ("intra-cell A-B hopping,   |R| = 0", "intra",
                            [0.0, 0.05, 0.10, 0.20, 0.40])]:
        print(f"\nSENSITIVITY to a symmetry-ALLOWED {lbl}")
        print(f"   {'amp':>7}{'spread|d|':>11}{'W(meV)':>9}{'Delta_iso':>11}"
              f"{'M_min':>9}{'2U ok?':>8}")
        for v in vals:
            r = metrics(**{key: v})
            print(f"   {v:>7.2f}{r['dn']:>11.4f}{r['W']*1000:>9.0f}"
                  f"{r['iso']*1000:>11.0f}{r['M']:>9.4f}"
                  f"{'yes' if r['iso'] >= 2*U else 'NO':>8}")
    print("""
READING.
  * M_min is essentially UNAFFECTED by either perturbation: 1.1016 -> 1.0900 at
    50% contamination of the |R| = 1 channel, 1.1016 -> 1.0606 at 40% of the
    intra-cell channel. C2 is robust to fine-tuning failure. The constant-|d|
    texture is not what carries the metric; the winding of dhat is, and that
    survives.
  * What degrades is the ISOLATION. Delta_iso falls 680 -> 405 meV by 20%
    contamination and crosses 2U = 480 meV at roughly 15%. Band width inflates
    in step (120 -> 238 meV).
  * So the fine-tuning requirement is "short-range A-B hopping below ~15% of
    the |R| = 2 hopping" -- a real chemical constraint, not a razor. AND it is
    recoverable: Delta_iso scales with t while 2U is fixed, so a larger t
    restores C3 at the same contamination ratio.
  * CONCLUSION: the vanishing short-range hopping is an assumption without a
    symmetry basis, but the construction is far more tolerant of violating it
    than a naive fine-tuning objection assumes. This weakens, and does not
    sustain, the C2' objection against this particular model.
""")
