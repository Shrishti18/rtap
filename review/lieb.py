"""
Is <tr g> of a TOUCHING flat band a number at all?

Lieb lattice, 3 sites/cell (corner A, edges B, C). PERIODIC GAUGE:

    H_AB = t(1 + e^{-i kx}),   H_AC = t(1 + e^{-i ky}),   H_BC = 0

so H(k + 2pi) = H(k) EXACTLY. (The textbook 2t cos(k/2) form is NOT 2pi
periodic -- H(0) = -H(2pi) -- and differencing across the BZ boundary in that
gauge injects a spurious discontinuity. That is the BZ-boundary term this
project already fixed once in harvest._trg_shifted; the first version of this
file walked into it and reported a divergence peaked at GAMMA, where the gap is
LARGEST. The tell was that the maximum sat at the biggest gap, not the smallest.)

Spectrum: E = 0 (flat) and E = +-sqrt(|H_AB|^2 + |H_AC|^2). At M = (pi,pi) both
off-diagonals vanish, so all three bands meet -- a threefold touching, and it is
protected: an on-site term on A leaves the flat state (which has zero amplitude
on A) exactly where it is.

THE QUESTION. At the touching the rank-1 projector is discontinuous, tr g ~
1/|k-M|^2, and int d^2k/|k|^2 is log divergent in 2D. If <tr g> has no limit,
then c_eff = D_s/(U*M) measured at one grid is not a suppressed constant -- it
is a finite, well-defined D_s divided by a quantity that does not converge.
"""
import numpy as np

T = 1.0


def hk(kx, ky, m=0.0):
    """Periodic gauge. m is an on-site energy on the A (corner) site."""
    H = np.zeros(kx.shape + (3, 3), complex)
    A = T * (1 + np.exp(-1j * kx))
    B = T * (1 + np.exp(-1j * ky))
    H[..., 0, 1] = A
    H[..., 1, 0] = np.conj(A)
    H[..., 0, 2] = B
    H[..., 2, 0] = np.conj(B)
    H[..., 0, 0] = m
    return H


def metric(nk, m=0.0, rcut=None):
    """<tr g> of the flat band. rcut excises a disc of that radius around M."""
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KX, KY = np.meshgrid(k, k, indexing='ij')
    dk = 2 * np.pi / nk
    E, V = np.linalg.eigh(hk(KX, KY, m))
    j = int(np.argmin([E[..., b].max() - E[..., b].min() for b in range(3)]))
    u = V[..., :, j]
    P = np.einsum('...i,...j->...ij', u, u.conj())
    trg = np.zeros(P.shape[:-2])
    for ax in (0, 1):
        dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * dk)
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    d = np.sqrt((KX - np.pi) ** 2 + (KY - np.pi) ** 2)
    mask = np.ones_like(trg, bool) if rcut is None else (d > rcut)
    others = [b for b in range(3) if b != j]
    gap = float(min(np.abs(E[..., b] - E[..., j]).min() for b in others))
    i0 = np.unravel_index(np.argmax(trg), trg.shape)
    return dict(mean=float(trg[mask].mean()), W=float(np.ptp(E[..., j])),
                gap=gap, kmax=(KX[i0] / np.pi, KY[i0] / np.pi),
                gap_at_max=float(np.abs(E[..., others[0]] - E[..., j])[i0]))


print(__doc__)
print("=" * 76)
print("1. Periodicity check (the thing the first version got wrong)")
k0 = np.array([[0.3]])
h1, h2 = hk(k0, k0), hk(k0 + 2 * np.pi, k0 + 2 * np.pi)
print(f"   max|H(k) - H(k+2pi)| = {np.abs(h1-h2).max():.2e}   (must be ~0)")
old = lambda k: 2 * T * np.cos(k / 2)
print(f"   textbook gauge: 2t cos(k/2) at k=0 vs k=2pi -> {old(0.0):+.1f} vs"
      f" {old(2*np.pi):+.1f}   NOT periodic")

print()
print("2. Is the touching protected against an on-site m on A?")
print(f"   {'m':>6}{'W_flat':>11}{'min gap':>11}")
for m in [0.0, 0.5, 2.0]:
    r = metric(192, m)
    print(f"   {m:>6.2f}{r['W']:>11.2e}{r['gap']:>11.2e}")
print("   flat band unmoved, gap stays at grid resolution -> protected.")

print()
print("3. Does <tr g> converge?   (where does the max sit?)")
print(f"   {'nk':>6}{'<tr g>':>11}{'ratio':>8}{'/log(nk)':>11}{'k_max/pi':>18}"
      f"{'gap there':>11}")
prev = None
for nk in [48, 96, 192, 384, 768]:
    r = metric(nk)
    ratio = r['mean'] / prev if prev else float('nan')
    print(f"   {nk:>6}{r['mean']:>11.4f}{ratio:>8.3f}{r['mean']/np.log(nk):>11.4f}"
          f"   ({r['kmax'][0]:.3f},{r['kmax'][1]:.3f}){r['gap_at_max']:>13.2e}")
    prev = r['mean']
print("   convergent -> ratio 1.0 ; log divergent -> /log(nk) constant.")
print("   The max must sit AT M=(1,1) where the gap closes, not at Gamma.")

print()
print("4. Excision control: same integral, disc of radius rcut around M removed.")
print("   If the divergence is the touching, excising it must restore convergence.")
print(f"   {'rcut':>7}{'nk=192':>11}{'nk=384':>11}{'nk=768':>11}{'ratio 768/384':>15}")
for rc in [0.30, 0.15, 0.05]:
    vals = [metric(nk, rcut=rc)['mean'] for nk in (192, 384, 768)]
    print(f"   {rc:>7.2f}{vals[0]:>11.4f}{vals[1]:>11.4f}{vals[2]:>11.4f}"
          f"{vals[2]/vals[1]:>15.3f}")
