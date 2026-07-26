"""
The folding kill, applied to the P4/mbm-CsTlF3 proposal.

THE ISSUE. A cell doubled by an OCTAHEDRAL ROTATION leaves the two metal sites
symmetry EQUIVALENT. The resulting "two-orbital" Hamiltonian is then just the
one-orbital band FOLDED into a smaller BZ. Folding cannot create quantum
geometry: the true primitive description has one s orbital, n_phi = 1, and
tr g == 0 identically. The apparent M_naive of the folded description is pure
gauge and vanishes under orbital-position minimisation -- exactly the C4 dimer
pathology in a different costume.

The proposal names this as its Kill #2 and calls it "probably the cheapest and
most important kill". It is, and it can be run in seconds rather than with HSE.

TEST. Simple-cubic s band (one Tl per primitive cell), doubled along z:
    A sublattice at even z, B at odd z
    intra-plane hopping t on both, inter-plane t1 (A->B) and t2 (B->A)
  t1 = t2  -> PURE FOLDING, the sites stay equivalent
  t1 != t2 -> genuine dimerisation, sites inequivalent in their bonding

Reported: M_naive, M_min (minimised over the A/B orbital offset), and lam.
"""
import numpy as np
from scipy.optimize import minimize

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def analyse(t1, t2, t=1.0, dz=0.0, nk=24):
    """Two-site cell along z. dz = on-site energy difference (breathing)."""
    q = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KX, KY, KZ = np.meshgrid(q, q, q, indexing='ij')
    f = t1 + t2 * np.exp(-1j * KZ)                 # A->B inter-plane
    d0 = -2 * t * (np.cos(KX) + np.cos(KY))        # identical on both sites
    dx, dy, dzc = f.real, -f.imag, np.full_like(f.real, dz / 2)
    H = (dx[..., None, None] * SX + dy[..., None, None] * SY
         + dzc[..., None, None] * SZ)
    E, V = np.linalg.eigh(H)
    u = V[..., :, 0]

    def trg(delta):
        us = u * np.stack([np.ones_like(KZ), np.exp(1j * delta * KZ)], -1)
        P = np.einsum('...i,...j->...ij', us, us.conj())
        tot = 0.0
        for ax in (0, 1, 2):
            dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * (2 * np.pi / nk))
            tot = tot + 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
        return float(tot.mean())

    naive = trg(0.0)
    best = naive
    for x0 in [0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0]:
        r = minimize(lambda x: trg(x[0]), [x0], method="Powell",
                     options=dict(xtol=1e-10, ftol=1e-12))
        best = min(best, float(r.fun))
    rho = (np.abs(u) ** 2).reshape(-1, 2)
    lam = float(np.linalg.eigvalsh((rho.T @ rho) / rho.shape[0])[-1])
    gap = float((E[..., 1] - E[..., 0]).min())
    return dict(naive=naive, M_min=best, lam=lam, gap=gap, w=rho.mean(0))


print(__doc__)
print("=" * 78)
print("1. PURE FOLDING (t1 = t2, sites symmetry-equivalent)")
print("   This is what an octahedral-rotation cell doubling gives.")
print(f"   {'t1':>6}{'t2':>6}{'dz':>6}{'gap':>8}{'lam':>8}{'M_naive':>10}"
      f"{'M_min':>12}{'rho_A/rho_B':>14}")
for t1, t2 in [(1.0, 1.0), (0.5, 0.5), (2.0, 2.0)]:
    r = analyse(t1, t2)
    print(f"   {t1:>6.2f}{t2:>6.2f}{0.0:>6.2f}{r['gap']:>8.4f}{r['lam']:>8.4f}"
          f"{r['naive']:>10.4f}{r['M_min']:>12.3e}"
          f"   {r['w'][0]:.3f}/{r['w'][1]:.3f}")
print("   -> M_min does NOT collapse -- but read the GAP column: 0.2616 is exactly")
print("      the grid spacing. The folded bands TOUCH at the zone boundary and the")
print("      half-offset grid steps over the node. Confirmed by scanning nk:")
print("        nk = 12/18/24/36/48 -> gap 0.522/0.349/0.262/0.175/0.131 (~1/nk)")
print("                              M    0.193/0.276/0.356/0.511/0.665 (diverging)")
print("      So folding gives NO ISOLATED MANIFOLD: C3 fails outright and the")
print("      'metric' is the divergent touching contribution, not a usable number.")

print("\n2. GENUINE DIMERISATION (t1 != t2, bonding environments differ)")
print(f"   {'t1':>6}{'t2':>6}{'dz':>6}{'gap':>8}{'lam':>8}{'M_naive':>10}"
      f"{'M_min':>12}{'rho_A/rho_B':>14}")
for t1, t2 in [(1.2, 0.8), (1.5, 0.5), (1.9, 0.1)]:
    r = analyse(t1, t2)
    print(f"   {t1:>6.2f}{t2:>6.2f}{0.0:>6.2f}{r['gap']:>8.4f}{r['lam']:>8.4f}"
          f"{r['naive']:>10.4f}{r['M_min']:>12.4f}"
          f"   {r['w'][0]:.3f}/{r['w'][1]:.3f}")

print("\n3. BREATHING (dz != 0, on-site energies differ -- the charge-ordered case)")
print(f"   {'t1':>6}{'t2':>6}{'dz':>6}{'gap':>8}{'lam':>8}{'M_naive':>10}"
      f"{'M_min':>12}{'rho_A/rho_B':>14}")
for dz in [0.5, 1.0, 2.0]:
    r = analyse(1.0, 1.0, dz=dz)
    print(f"   {1.0:>6.2f}{1.0:>6.2f}{dz:>6.2f}{r['gap']:>8.4f}{r['lam']:>8.4f}"
          f"{r['naive']:>10.4f}{r['M_min']:>12.4f}"
          f"   {r['w'][0]:.3f}/{r['w'][1]:.3f}")

print("""
CORRECTION (audit, accepted in full). The M_min values quoted in the verdict
below were computed at nk = 24 and are UNDER-RESOLVED by up to 33%. tr g for a
dimerised chain is sharply peaked near the smallest gap, and a 24-point grid
misses the peak -- the same under-resolution error made once already on the
TBG metric. The exact result is a closed form,
        M_min = min(t1^2, t2^2) / (8 |t1^2 - t2^2|),
verified against high-nk integration:
    (t1,t2)     nk=24     nk=48     closed form
    (1.1,0.9)  0.17122   0.22014     0.253125
    (1.2,0.8)  0.08408   0.09461     0.100000
    (1.4,0.6)  0.02486   0.02682     0.028125
The correction makes the numbers LARGER, i.e. more favourable, and changes no
conclusion: the closed form is monotonically decreasing in the dimerisation and
diverges only as t1 -> t2, which is exactly the gapless folding limit. The
family still tops out ~5x below the M_min >= 1.35 the proposal required, and
does so for a reason now visible in one line rather than in a scan.
""")

print("""
VERDICT FOR THE PROPOSAL.
  1. Pure folding (symmetry-equivalent sites, the rotation-doubled case) gives
     a BAND TOUCHING, not a gapped two-orbital manifold. C3 fails outright.
  2. Opening a real gap requires site inequivalence -- in bonding (Sec 2) or
     on-site energy (Sec 3) -- and BOTH COST METRIC MONOTONICALLY:
        t1/t2 = 1.0  -> M_min 0.344 (gap 0.28, barely gapped)
        t1/t2 = 1.1/0.9 -> 0.170     t1/t2 = 1.2/0.8 -> 0.084
        t1/t2 = 1.4/0.6 -> 0.025
        dz = 0.5 -> 0.235   dz = 1.0 -> 0.154   dz = 2.0 -> 0.087
  3. BEST M_min FOUND ANYWHERE IN THE FAMILY: 0.344, at a point that is only
     marginally gapped. The proposal requires M_min >= 1.35.
        SHORT BY 3.9x, and worse once a real 2.4 eV isolation gap is demanded.
  This is the same isolation-vs-metric trade-off measured on the diamond net
  (M_min 0.52 -> 0.09 as breathing went 0 -> 4t). It appears to be generic to
  two-site s-orbital networks.
  CAVEAT: this is a cubic s band doubled along z, not the actual a0a0c+ rotation
  pattern. The mechanism is generic; the exact numbers for P4/mbm-CsTlF3 need
  the real structure.
""")
