"""
LiPbAu2 bands 24/25: is the metric from WINDING, and is 185 K right?

Two questions, one computation.

(1) CONSISTENCY. The stored Tc = 185 K pairs M_trg and lam from the SINGLE-BAND
    descriptors of band 24 with d_iso = 0.390 eV, which is the isolation gap of
    the two-band GROUP. Those are different manifolds. Band 24's own nearest
    neighbour is band 25 at 0.024 eV. Either
      (a) pair in the rank-2 manifold -> U <= 0.195 eV, use M_manifold, lam_man
      (b) pair in band 24 alone       -> U <= 0.012 eV, use M_trg, lam
    but not one from each. This is the same rank-1-descriptor-with-rank-2-gap
    error as D1, in its fifth costume.

(2) WINDING. For ANY two-band model H = d0*I + d.sigma the rank-1 metric is
    EXACTLY tr g = (1/4) sum_i |d_i dhat|^2. So if band 24's metric came from a
    winding of dhat inside the 24/25 pair, a two-band effective model built from
    that pair would reproduce M_trg = 0.547. It can only do so if the 24<->25
    channel dominates. The manifold metric already argues it does not:
        M_manifold = 1.0531   vs   M_24 + M_25 = 1.0869
        cross term = -0.034 = 3% of the total
    A rank-2 projector removes exactly the intra-pair channel, so a 3% cross
    term means 97% of the metric comes from bands OUTSIDE the group. This file
    tests that directly by building the effective two-band model.
"""
import os, sys, zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastops

ZIPS = "/home/user/rtap/data/priority1_wannier_jarvis/zips"
JID, BANDS, DIM, NK = "JVASP-81597", (24, 25), 3, 24
KB = 0.086173e-3
C3 = 0.55 * 0.67

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)

print(__doc__)
print("=" * 78)

with zipfile.ZipFile(os.path.join(ZIPS, "%s.zip" % JID)) as z:
    name = next(n for n in z.namelist() if n.endswith("_hr.dat"))
    R, H, deg = fastops.read_hr_bytes(z.read(name))
norb = H.shape[1]
Hk, dks = fastops.hk_grid_fast(R, H, deg, NK, dim=DIM)
E, V = np.linalg.eigh(Hk)
print(f"{JID}  norb = {norb}  nk = {NK}^3")

b0, b1 = BANDS
gap_pair = (E[..., b1] - E[..., b0])
lo = E[..., b0] - E[..., b0 - 1]
hi = E[..., b1 + 1] - E[..., b1]
iso = np.minimum(lo, hi)
print(f"\nintra-pair splitting  E25-E24 : min {gap_pair.min():.4f}  "
      f"mean {gap_pair.mean():.4f}  max {gap_pair.max():.4f} eV")
print(f"group isolation gap           : min {iso.min():.4f}  "
      f"mean {iso.mean():.4f}  max {iso.max():.4f} eV")

# ---------------------------------------------------------------- (1)
print("\n" + "-" * 78)
print("(1) CONSISTENCY -- the two self-consistent treatments\n")
d = fastops.descriptors_fast(E, V, dks, b0)
dm = fastops.descriptors_manifold_fast(E, V, dks, list(BANDS))
M1, lam1, nphi1 = d['M_trg'], d['lam'], d['nphi']
Mm, lamm = dm['M_trg'], dm['lam']
print(f"   single band 24 : M_trg = {M1:.4f}  lam = {lam1:.4f}  n_phi = {nphi1:.2f}")
print(f"   rank-2 manifold: M_trg = {Mm:.4f}  lam = {lamm:.4f}  n_phi = {1/lamm:.2f}")
print(f"   cross term = M_manifold - (M24 + M25) = "
      f"{Mm - (M1 + fastops.descriptors_fast(E, V, dks, b1)['M_trg']):+.4f}")
print()
print(f"   {'treatment':<34}{'U (eV)':>9}{'M':>9}{'lam':>8}"
      f"{'amp(K)':>9}{'geo(K)':>9}{'Tc(K)':>9}")
for lab, U, M, lm in [
        ("(a) rank-2 manifold, U<=iso/2", float(iso.min()) / 2, Mm, lamm),
        ("(b) single band, U<=pair_gap/2", float(gap_pair.min()) / 2, M1, lam1),
        ("    STORED (mixes a's U, b's M)", float(iso.min()) / 2, M1, lam1)]:
    amp = U * lm / 4 / KB
    geo = C3 * U * M / KB
    print(f"   {lab:<34}{U:>9.4f}{M:>9.4f}{lm:>8.4f}{amp:>9.0f}{geo:>9.0f}"
          f"{min(amp, geo):>9.0f}")

# ---------------------------------------------------------------- (2)
print("\n" + "-" * 78)
print("(2) WINDING -- build the effective two-band model and see if it")
print("    reproduces the metric.\n")
# smooth gauge: project the manifold onto the 2 Wannier orbitals carrying the
# most weight, then Loewdin-orthonormalise.
Um = V[..., :, list(BANDS)]                       # (...,norb,2)
wt = (np.abs(Um) ** 2).sum(axis=(0, 1, 2, 4))     # weight per orbital
pick = list(np.argsort(wt)[::-1][:2])
print(f"   dominant Wannier orbitals: {pick}  "
      f"(weights {wt[pick[0]]:.1f}, {wt[pick[1]]:.1f} of {wt.sum():.1f} total)")
Pi = Um[..., pick, :]                             # (...,2,2) overlap Pi = <a|n>
S = np.einsum('...an,...am->...nm', Pi.conj(), Pi)
ev, evec = np.linalg.eigh(S)
bad = float((ev[..., 0] < 1e-6).mean())
Sm12 = np.einsum('...ij,...j,...kj->...ik', evec, 1 / np.sqrt(np.maximum(ev, 1e-12)), evec.conj())
Eman = np.einsum('...n,...nm->...nm', E[..., list(BANDS)], np.eye(2)[None, None, None])
Heff = np.einsum('...na,...nm,...mb->...ab', Pi.conj(), Eman, Pi)
Heff = np.einsum('...ai,...ij,...jb->...ab', Sm12, Heff, Sm12)
print(f"   fraction of k with near-singular projection: {bad:.4f}")

d0 = 0.5 * np.real(np.einsum('...aa->...', Heff))
dv = np.stack([0.5 * np.real(np.einsum('...ab,ba->...', Heff, SX)),
               0.5 * np.real(np.einsum('...ab,ba->...', Heff, SY)),
               0.5 * np.real(np.einsum('...ab,ba->...', Heff, SZ))], -1)
dn = np.linalg.norm(dv, axis=-1)
print(f"\n   |d(k)|: min {dn.min():.4f}  mean {dn.mean():.4f}  max {dn.max():.4f}"
      f"   uniformity min/max = {dn.min()/dn.max():.4f}")
print(f"   (a true winding structure has |d| EXACTLY constant, uniformity 1.000)")

dhat = dv / np.maximum(dn, 1e-12)[..., None]
trg = np.zeros(dhat.shape[:-1])
for i in range(3):
    dd = (np.roll(dhat, -1, i) - np.roll(dhat, 1, i)) / (2 * dks[i])
    trg += 0.25 * np.sum(dd ** 2, -1)
print(f"\n   M_trg from the effective 2-band model : {trg.mean():.4f}")
print(f"   M_trg of band 24, full Hamiltonian    : {M1:.4f}")
print(f"   ratio = {trg.mean()/M1:.3f}")
print()
if trg.mean() < 0.5 * M1:
    print("   => the 24<->25 channel does NOT carry the metric. The large M_trg")
    print("      comes from mixing with bands OUTSIDE the group, which are")
    print("      >= 0.39 eV away. This is NOT a winding structure.")
else:
    print("   => the effective two-band model accounts for the metric; the")
    print("      winding interpretation survives and the winding number should")
    print("      be computed next.")
