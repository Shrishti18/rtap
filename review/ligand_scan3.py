"""Loophole test for the 8 K ceiling of ligand_scan2/3.

The cubic model isolates the j=1/2 doublet by SOC alone, so isolation dies at
W = 3*lam/2 and the ceiling follows. Real j=1/2 materials (Sr2IrO4) also get a
TETRAGONAL crystal field, which is why they are isolated at t = 0.26 eV where
the cubic model says iso = -1.44 eV. If Delta_tet can hold the manifold apart at
large t, the metric could be large AND isolated, and the 8 K ceiling is an
artefact of the cubic model rather than a real bound.

Test: add Delta_tet to the xy orbital and scan (t, Delta_tet).
"""
import numpy as np
import jeff

KB = 0.086173e-3
LAM = 0.417
C3 = 0.55 * 0.67
NK = 32


def hk_tet(k, t, lam, dtet, tp=0.0):
    H = jeff.hk(k, t, lam, tp)
    cf = np.diag([dtet, 0.0, 0.0])                 # xy shifted by dtet
    H = H + np.kron(cf, np.eye(2))
    return H


def analyse_tet(t, lam, dtet, nk=NK, bands=(4, 5)):
    kk = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), axis=-1)
    dk = 2 * np.pi / nk
    E, V = np.linalg.eigh(hk_tet(K, t, lam, dtet))
    bl = list(bands)
    Um = V[..., :, bl]
    P = np.einsum('...ia,...ja->...ij', Um, Um.conj())
    trg = np.zeros(P.shape[:-2])
    for i in range(3):
        dP = (np.roll(P, -1, i) - np.roll(P, 1, i)) / (2 * dk)
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    r = len(bl)
    rho = np.real(np.einsum('...ii->...i', P)).reshape(-1, 6)
    orb = rho.reshape(-1, 3, 2).sum(-1)
    Am = (orb.T @ orb) / orb.shape[0]
    w = orb.mean(0)
    iso = float(E[..., bl].min() - E[..., bl[0] - 1].max())
    return dict(iso=iso, M_trg=float(trg.mean()),
                lam_pair=float(np.linalg.eigvalsh(Am)[-1]) / r ** 2,
                nphi=float(r ** 2 / np.sum(w ** 2)),
                W=float(max(E[..., j].max() - E[..., j].min() for j in bl)))


print(f"lam = {LAM} eV, nk = {NK}, U capped by the project rule U <= iso/2, R = 1\n")
print(f"{'dtet':>7}{'t':>7}{'iso':>9}{'M_trg':>9}{'nphi':>7}{'lam_p':>8}"
      f"{'U=iso/2':>9}{'Tc_geo':>9}")
print('-' * 65)
best = (0.0,)
for dtet in [0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0]:
    for t in [0.05, 0.10, 0.15, 0.20, 0.26, 0.35, 0.50]:
        d = analyse_tet(t, LAM, dtet)
        U = max(d['iso'], 0.0) / 2
        tc = C3 * U * d['M_trg'] / KB
        if tc > best[0]:
            best = (tc, dtet, t, d, U)
        if d['iso'] > 0:
            print(f"{dtet:>7.2f}{t:>7.2f}{d['iso']:>9.3f}{d['M_trg']:>9.4f}"
                  f"{d['nphi']:>7.3f}{d['lam_pair']:>8.4f}{U:>9.3f}{tc:>9.1f}")
tc, dtet, t, d, U = best
print(f"\nBEST: Tc = {tc:.1f} K at dtet = {dtet:.2f} eV, t = {t:.2f} eV")
print(f"      iso = {d['iso']:.3f}  M_trg = {d['M_trg']:.4f}  nphi = {d['nphi']:.3f}"
      f"  lam_pair = {d['lam_pair']:.4f}  U = {U:.3f} eV")
