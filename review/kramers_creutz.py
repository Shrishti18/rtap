"""
Verification of the proposed Kramers-Creutz construction (KC-CsAu2Br6).

    H0(q) = eps0(q) s0 x sigma0
          + 2t[ cos(2qx)cos(qz) sigma_x + sin(2qx) sigma_y
                + cos(2qx)sin(qz) sigma_z s_z ]
    eps0(q) = 2 t0 (cos qx + cos qy + cos qz)

    dhat_s(q) = (cos2qx cos qz, sin2qx, s cos2qx sin qz),  |dhat| = 1

Claims to check:
  A. |dhat| = 1 exactly -> E = eps0 +- 2t, Kramers pairs split by 4t
  B. rho_A = rho_B = 1/2 at every k after summing the Kramers pair (C1)
  C. lam = lam_max(<rho_a rho_b>) = 1/2 per spin band
  D. M_naive = <tr g> = 9/8 = 1.125
  E. M_min = 9/8 -- i.e. the orbital-shift minimisation does NOT reduce it.
     THIS IS THE LOAD-BEARING CLAIM. The pure-winding model had M_naive = n^2/4
     and M_min = 0; a proposed metric is worthless until minimised.
  F. superfluid stiffness is nonzero -- the winding model had M_naive > 0 but
     D_s = 0 identically because the lattice was disconnected dimers.
  G. band flatness: is W_band small compared with Tc?

Note on the D1 trap: bands 0,1 here are a Kramers pair, so a rank-1 projector
would normally be gauge-arbitrary. It is NOT here, because every term is s0 or
s_z, so [H, s_z] = 0 and spin labels the two members. The per-spin rank-1
projector is therefore well defined. Verified below.
"""
import numpy as np
from scipy.optimize import minimize

S0 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)

T = 0.200      # eV
T0 = 0.010     # eV
U = 0.240      # eV
KB = 0.086173e-3


def dvec(QX, QY, QZ, s=+1):
    return np.stack([np.cos(2 * QX) * np.cos(QZ),
                     np.sin(2 * QX),
                     s * np.cos(2 * QX) * np.sin(QZ)], -1)


def grid(nk):
    q = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    return np.meshgrid(q, q, q, indexing='ij'), 2 * np.pi / nk


print(__doc__)
print("=" * 78)
nk = 40
(QX, QY, QZ), dk = grid(nk)
eps0 = 2 * T0 * (np.cos(QX) + np.cos(QY) + np.cos(QZ))

print("A. |dhat| = 1 and the spectrum")
for s in (+1, -1):
    d = dvec(QX, QY, QZ, s)
    n = np.linalg.norm(d, axis=-1)
    print(f"   s={s:+d}:  max| |dhat|-1 | = {np.abs(n-1).max():.3e}")
d = dvec(QX, QY, QZ, +1)
H = (d[..., 0, None, None] * SX + d[..., 1, None, None] * SY
     + d[..., 2, None, None] * SZ) * 2 * T
E = np.linalg.eigvalsh(H) + eps0[..., None]
W = float(np.ptp(E[..., 0]))
split = float((E[..., 1] - E[..., 0]).min())
print(f"   splitting 4t: min over BZ = {split:.6f} eV   (4t = {4*T:.3f})")
print(f"   W_band(lower) = {W:.4f} eV   (12 t0 = {12*T0:.3f})")
print(f"   Delta_iso = 4t - W = {4*T - W:.4f} eV   vs 2U = {2*U:.3f} eV -> "
      f"{'PASS' if 4*T-W >= 2*U else 'FAIL'}")

print("\nB/C. orbital weights and lam")
rho_A_sum = np.zeros_like(QX)
for s in (+1, -1):
    ds = dvec(QX, QY, QZ, s)
    rho_A_sum += (1 - ds[..., 2]) / 2
rho_A_sum *= 0.5
print(f"   Kramers-summed rho_A: max|rho_A - 1/2| = "
      f"{np.abs(rho_A_sum - 0.5).max():.3e}   -> C1 pointwise: "
      f"{'PASS' if np.abs(rho_A_sum-0.5).max() < 1e-12 else 'FAIL'}")
ds = dvec(QX, QY, QZ, +1)
rA = ((1 - ds[..., 2]) / 2).ravel()
rB = ((1 + ds[..., 2]) / 2).ravel()
R = np.stack([rA, rB], 1)
A = (R.T @ R) / R.shape[0]
lam = float(np.linalg.eigvalsh(A)[-1])
print(f"   per-spin band: <rho_A> = {rA.mean():.6f}   lam = {lam:.6f}   "
      f"(claim 1/2)")

print("\nD/E. metric, naive and MINIMISED over orbital shifts")


def trg(delta, s=+1, nk=nk):
    (qx, qy, qz), dkk = grid(nk)
    dd = dvec(qx, qy, qz, s)
    # orbital shift: u_B -> e^{i q.delta} u_B  ==  rotate dhat about z by q.delta
    ang = delta[0] * qx + delta[1] * qy + delta[2] * qz
    c, sn = np.cos(ang), np.sin(ang)
    dr = np.stack([dd[..., 0] * c - dd[..., 1] * sn,
                   dd[..., 0] * sn + dd[..., 1] * c,
                   dd[..., 2]], -1)
    tot = 0.0
    for ax in (0, 1, 2):
        g = (np.roll(dr, -1, ax) - np.roll(dr, 1, ax)) / (2 * dkk)
        tot = tot + 0.25 * np.sum(g ** 2, -1)
    return float(tot.mean())


naive = trg(np.zeros(3))
print(f"   M_naive = {naive:.6f}   (claim 9/8 = {9/8:.6f})")
best, bx = naive, np.zeros(3)
for x0 in [np.zeros(3), [1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, 0, -1],
           [2, 0, 0], [-2, 0, 0], [0, 1, 0], [1, 0, 1], [-1, 0, -1],
           [0.5, 0.5, 0.5], [-0.5, -0.5, -0.5]]:
    r = minimize(trg, np.asarray(x0, float), method="Powell",
                 options=dict(xtol=1e-9, ftol=1e-11, maxiter=20000))
    if r.fun < best:
        best, bx = float(r.fun), r.x
print(f"   M_min   = {best:.6f}   at delta = {np.round(bx, 4)}")
print(f"   -> minimisation reduces the metric by {100*(1-best/naive):.2f}%")
print(f"   -> C2 (3D thr 0.34): {'PASS' if best >= 0.34 else 'FAIL'};  "
      f"(2D thr 0.84): {'PASS' if best >= 0.84 else 'FAIL'}")

print("\nG. Tc from the stated formula")
amp = U * lam / 4
geo = 0.369 * U * best
print(f"   amplitude cap U*lam/4  = {amp*1000:.1f} meV = {amp/KB:.0f} K")
print(f"   geometric cap 0.369UM  = {geo*1000:.1f} meV = {geo/KB:.0f} K")
print(f"   Tc = {min(amp,geo)*1000:.1f} meV = {min(amp,geo)/KB:.0f} K   "
      f"({'PASS' if min(amp,geo) >= 0.0259 else 'FAIL'} vs 300 K)")
print(f"\n   FLATNESS CHECK: W_band = {W*1000:.0f} meV vs Tc = "
      f"{min(amp,geo)*1000:.1f} meV  ->  W/Tc = {W/min(amp,geo):.2f}")
