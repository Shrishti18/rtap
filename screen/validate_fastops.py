#!/usr/bin/env python3
"""Prove fastops reproduces the shipped harvest v4 routines exactly."""
import numpy as np, harvest, fastops

print("=== 1. hk_grid: shipped vs fast ===")
for tag in ["M0.0", "M0.9"]:
    R, H, d = harvest.read_hr("cb_%s_hr.dat" % tag)
    A, dka = harvest.hk_grid(R, H, d, 60, dim=2)
    B, dkb = fastops.hk_grid_fast(R, H, d, 60, dim=2)
    print("  %-6s max|diff| = %.3e   dks equal: %s"
          % (tag, np.abs(A - B).max(), dka == dkb))
R, H, d = harvest.read_hr("hr/JVASP-13733_hr.dat")
A, _ = harvest.hk_grid(R, H, d, 8, dim=3)
B, _ = fastops.hk_grid_fast(R, H, d, 8, dim=3)
print("  real 60-orb  max|diff| = %.3e" % np.abs(A - B).max())

print("\n=== 2. unshifted <tr g>: metric_trace vs rank-1 identity ===")
for tag in ["M0.0", "M0.9"]:
    R, H, d = harvest.read_hr("cb_%s_hr.dat" % tag)
    Hk, dks = harvest.hk_grid(R, H, d, 120, dim=2)
    for b in (0, 1):
        P, _ = harvest.projector(Hk, [b])
        ref = harvest.metric_trace(P, dks).mean()
        got = fastops.trg_from_u(fastops.band_vectors(Hk, b), dks)
        print("  %-6s band %d  shipped %.12f  fast %.12f  diff %.2e"
              % (tag, b, ref, got, abs(ref - got)))

print("\n=== 3. SHIFTED <tr g>, v4 commutator form (the Powell objective) ===")
R, H, d = harvest.read_hr("cb_M0.9_hr.dat")
Hk, dks = harvest.hk_grid(R, H, d, 60, dim=2)
P, _ = harvest.projector(Hk, [0])
u = fastops.band_vectors(Hk, 0)
norb, dim = P.shape[-1], 2
rng = np.random.default_rng(1)
worst = 0.0
for trial in range(6):
    s = rng.uniform(-1.5, 1.5, (norb, dim))
    s[0] = 0
    ref = harvest._trg_shifted(P, dks, None, s)   # kgrids ignored in v4
    got = fastops.trg_from_u(u, dks, s)
    worst = max(worst, abs(ref - got))
    print("  shift %d  shipped %.12f  fast %.12f  diff %.2e"
          % (trial, ref, got, abs(ref - got)))
print("  worst diff over 6 random shifts: %.2e" % worst)

print("\n=== 4. real 3D material, 60 orbitals, random shifts ===")
R, H, d = harvest.read_hr("hr/JVASP-13733_hr.dat")
Hk, dks = fastops.hk_grid_fast(R, H, d, 10, dim=3)
P, _ = harvest.projector(Hk, [0])
u = fastops.band_vectors(Hk, 0)
norb = P.shape[-1]
for trial in range(3):
    s = rng.uniform(-1.0, 1.0, (norb, 3))
    s[0] = 0
    ref = harvest._trg_shifted(P, dks, None, s)
    got = fastops.trg_from_u(u, dks, s)
    print("  shift %d  shipped %.12f  fast %.12f  diff %.2e"
          % (trial, ref, got, abs(ref - got)))

print("\n=== 5. minimal_metric: shipped vs fast ===")
R, H, d = harvest.read_hr("cb_M0.9_hr.dat")
Hk, dks = harvest.hk_grid(R, H, d, 60, dim=2)
a = harvest.minimal_metric(Hk, dks, 0, 60, 2, restarts=2, seed=0)
b = fastops.minimal_metric_fast(Hk, dks, 0, 2, restarts=2, seed=0)
print("  shipped  M_naive %.9f  M_min %.9f" % (a["M_naive"], a["M_min"]))
print("  fast     M_naive %.9f  M_min %.9f" % (b["M_naive"], b["M_min"]))
print("  diff     M_naive %.2e  M_min %.2e"
      % (abs(a["M_naive"] - b["M_naive"]), abs(a["M_min"] - b["M_min"])))

print("\n=== 6. the two properties v4 claims (checked on fastops) ===")


def hk_offset(R, H, deg, nk, dim, off):
    ks = [2 * np.pi * (np.arange(nk) + off) / nk for _ in range(dim)]
    mesh = np.meshgrid(*ks, indexing="ij")
    kf = np.stack([m.ravel() for m in mesh], 1)
    Hf = (H / deg[:, None, None]).reshape(H.shape[0], -1)
    out = (np.exp(1j * (kf @ R[:, :dim].T.astype(float))) @ Hf)
    return out.reshape(mesh[0].shape + H.shape[1:]), [2 * np.pi / nk] * dim


print("  offset-independence of M_naive (checkerboard M0.9, nk=120):")
for off in (0.0, 0.25, 0.5):
    Hk2, dks2 = hk_offset(R, H, d, 120, 2, off)
    v = fastops.trg_from_u(fastops.band_vectors(Hk2, 0), dks2)
    print("     offset %.2f -> M_naive %.6f" % (off, v))

# decoupled dimers: A in cell i binds B in cell i-1; shift=1 must give exactly 0
n = 2
Rd = np.array([[0, 0, 0], [-1, 0, 0]])
Hd = np.zeros((2, 2, 2), complex)
Hd[1, 0, 1] = 1.0        # A(0) - B(-1)
Hd[0, 1, 0] = 0.0
Hd = np.array([np.array([[0, 0], [0, 0]], complex),
               np.array([[0, 1.0], [0, 0]], complex)])
Hfull = np.concatenate([Hd, np.conj(np.transpose(Hd[::-1], (0, 2, 1)))])
Rfull = np.concatenate([Rd, -Rd[::-1]])
Hk3, dks3 = fastops.hk_grid_fast(Rfull, Hfull, np.ones(len(Rfull)), 200, dim=1)
u3 = fastops.band_vectors(Hk3, 0)
print("  decoupled dimers, <tr g> vs assumed position of orbital B:")
for s in (0.0, 0.5, 1.0):
    sh = np.array([[0.0], [s]])
    print("     B at %.2f a -> %.9f" % (s, fastops.trg_from_u(u3, dks3, sh)))
