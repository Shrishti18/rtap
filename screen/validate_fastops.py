#!/usr/bin/env python3
"""Prove fastops reproduces the shipped harvest routines exactly."""
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

print("\n=== 2. <tr g>: shipped metric_trace vs rank-1 identity ===")
for tag in ["M0.0", "M0.9"]:
    R, H, d = harvest.read_hr("cb_%s_hr.dat" % tag)
    Hk, dks = harvest.hk_grid(R, H, d, 120, dim=2)
    for b in (0, 1):
        P, _ = harvest.projector(Hk, [b])
        ref = harvest.metric_trace(P, dks).mean()
        u = fastops.band_vectors(Hk, b)
        got = fastops.trg_from_u(u, dks)
        print("  %-6s band %d   shipped %.12f   fast %.12f   diff %.2e"
              % (tag, b, ref, got, abs(ref - got)))

print("\n=== 3. shifted <tr g> (the minimal_metric objective) ===")
R, H, d = harvest.read_hr("cb_M0.9_hr.dat")
Hk, dks = harvest.hk_grid(R, H, d, 60, dim=2)
P, _ = harvest.projector(Hk, [0])
u = fastops.band_vectors(Hk, 0)
nk, norb, dim = 60, P.shape[-1], 2
ks = [2 * np.pi * (np.arange(nk) + 0.5) / nk for _ in range(dim)]
kg = np.meshgrid(*ks, indexing="ij")
rng = np.random.default_rng(1)
for trial in range(4):
    s = rng.uniform(-1.5, 1.5, (norb, dim))
    s[0] = 0
    ref = harvest._trg_shifted(P, dks, kg, s)
    got = fastops.trg_from_u(u, dks, s)
    print("  shift %d   shipped %.12f   fast %.12f   diff %.2e"
          % (trial, ref, got, abs(ref - got)))

print("\n=== 4. minimal_metric: shipped vs fast (small case) ===")
a = harvest.minimal_metric(Hk, dks, 0, nk, dim, restarts=1, seed=0)
b = fastops.minimal_metric_fast(Hk, dks, 0, dim, restarts=1, seed=0)
print("  shipped  M_naive %.9f  M_min %.9f" % (a["M_naive"], a["M_min"]))
print("  fast     M_naive %.9f  M_min %.9f" % (b["M_naive"], b["M_min"]))
print("  diff     M_naive %.2e  M_min %.2e"
      % (abs(a["M_naive"] - b["M_naive"]), abs(a["M_min"] - b["M_min"])))

print("\n=== 5. timing (real material, nk=32, 60 orbitals) ===")
import time
R, H, d = harvest.read_hr("hr/JVASP-13733_hr.dat")
t = time.time(); Hk, dks = fastops.hk_grid_fast(R, H, d, 32, dim=3)
print("  hk_grid_fast      %6.1f s  (%.2f GB)" % (time.time() - t, Hk.nbytes / 1e9))
t = time.time(); u = fastops.band_vectors(Hk, 0)
print("  band_vectors      %6.1f s" % (time.time() - t))
t = time.time(); v = fastops.trg_from_u(u, dks)
print("  trg_from_u        %6.3f s  -> Powell eval cost" % (time.time() - t))
