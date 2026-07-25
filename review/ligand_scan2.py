"""Where does the j=1/2 doublet actually become isolated, and what is its metric
there? Small-t (molecular / weakly-coupled-octahedra) branch of the ligand scan.
"""
import numpy as np
import jeff

KB = 0.086173e-3
LAM = 0.417
COEF3 = 0.55 * 0.67
NK = 32

print(f"lam_SOC = {LAM} eV fixed, nk = {NK}, rank-2 projector\n")
print(f"{'t (eV)':>8}{'lam/t':>8}{'W=8t':>8}{'iso (eV)':>10}{'M_trg':>9}"
      f"{'U<=iso/2':>10}{'Tc_geo(K)':>11}{'Tc_pair(K)':>11}")
print('-' * 76)
for t in [0.01, 0.02, 0.03, 0.05, 0.078, 0.10, 0.15, 0.20, 0.26]:
    out, _ = jeff.analyse(t=t, lam=LAM, tp=0.0, nk=NK, bands=(4, 5))
    d = out[0]
    iso = d['iso']
    Ucap = max(iso, 0.0) / 2          # the project's own projection rule
    tc_geo = COEF3 * Ucap * d['M_trg'] / KB
    tc_pair = Ucap * d['lam_pair'] / 4 / KB
    print(f"{t:>8.3f}{LAM/t:>8.2f}{8*t:>8.3f}{iso:>10.3f}{d['M_trg']:>9.4f}"
          f"{Ucap:>10.3f}{tc_geo:>11.0f}{tc_pair:>11.0f}")

print("\nM_trg threshold for 3D is 0.339. Isolation (iso>0) and M_trg>0.339")
print("are satisfied on disjoint ranges of t -- that is D2, in one table.")
