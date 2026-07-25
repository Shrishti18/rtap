import numpy as np
KB=0.086173e-3  # eV/K
# EXACT sum rule:  M = sum_m |v_nm|^2/(E_n-E_m)^2  <=  <v^2>/Delta^2
# Constraints:     Tc = min(U/(4 n_phi), coef*U*M),   U <= Delta/2
# => Tc <= min( Delta/(8 n_phi) , coef*<v^2>/(2 Delta) )
# The two branches cross at Delta* = 2 sqrt(coef * n_phi * <v^2>), giving
#        Tc_max = sqrt(coef*<v^2>/n_phi)/4
# LARGE METRIC NEEDS SMALL GAP; PROJECTION NEEDS LARGE GAP. The optimum is
# where they meet, and it is a CEILING no geometry can exceed.
coef3=0.55*0.67; coef2=0.2225*0.67
def ceiling(vrms,nphi,coef): return np.sqrt(coef*vrms**2/nphi)/4
print("Tc_max = sqrt(coef*<v^2>/n_phi)/4   [v_rms in eV, lattice units a=1]\n")
print(f"{'n_phi':>6} {'3D: Tc/v_rms':>13} {'2D: Tc/v_rms':>13}")
for n in [1,2,3,4,8]:
    print(f"{n:>6} {ceiling(1,n,coef3):13.4f} {ceiling(1,n,coef2):13.4f}")
print("\nn_phi=1 is excluded (g=0 identically), so n_phi=2 is the optimum: Tc <= 0.107 v_rms (3D)")
print()
# convert v_rms to a hopping scale: 1D cosine band dH/dk = -2t sin k -> <v^2> = 2t^2
print("For a cosine band, <v^2> = 2 t^2 => v_rms = 1.414 t, so")
for d,c,lab in [(3,coef3,'3D'),(2,coef2,'2D')]:
    r=ceiling(np.sqrt(2),2,c)
    print(f"   {lab}: Tc <= {r:.3f} t     -> 300 K needs t >= {25.85e-3/r*1000:.0f} meV")
print()
print("CROSS-DOMAIN CHECK -- attractive Hubbard model, QMC:")
print("   Tc/t peaks at ~0.10-0.15 near U/t ~ 5-8 (2D and 3D).")
print("   My bound, derived from quantum geometry + dilution + the gap constraint:")
print("   Tc/t <= 0.152 (3D).  SAME NUMBER, INDEPENDENT DERIVATION.")
print("   => geometry REACHES the attractive-Hubbard ceiling; it does not raise it.")
print()
print("APPLY TO BKBO (the one material with a measured U_eff of the right size):")
U=1.9; W=1.6; t=W/12.0   # 3D simple-cubic-like: W = 12t
print(f"   U_eff ~ {U} eV, W ~ {W} eV -> t ~ {t:.3f} eV, U/t ~ {U/t:.1f}")
print(f"   optimum is U/t ~ 6, so BKBO is ~{U/t/6:.1f}x TOO STRONGLY COUPLED")
print(f"   ceiling at its own t: Tc <= {0.152*t*1000:.1f} meV = {0.152*t/KB:.0f} K")
print(f"   measured Tc = 30 K -> {30/(0.152*t/KB)*100:.0f}% of its own ceiling")
print()
print("   and BKBO's carrier is Bi 6s: n_phi = 1 -> tr g = 0 IDENTICALLY.")
print("   It receives NO geometric enhancement at all. That is the gap in the")
print("   literature: the one family with the right U has the wrong orbital count.")
