"""
Ligand-series test of the C8 inversion scenario.

QUESTION (user): is cRPA U for Ir4+ in a telluride or iodide host below 1.25 eV?
If yes, C8's window (lam_SOC > U_rep/3) opens.

The literature answer is developed in LIGAND_U.md. This script tests what the
project's OWN model says once U is granted at the telluride value, because the
same ligand polarizability that lowers U raises the hopping t, and D2 showed
isolation and metric trade off against each other along exactly that axis.

lam_SOC(Ir4+) = 0.417 eV is fixed by the ion (RIXS) and does not move with the
ligand. t does. So the ligand series is a scan in lam/t at fixed lam.
"""
import numpy as np
import jeff

KB = 0.086173e-3          # eV/K
LAM = 0.417               # RIXS, Sr2IrO4
COEF3 = 0.55 * 0.67       # 3D prefactor in Tc <= coef * R * U * M_trg
NK = 32

# t2g nearest-neighbour hopping by ligand. Oxide anchor is the Sr2IrO4 j=1/2
# band (W ~ 2 eV, W = 8t). Sulfide/telluride scaled by the p-d covalency trend:
# the RuX3 cRPA paper (npj Quant Mater 7, 75) gives NN |t| roughly doubling
# Cl -> I; the MX2 cRPA paper (PRB 109, 155107) shows the same direction via
# Wannier spread. Telluride t is bracketed rather than claimed to one figure.
LIGANDS = [('oxide   (O 2p)',    0.26),
           ('sulfide (S 3p)',    0.36),
           ('selenide(Se 4p)',   0.44),
           ('telluride(Te 5p)',  0.52),
           ('telluride, hi-end', 0.65)]

# U from the ligand series, anchored on cRPA Sr2IrO4 = 1.9 eV and reduced by the
# measured chalcogen factors (see LIGAND_U.md). Deliberately GENEROUS to the
# inversion scenario: the full factor-2 reduction is applied by Te.
U_LIG = {'oxide   (O 2p)': 1.90, 'sulfide (S 3p)': 1.35,
         'selenide(Se 4p)': 1.15, 'telluride(Te 5p)': 0.95,
         'telluride, hi-end': 0.80}

print(f"lam_SOC = {LAM} eV (fixed by the ion), nk = {NK}, rank-2 projector\n")
print(f"{'ligand':<20}{'t':>6}{'lam/t':>7}{'W=8t':>7}{'iso':>8}{'M_trg':>8}"
      f"{'U_cRPA':>8}{'U<=3lam':>9}{'U<=3lam/4':>11}{'Tc_geo(K)':>11}")
print('-' * 96)
rows = []
for name, t in LIGANDS:
    out, gapm = jeff.analyse(t=t, lam=LAM, tp=0.0, nk=NK, bands=(4, 5))
    d = out[0]
    U = U_LIG[name]
    # C8 as stated (loose):  window open iff U < 3 lam
    ok_c8 = U < 3 * LAM
    # the project's OWN projection rule U <= Delta_iso/2 with Delta_iso = 3lam/2
    ok_own = U <= 3 * LAM / 4
    # geometric cap, granting U in full and R = 1 (no mediator penalty)
    tc_geo = COEF3 * U * d['M_trg'] / KB
    rows.append((name, t, d, U, ok_c8, ok_own, tc_geo))
    print(f"{name:<20}{t:>6.2f}{LAM/t:>7.2f}{8*t:>7.2f}{d['iso']:>8.3f}"
          f"{d['M_trg']:>8.4f}{U:>8.2f}{str(ok_c8):>9}{str(ok_own):>11}"
          f"{tc_geo:>11.0f}")

print("\nlam_pair and nphi are ligand-independent (orbitally isotropic doublet):")
print(f"   lam_pair = {rows[0][2]['lam_pair']:.4f}   nphi = {rows[0][2]['nphi']:.4f}")

print("\nD2 isolation criterion 3 lam/2 > W = 8t  ->  t <= 3*lam/16 = "
      f"{3*LAM/16*1000:.0f} meV  (lam/t >= 5.33)")
print("   none of the ligands above reaches it; the oxide misses by "
      f"{0.26/(3*LAM/16):.1f}x, the telluride by {0.52/(3*LAM/16):.1f}x")

print("\nBoth caps, at the telluride point, granting the U inversion in full:")
name, t, d, U, _, _, tc_geo = rows[3]
tc_lam = U * d['lam_pair'] / 4 / KB
print(f"   Tc <= U*lam/4        = {tc_lam:8.0f} K   (pairing cap, U={U} lam=1/3)")
print(f"   Tc <= coef*R*U*M_trg = {tc_geo:8.0f} K   (geometric cap, R=1, M={d['M_trg']:.4f})")
print(f"   binding cap = {min(tc_lam, tc_geo):.0f} K")
