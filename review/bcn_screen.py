"""
B2C5N (ordered diamondoid B8C20N4) -- the cheap kills, in the order specified.

The candidate's own procedure says: check the electron count, check the
structure is possible, find the Fermi-level orbitals, CALCULATE THE PAIRING
SCALE REQUIRED, and reject at the first hard contradiction. Steps 2-3 pass.
Step 5 is where it ends, and it ends against published ab initio results on
the two nearest neighbours of this exact composition -- one of which has been
synthesised.

LITERATURE ANCHORS (not my calculations):
  Calandra & Mauri, PRL 101, 016401 (2008): diamond-like BC5, SYNTHESISED at
    high pressure and recoverable, lam = 0.89, omega_log = 67.4 meV, and
    "superconductivity is mostly sustained by concerted vibrations of the B
    atom and its C neighbors".
  Smith, Vohra & Chen, arXiv:2410.02104 / PRB 113, 214509: first-principles
    electron-phonon for the superhard B-C-N metals B2C3N and B4C5N3, Table 1:

     material   N_F     lam    omega_log(meV)  Tc_AllenDynes  Tc_anisoEliashberg
     BC5       0.133   0.764       70.09           34.6 K          41.2 K
     B2C3N     0.146   0.693       74.32           29.3 K          38.2 K
     B4C5N3    0.081   0.570       74.98           16.9 K          22.9 K
     MgB2      0.232   0.694       56.97           22.5 K          40.7 K
   (mu* = 0.1 throughout)

The acceptance window asked for lam >= 2.7 and omega_log >= 185 meV.
"""
import math
from scipy.optimize import brentq

KB = 8.617333e-5


def allen_dynes(lam, wlog_meV, mu):
    if lam <= mu * (1 + 0.62 * lam):
        return 0.0
    e = -1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))
    return (wlog_meV / 1.2) * math.exp(e) / 1000.0 / KB


ANCH = [("BC5", 0.133, 0.764, 70.09, 34.6, 41.2, 1 / 6),
        ("B2C3N", 0.146, 0.693, 74.32, 29.3, 38.2, 1 / 6),
        ("B4C5N3", 0.081, 0.570, 74.98, 16.9, 22.9, 1 / 12),
        ("MgB2", 0.232, 0.694, 56.97, 22.5, 40.7, float('nan'))]

print(__doc__)
print("=" * 78)

# ---------------------------------------------------------------------------
print("STEP 2. ELECTRON COUNTING -- passes, exactly as claimed")
for name, nB, nC, nN in (("B2C5N", 2, 5, 1), ("B8C20N4", 8, 20, 4),
                         ("B2C3N", 2, 3, 1), ("B4C5N3", 4, 5, 3),
                         ("BC5", 1, 5, 0)):
    n = nB + nC + nN
    ve = 3 * nB + 4 * nC + 5 * nN
    print(f"   {name:<9} {n:>3} atoms, {ve:>4} valence e- vs {4*n:>4} for a"
          f" filled sp3 net -> {4*n-ve} holes, {(4*n-ve)/n:.4f} holes/atom")

print("""
STEP 3. STRUCTURE -- passes, but note the constraints are automatic.
   Diamond is BIPARTITE: every site on one sublattice bonds only to the other.
   So "all B and N on one sublattice" ALREADY implies no B-B bonds, no B-N
   bonds, and B/N bonded only to C. Those three stated rules add nothing.
   Counting: a 32-atom cell has 16 sites per sublattice, and B8+N4 = 12, so
   4 carbons must also sit on the B/N sublattice. The rules do not say where.
   No contradiction -- but the ordering is underspecified, not determined.""")

# ---------------------------------------------------------------------------
print("STEP 5. THE PAIRING SCALE. This is where it ends.")
print(f"   {'material':<9}{'holes/atom':>12}{'N_F':>8}{'lam':>8}"
      f"{'w_log(meV)':>12}{'Tc_AD':>9}{'Tc_aME':>9}{'my AD check':>13}")
for name, nf, lam, wl, tad, tame, hpa in ANCH:
    print(f"   {name:<9}{hpa:>12.4f}{nf:>8.3f}{lam:>8.3f}{wl:>12.2f}"
          f"{tad:>8.1f}K{tame:>8.1f}K{allen_dynes(lam, wl, 0.1):>12.1f}K")
print("   (my Allen-Dynes reproduces their Tc_AD column, so the anchors are"
      " being read correctly)")

print("""
   B2C5N IS B2C3N DILUTED. Both have B:N = 2:1, so both carry exactly one hole
   per formula unit -- but B2C5N spreads it over 8 atoms instead of 6. The
   paper names the lever explicitly: "a larger difference in the B-to-N ratio
   will result in stronger hole doping. A higher N_F can lead to stronger
   electron-phonon coupling." Adding two carbons per formula unit moves the
   candidate the WRONG WAY along that lever.""")
# N_F is close to linear in hole density across the two B-C-N points
h1, f1 = 1 / 12, 0.081
h2, f2 = 1 / 6, 0.146
hx = 1 / 8
nfx = f1 + (f2 - f1) * (hx - h1) / (h2 - h1)
lx = 0.570 + (0.693 - 0.570) * (nfx - f1) / (f2 - f1)
wlx = 74.5
print(f"   linear in hole density between B4C5N3 and B2C3N:")
print(f"      N_F(B2C5N) ~ {nfx:.4f}     lam(B2C5N) ~ {lx:.3f}"
      f"     omega_log ~ {wlx:.1f} meV")
tad = allen_dynes(lx, wlx, 0.1)
print(f"      Tc(Allen-Dynes)          ~ {tad:>6.1f} K")
print(f"      Tc(x1.30, the aME/AD ratio in their table) ~ {1.30*tad:>5.1f} K")
print(f"      TARGET                     {300:>6} K"
      f"   -> SHORT BY {300/(1.30*tad):.0f}x")

# ---------------------------------------------------------------------------
print("\nWHY THE PREMISE FAILS: light atoms give a high DEBYE temperature,")
print("not a high omega_log.")
print("""   Diamond's optical mode is ~165 meV and these are superhard networks
   with Debye temperatures near 2000 K -- yet omega_log across the whole
   family is 70.09, 74.32, 74.98 meV. It barely moves while the composition
   changes a lot. The reason is in Calandra & Mauri: the coupling is carried
   by the SOFT concerted B-centred vibrations, not the hard C-C stretches.
   omega_log is the coupling-weighted average, so stiffening the lattice
   raises the modes that do not pair and leaves omega_log where it was.
   The anticorrelation is visible in the anchor table itself:""")
print(f"   {'material':<9}{'w_log':>9}{'N_F':>8}{'lam':>8}"
      f"{'lam*w_log^2 (Hopfield eta/M)':>32}")
for name, nf, lam, wl, *_ in ANCH:
    print(f"   {name:<9}{wl:>9.2f}{nf:>8.3f}{lam:>8.3f}{lam*wl**2:>32.0f}")
print("   -> the highest omega_log in the family (B4C5N3) has the LOWEST N_F")
print("      and the LOWEST lambda. They are not independent knobs.")

eta_bc5 = 0.764 * 70.09 ** 2
need = 2.7 * 185.0 ** 2
print(f"\n   To hold lam = 2.7 AT omega_log = 185 meV, the Hopfield parameter")
print(f"   eta = lam*M*<w^2> must rise from {eta_bc5:.0f} (BC5) to "
      f"{need:.0f}:  {need/eta_bc5:.1f}x")
print("   while the same acceptance window demands W_active >= 5 eV, which")
print("   pushes N(0) -- and therefore eta -- the other way.")

# ---------------------------------------------------------------------------
print("\nTHE FAMILY-LEVEL RESULT: omega_log, NOT lambda, IS THE BINDING WALL")
print("""   Allen-Dynes has a hard ceiling as lam -> infinity:
       Tc -> (omega_log/1.2) * exp[-1.04/(1 - 0.62*mu*)]
   because the exponent saturates. So for a given omega_log there is a maximum
   Tc that NO coupling strength can beat.""")
print(f"   {'mu*':>6}{'Tc_max/omega_log':>20}{'w_log for 300 K':>18}"
      f"{'w_log for 350 K':>18}")
for mu in (0.08, 0.10, 0.12):
    c = math.exp(-1.04 / (1 - 0.62 * mu)) / 1.2
    print(f"   {mu:>6.2f}{c:>20.4f}{300*KB*1000/c:>16.1f} meV"
          f"{350*KB*1000/c:>16.1f} meV")
print(f"""
   The whole superhard B-C-N family sits at omega_log = 70-75 meV. At
   mu* = 0.08 that caps Tc at {allen_dynes(1e6, 74.5, 0.08):.0f} K AT INFINITE COUPLING.
   300 K is not merely hard for this family -- it is UNREACHABLE at any
   lambda, by the candidate's own screening formula, unless omega_log roughly
   doubles. And omega_log is the one quantity that barely moved across BC5,
   B2C3N and B4C5N3.""")

for tgt in (300, 350):
    f = brentq(lambda w: allen_dynes(1e6, w, 0.08) - tgt, 10, 500)
    print(f"   minimum omega_log for {tgt} K at infinite lambda: {f:.1f} meV"
          f"   (family has 70-75)")

# ---------------------------------------------------------------------------
print("\nAND lam = 2.7 AT 1 BAR HAS NEVER BEEN OBSERVED")
print(f"   {'material':<22}{'lambda':>9}{'pressure':>14}")
for m, l, p in (("Pb", 1.55, "ambient"), ("Nb3Ge", 1.7, "ambient"),
                ("MgB2", 0.69, "ambient"), ("BC5 (superhard)", 0.76, "ambient"),
                ("H3S", 2.0, "155 GPa"), ("LaH10", 2.2, "170 GPa")):
    print(f"   {m:<22}{l:>9.2f}{p:>14}")
print("""   Coupling of order 2-3 appears only in the hydrides, and the megabar
   pressure is precisely what holds the lattice against the softening that
   such coupling implies. This design forbids pressure by construction, so it
   asks for hydride-scale coupling in a diamond lattice at 1 bar.""")

print("""
================================================================================
VERDICT ON B2C5N
================================================================================
REJECT, at step 5 of its own procedure, on published ab initio results for its
own composition class rather than on any argument of mine.

  Passes: exact intrinsic hole count (1 per 8 atoms, no dopant); a bipartite
  diamond net that makes the B/B, B/N and coordination rules automatic; no
  mixed valence; no mobile anions; no invented negative U; no flat band.
  Every previously fatal pattern is genuinely avoided. It is a much better
  guess than the last three.

  Fails on the number it was built to deliver:
    asked for  lam >= 2.7,  omega_log >= 185 meV,  Tc ~ 380 K
    the class gives lam = 0.57-0.76, omega_log = 70-75 meV, Tc = 23-41 K
    and B2C5N is B2C3N (38 K) diluted with two extra carbons per formula
    unit, moving it DOWN the paper's own N_F lever to ~23-30 K.

  The premise itself is wrong, and this is the transferable part: light atoms
  buy a high Debye temperature, not a high omega_log. omega_log is weighted by
  what couples, and what couples in a covalent hole-doped network are the soft
  B-centred modes. Across BC5, B2C3N and B4C5N3 the composition changes
  drastically and omega_log moves by 7%.

  FAMILY-LEVEL, not compound-level: at omega_log = 74.5 meV the Allen-Dynes
  ceiling at infinite coupling is ~241 K. No sp3 B-C-N metal reaches 300 K for
  ANY lambda. The search should leave this family rather than reorder it.

  WHAT THE NEXT GUESS MUST DO DIFFERENTLY. The binding constraint is now
  named and it is omega_log >= ~93 meV at mu* = 0.08 as a floor, with enough
  N(0) left over to give lam >~ 1.5-2 at that frequency. That means the
  COUPLING modes must be the light ones -- hydrogen or a light-atom mode that
  sits at the Fermi-level orbital, not a heavy-atom mode decorated by light
  neighbours. Which is exactly why the hydrides work and why they need
  pressure. Any ambient-pressure route has to find that combination without
  the pressure, and no amount of network stiffness substitutes for it.
""")
