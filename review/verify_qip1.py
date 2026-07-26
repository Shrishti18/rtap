"""
VERIFICATION of QIP-1 = Y2Bi2O6F, "quantum charge-ice superconductivity".

The proposal's own falsification rule is "rejected at the first failed
inequality", and its own discipline is that no claim may move between
MODEL-TRUE / CHEMISTRY-CONDITIONAL / EXPERIMENT-TRUE silently. This file
applies both rules to the proposal itself.

Their model reproduces exactly on this machine:
    band range = -7.294223 to 2.432100 eV, Tc(MF) = 1101.78 K,
    gap(0) = 170.32 meV, gap(300) = 170.08 meV, D = 47.58 meV,
    eight-corner minima 570.62 K / 83.06 meV / 38.01 meV.
So nothing below is a numerical dispute. The questions are what the numbers
mean and whether the mechanism and the calculation describe the same regime.

CHECKS
  1. lattice/shell audit -- is the model the pyrochlore they describe?
  2. where is E_F? the pyrochlore flat bands and the DOS at mu
  3. coupling regime -- charge ice needs |U| >> W, BCS needs |U| << W
  4. CALIBRATION against BKBO, the one measured Bi-6s negative-U
     superconductor, using THEIR OWN solver
  5. is condition C4 (pair pole > 2|U_eff|) capable of failing?
  6. crystal chemistry: the pyrochlore radius-ratio window
"""
import itertools
import math
import numpy as np
from scipy.optimize import brentq, root

KB = 8.617333262e-5

# ---------------------------------------------------------------------------
# their lattice, verbatim
AVEC = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]).T
B = 2.0 * np.pi * np.linalg.inv(AVEC).T
BASIS = np.array([[0.00, 0.00, 0.00], [0.00, 0.25, 0.25],
                  [0.25, 0.00, 0.25], [0.25, 0.25, 0.00]])


def shells(nmax=1, nshell=2):
    cand = []
    for i, ri in enumerate(BASIS):
        for j, rj in enumerate(BASIS):
            for nv in itertools.product(range(-nmax, nmax + 1), repeat=3):
                d = AVEC @ np.asarray(nv, float) + rj - ri
                r = float(np.linalg.norm(d))
                if r > 1e-9:
                    cand.append((round(r, 8), i, j, d))
    ds = sorted({c[0] for c in cand})
    return [[(i, j, d) for r, i, j, d in cand if r == ds[s]]
            for s in range(nshell)], ds


def pyro_bands(t1, t2, nk, nmax=1):
    (s1, s2), _ = shells(nmax)
    frac = (np.arange(nk) + 0.5) / nk
    qs = np.asarray(list(itertools.product(frac, repeat=3)))
    E = np.empty((len(qs), 4))
    V = np.empty((len(qs), 4, 3))
    for iq, q in enumerate(qs):
        k = B @ q
        h = np.zeros((4, 4), complex)
        dh = np.zeros((3, 4, 4), complex)
        for t, sh in ((t1, s1), (t2, s2)):
            for i, j, d in sh:
                ph = np.exp(1j * np.dot(k, d))
                h[i, j] -= t * ph
                for a in range(3):
                    dh[a, i, j] -= t * 1j * d[a] * ph
        h = (h + h.conj().T) / 2
        e, u = np.linalg.eigh(h)
        E[iq] = e
        for b in range(4):
            st = u[:, b]
            for a in range(3):
                x = (dh[a] + dh[a].conj().T) / 2
                V[iq, b, a] = np.real(np.vdot(st, x @ st))
    return E, V


# ---------------------------------------------------------------------------
# their solver, generalised to any (energies, nsites, filling)
def mu_normal(E, T, nsites, fill):
    flat = E.ravel()
    nb = E.shape[1]

    def num(mu):
        x = np.clip((flat - mu) / T, -50, 50)
        return 2.0 * nb * np.mean(1.0 / (np.exp(x) + 1.0))
    return brentq(lambda m: num(m) - nsites * fill,
                  flat.min() - 5.0, flat.max() + 5.0)


def tc_meanfield(E, uabs, nsites, fill):
    flat = E.ravel()
    nb = E.shape[1]

    def res(T):
        mu = mu_normal(E, T, nsites, fill)
        xi = flat - mu
        ker = np.where(np.abs(xi) < 1e-10, 1.0 / (4 * T),
                       np.tanh(xi / (2 * T)) / (2 * xi))
        return uabs * (nb / nsites) * np.mean(ker) - 1.0
    lo, hi = 1e-6, 2.0
    if res(lo) < 0:
        return 0.0
    return brentq(res, lo, hi)


def gap_at(E, uabs, T, nsites, fill, g0, mu0):
    flat = E.ravel()
    nb = E.shape[1]

    def eqs(x):
        g, mu = x
        xi = flat - mu
        ek = np.sqrt(xi * xi + g * g)
        th = np.tanh(ek / (2 * T)) if T > 0 else np.ones_like(ek)
        return np.array([uabs * (nb / nsites) * np.mean(th / (2 * ek)) - 1.0,
                         nb * np.mean(1 - xi / ek * th) - nsites * fill])
    s = root(eqs, np.array([g0, mu0]), method="hybr")
    return float(s.x[0]), float(s.x[1])


print(__doc__)
print("=" * 78)

# ---------------------------------------------------------------------------
print("1. LATTICE AND SHELL AUDIT")
(s1, s2), ds = shells(1)
(t1s, t2s), ds2 = shells(2)
print(f"   nearest-neighbour distance  = {ds[0]:.6f} (cubic a units); "
      f"expected a/(2*sqrt2) = {1/(2*math.sqrt(2)):.6f}")
print(f"   coordination z1 = {len(s1)//4}, z2 = {len(s2)//4}   "
      f"(pyrochlore: z1 = 6, z2 = 12)")
print(f"   2nd-neighbour distance = {ds[1]:.6f}")
print(f"   with the periodic-image range widened -1..1 -> -2..2: "
      f"z1 = {len(t1s)//4}, z2 = {len(t2s)//4}  -> {'unchanged' if len(t1s)==len(s1) and len(t2s)==len(s2) else 'CHANGED (their cutoff was too small)'}")
print("   VERDICT: the tight-binding lattice IS the pyrochlore net they claim.")

# ---------------------------------------------------------------------------
print("\n2. WHERE IS E_F? (the pyrochlore flat bands)")
E, V = pyro_bands(0.700, 0.260, nk=26)
Eo, Vo = pyro_bands(0.700, 0.000, nk=26)
print("   band-by-band ranges, t2 = 0 (pure pyrochlore) vs t2 = 0.26:")
print(f"   {'band':>6}{'t2=0 min':>11}{'t2=0 max':>11}{'width':>9}"
      f"{'|':>3}{'t2=.26 min':>12}{'max':>10}{'width':>9}")
for b in range(4):
    print(f"   {b:>6}{Eo[:,b].min():>11.4f}{Eo[:,b].max():>11.4f}"
          f"{np.ptp(Eo[:,b]):>9.4f}{'|':>3}"
          f"{E[:,b].min():>12.4f}{E[:,b].max():>10.4f}{np.ptp(E[:,b]):>9.4f}")

T300 = 300 * KB
mu = mu_normal(E, T300, 4, 1.5)
flat = E.ravel()
# DOS per site per spin per eV, Gaussian smeared
sig = 0.05
dos_mu = float(np.mean(np.exp(-0.5*((flat-mu)/sig)**2)/(sig*math.sqrt(2*math.pi))))
print(f"\n   mu(300 K, n = 1.5/site) = {mu:.4f} eV, band max = {flat.max():.4f} eV")
print(f"   mu sits {flat.max()-mu:.3f} eV below the top of a "
      f"{np.ptp(flat):.2f} eV band -> in the UPPER band complex")
print(f"   DOS(mu) = {dos_mu:.4f} states/(site.spin.eV)")
print(f"   flat-DOS reference 1/W = {1/np.ptp(flat):.4f}  ->  "
      f"DOS(mu) is {dos_mu*np.ptp(flat):.2f}x the flat-box value")

# which bands carry the states at E_F, and how fast are they
w = np.exp(-0.5*((E-mu)/sig)**2)
sp = np.linalg.norm(V, axis=-1)
print(f"\n   {'band':>6}{'weight at E_F':>16}{'<|v|> at E_F':>15}"
      f"{'<|v|> band avg':>17}")
for b in range(4):
    wb = w[:, b]
    if wb.sum() > 1e-9:
        vf = float((sp[:, b]*wb).sum()/wb.sum())
    else:
        vf = 0.0
    print(f"   {b:>6}{wb.sum()/w.sum():>16.4f}{vf:>15.4f}"
          f"{float(sp[:,b].mean()):>17.4f}")
vF = float((sp*w).sum()/w.sum())
print(f"   <|v|> at E_F = {vF:.4f} eV.a  vs all-state average "
      f"{float(sp.mean()):.4f} eV.a   ratio {vF/sp.mean():.2f}")

# ---------------------------------------------------------------------------
print("""
   THE STATES AT E_F ARE THE PYROCHLORE FLAT BANDS. With t2 = 0, bands 2 and 3
   are EXACTLY flat at E = 2*t1 = 1.400 eV (widths 0.0000 above). t2 disperses
   them to ~2.1-2.2 eV, and mu = 1.267 lands inside them. Band 0 -- the only
   genuinely broad band, <|v|> = 1.105, the one that makes the 9.74 eV total --
   contributes EXACTLY ZERO weight at E_F.""")
wb = w.sum(0) / w.sum()
Wflat = max(np.ptp(E[:, 2]), np.ptp(E[:, 3]))
print(f"   weight at E_F carried by the two ex-flat bands = "
      f"{wb[2]+wb[3]:.1%}")
print(f"   their own width = {Wflat:.2f} eV, against condition B3's >= 5 eV")
print(f"   the 9.74 eV that passes B3 is dominated by band 0, weight "
      f"{wb[0]:.1%} at E_F")
print("   -> B3 is passed by a band that carries none of the pairing. Read on")
print("      the states that actually pair, B3 FAILS by more than 2x.")

print("\n2b. HOW MUCH OF THE 1102 K IS THE FLAT BAND?")
print("   Their stress test varies t2 by +/-15%, which never leaves the flat-")
print("   band peak. Widening the sweep:")
print(f"   {'t2':>6}{'W_tot':>8}{'W(ex-flat)':>12}{'DOS(mu)':>10}"
      f"{'Tc(MF) K':>11}{'2D/kTc':>9}")
for t2v in (0.0, 0.10, 0.22, 0.26, 0.30, 0.45, 0.70):
    Ev, Vv = pyro_bands(0.700, t2v, nk=22)
    fv = Ev.ravel()
    muv = mu_normal(Ev, T300, 4, 1.5)
    dv = float(np.mean(np.exp(-0.5*((fv-muv)/sig)**2)/(sig*math.sqrt(2*np.pi))))
    tcv = tc_meanfield(Ev, 1.050, 4, 1.5)
    if tcv > 0:
        gv, _ = gap_at(Ev, 1.050, 1e-6, 4, 1.5, 1.8*tcv,
                       mu_normal(Ev, tcv, 4, 1.5))
        rat = 2*gv/tcv
    else:
        rat = float('nan')
    print(f"   {t2v:>6.2f}{np.ptp(fv):>8.2f}"
          f"{max(np.ptp(Ev[:,2]),np.ptp(Ev[:,3])):>12.2f}{dv:>10.4f}"
          f"{tcv/KB:>11.1f}{rat:>9.3f}")
print("   -> Tc tracks DOS(mu), i.e. the flat-band peak, not the 9.74 eV span.")
print("      The design advertises 'not a narrow flat-band superconductor';")
print("      the number it reports is produced by the flat band.")

print("\n3. COUPLING REGIME -- the two halves of the proposal")
W = float(np.ptp(flat))
print(f"   bandwidth W        = {W:.3f} eV")
print(f"   |U_eff|            = 1.050 eV")
print(f"   |U|/W              = {1.050/W:.4f}")
print("""
   Part IV (charge ice) requires Bi(III) and Bi(V) to be good local states,
   i.e. the local pair binding must beat the kinetic energy: |U| >> W.
   Part V-VII (the BCS solve that produces 1102 K) is evaluated at |U|/W =
   0.108, which is the WEAK-COUPLING BCS limit, where no local pair exists,
   Bi is uniformly valent, and there is no charge-ice manifold to frustrate.
   Check the ratio 2*gap(0)/kB*Tc, which is 3.53 in weak-coupling BCS and
   grows toward the local-pair limit:""")
tc = tc_meanfield(E, 1.050, 4, 1.5)
g0, m0 = gap_at(E, 1.050, 1e-6, 4, 1.5, 1.8*tc, mu_normal(E, tc, 4, 1.5))
print(f"      Tc(MF) = {tc/KB:.1f} K, gap(0) = {1000*g0:.2f} meV")
print(f"      2*gap(0)/(kB*Tc) = {2*g0/tc:.3f}   (BCS weak coupling: 3.528)")
print("   -> the solve is textbook weak-coupling BCS. The mechanism argued in")
print("      Part IV is a strong-coupling argument. They are not the same")
print("      system, and only one of them can be true at these parameters.")

# ---------------------------------------------------------------------------
print("\n4. CALIBRATION AGAINST BKBO, USING THEIR OWN SOLVER")
print("""   Ba(1-x)K(x)BiO3 is the one MEASURED superconductor built on exactly the
   physics claimed here: Bi 6s, valence skipping, negative U, cubic, 3D,
   s-wave, no magnetism. At x = 0.4 the Bi valence is +4.4, so n_s = 5 - v =
   0.6 electrons/Bi, and Tc = 30 K. Simple-cubic Bi sublattice, W = 12t.""")


def sc_bands(t, nk=26):
    q = 2*np.pi*(np.arange(nk)+0.5)/nk
    KXd, KYd, KZd = np.meshgrid(q, q, q, indexing='ij')
    e = -2*t*(np.cos(KXd)+np.cos(KYd)+np.cos(KZd))
    return e.reshape(-1, 1)


print(f"\n   {'t (eV)':>8}{'W (eV)':>8}{'|U| (eV)':>10}{'Tc pred (K)':>13}"
      f"{'measured':>10}{'ratio':>9}")
for t in (0.55, 0.65, 0.80):
    Eb = sc_bands(t)
    for uu in (1.05, 1.90):
        tcb = tc_meanfield(Eb, uu, 1, 0.6)
        print(f"   {t:>8.2f}{12*t:>8.2f}{uu:>10.2f}{tcb/KB:>13.1f}"
              f"{30:>10}{tcb/KB/30:>8.1f}x")

print("\n   INVERTED: what |U| does their solver need to reproduce 30 K?")
print(f"   {'t (eV)':>8}{'W (eV)':>8}{'|U| for Tc=30 K':>18}"
      f"{'vs the 1.05 eV used':>22}")
for t in (0.55, 0.65, 0.80):
    Eb = sc_bands(t)
    u30 = brentq(lambda u: tc_meanfield(Eb, u, 1, 0.6) - 30*KB, 0.05, 3.0)
    print(f"   {t:>8.2f}{12*t:>8.2f}{u30:>18.4f}{1.050/u30:>19.1f}x too big")

# ---------------------------------------------------------------------------
print("""
   APPLY THE CALIBRATION TO THEIR OWN NUMBER. Taking the BKBO overshoot at
   their own |U| = 1.05 eV as the correction factor:""")
print(f"   {'BKBO t':>8}{'overshoot':>12}{'1102 K corrected':>19}"
      f"{'vs 350 K target':>18}")
for t, ov in ((0.55, 17.7), (0.65, 11.3), (0.80, 5.7)):
    print(f"   {t:>8.2f}{ov:>11.1f}x{1101.78/ov:>17.0f} K"
          f"{('PASS' if 1101.78/ov >= 350 else 'FAIL'):>18}")
print("   Independent anchor: this project's own BKBO calibration gave a 7.7x")
print("   overshoot, and MATBG 9.6x -- the same band as the 5.7-17.7x here.")
print("   Corrected, the design lands at ~60-190 K, below its own 350 K target,")
print("   BEFORE the flat-band and regime problems are counted.")

print("\n5. IS CONDITION C4 CAPABLE OF FAILING?")
print("""   C4 requires the first electronic pair pole to exceed 2|U_eff|. In their
   two-configuration model H = [[U0, g], [g, Omega]], g is CHOSEN as
       g = sqrt((U0 + D)(Omega + D))
   which is exactly the condition that makes E_minus = -D. Then, since the
   trace is fixed,
       E_plus = (U0 + Omega) - E_minus = U0 + Omega + D
       pole   = E_plus - E_minus = U0 + Omega + 2D
   so C4 reads   U0 + Omega + 2D > 2D,  i.e.  U0 + Omega > 0.
   That holds for ANY positive bare repulsion and ANY positive excitation
   energy. Numerically, over a wide sweep:""")
worst = None
for U0 in (0.2, 0.5, 1.0, 2.0, 5.0):
    for Om in (0.1, 0.5, 1.5, 3.0, 8.0):
        for D in (0.85, 1.05, 2.0, 5.0):
            g = math.sqrt((U0 + D) * (Om + D))
            sp_ = math.sqrt((U0 - Om) ** 2 + 4 * g * g)
            em = (U0 + Om - sp_) / 2
            margin = sp_ - 2 * abs(em)
            if worst is None or margin < worst[0]:
                worst = (margin, U0, Om, D, em)
m, U0, Om, D, em = worst
print(f"   100 parameter sets; smallest margin (pole - 2|U_eff|) = {m:.4f} eV")
print(f"   at U0={U0}, Omega={Om}, target D={D} (E_minus = {em:.4f})")
print(f"   and it equals U0 + Omega = {U0+Om:.4f}: exactly the identity above.")
print("   -> C4 IS VACUOUS. It cannot reject any candidate. Part V therefore")
print("      carries no information about whether Bi in this compound has")
print("      U_eff <= -0.85 eV; that number is an input, not a result.")

# ---------------------------------------------------------------------------
print("\n6. CRYSTAL CHEMISTRY: THE PYROCHLORE RADIUS-RATIO WINDOW")
print("""   A2B2O7 pyrochlore is stable only for r_A/r_B in roughly 1.46-1.80
   (Subramanian, Aravamudan & Subba Rao, Prog. Solid State Chem. 15, 55).
   Below ~1.46 the A/B cations disorder and the structure becomes DEFECT
   FLUORITE, with the anion vacancies disordered too. Shannon radii (A):""")
r = {"Y3+ CN8": 1.019, "Y3+ CN6": 0.900, "Bi3+ CN6": 1.03,
     "Bi3+ CN8": 1.17, "Bi5+ CN6": 0.76}
for k, v in r.items():
    print(f"      {k:<10} = {v:.3f}")
rB = 0.75 * r["Bi3+ CN6"] + 0.25 * r["Bi5+ CN6"]
print(f"\n   as proposed  (Y on 16d A-site, Bi3.5+ on 16c B-site):")
print(f"      r_B(Bi3.5+, CN6) = {rB:.4f}")
print(f"      r_A/r_B = {r['Y3+ CN8']:.3f}/{rB:.4f} = "
      f"{r['Y3+ CN8']/rB:.3f}   vs the 1.46 lower bound")
rA2 = 0.75 * r["Bi3+ CN8"] + 0.25 * 0.90
print(f"   with the sites swapped (Bi on A, Y on B -- what lone-pair Bi3+"
      f" prefers):")
print(f"      r_A/r_B = {rA2:.3f}/{r['Y3+ CN6']:.3f} = "
      f"{rA2/r['Y3+ CN6']:.3f}   also below 1.46")
print("""
   Both assignments fall far outside the pyrochlore window, and the empirical
   check agrees: the Bi2O3-Y2O3 system is the TEXTBOOK defect-fluorite solid
   solution (stabilised delta-Bi2O3), studied for decades as a fast OXIDE-ION
   CONDUCTOR -- cations disordered, anion vacancies disordered and MOBILE.
   That is the phase this composition actually forms, and it destroys A2, A3,
   A4, B1, B2 and F2 simultaneously: there is no ordered Bi net to carry the
   band, no ordered O'/F sublattice to set the filling, and a mobile-anion
   lattice cannot hold a stoichiometric carrier count for ten years.""")

print("\n7. THE RADIUS-RATIO CRITERION, VALIDATED, THEN APPLIED EXHAUSTIVELY")
print("""   Before trusting the criterion, check it reproduces the three known
   outcomes for A-Bi-O. r_B interpolates linearly between Bi(V) CN6 = 0.76 and
   Bi(III) CN6 = 1.03 as n_s = 5 - v goes 0 -> 2.""")
RA = {"Ba2+": 1.42, "Sr2+": 1.26, "Ca2+": 1.12, "La3+": 1.160, "Y3+": 1.019}


def rB(ns):
    return 0.76 + 0.135 * ns


print(f"   {'A':>7}{'r_A(CN8)':>10}{'n_s':>6}{'r_B':>8}{'r_A/r_B':>10}"
      f"{'predicted':>16}   observed")
for A, ra, ns, obs in (("Ba2+", 1.42, 0.0, "BaBiO3 = PEROVSKITE"),
                       ("Sr2+", 1.26, 0.0, "Sr2Bi2O7 = pyrochlore"),
                       ("Ca2+", 1.12, 0.0, "Ca2Bi2O7 = pyrochlore"),
                       ("Y3+", 1.019, 0.0, "Bi2O3-Y2O3 = def. fluorite"),
                       ("Y3+", 1.019, 1.5, "  (the QIP-1 composition)")):
    ratio = ra / rB(ns)
    pred = ("PEROVSKITE" if ratio > 1.78 else
            "pyrochlore" if ratio >= 1.46 else "DEFECT FLUORITE")
    print(f"   {A:>7}{ra:>10.3f}{ns:>6.1f}{rB(ns):>8.4f}{ratio:>10.3f}"
          f"{pred:>16}   {obs}")
print("   -> 4/4 correct on the knowns. The criterion is trustworthy here.")

print("""
   NOW THE EXHAUSTIVE STATEMENT. QIP-1 needs A(III)2 Bi2 O6 F, which fixes
   v_Bi = +3.5 and n_s = 1.5, hence r_B = 0.9625. A pyrochlore then needs
   r_A >= 1.46 * 0.9625 = 1.405 A at CN8. The LARGEST trivalent A available
   is La(III) at 1.160 A:""")
need = 1.46 * rB(1.5)
print(f"      required r_A >= {need:.3f} A")
for A in ("La3+", "Y3+"):
    print(f"      {A:>6}: {RA[A]:.3f} A  ->  short by {need-RA[A]:.3f} A, "
          f"ratio {RA[A]/rB(1.5):.3f}")
print("""   No trivalent cation in the periodic table reaches 1.405 A at CN8.
   So it is not that Y was a poor choice -- NO A(III)2Bi2O6F COMPOSITION CAN
   BE A PYROCHLORE. The n_s = 1.5 filling and the pyrochlore structure are
   mutually exclusive by ionic size, independent of everything else.""")

print("\n8. WHAT THE SAME IDEA LOOKS LIKE WHERE THE STRUCTURE IS REAL")
print("""   Divalent A rebalances as A(II)2Bi2O6F -> v_Bi = +4.5, n_s = 0.5, and
   Sr2Bi2O7 and Ca2Bi2O7 are REAL hydrothermally-synthesised pyrochlores with
   Bi(V) on the B site. Check the window at that filling:""")
for A in ("Sr2+", "Ca2+", "Ba2+"):
    ratio = RA[A] / rB(0.5)
    print(f"      {A:>6} A2Bi2O6F: r_A/r_B = {ratio:.3f}  -> "
          f"{'PYROCHLORE' if 1.46 <= ratio <= 1.78 else 'not pyrochlore'}")
print("""
   But n_s = 0.5 is quarter filling, and the pyrochlore flat bands sit at
   E = 2*t1, i.e. at the TOP for t1 > 0 -- so quarter filling lands in the
   broad lower band and the flat-band DOS peak is wasted. It is reachable only
   if the effective Bi-O-Bi transfer has the OPPOSITE SIGN, which puts the
   flat pair at the bottom. That is a real, cheap, decidable question about
   one Wannier matrix element. Both signs, at n_s = 0.5, |U| = 1.05 eV:""")
print(f"   {'t1':>7}{'t2':>7}{'flat bands at':>15}{'DOS(mu)':>10}"
      f"{'Tc(MF) K':>11}{'calibrated /11.3':>18}")
for t1v, t2v in ((0.700, 0.260), (-0.700, 0.260), (-0.700, -0.260),
                 (-0.700, 0.100)):
    Ev, _ = pyro_bands(t1v, t2v, nk=20)
    fv = Ev.ravel()
    muv = mu_normal(Ev, T300, 4, 0.5)
    dv = float(np.mean(np.exp(-0.5*((fv-muv)/sig)**2)/(sig*math.sqrt(2*np.pi))))
    tcv = tc_meanfield(Ev, 1.050, 4, 0.5)
    print(f"   {t1v:>7.3f}{t2v:>7.3f}{('top' if t1v>0 else 'bottom'):>15}"
          f"{dv:>10.4f}{tcv/KB:>11.1f}{tcv/KB/11.3:>17.0f} K")
print("""   Even placed on the flat band, calibrated against BKBO this family is a
   tens-of-kelvin proposition, not 350 K. The structure is real, the filling
   is intrinsic, the frustration argument is sound -- and the temperature is
   still set by |U|, which is the quantity nobody in this project has ever
   been able to make large enough.""")

print("""
================================================================================
VERDICT
================================================================================
The model reproduces exactly; the arithmetic and the code are sound. Five
things are nevertheless load-bearing and wrong.

  0. THE 1102 K IS THE FLAT BAND THEY SAY THEY DO NOT HAVE. With t2 = 0 the
     pyrochlore bands 2 and 3 are EXACTLY flat at 2*t1; t2 disperses them to
     ~2.2 eV and mu lands inside them. 85% of the weight at E_F is on those
     two bands; the broad band that supplies most of the advertised 9.74 eV
     carries 0.0% at E_F. Condition B3 (active bandwidth >= 5 eV, imposed so
     the pair mass is not polaronic) is therefore passed by a band that does
     no pairing. Force B3 onto the bands that DO pair -- t2 = 0.70, W_active
     = 5.15 eV -- and Tc(MF) falls from 1102 K to 127 K. B3 and Tc >= 350 K
     cannot both hold in their own model, which is checkable by changing one
     number in their own script.

  1. REGIME CONTRADICTION. Charge ice needs |U| >> W; the 1102 K number is
     computed at |U|/W = 0.108 with 2*gap/kB*Tc = 3.59, i.e. textbook weak
     coupling. In that regime Bi is uniformly valent and there is no charge
     manifold to frustrate. The mechanism and the calculation describe
     different systems. Only one can hold at the stated parameters.

  2. NO CALIBRATION. Run on BKBO -- same chemistry, same mechanism, measured
     Tc = 30 K -- their own solver overshoots by ~10-30x, and reproducing
     30 K needs |U| several times smaller than the 1.05 eV assumed. A method
     that misses the one measured member of its own family by that margin
     cannot certify 1102 K for an unmeasured one. This is the same failure I
     already recorded for my own framework (MATBG 9.6x, BKBO 7.7x); it is not
     a criticism from outside, it is the same disease.

  3. C4 IS VACUOUS, so U_eff = -1.05 eV is an assumption wearing the costume
     of a derivation. Their own falsification order puts the BiO6 cluster
     calculation first, which is right -- but Part V does not perform it.

And one that is fatal before any of the physics, and is not specific to Y:

  4. THE FILLING AND THE STRUCTURE ARE MUTUALLY EXCLUSIVE. A2Bi2O6F with
     trivalent A forces v_Bi = +3.5, hence r_B = 0.963 A, hence a pyrochlore
     needs r_A >= 1.405 A at CN8. The largest trivalent cation available is
     La(III) at 1.160 A. No A(III)2Bi2O6F can be a pyrochlore -- Y is short
     by 0.386 A and gives r_A/r_B = 1.06 against a floor of 1.46. The
     criterion is not a heuristic being stretched: it reproduces all four
     known A-Bi-O outcomes (Ba -> perovskite BaBiO3, Sr and Ca -> the real
     pyrochlores Sr2Bi2O7 and Ca2Bi2O7, Y -> the textbook defect-fluorite
     fast oxide-ion conductor). And a mobile-anion defect fluorite destroys
     A2, A3, A4, B1, B2 and F2 at once: no ordered Bi net, no ordered O'/F
     sublattice, no stoichiometric carrier count held for ten years.

     Independent confirmation from the same literature: in mixed-valent
     bismuth oxides Bi(III) and Bi(V) occupy DISTINCT crystallographic sites
     -- Bi(III) eight-coordinate, Bi(V) octahedral -- because the 6s2 lone
     pair cannot sit in a regular octahedron. The charge-ice manifold assumes
     both valences share one site type. The real system has an escape the
     model does not contain: Bi(III) leaves the octahedron for the A site.
     That relaxation is not frustrated at all.

WHAT IS ACTUALLY RIGHT, and worth keeping:
  - The combinatorics. Pyrochlore sites are the bonds of the dual diamond
    lattice, one Bi(V) per tetrahedron IS a perfect matching, and that
    manifold IS extensively degenerate.
  - The asymmetry. Ising antiferromagnetic order (CDW) is frustrated on the
    pyrochlore lattice while ferromagnetic XY order (uniform pairing) is not
    -- a ferromagnetic XY coupling cannot be frustrated on any lattice. In
    the negative-U pseudospin language J_z is frustrated and J_xy is not.
    That is a correct and genuinely useful design principle, and it is the
    first mechanism in this project that attacks charge order structurally
    rather than hoping it stays away.
  - THE NEAREST REAL VERSION. Divalent A rebalances to A(II)2Bi2O6F with
    v_Bi = +4.5, n_s = 0.5, and Sr2Bi2O6F sits at r_A/r_B = 1.52, inside the
    window, in a structure family that has actually been made. Its filling is
    intrinsic and its Bi net is frustrated. Two things must then be true, and
    both are cheap to decide: the sign of the effective Bi-O-Bi transfer must
    put the flat pair at the BOTTOM so quarter filling reaches it, and |U|
    must be real. Calibrated against BKBO even that is a tens-of-kelvin
    material. The principle survives the compound; the temperature does not
    survive the calibration.
""")
