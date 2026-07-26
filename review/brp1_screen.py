"""
BRP-1 = Ba3LaSn3BiO12, "stoichiometric Bi-Sn electronic pair resonance".

Screened by the candidate's own procedure. Steps 2 and 3 pass cleanly -- this
is the best-constructed guess so far, and the redox identity, the 1:3 B-site
topology and the tolerance factor are all correct. It fails at step 7 (the
dominant instability) and, decisively, against a measured member of its own
structure family.

WHAT THIS FILE CHECKS
  1. the charge arithmetic and the tolerance factor (both endpoints)
  2. THE MATERIAL-INDEPENDENCE TEST: does Tc depend on the Sn network at all?
  3. what the resonant boson is actually doing to the coupling
  4. the empirical anchor: this composition IS BaPb(0.75)Bi(0.25)O3 with Pb
     replaced by Sn, and that compound is measured
  5. delta0 from real band alignment, against their 1.5-3.0 eV window
  6. the size/breathing tension their own ordering argument creates
"""
import math
import sys
import numpy as np

sys.path.insert(0, "/tmp/brp")
import brp1_model as M                                    # noqa: E402

KB = 8.617333262e-5
print(__doc__)
print("=" * 78)

# ---------------------------------------------------------------------------
print("1. ARITHMETIC AND GEOMETRY -- both pass, and the tolerance factor is")
print("   right at BOTH redox endpoints, which the pyrochlore candidate never")
print("   managed.")
print(f"   A-site 3Ba2+ + La3+ = +{3*2+3};  O12 = -24;  B-sites must be "
      f"+{24-9}")
for lbl, tot in (("Bi3+ + 3Sn4+", 3 + 3 * 4), ("Bi5+ + 2Sn3+ + Sn4+",
                                               5 + 2 * 3 + 4)):
    print(f"      {lbl:<22} = +{tot}   {'OK' if tot == 15 else 'MISMATCH'}")
rA = (3 * 1.61 + 1.36) / 4          # Ba2+, La3+ at CN12
rO = 1.40
for lbl, rB in (("Bi3+ endpoint", (3 * 0.69 + 1.03) / 4),
                ("Bi5+ endpoint", (3 * 0.69 + 0.76) / 4)):
    t = (rA + rO) / (math.sqrt(2) * (rB + rO))
    print(f"      {lbl:<16} rB = {rB:.4f} A   Goldschmidt t = {t:.4f}"
          f"   {'perovskite OK' if 0.85 <= t <= 1.06 else 'FAILS'}")

# ---------------------------------------------------------------------------
print("\n2. THE MATERIAL-INDEPENDENCE TEST")
print("""   The design's stated virtue is a broad dispersive Sn network carrying
   both the pairing and the stiffness. If that is what produces 919 K, then
   Tc must depend strongly on the Sn hopping. Vary t and re-run THEIR solver
   with everything else fixed:""")
print(f"   {'t (eV)':>8}{'Sn span (eV)':>14}{'Tc_MF (K)':>12}{'Tc*t':>10}")
rows = []
for t in (0.20, 0.35, 0.70, 1.40):
    M.T_HOP = t
    M.HOPS = M.nearest_hops()
    E, V = M.precompute()
    tc = M.mean_field_tc(E) / KB
    rows.append((t, float(E.max()), tc))
    print(f"   {t:>8.2f}{float(E.max()):>14.3f}{tc:>12.1f}{tc*t:>10.0f}")
print("""   Tc is NOT monotonic in t, and it is not material-independent either. It
   PEAKS AT t = 0.70 eV -- the exact value chosen -- and falls off hard on
   both sides: -22% at half the hopping, -57% at double, and to ZERO at
   t = 0.20. That is worse than material-independence. delta0 is fixed at
   2.00 eV while mu depends on t, so the near-resonance delta0 = 2*mu is a
   COINCIDENCE BETWEEN TWO INDEPENDENTLY DETERMINED QUANTITIES: a redox
   energy set by Bi/Sn chemistry and a chemical potential set by orbital
   overlap. Their published scan varies delta0 at fixed t and finds 765-1019
   K, which looks robust; it is robust along one axis of a two-axis
   resonance. The joint requirement is the real one, and the design point
   sits on its maximum.""")
M.T_HOP = 0.700
M.HOPS = M.nearest_hops()
E, V = M.precompute()

# ---------------------------------------------------------------------------
print("\n3. WHAT THE RESONANT BOSON IS DOING")
print("""   The mediated interaction in their linearised equation is
       v_pair = 2 g^2 * tanh(d/2T)/d  -  U_c ,      d = delta0 - 2*mu
   and charge conservation pins d near zero, because a partially occupied
   reservoir IS a chemical-potential pin. tanh(d/2T)/d grows as T falls,
   toward its T -> 0 limit of 1/|d| -- so the attraction is a propagator
   evaluated near its own pole, and |d| is the only thing bounding it.
   Evaluated at their own self-consistent detuning d = -59.5 meV:""")
print(f"   {'T (K)':>8}{'tanh(d/2T)/d':>16}{'v_pair (eV)':>14}"
      f"{'v_pair/g':>11}{'v_pair/W':>11}")
d = -0.059519
W = float(E.max())
for TK in (2000, 919, 300, 100, 30):
    T = TK * KB
    chib = math.tanh(d / (2 * T)) / d
    vp = 2 * 0.45 ** 2 * chib - 1.0
    print(f"   {TK:>8}{chib:>16.3f}{vp:>14.3f}{vp/0.45:>11.2f}{vp/W:>11.3f}")
print(f"""   The attraction is not a material constant -- it is a propagator sitting
   near its pole. It grows as T falls until T ~ |d|/2 ~ 345 K, then saturates
   at 2g^2/|d| - U_c = {2*0.45**2/abs(d)-1.0:.2f} eV. Compare that with the entire Sn
   manifold, W = {W:.2f} eV:
       |V| / W  =  {(2*0.45**2/abs(d)-1.0)/W:.2f}
   The effective attraction is comparable to the FULL BANDWIDTH. At that
   coupling the attractive-Hubbard problem is past the BCS-BEC crossover:
   pairs are local, Tc is set by pair HOPPING rather than by the gap equation,
   and a BCS mean-field Tc is not a valid estimate -- it is exactly the regime
   where mean field overestimates most. The 919 K is computed with a
   weak-coupling formula at |V| ~ W.

   The physical content that is missing is RETARDATION. A real redox
   intermediate has dynamics; the interaction is 2g^2/(w - d), and the scale
   that cuts off the divergence is the boson's own dispersion, which this
   model sets to zero. The design explicitly dropped the instantaneous-
   interaction requirement for the phonon route -- but the electronic route
   then reinstates it in a stronger form: a boson with NO dispersion at all.""")

# ---------------------------------------------------------------------------
print("\n4. THE ANCHOR: THIS IS BaPb(0.75)Bi(0.25)O3 WITH Pb REPLACED BY Sn")
print("""   Strip the labels. BRP-1 is a cubic perovskite whose B site carries a
   group-14 s-band former plus Bi at exactly 1/4 occupancy, with Ba on A,
   at ambient pressure, pairing by Bi(III)/Bi(V) two-electron transfer into
   the group-14 band. That compound was made in 1975.""")
print(f"   {'compound':<28}{'Bi fraction':>13}{'Tc measured':>14}")
for c, x, tc in (("BaPb(1-x)Bi(x)O3, best x", "0.05-0.3", "13 K"),
                 ("Ba(1-x)K(x)BiO3, x~0.4", "1.0", "30 K"),
                 ("BRP-1 Ba3LaSn3BiO12", "0.25", "919 K claimed")):
    print(f"   {c:<28}{x:>13}{tc:>14}")
print("""   Sleight's BaPb(1-x)Bi(x)O3 peaks at 13 K where "the 6s band is about
   25% filled" -- the same 1:3 ratio BRP-1 selects. Ordering the Bi rather
   than leaving it random is a real difference, but it is a factor-of-two
   kind of difference, not a factor of 70.

   AND THE SUBSTITUTION RUNS THE WRONG WAY. BaPbO3 is a SEMIMETAL: the Pb 6s
   hybridises so strongly with O 2p that the bands inverted and there is no
   gap. BaSnO3 is a 3.1 eV BAND INSULATOR with the Sn 5s conduction band well
   above O 2p. Alloying Sn into BaPbO3 is a documented route to DESTROY the
   metallicity -- Ba(Pb,Sn)O3 is studied precisely as a metal-to-semiconductor
   series. BRP-1 takes the one working member of this family and swaps out the
   element responsible for it working.""")

# ---------------------------------------------------------------------------
print("\n5. delta0 FROM BAND ALIGNMENT, AGAINST THEIR 1.5-3.0 eV WINDOW")
print("""   delta0 is the cost of moving the Bi 6s^2 pair to the bottom of the
   group-14 s band -- two electrons, so twice the level separation. The Bi(III)
   6s^2 lone pair sits at the top of the O 2p valence band (that is the
   standard lone-pair picture, and why Bi(III) perovskites are studied as
   p-type oxides). So delta0 ~ 2 * (CBM - VBM) = 2 * Eg of the host stannate
   or plumbate:""")
print(f"   {'host':<12}{'character':<22}{'Eg (eV)':>9}{'delta0 ~ 2Eg':>15}"
      f"{'vs 1.5-3.0 window':>20}")
for h, ch, eg in (("BaPbO3", "semimetal, inverted", 0.0),
                  ("BaSnO3", "band insulator", 3.1)):
    d0 = 2 * eg
    v = ("BELOW" if d0 < 1.5 else "inside" if d0 <= 3.0 else "ABOVE")
    print(f"   {h:<12}{ch:<22}{eg:>9.2f}{d0:>15.2f}{v:>20}")
print("""   The window is bracketed but never hit: Pb gives ~0 (the pair is already
   delocalised -- which is why BaPbO3 is a metal and why BPBO superconducts at
   all), Sn gives ~6 eV (the pair stays on Bi and the compound is an
   insulator). Their design point delta0 = 2.0 eV corresponds to a host gap of
   ~1 eV, i.e. a group-14 perovskite that is neither BaPbO3 nor BaSnO3. Such a
   host would have to be engineered, and the obvious way to get it is a
   Ba(Pb,Sn)O3 alloy -- which returns the candidate to the BPBO family and its
   13 K.""")

# ---------------------------------------------------------------------------
print("\n6. THE ORDERING ARGUMENT AND THE BREATHING MODE ARE THE SAME FACT")
r = {"Bi3+": 1.03, "Bi5+": 0.76, "Sn4+": 0.69}
print(f"   Bi(III) CN6 = {r['Bi3+']:.2f} A, Bi(V) CN6 = {r['Bi5+']:.2f} A, "
      f"Sn(IV) CN6 = {r['Sn4+']:.2f} A")
print(f"   Bi/Sn size contrast at the Bi(III) endpoint: "
      f"{(r['Bi3+']-r['Sn4+'])/r['Sn4+']:.0%}   -> strongly drives the 1:3"
      f" B-site order (a genuine strength)")
print(f"   Bi(III) -> Bi(V) radius change on the SAME site:  "
      f"{(r['Bi3+']-r['Bi5+'])/r['Bi3+']:.0%}")
print("""   Those are one fact, not two. A 26% radius change accompanying the
   two-electron transfer IS a large breathing-mode coupling on the Bi site,
   and a large breathing coupling is exactly what makes BaBiO3 a
   charge-ordered insulator rather than a metal. The size contrast that makes
   the ordered structure favourable is the same contrast that localises the
   transferred pair. Their own hard condition 7 ("no static Bi(III)/Bi(V)
   order, phase separation, or breathing insulator") is therefore not an
   independent risk to be checked later -- it is in tension with condition 1
   (the ordering forms) by construction, and the model contains no phonon at
   all with which to test it.""")

print("""
================================================================================
VERDICT ON BRP-1
================================================================================
This is the best-built candidate of the series and it should be said plainly.
The redox identity is exact and derived rather than assigned; g_pair comes
from eliminating a named physical intermediate rather than from inverting a
matrix; the 1:3 topology genuinely leaves a 4-connected 3-D Sn network with no
Bi-Bi bonds; the Fermi-level bands are dispersive (1.4-2.5 eV) rather than
flat; the tolerance factor is ~0.96 at BOTH redox endpoints; and both Bi
valences stay octahedral. Every specific failure of the previous four
candidates is genuinely avoided.

It fails on three things instead.

  1. THE RESONANCE IS A TWO-PARAMETER COINCIDENCE, AND THE DESIGN POINT SITS
     ON ITS PEAK. Tc(t) at fixed delta0 = 2.0 eV runs 0 / 720 / 919 / 397 K
     at t = 0.20 / 0.35 / 0.70 / 1.40. The maximum is at t = 0.70, the value
     chosen. delta0 comes from Bi/Sn redox chemistry and mu comes from
     orbital overlap; nothing links them, so delta0 = 2*mu is a coincidence
     that must be arranged. The published robustness scan varies delta0 at
     fixed t and therefore moves ALONG the resonance ridge rather than across
     it.

  2. THE GAP EQUATION IS SOLVED OUTSIDE ITS VALIDITY. Charge conservation
     pins d = delta0 - 2*mu near zero -- correctly, that is what a partially
     filled reservoir does -- and the mediated attraction then saturates at
     2g^2/|d| - U_c = 5.81 eV against a total Sn bandwidth of 5.57 eV. At
     |V|/W ~ 1 the pairs are local, Tc is set by pair hopping, and BCS mean
     field is at its least reliable. A 919 K weak-coupling number computed at
     strong coupling is not a screen, it is a category error. The static
     treatment is the root: a real redox intermediate has dispersion, and
     that dispersion is what would cut off 2g^2/|d|.

  3. THE ANCHOR IS MEASURED AND IT IS 13 K. BRP-1 is BaPb(1-x)Bi(x)O3 at the
     x ~ 0.25 optimum, ordered rather than random, with Pb replaced by Sn.
     Sleight 1975: 13 K. BKBO: 30 K. And the Pb -> Sn swap is the documented
     way to turn BaPbO3 from a semimetal into a semiconductor, because
     BaSnO3's Sn 5s band sits 3.1 eV up while BaPbO3's Pb 6s is inverted into
     O 2p. Band alignment then puts delta0 ~ 6 eV for Sn against their
     1.5-3.0 eV window, and ~0 eV for Pb. The window is bracketed, not hit.

WHAT TO KEEP. The 1:3 ordered pair-centre topology is a good idea and is new
here: periodic rather than random negative-U centres, no centre-centre bonds,
conductor network still 3-D and 4-connected. If it is worth anything it is
worth trying on the host that already works -- ORDERED Ba(Pb,Bi)O3 at x = 1/4,
where the measured disordered version gives 13 K and the question "how much
does ordering the Bi buy?" is answerable and has never been answered. That is
a real, cheap, publishable question. It is not a 300 K question.

THE PATTERN, FIVE CANDIDATES IN. Every design has now died at the same place:
a coupling constant asserted at a value no measured member of its family
reaches. Bi-Sn g_pair = 0.45 eV here, lam = 2.7 there, U_eff = -1.05 eV
before that, t1 = 0.70 eV before that. The next candidate should be chosen so
that its central coupling is one that has ALREADY BEEN MEASURED in some
compound, and the design's job is to arrange geometry around a known number
rather than to require a new one.
""")
