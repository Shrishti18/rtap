"""
The self-consistency selector, as a function. STANDING RULE.

History: six appearances of the same error family -- a quantity computed on a
projector whose rank does not match the degeneracy structure of what it
projects onto. The sixth was in a selector built to prevent the fifth, and the
seventh candidate was in the fix for the sixth. The recurring cause, stated by
the author and sharpened here:

    THE VALIDITY TESTS FOR THE TWO RANKS HAVE UNEQUAL STRENGTH. The
    single-band test is VACUOUS: its own U-cap enforces
        kT(Tc_single) <= U_s*lam/4 <= (d_near/2)/4 = d_near/8 <= s/8
    (lam = lam_max(A) <= tr A = <sum rho^2> <= 1), so "s >> Tc_single" holds
    IDENTICALLY. Any rule that consults it -- "which test passes", or a min/max
    over both predictions -- inherits a structural bias toward single-band,
    which is precisely the systematic the rule exists to prevent. At the
    degenerate limit s -> 0 a min-rule returns NEITHER for an exact Kramers
    pair, which is wrong by inspection.

Therefore the decision uses the ONE non-vacuous scale:

    kT_ref = U_m * lam_m / 4        (the manifold treatment's AMPLITUDE cap)

the pairing scale, which is what coherence across the splitting compares to,
and which -- unlike kT(Tc_manifold) = kT(min(amp, geo)) -- survives
M_manifold -> 0. Without that, structural zeros (bands spanning the same
orbitals, M == 0 by C2) drive Tc_manifold to 0 through the stiffness cap and
hand the decision back to single-band: the bias re-entering through the
structural-zero door.

The rule:
    s <  kT_ref          -> MANIFOLD   (pairs coherent at the pairing scale)
    s >  margin * kT_ref -> SINGLE     (partner spectrally inert)
    otherwise            -> NEITHER    (the answer depends on which Tc is
                                        right, which is what is being decided:
                                        circular region, report conservatively)
plus:
    M_manifold < mzero inside the MANIFOLD branch -> DEAD_BY_C2. Tc = 0 is a
    THEOREM (same-orbital manifold, k-independent projector, tr g == 0), not a
    failed treatment. It must never fall through to single-band.
"""
import numpy as np

KB = 0.086173e-3   # eV/K


def select(s, tc_single, tc_manifold, U_m, lam_m, M_manifold,
           rank=2, margin=3.0, mzero=1e-6):
    """Choose the self-consistent treatment for one band.

    s            intra-group splitting, eV (min over BZ)
    tc_single    single-band Tc, K  (U capped by the band's OWN nearest gap)
    tc_manifold  manifold Tc, K    (U capped by the GROUP's isolation gap)
    U_m, lam_m   the manifold treatment's U (eV) and pairing eigenvalue
    M_manifold   the manifold metric (structural-zero detection)
    rank         group size; rank 1 is single by definition

    Returns dict(verdict, Tc, kT_ref, note).
      verdict in {SINGLE, MANIFOLD, DEAD_BY_C2, NEITHER}
      NEITHER reports Tc = min(both) flagged as an unresolved upper bound.
    """
    if rank == 1:
        return dict(verdict="SINGLE", Tc=tc_single, kT_ref=np.nan,
                    note="rank-1 group")
    kT_ref = U_m * lam_m / 4.0                      # eV; amplitude scale
    if s < kT_ref:
        if M_manifold < mzero:
            return dict(verdict="DEAD_BY_C2", Tc=0.0, kT_ref=kT_ref,
                        note="same-orbital manifold: tr g == 0 by theorem")
        return dict(verdict="MANIFOLD", Tc=tc_manifold, kT_ref=kT_ref,
                    note="s < pairing scale: coherent")
    if s > margin * kT_ref:
        return dict(verdict="SINGLE", Tc=tc_single, kT_ref=kT_ref,
                    note="s > %g x pairing scale: partner inert" % margin)
    return dict(verdict="NEITHER", Tc=min(tc_single, tc_manifold),
                kT_ref=kT_ref,
                note="s inside the margin band: circular region, "
                     "Tc is an unresolved upper bound")


# ------------------------------------------------------------------ unit tests
def _tests():
    # 1. exact Kramers pair: s = 0 must select MANIFOLD (the min-rule fails
    #    this: tc_single = 0 makes s < kT(min) unsatisfiable).
    r = select(s=0.0, tc_single=0.0, tc_manifold=120.0, U_m=0.2, lam_m=0.33,
               M_manifold=0.8)
    assert r["verdict"] == "MANIFOLD" and r["Tc"] == 120.0, r

    # 2. same-orbital degenerate manifold: MUST NOT fall through to single.
    r = select(s=0.0, tc_single=855.0, tc_manifold=0.0, U_m=0.4, lam_m=0.4,
               M_manifold=1e-12)
    assert r["verdict"] == "DEAD_BY_C2" and r["Tc"] == 0.0, r

    # 3. wide splitting (the 855 K fixture): s = 1.6 eV >> pairing scale.
    r = select(s=1.6, tc_single=855.0, tc_manifold=300.0, U_m=0.4, lam_m=0.4,
               M_manifold=0.5)
    assert r["verdict"] == "SINGLE" and r["Tc"] == 855.0, r

    # 4. LiPbAu2 24/25: s = 0.0194, U_m = 0.195, lam_m = 0.193
    #    -> kT_ref = 0.00941, margin band [0.0094, 0.0282]: NEITHER, Tc = 9 K.
    r = select(s=0.0194, tc_single=9.0, tc_manifold=109.0, U_m=0.1952,
               lam_m=0.1931, M_manifold=1.0531)
    assert r["verdict"] == "NEITHER" and abs(r["Tc"] - 9.0) < 1e-9, r

    # 5. vacuity guard: the decision must not consult tc_single's scale.
    #    Same inputs, tc_single scaled 100x -> same verdict.
    a = select(s=0.05, tc_single=1.0, tc_manifold=200.0, U_m=0.3, lam_m=0.3,
               M_manifold=0.7)
    b = select(s=0.05, tc_single=100.0, tc_manifold=200.0, U_m=0.3, lam_m=0.3,
               M_manifold=0.7)
    assert a["verdict"] == b["verdict"], (a, b)
    print("selector: 5/5 tests pass")


if __name__ == "__main__":
    _tests()
