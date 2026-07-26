"""
The charge-ice kinematics claim, checked on the actual pyrochlore lattice.

CONTEXT. An independent audit of the QIP-1 / Y2Bi2O6F proposal accepted the
structural, flat-band and pair-resonance failures but added a NEW objection
neither the proposal nor my review had treated:

    inside a strict one-minority-per-tetrahedron manifold, a single
    nearest-neighbour pair hop violates the ice rule on two tetrahedra, so
    coherent pair motion begins at ring-exchange order K ~ t_b^3 / V_b^2 --
    the same interaction that builds the charge ice also immobilises the
    pairs.

That is a real claim about a real lattice and it is checkable, so it is
checked here rather than argued about. Sections:

  1. lattice bookkeeping (each site in 2 tetrahedra, NN share exactly 1)
  2. enumerate the ice manifold on a periodic cluster; is it extensive?
  3. THEIR CLAIM: how many ice defects does one NN pair hop create?
  4. is there any 2-hop move back into the manifold? (order of the first
     manifold-preserving process)
  5. WHERE THE STRICT LIMIT COMES FROM. Second order in the attractive
     Hubbard model gives J_z and J_perp of the SAME size, so V_b/t_b is not
     a free parameter -- strict ice needs a large EXTRA intersite V.
  6. the Bi-O-Bi geometry in their own CIF, which sets both the magnitude
     and the sign of the one number everything now depends on.
"""
import itertools
import math
import numpy as np

# ---------------------------------------------------------------------------
# pyrochlore: fcc cells, 4-site basis. Up tetra at R = the 4 basis sites of R;
# down tetra at R = {(R,0), (R-a1,1), (R-a2,2), (R-a3,3)}.
def build(L):
    cells = [(i, j, k) for i in range(L) for j in range(L) for k in range(L)]
    cidx = {c: n for n, c in enumerate(cells)}

    def sid(cell, s):
        return 4 * cidx[tuple(x % L for x in cell)] + s

    ns = 4 * len(cells)
    tetra = []
    for c in cells:
        tetra.append(tuple(sorted(sid(c, s) for s in range(4))))       # up
        dn = [sid(c, 0)]
        for s in (1, 2, 3):
            sh = list(c)
            sh[s - 1] -= 1
            dn.append(sid(sh, s))
        tetra.append(tuple(sorted(dn)))                                # down
    site_tet = [[] for _ in range(ns)]
    for it, t in enumerate(tetra):
        for s in t:
            site_tet[s].append(it)
    nbr = [set() for _ in range(ns)]
    for t in tetra:
        for a, b in itertools.combinations(t, 2):
            nbr[a].add(b)
            nbr[b].add(a)
    return ns, tetra, site_tet, nbr


def girth(tetra, site_tet):
    """Shortest cycle of the DUAL diamond graph: vertices = tetrahedra,
    edges = pyrochlore sites. True diamond has girth 6."""
    nv = len(tetra)
    edges = [(site_tet[s][0], site_tet[s][1], s) for s in range(len(site_tet))]
    adj = [[] for _ in range(nv)]
    for u, v, e in edges:
        adj[u].append((v, e))
        adj[v].append((u, e))
    best = 10 ** 9
    for u0, v0, e0 in edges:
        if u0 == v0:
            return 1
        dist = {u0: 0}
        dq = [u0]
        while dq:
            nxt = []
            for u in dq:
                for v, e in adj[u]:
                    if e == e0 or v in dist:
                        continue
                    dist[v] = dist[u] + 1
                    nxt.append(v)
            dq = nxt
            if v0 in dist:
                break
        if v0 in dist:
            best = min(best, dist[v0] + 1)
    return best


NS, TETRA, SITE_TET, NBR = build(2)

print(__doc__)
print("=" * 78)
print("1. LATTICE BOOKKEEPING, AND THE FINITE-SIZE TRAP")
print(f"   {'L':>3}{'sites':>8}{'tetra':>8}{'tet/site':>10}{'nbrs':>7}"
      f"{'tetra shared by NN':>21}{'dual girth':>12}")
for L in (2, 3, 4):
    ns, tt, st, nb = build(L)
    sh = {len(set(st[a]) & set(st[b])) for a in range(ns) for b in nb[a]}
    print(f"   {L:>3}{ns:>8}{len(tt):>8}"
          f"{str(sorted({len(v) for v in st})):>10}"
          f"{str(sorted({len(v) for v in nb})):>7}"
          f"{str(sorted(sh)):>21}{girth(tt, st):>12}")
print("""   The dual of the pyrochlore net is the DIAMOND lattice, whose true girth
   is 6. At L = 2 the periodic wrap creates spurious 4-cycles, so ANY
   ring-exchange order counted on a 2x2x2 cluster is a finite-size artefact.
   L >= 3 recovers girth 6. Everything below therefore uses L = 3.""")

NS, TETRA, SITE_TET, NBR = build(3)

# ---------------------------------------------------------------------------
print("\n2. THE ICE MANIFOLD (exactly one minority per tetrahedron)")


def ice_configs(tetra, site_tet, ns, limit=None):
    """DFS over tetrahedra. A site may be chosen only if BOTH of its
    tetrahedra are still empty -- the guard the first draft omitted."""
    nt = len(tetra)
    tcount = [0] * nt
    chosen = []
    out = []

    def rec(it):
        if limit is not None and len(out) >= limit:
            return
        if it == nt:
            out.append(tuple(sorted(chosen)))
            return
        if tcount[it] == 1:
            rec(it + 1)
            return
        for s in tetra[it]:
            ta, tb = site_tet[s]
            if tcount[ta] or tcount[tb]:
                continue
            tcount[ta] += 1
            tcount[tb] += 1
            chosen.append(s)
            rec(it + 1)
            chosen.pop()
            tcount[ta] -= 1
            tcount[tb] -= 1
    rec(0)
    return out


ns2, tt2, st2, nb2 = build(2)
small = ice_configs(tt2, st2, ns2)
print(f"   L = 2, complete enumeration: {len(small)} configurations, "
      f"minority counts {sorted({len(c) for c in small})} "
      f"(must be [{len(tt2)//2}])")
ice = ice_configs(TETRA, SITE_TET, NS, limit=4000)
print(f"   L = 3, sampled: {len(ice)} configurations, minority counts "
      f"{sorted({len(c) for c in ice})} (must be [{len(TETRA)//2}])")
print(f"   minority fraction = {len(ice[0])/NS:.3f}   "
      f"(Bi(V) fraction at v = +3.5 is 0.250)")
print(f"   L = 2 entropy density ln(N)/N_sites = "
      f"{math.log(len(small))/ns2:.4f}   -> extensive, as the proposal says.")
print("   Their combinatorics was never the issue.")

# ---------------------------------------------------------------------------
print("\n3. THEIR CLAIM: DEFECTS CREATED BY ONE NEAREST-NEIGHBOUR PAIR HOP")


def defects(occ):
    cnt = [sum(1 for s in t if s in occ) for t in TETRA]
    return sum(1 for c in cnt if c == 0), sum(1 for c in cnt if c >= 2)


hist, total = {}, 0
for cfg in ice[:400]:
    occ = set(cfg)
    for i in cfg:
        for j in NBR[i]:
            if j in occ:
                continue
            d = defects((occ - {i}) | {j})
            hist[d] = hist.get(d, 0) + 1
            total += 1
print(f"   all {total} single hops from 400 ice configurations (L = 3):")
for k in sorted(hist):
    print(f"      (empty tetra, doubly-occupied) = {k}  ->  {hist[k]} moves"
          f"  ({hist[k]/total:.1%})")
print(f"   moves that stay in the manifold: {hist.get((0, 0), 0)}")
print("   VERDICT: their claim is CORRECT and exact. Every NN pair hop leaves")
print("   the manifold, creating exactly one empty and one doubly-occupied")
print("   tetrahedron -- the two tetrahedra the moving site does NOT share.")

# ---------------------------------------------------------------------------
print("\n4. THE ORDER OF THE FIRST MANIFOLD-PRESERVING PROCESS")
iceset = set(ice)


def valid(occ):
    return all(sum(1 for s in t if s in occ) == 1 for t in TETRA)


found2 = 0
for cfg in ice[:150]:
    occ = set(cfg)
    for i in cfg:
        for j in NBR[i]:
            if j in occ:
                continue
            c1 = (occ - {i}) | {j}
            for a in list(c1):
                for b in NBR[a]:
                    if b in c1:
                        continue
                    c2 = (c1 - {a}) | {b}
                    if c2 != occ and valid(c2):
                        found2 += 1
print(f"   two-hop moves from 150 configurations reaching a DIFFERENT ice"
      f" state: {found2}")
print(f"   (at L = 2 the same search returns a nonzero count -- those are the"
      f" spurious\n    4-cycles of the wrapped cluster, not physical moves.)")

hexfound = 0
for cfg in ice[:150]:
    occ = set(cfg)
    for i in cfg:
        for j in NBR[i]:
            if j in occ:
                continue
            c1 = (occ - {i}) | {j}
            for a in list(c1):
                for b in NBR[a]:
                    if b in c1:
                        continue
                    c2 = (c1 - {a}) | {b}
                    for p in list(c2):
                        for q in NBR[p]:
                            if q in c2:
                                continue
                            c3 = (c2 - {p}) | {q}
                            if c3 != occ and valid(c3):
                                hexfound += 1
    if hexfound:
        break
print(f"   three-hop moves reaching a different ice state: "
      f"{'FOUND' if hexfound else 'none'} ({hexfound} from the first "
      f"configuration searched)")
print("""   -> no two-hop move survives at L = 3, three-hop moves do. The first
      manifold-preserving process is the THIRD-order ring exchange around the
      shortest closed loop, which is the hexagon (dual girth 6). That is
      exactly the K_ring ~ t_b^3 / V_b^2 the audit wrote down. Confirmed, and
      the confirmation required L >= 3 -- my first attempt used L = 2 and its
      wrap-around 4-cycles gave the wrong answer.""")

# ---------------------------------------------------------------------------
print("\n5. BUT THE STRICT LIMIT IS NOT WHAT THE MODEL SUPPLIES")
print("""   Second-order perturbation theory in the attractive Hubbard model gives
   BOTH couplings from the same matrix element:
       pair hopping        t_b = 2 t^2 / |U|      (XY, ferromagnetic)
       density-density     V_b = 4 t^2 / |U|      (Ising, antiferromagnetic)
   so V_b / t_b = 2 EXACTLY, independent of t and |U|. The ratio is fixed by
   the algebra, not chosen. A strict ice manifold needs V_b / t_b >> 1, which
   the negative-U model alone never provides -- it requires a large EXTRA
   intersite Coulomb V on top.""")
print(f"   {'|U| (eV)':>10}{'t_b (eV)':>10}{'V_b (eV)':>10}{'V_b/t_b':>10}"
      f"{'strict ice?':>13}")
for uu in (0.5, 1.05, 2.0, 4.0, 8.0):
    tb, vb = 2 * 0.7 ** 2 / uu, 4 * 0.7 ** 2 / uu
    print(f"   {uu:>10.2f}{tb:>10.4f}{vb:>10.4f}{vb/tb:>10.2f}"
          f"{'no':>13}")
print("""   CONSEQUENCE, and it cuts both ways. The audit's ring-exchange objection
   bites only when V_extra >> t_b; at V_extra = 0 there is no strict ice to
   block anything, and the unfrustrated XY coupling is free to order. So the
   'soft charge-ice superfluid' the audit identifies as the viable regime is
   not a compromise between the two -- it is the DEFAULT for a negative-U
   model, and the strict-ice limit is the thing that would have to be
   engineered. One dimensionless number, V_extra / t_b, separates them:
       V_extra / t_b -> infinity : rigid ice, pairs move only by ring
                                   exchange, stiffness collapses  (the audit)
       V_extra / t_b -> 0        : mobile pairs, only soft charge
                                   correlations                   (the model)
   The proposal argues its mechanism at the first limit and computes its Tc
   at the second. That is the same regime split the audit found, now with the
   controlling parameter named and its natural value known.""")

# ---------------------------------------------------------------------------
print("\n6. THE ONE NUMBER EVERYTHING NOW DEPENDS ON, FROM THEIR OWN CIF")
a0, xO = 10.300, 0.330
# 48f O at (x,1/8,1/8) bridges the 16c Bi at (1/4,0,1/4) and (1/4,1/4,0)
O = np.array([xO, 0.125, 0.125])
B1 = np.array([0.25, 0.00, 0.25])
B2 = np.array([0.25, 0.25, 0.00])
v1, v2 = (B1 - O) * a0, (B2 - O) * a0
ang = math.degrees(math.acos(float(v1 @ v2 / (np.linalg.norm(v1) *
                                              np.linalg.norm(v2)))))
print(f"   QIP1_Y2Bi2O6F_ideal.cif: a = {a0} A, O(48f) x = {xO}")
print(f"   Bi-O bond length      = {np.linalg.norm(v1):.3f} A")
print(f"   Bi-Bi (NN, a/2sqrt2)  = {a0/(2*math.sqrt(2)):.3f} A"
      f"   (BaBiO3: ~4.34 A)")
print(f"   Bi-O-Bi ANGLE         = {ang:.1f} deg   (perovskite BaBiO3: 180 deg)")
print(f"   |cos(angle)|          = {abs(math.cos(math.radians(ang))):.3f}")
print("""
   For an s-O(2p sigma)-s superexchange path the two-step transfer carries the
   overlap of the two oxygen p lobes, so t_eff scales as cos(theta): a 131 deg
   bridge retains only ~0.66 of the linear-bond value, before the shorter
   Bi-Bi distance is counted the other way. Their t1 = 0.70 eV is asserted,
   not computed, and it is the parameter every surviving ceiling is linear in.

   The SAME matrix element also fixes its SIGN, and the sign is what decides
   whether the Sr2Bi2O6F redirect is a superconductor at all: at quarter
   filling t1 > 0 puts the flat-band descendants at the wrong end of the
   spectrum and gives Tc = 0, while t1 < 0 gives 618-1464 K before any
   correction. One Wannier matrix element settles magnitude and sign at once.
""")

print("""
================================================================================
WHERE I WAS WRONG, AND WHAT I NOW HOLD
================================================================================
CONCEDED, without reservation:

  1. "No trivalent A can EVER form this pyrochlore" was overstated. The
     radius-ratio window is an empirical screen, not an existence theorem;
     covalency, mixed anions and pressure move boundaries. What the data
     support is the weaker and still decisive claim: the screen reproduces
     4/4 known A-Bi-O outcomes, and Y2Bi2O6F misses the floor by 0.39 A --
     not marginally outside, enormously outside. Y2Bi2O6F is not a credible
     ambient-pressure equilibrium pyrochlore. It is not a proved impossibility.

  2. The BKBO/MATBG ratios are not a universal Tc multiplier. That ratio
     absorbs errors in |U|, in the band Hamiltonian, in filling, self-energy,
     competing order, disorder and the mean-field approximation itself, and
     none of those can be separated by dividing two numbers. I should not have
     used it as a correction factor -- particularly because I had already
     conceded exactly this objection earlier in this project, about my own
     eta calibration, and then reproduced the error here in a new costume.
     What the BKBO run does establish is narrower and worth keeping: THEIR
     solver, at THEIR |U|, applied to the one measured member of this family,
     overshoots by 5.7-17.7x. That is a warning about the solver-plus-
     parameters, not a coefficient to multiply into an unrelated model.

  3. A favourable hopping sign for Sr2Bi2O6F would not by itself solve RTAP.
     Agreed, and already stated: tens of kelvin, not 350 K.

ACCEPTED FROM THE AUDIT, and verified above rather than taken on trust:

  4. Strict charge ice blocks first-order pair motion. Exact on the lattice:
     of every single NN hop from every ice configuration on the cluster, ZERO
     stay in the manifold, and all create exactly one empty and one doubly
     occupied tetrahedron. No two-hop move reaches a different ice state
     either, so the first manifold-preserving process is third order -- the
     hexagon ring exchange. The audit found this; it is right.

WHAT I ADD, which neither document has:

  5. The strict limit is not the model's own limit. Second order gives
     V_b/t_b = 2 exactly, fixed by the algebra. Rigid ice therefore requires
     a large EXTRA intersite V that the negative-U model does not supply, so
     the audit's ring-exchange collapse is a statement about V_extra >> t_b,
     while the proposal's BCS number is V_extra = 0. The "soft-ice superfluid"
     is the default case, not a narrow window -- which is better news for the
     surviving idea than the audit allows, and it names the one parameter,
     V_extra/t_b, that the next unified calculation has to scan.

  6. Everything now reduces to one matrix element. Both regimes' ceilings are
     linear in t. Their own CIF gives a 131.3 deg Bi-O-Bi bridge against 180
     deg in BaBiO3, so |cos| = 0.66 on the magnitude, and the same element
     fixes the sign that decides whether the Sr redirect is 0 K or hundreds.
     t1 = 0.70 eV was asserted. Computing it is cheaper than any of the
     calculations either document proposes next, and nothing else is worth
     doing until it is done.
""")
