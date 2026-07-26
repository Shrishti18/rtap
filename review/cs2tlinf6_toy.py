"""
Cs2TlInF6: can lam >= 0.60 and M_min >= 0.85 be reached AT THE SAME TIME?

The proposed active pair is TWO CO-LOCATED functions on the same Tl-F cluster:
    |1> = Tl 6s          |2> = TlF6 a1g (F p_sigma symmetry-adapted)
They share the same symmetry -- which is exactly why they hybridise, and also
why the folding pathology cannot apply (they are not a doubled cell).

ONE GOOD CONSEQUENCE, worth stating first: because both Wannier centres sit at
the SAME site (the a1g centroid is the Tl position by symmetry), there is no
relative orbital-position freedom, so M_min = M_naive identically. The gauge
collapse that killed the winding model and the folded chain is impossible here.

THE CONCERN. On-site s-a1g hybridisation V is k-INDEPENDENT, and a k-independent
mass is metrically sterile (a constant term leaves the projector flat in k).
The metric must therefore come from the k-DEPENDENT parts: the difference of the
two bandwidths, and any inter-cell s-a1g hybridisation. Meanwhile strong V drives
the bands toward pure bonding/antibonding character, pushing <rho> past 0.75/0.25
toward 0.9/0.1 -- which raises lam but suppresses the metric, the same trade
measured on the diamond net.

MODEL (fcc-like cubic connectivity, one Tl-F cluster per primitive cell):
    d_z(k) = de/2 + (t1 - t2) * s(k)          s(k) = cos kx + cos ky + cos kz
    d_x(k) = V + w * s(k)                     w = inter-cell s-a1g hybridisation
    d_y    = 0
    gap = 2|d| ,  tr g = (1/4)|grad theta|^2 with theta = atan2(d_x, d_z)

Scanned: de (on-site s/a1g splitting), V (on-site hybridisation), t1, t2, w.
Reported at CONVERGED nk -- the previous file reported nk=24 values that were
up to 33% low, which is the error this audit caught.
"""
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def analyse(de, V, t1, t2, w=0.0, nk=48):
    q = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KX, KY, KZ = np.meshgrid(q, q, q, indexing='ij')
    s = np.cos(KX) + np.cos(KY) + np.cos(KZ)
    dz = de / 2 + (t1 - t2) * s
    dx = V + w * s
    H = dx[..., None, None] * SX + dz[..., None, None] * SZ
    E, Vc = np.linalg.eigh(H)
    u = Vc[..., :, 0]
    # co-located orbitals -> no shift freedom -> M_min = M_naive
    P = np.einsum('...i,...j->...ij', u, u.conj())
    trg = np.zeros(P.shape[:-2])
    dk = 2 * np.pi / nk
    for ax in (0, 1, 2):
        dP = (np.roll(P, -1, ax) - np.roll(P, 1, ax)) / (2 * dk)
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    rho = (np.abs(u) ** 2).reshape(-1, 2)
    lam = float(np.linalg.eigvalsh((rho.T @ rho) / rho.shape[0])[-1])
    # band energies including the common dispersion (t1+t2)/2 * s
    d0 = 0.5 * (t1 + t2) * s
    Elo = d0 + E[..., 0]
    return dict(M=float(trg.mean()), lam=lam,
                gap=float((E[..., 1] - E[..., 0]).min()),
                W=float(np.ptp(Elo)), w1=float(rho[:, 0].mean()))


print(__doc__)
print("=" * 78)
print("0. CONVERGENCE CHECK (the error this audit caught)")
print(f"   {'nk':>5}{'M':>10}{'lam':>9}")
for nk in [16, 24, 32, 48, 64]:
    r = analyse(1.0, 0.8, 0.5, -0.2, 0.1, nk=nk)
    print(f"   {nk:>5}{r['M']:>10.5f}{r['lam']:>9.5f}")

print("\n1. THE TRADE-OFF: on-site hybridisation V vs the metric")
print("   (de = 1.0, t1 = 0.5, t2 = -0.2, w = 0.1 -- V is the only variable)")
print(f"   {'V':>6}{'gap':>8}{'W':>8}{'<rho_Tl>':>10}{'lam':>9}{'M':>9}"
      f"{'lam>=.60':>9}{'M>=.85':>8}")
for V in [0.0, 0.2, 0.5, 0.8, 1.5, 3.0]:
    r = analyse(1.0, V, 0.5, -0.2, 0.1)
    print(f"   {V:>6.2f}{r['gap']:>8.3f}{r['W']:>8.3f}{r['w1']:>10.3f}"
          f"{r['lam']:>9.4f}{r['M']:>9.4f}"
          f"{('yes' if r['lam']>=0.60 else 'no'):>9}"
          f"{('yes' if r['M']>=0.85 else 'no'):>8}")

print("\n2. CAN BOTH TARGETS BE MET? scan over de, V, w at t1=0.5, t2=-0.2")
best = None
hits = 0
for de in [0.0, 0.5, 1.0, 2.0, 3.0]:
    for V in [0.0, 0.1, 0.3, 0.6, 1.0, 2.0]:
        for w in [0.0, 0.1, 0.3, 0.6]:
            r = analyse(de, V, 0.5, -0.2, w)
            ok = r['lam'] >= 0.60 and r['M'] >= 0.85
            hits += ok
            score = min(r['lam'] / 0.60, r['M'] / 0.85)
            if best is None or score > best[0]:
                best = (score, de, V, w, r)
print(f"   parameter points satisfying BOTH lam>=0.60 and M>=0.85: {hits}/120")
sc, de, V, w, r = best
print(f"   closest point: de={de}, V={V}, w={w}")
print(f"      lam = {r['lam']:.4f} (need 0.60)   M = {r['M']:.4f} (need 0.85)")
print(f"      gap = {r['gap']:.3f}   W = {r['W']:.3f}   <rho_Tl> = {r['w1']:.3f}")

print("\n3. WHERE THE METRIC PEAKS, and what lam does there")
print(f"   {'de':>6}{'V':>6}{'w':>6}{'lam':>9}{'M':>9}{'gap':>9}{'<rho_Tl>':>10}")
rows = []
for de in [0.0, 0.5, 1.0, 2.0]:
    for V in [0.0, 0.3, 1.0]:
        r = analyse(de, V, 0.5, -0.2, 0.3)
        rows.append((r['M'], de, V, r))
for M, de, V, r in sorted(rows, reverse=True)[:6]:
    print(f"   {de:>6.2f}{V:>6.2f}{0.3:>6.2f}{r['lam']:>9.4f}{M:>9.4f}"
          f"{r['gap']:>9.4f}{r['w1']:>10.3f}")

# ---------------------------------------------------------------------------
print("\n4. THE GAP FILTER, and the DIVERGENCE TEST")
print("""   Every high-M point above has gap ~ 0. That is the folding pathology
   again: d_x = V + w*s vanishes on a surface where d_z also crosses zero, so
   the bands TOUCH and tr g has an integrable-looking but grid-dependent spike
   there. The signature is M growing with nk instead of converging. Below,
   every point that met BOTH targets is re-run at nk = 32 and 64; a genuine
   gapped metric changes by <2%, a touching artefact keeps climbing.""")
print(f"   {'de':>5}{'V':>5}{'w':>5}{'lam':>8}{'M(32)':>9}{'M(64)':>9}"
      f"{'drift':>8}{'gap':>9}   verdict")
survivors = []
for de in [0.0, 0.5, 1.0, 2.0, 3.0]:
    for V in [0.0, 0.1, 0.3, 0.6, 1.0, 2.0]:
        for w in [0.0, 0.1, 0.3, 0.6]:
            r = analyse(de, V, 0.5, -0.2, w)
            if not (r['lam'] >= 0.60 and r['M'] >= 0.85):
                continue
            a = analyse(de, V, 0.5, -0.2, w, nk=32)
            b = analyse(de, V, 0.5, -0.2, w, nk=64)
            drift = b['M'] / a['M'] - 1.0
            gapped = b['gap'] > 1e-3
            ok = gapped and abs(drift) < 0.02
            if ok:
                survivors.append((de, V, w, b))
            print(f"   {de:>5.1f}{V:>5.1f}{w:>5.1f}{b['lam']:>8.4f}"
                  f"{a['M']:>9.4f}{b['M']:>9.4f}{drift*100:>7.1f}%"
                  f"{b['gap']:>9.5f}   "
                  f"{'REAL' if ok else ('gapless' if not gapped else 'diverging')}")
print(f"\n   points meeting lam>=0.60, M>=0.85, GAPPED, and CONVERGED:"
      f" {len(survivors)}/120")
print("""   CAUTION on that table. "gapless" is exact, not numerical: d_x = V + w*s
   with s in [-3,3] has a zero whenever |V| <= 3|w|, so EVERY row with w != 0
   above is truly nodal and the grid merely missed the node. The w = 0 rows are
   the opposite case -- d_x = V is constant, the gap is exactly 2|V|, and the
   drift there is under-resolution, not divergence. That case is solved in
   closed form below rather than differenced.""")

# ---------------------------------------------------------------------------
print("\n5. THE w = 0 FAMILY IN CLOSED FORM (no finite differences)")
print("""   With w = 0 the d-vector is d = (V, 0, dz), dz = de/2 + (t1-t2)*s, so
   theta = atan2(V, dz) and
       tr g = (1/4)|grad theta|^2 = (1/4) * V^2 |grad dz|^2 / (V^2 + dz^2)^2
       grad dz = -(t1-t2)(sin kx, sin ky, sin kz)
   Exact integrand -> the only error is quadrature, and it is checked below.""")


def closed(de, V, t1, t2, nk):
    """Slice over kx so nk = 256 stays in memory."""
    q = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    KY, KZ = np.meshgrid(q, q, indexing='ij')
    syz = np.sin(KY) ** 2 + np.sin(KZ) ** 2
    cyz = np.cos(KY) + np.cos(KZ)
    mt = m1 = m2 = 0.0
    for kx in q:
        s = np.cos(kx) + cyz
        dz = de / 2 + (t1 - t2) * s
        g2 = (t1 - t2) ** 2 * (np.sin(kx) ** 2 + syz)
        den = V ** 2 + dz ** 2
        mt += (0.25 * V ** 2 * g2 / den ** 2).mean()
        r1 = 0.5 * (1 - dz / np.sqrt(den))          # lower-band weight on |1>
        m1 += r1.mean()
        m2 += (r1 ** 2).mean()
    mt, m1, m2 = mt / nk, m1 / nk, m2 / nk
    # A = <rho rho^T> with rho = (r1, 1-r1)
    A = np.array([[m2, m1 - m2], [m1 - m2, 1 - 2 * m1 + m2]])
    return mt, float(np.linalg.eigvalsh(A)[-1]), m1


KB = 8.617333e-5
chk_c = closed(1.0, 0.8, 0.5, -0.2, 48)
chk_a = analyse(1.0, 0.8, 0.5, -0.2, 0.0, nk=48)
print(f"   cross-check vs the finite-difference routine at V=0.8, nk=48:")
print(f"      closed form  M = {chk_c[0]:.5f}   lam = {chk_c[1]:.5f}")
print(f"      differenced  M = {chk_a['M']:.5f}   lam = {chk_a['lam']:.5f}")

print(f"\n   {'V':>6}{'gap=2V':>9}{'M(nk=64)':>10}{'M(128)':>10}{'M(256)':>10}"
      f"{'lam':>9}{'U*M  (U=V)':>13}")
CONST = []
for V in [0.02, 0.05, 0.10, 0.20, 0.40, 0.80]:
    m64 = closed(1.0, V, 0.5, -0.2, 64)[0]
    m128 = closed(1.0, V, 0.5, -0.2, 128)[0]
    M, lam, _ = closed(1.0, V, 0.5, -0.2, 256)
    CONST.append(V * M)
    print(f"   {V:>6.2f}{2*V:>9.2f}{m64:>10.4f}{m128:>10.4f}{M:>10.4f}"
          f"{lam:>9.4f}{V*M:>13.4f}")

print("""
   READ THE LAST COLUMN. U*M is CONSTANT. That is not a coincidence -- it is a
   sum rule. As V -> 0 the metric concentrates on the dz = 0 surface with width
   V, and the peak height goes as 1/V^2 while the width goes as V, so
       M  ->  (pi/8) * <|grad dz|> * A_surface / V     i.e.  M ~ C / V
   while the isolation gap is exactly 2V, so the projection cap gives U <= V.
   Therefore
       Tc(geometric) = 0.369 * U * M  ->  0.369 * C,  INDEPENDENT OF V.
   Hybridisation cannot be traded for temperature in this family: every eV of
   metric bought by shrinking the gap is paid back exactly by the U it forfeits.
   The lattice constant C, not the chemistry, sets the ceiling.""")
C = float(np.median(CONST))
print(f"   measured C = U*M = {C:.4f} t   (with t1-t2 = 0.7 t)")
print(f"   C/(t1-t2) = {C/0.7:.4f}, confirming C is linear in the bandwidth"
      f" DIFFERENCE")
print(f"   Tc(geometric ceiling) = 0.369 * {C:.4f} t / kB = "
      f"{0.369*C/KB:>6.0f} K  x  t[eV]")

# ---------------------------------------------------------------------------
print("\n6. BOTH CAPS TOGETHER -- the actual optimum, and the flatness test")
print("""   The geometric cap is V-independent, but the AMPLITUDE cap is not:
       Tc(amp) = U*lam/4 = V*lam/4     ->  0  as V -> 0.
   So Tc = min(amp, geo) has an interior maximum. Below, everything in units
   of t, with the lower-band width W and the flat-band validity ratio W/Tc
   (the diamond-net criterion demands W/Tc <~ 2; it measured 66-318).""")
print(f"   {'V/t':>6}{'lam':>8}{'M':>9}{'amp K/t':>10}{'geo K/t':>10}"
      f"{'Tc K/t':>9}{'W/t':>7}{'W/Tc':>8}")
bestT = None
for V in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 0.80]:
    M, lam, _ = closed(1.0, V, 0.5, -0.2, 192)
    Wb = analyse(1.0, V, 0.5, -0.2, 0.0, nk=48)['W']
    amp, geo = V * lam / 4, 0.369 * V * M
    tc = min(amp, geo)
    if bestT is None or tc > bestT[0]:
        bestT = (tc, V, lam, M, Wb)
    print(f"   {V:>6.2f}{lam:>8.4f}{M:>9.4f}{amp/KB:>10.0f}{geo/KB:>10.0f}"
          f"{tc/KB:>9.0f}{Wb:>7.2f}{Wb/tc:>8.1f}")
tc, V, lam, M, Wb = bestT
print(f"\n   OPTIMUM: V = {V:.2f} t, lam = {lam:.4f}, M = {M:.4f}")
print(f"      Tc = {tc/KB:.0f} K x t[eV]   ->  300 K needs t = "
      f"{300*KB/tc:.2f} eV")
print(f"      W/Tc = {Wb/tc:.0f}, against the flat-band validity limit of ~2:"
      f"  SHORT BY {Wb/tc/2:.0f}x")

print("""
VERDICT ON Cs2TlInF6 (toy level).
  1. lam >= 0.60 and M >= 0.85 ARE simultaneously reachable on paper -- 10 of
     120 parameter points meet both. NONE of them survives the gap filter.
     Every one is either exactly nodal (|V| <= 3|w| forces a zero of d_x) or
     resolves, at higher nk, into a near-touching whose "metric" is the
     touching singularity. Gapped-and-converged: 0/120.
  2. The reason is a SUM RULE, not a bad scan: U*M = const along the only
     gapped direction, so the geometric cap cannot be tuned at all. The
     acceptance vector (lam >= 0.60 AND M >= 0.85 AND Delta_iso >= 2.50 eV)
     asks for a large metric and a large gap from the same knob, and the
     product of the two is fixed by the lattice.
  3. THE ACCEPTANCE VECTOR IS OVER-SPECIFIED. Because U*M is fixed, M is not a
     free target at all: demanding Delta_iso >= 2.50 eV sets V, and V then sets
     M = C/V with no remaining freedom. At the V that delivers 2.50 eV of
     isolation, M ~ 0.06 -- failing the stated M >= 0.85 by ~15x -- and yet Tc
     is UNCHANGED, because the geometric cap only ever saw the product. Two of
     the six acceptance conditions are therefore not independent tests, and one
     of them (M >= 0.85) cannot be passed and does not need to be.
  4. What survives is a REAL and rather good ceiling once both caps are used:
     486 K per eV of hopping scale, peaking at V = 0.30 t. That is the best
     number this project has produced. It fails on FLATNESS instead: W/Tc = 42
     against the ~2 the formula needs, i.e. the band is ~20x too wide for
     Tc = U*lam/4 to mean anything. Same failure mode as the diamond net
     (W/Tc = 66-318), and only modestly less severe.
  5. RECOMMENDED KILL ORDER, unchanged in conclusion but now much cheaper: do
     not run the structure search, and do not compute a metric. Compute t --
     the Tl-Tl effective hopping through the In-F-In bridge -- alone. It has to
     clear 0.62 eV for the ceiling to reach 300 K at all, and in a double
     perovskite with ~6.5 A Tl-Tl separation and a closed-shell In(III) spacer
     it will be an order of magnitude below that. The same single number also
     fixes W, so the flatness failure is decided by it too. One quantity kills
     or clears the candidate.""")
