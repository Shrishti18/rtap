"""
Is the spin-triplet channel a door?

CLAIM UNDER TEST (the author's): "Triplet pairing evades the n_phi dilution
entirely because the pair wavefunction is antisymmetric in orbital rather than
spin."

The singlet dilution is lam = lam_max(A), A_ab = <rho_a rho_b>, = 1/n_phi for
uniform weights. The naive triplet analogue is sum_{a!=b} rho_a rho_b
= 1 - 1/n_phi, which GROWS with n_phi. If that were the whole story the claim
would be right.

Three things are tested here, in increasing order of severity.

T1. THE CHANNEL IS EMPTY WITHOUT SPIN-ORBIT COUPLING.
    The gauge-covariant on-site pair amplitude in orbital channel (a,b) is
        G_ab(k) = u_a(k) u_b(-k) - u_b(k) u_a(-k)
    (built with the -k EIGENVECTOR, not with a conjugate -- |G| is then
    invariant under the per-k phase, and under SU(2) mixing inside a doublet).
    Symmetric part: G_aa = u_a(k)u_a(-k), which reduces to rho_a >= 0, the
    singlet kernel of C1.

    For orbitals on the SAME inversion-symmetric site with the SAME parity --
    which is exactly the d-orbital case -- time reversal gives u(-k) = u(k)* and
    inversion gives u(-k) = u(k), so u is REAL and G_ab == 0 identically for
    a != b. The triplet channel is not diluted, it is EMPTY.

    SOC is what opens it: with spin, u is genuinely complex even with inversion.
    So the mechanism the author wants exists only where SOC is strong -- which
    is where the project already is, and where Delta_iso is already spent.

T2. THE ON-SITE TRIPLET CHANNEL IS REPULSIVE IN EVERY REAL d-ELECTRON SYSTEM.
    Kanamori two-electron energies on one site:
        intra-orbital singlet   U
        inter-orbital singlet   U' + J = U - J
        inter-orbital TRIPLET   U' - J = U - 3J
    Hund's J makes the triplet the LEAST repulsive configuration; it makes it
    attractive only if U < 3J, i.e. J/U > 1/3. Measured J/U is 0.15-0.25 (3d)
    and ~0.1 (5d). So U - 3J > 0 always. Hund's coupling aligns spins; it does
    not supply a pair attraction.

T3. EVEN GRANTING BOTH, THE CAP THAT BINDS IS NOT THE ONE TRIPLET MOVES.
    lam enters the AMPLITUDE cap Tc <= U*lam/4. The terminus is set by the
    STIFFNESS cap Tc <= coef*R*U*M_trg, which is 12.6x smaller and depends on
    the band's quantum metric, not on the pairing channel.
"""
import numpy as np
import jeff


KB = 0.086173e-3
C3 = 0.55 * 0.67


# ------------------------------------------------- T1/T3, j=1/2 (complex, SOC)
def channels_jeff(t=0.052, lam=0.417, nk=24):
    """j=1/2 Kramers doublet. SOC makes the Bloch functions genuinely complex,
    so the antisymmetric channel is NOT empty here. Pair amplitude inside a
    Kramers doublet is the pseudospin singlet Delta_nm = Delta * eps_nm, and

        F_{alpha beta}(k) = u1_alpha(k) u2_beta(-k) - u2_alpha(k) u1_beta(-k)

    is the epsilon-contraction, which is invariant under SU(2) mixing inside the
    doublet (it picks up only det U = 1). alpha = (orbital, spin), 6 components,
    orbital-major: index = 2*orb + spin.
    """
    kk = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), axis=-1)
    E, V = np.linalg.eigh(jeff.hk(K, t, lam, 0.0))
    u1 = V[..., :, 4].reshape(-1, 6)
    u2 = V[..., :, 5].reshape(-1, 6)
    # -k on the same grid: k -> -k maps index i -> (nk-1-i) for a half-offset grid
    def flip(a):
        a = a.reshape(nk, nk, nk, 6)
        return a[::-1, ::-1, ::-1].reshape(-1, 6)
    u1m, u2m = flip(u1), flip(u2)
    # F_{alpha beta}(k), antisymmetric under (alpha,k)<->(beta,-k)
    F = (np.einsum('ka,kb->kab', u1, u2m) - np.einsum('ka,kb->kab', u2, u1m))
    # channel decomposition on the physical on-site index alpha=(orb,spin)
    idx = lambda o, s: 2 * o + s
    intra = 0.0        # same orbital, opposite spin  -> the singlet channel
    for o in range(3):
        intra += np.abs(F[:, idx(o, 0), idx(o, 1)]) ** 2
    inter_eq = 0.0     # different orbital, SAME spin -> the equal-spin triplet
    for o in range(3):
        for p in range(o + 1, 3):
            for s in range(2):
                inter_eq += np.abs(F[:, idx(o, s), idx(p, s)]) ** 2
    inter_op = 0.0     # different orbital, opposite spin
    for o in range(3):
        for p in range(o + 1, 3):
            for s in range(2):
                inter_op += np.abs(F[:, idx(o, s), idx(p, 1 - s)]) ** 2
    tot = np.abs(F) ** 2
    tot = tot.sum(axis=(1, 2)) / 2      # each unordered pair counted twice
    return dict(intra=float(intra.mean()), inter_eq=float(inter_eq.mean()),
                inter_op=float(inter_op.mean()), total=float(tot.mean()))


print(__doc__)
print("=" * 72)
print("T1. On-site channel weights vs SOC. Same 3 t2g orbitals on one site, so")
print("    this is a genuine ON-SITE pair; lam=0 is the no-SOC control.")
print(f"    {'lam_SOC':>9}{'intra (singlet)':>17}{'inter, SAME spin':>18}"
      f"{'inter, opp spin':>17}{'triplet %':>11}")
for lam in [0.0, 0.001, 0.01, 0.05, 0.15, 0.417]:
    c = channels_jeff(lam=lam)
    pct = 100 * c['inter_eq'] / c['total'] if c['total'] > 1e-14 else 0.0
    print(f"    {lam:>9.3f}{c['intra']:>17.3e}{c['inter_eq']:>18.3e}"
          f"{c['inter_op']:>17.3e}{pct:>11.1f}")
c = channels_jeff()
print("    -> at lam=0 the equal-spin channel vanishes to machine precision.")
print("       SOC opens it, and at the physical Ir value it carries ~1/3 of the")
print("       on-site pair weight. This is the ONE place the mechanism is real.")

print()
print("T2. Kanamori sign check: is the inter-orbital triplet attractive?")
print(f"     {'ion':<10}{'U (eV)':>8}{'J (eV)':>8}{'J/U':>7}{'U-3J (eV)':>11}  channel")
for name, U, J in [('3d (Mn2+)', 4.0, 0.90), ('3d (Ni2+)', 4.5, 0.85),
                   ('4d (Ru)', 2.5, 0.40), ('5d (Ir4+)', 1.9, 0.20),
                   ('5d, Te host', 0.95, 0.15)]:
    s = 'ATTRACTIVE' if U - 3 * J < 0 else 'repulsive'
    print(f"     {name:<10}{U:>8.2f}{J:>8.2f}{J/U:>7.2f}{U-3*J:>11.2f}  {s}")
print("     attraction needs J/U > 1/3; measured J/U is 0.10-0.25. None qualify.")

print()
print("T3. Grant the channel AND the attraction. Does the terminus move?")
out, _ = jeff.analyse(t=0.052, lam=0.417, tp=0.0, nk=32, bands=(4, 5))
d = out[0]
iso, M, nphi = d['iso'], d['M_trg'], d['nphi']
U = iso / 2
for label, lm in [('singlet  lam = 1/n_phi', 1 / nphi),
                  ('triplet  lam = 1-1/n_phi', 1 - 1 / nphi),
                  ('triplet, best case lam = 1', 1.0)]:
    amp = U * lm / 4 / KB
    geo = C3 * U * M / KB
    print(f"     {label:<28} lam={lm:.3f}  amp={amp:7.1f} K  geo={geo:6.1f} K"
          f"   Tc = {min(amp, geo):6.1f} K")
print(f"     M_trg = {M:.4f} and U = {U:.4f} eV are untouched by the pairing")
print(f"     channel. The stiffness cap does not move. Terminus unchanged.")
