"""
Does the winding model carry current?

C4 (winding.py) is load-bearing: M_trg = n^2/4 independent of the gap, which is
what turned the 8.1 K terminus into 334 K and left the winding class as "the one
open door". The model is

    d(k) = D[cos(nk) sx + sin(nk) sy]   =>   H_12(k) = d_x - i d_y = D e^{-ink}

Fourier transforming: H_12(R) = D delta_{R,n}. A SINGLE hop, at range n, and
nothing else -- no A-A, no B-B, no A-B at any other range. So the lattice does
not merely have "suppressed short-range hopping": it decomposes into DECOUPLED
DIMERS, each pairing A_R with B_{R+n}.

That is the pathology the project already met once. From sawtooth.py's own
docstring, on why two earlier test lattices were discarded:

    "H(R=1) = [[t,0],[0,0]] is NOT nilpotent -> a pair CAN circulate the ring,
     unlike the two earlier test lattices where H(R=1)^2 = 0 forced zero
     transport and made the ED/BdG comparison meaningless"

If the winding model is decoupled dimers, a flux threading the ring can be
gauged away on each dimer independently, so the spectrum is flux-independent and
the superfluid stiffness is ZERO -- while tr g is n^2/4. The metric would be the
spread of a molecular orbital across a dimer of length n (exactly n^2/4 for
weight 1/2 on sites 0 and n), not a transport property.

Three checks, each decisive on its own.
"""
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)


def hk(k, n, D=1.0):
    k = np.atleast_1d(k)
    return D * (np.cos(n * k)[:, None, None] * SX
                + np.sin(n * k)[:, None, None] * SY)


print(__doc__)
print("=" * 78)

# ---------------------------------------------------------------- 1. real space
print("1. REAL-SPACE HOPPING  H(R) = (1/N) sum_k H(k) e^{+ikR}")
for n in [1, 2, 3]:
    nk = 512
    k = 2 * np.pi * np.arange(nk) / nk
    H = hk(k, n)
    print(f"   n = {n}:", end="")
    for R in range(0, 5):
        HR = (H * np.exp(1j * k * R)[:, None, None]).mean(0)
        if np.abs(HR).max() > 1e-10:
            print(f"   H(R={R})_12 = {HR[0,1]:+.3f}", end="")
    print()
print("   -> exactly one nonzero hop, at R = n. The A sublattice touches the B")
print("      sublattice ONLY at range n. Decoupled dimers {A_R, B_{R+n}}.")

# ------------------------------------------------------- 2. flux / stiffness
print("\n2. SUPERFLUID STIFFNESS BY DIRECT FLUX THREADING")
print("   N-site ring, flux theta, every hop of range n picks up e^{i n theta/N}.")


def ring_spectrum(N, n, theta, D=1.0):
    """2N x 2N single-particle matrix. site index = 2*R + (0=A, 1=B)."""
    M = np.zeros((2 * N, 2 * N), complex)
    ph = np.exp(1j * n * theta / N)
    for R in range(N):
        a, b = 2 * R, 2 * ((R + n) % N) + 1
        M[b, a] += D * ph
        M[a, b] += D * np.conj(ph)
    return np.linalg.eigvalsh(M)


N = 24
for n in [1, 2, 3]:
    es = [ring_spectrum(N, n, th).sum() for th in (0.0, 0.3, 0.7, 1.1)]
    spread = max(es) - min(es)
    print(f"   n = {n}:  sum of all eigenvalues vs theta -> spread = {spread:.3e}")
    occ = [np.sort(ring_spectrum(N, n, th))[:N].sum() for th in (0.0, 0.3, 0.7)]
    print(f"            half-filled band energy vs theta -> spread = "
          f"{max(occ)-min(occ):.3e}")
print("   -> flux-INDEPENDENT to machine precision. D_s = 0 EXACTLY.")
print("      Each dimer is a two-site problem with eigenvalues +-D|e^{i phi}| =")
print("      +-D; the phase is removable by a local gauge on each dimer.")

# ------------------------------------------------- 3. the metric is a spread
print("\n3. WHAT tr g IS MEASURING")
nk = 4000
k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
for n in [1, 2, 3, 4]:
    E, V = np.linalg.eigh(hk(k, n))
    u = V[:, :, 0]
    P = np.einsum('ki,kj->kij', u, u.conj())
    dP = (np.roll(P, -1, 0) - np.roll(P, 1, 0)) / (2 * (2 * np.pi / nk))
    trg = 0.5 * np.real(np.einsum('kij,kji->k', dP, dP)).mean()
    # Wannier spread of (|A_0> + |B_n>)/sqrt2 : weight 1/2 at 0 and at n
    spread = 0.5 * 0 ** 2 + 0.5 * n ** 2 - (0.5 * 0 + 0.5 * n) ** 2
    print(f"   n = {n}:  <tr g> = {trg:.5f}   n^2/4 = {n**2/4:.5f}   "
          f"dimer spread <r^2>-<r>^2 = {spread:.5f}")
print("   -> tr g is EXACTLY the spread of a molecular orbital across a dimer")
print("      of length n. Intra-dimer geometry, not inter-cell coherence.")

# --------------------------------------------- 4. minimal metric = the fix
print("\n4. THE MINIMAL METRIC (the project's own v4 correction) SEES THIS")
print("   Shift convention u_a -> exp(+i k d_a) u_a, and seeding from a Wannier")
print("   centre r_a needs d_a = -r_a (REVIEW_README, 'conventions that have")
print("   bitten before'). B sits at r_B = n, so d_B = -n:")
for n in [1, 2, 3]:
    E, V = np.linalg.eigh(hk(k, n))
    u = V[:, :, 0].copy()
    u[:, 1] *= np.exp(-1j * k * n)         # d_B = -r_B, r_B = n
    P = np.einsum('ki,kj->kij', u, u.conj())
    dP = (np.roll(P, -1, 0) - np.roll(P, 1, 0)) / (2 * (2 * np.pi / nk))
    trg = 0.5 * np.real(np.einsum('kij,kji->k', dP, dP)).mean()
    print(f"   n = {n}:  M_naive = {n**2/4:.5f}   M_shifted = {trg:.3e}")
print("   -> zero to machine precision. The whole metric was gauge, removable")
print("      by placing the orbitals at their actual positions -- which is")
print("      exactly the non-gauge-invariance the author corrected in harvest")
print("      v4 and then never applied back to C4.")
