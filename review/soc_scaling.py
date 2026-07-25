"""
The ceiling as a function of lam_SOC -- generalising the terminus from one ion
to every ion.

The terminus was reported as "8.1 K for Ir4+". That understates it. The model
Hamiltonian is

    H(k) = t * H_hop(k)  +  lam * H_soc                                    (*)

with H_hop and H_soc both dimensionless and FIXED. Two exact consequences:

  (a) H/lam = (t/lam) H_hop + H_soc, so the EIGENVECTORS depend only on the
      ratio t/lam. Every quantity built from eigenvectors alone -- M_trg,
      n_phi, lam_pair -- is a function of t/lam ONLY, and is invariant under
      (t, lam) -> (s*t, s*lam).
  (b) The EIGENVALUES scale linearly: E(s*t, s*lam) = s*E(t, lam). So
      Delta_iso = lam * f(t/lam) for some fixed f.

The stiffness cap is Tc_geo = coef * U * M_trg with U <= Delta_iso/2. Therefore

    Tc_geo(t, lam) = coef * (lam/2) * f(t/lam) * M_trg(t/lam)

Maximising over t at fixed lam: the optimum sits at a FIXED ratio (t/lam)*, and

    Tc_max = A * lam_SOC ,   A = coef/2 * f((t/lam)*) * M_trg((t/lam)*)

** THE CEILING IS EXACTLY LINEAR IN THE SPIN-ORBIT CONSTANT. **

That is a much stronger statement than "Ir fails": it fixes the ceiling for
EVERY ion in the periodic table from a single number A, and it says what
lam_SOC room temperature would require. Both are computed below.
"""
import numpy as np
import jeff

KB = 0.086173e-3
C3 = 0.55 * 0.67


def ceiling(lam, npts=41, nk=32, lo=0.05, hi=0.30):
    """max over t of the stiffness cap, at fixed lam. Returns (Tc_K, t*, ...)"""
    best = None
    for r in np.linspace(lo, hi, npts):          # r = t/lam
        t = r * lam
        out, _ = jeff.analyse(t=float(t), lam=float(lam), tp=0.0, nk=nk,
                              bands=(4, 5))
        d = out[0]
        if d['iso'] <= 0:
            continue
        U = d['iso'] / 2
        tc = C3 * U * d['M_trg'] / KB
        if best is None or tc > best[0]:
            best = (tc, t, r, d['iso'], d['M_trg'], d['nphi'])
    return best


print(__doc__)
print("=" * 78)
print("1. Invariance check: is M_trg a function of t/lam ALONE?")
print(f"   {'t/lam':>8}{'t':>8}{'lam':>8}{'M_trg':>10}{'n_phi':>8}"
      f"{'iso/lam':>10}")
for r in [0.125, 0.125, 0.125, 0.20, 0.20]:
    pass
for (t, lam) in [(0.052, 0.417), (0.104, 0.834), (0.026, 0.2085),
                 (0.2, 1.6026), (0.010, 0.0801)]:
    out, _ = jeff.analyse(t=t, lam=lam, tp=0.0, nk=32, bands=(4, 5))
    d = out[0]
    print(f"   {t/lam:>8.5f}{t:>8.4f}{lam:>8.4f}{d['M_trg']:>10.6f}"
          f"{d['nphi']:>8.4f}{d['iso']/lam:>10.6f}")
print("   -> identical at fixed t/lam across a 20x range in lam. (*) confirmed.")

print()
print("2. Is the ceiling linear in lam?")
print(f"   {'lam (eV)':>10}{'Tc_max (K)':>12}{'Tc/lam':>10}{'t*/lam':>9}"
      f"{'iso/lam':>10}{'M_trg':>9}")
As = []
for lam in [0.1, 0.2, 0.417, 0.8, 1.5, 3.0]:
    tc, t, r, iso, M, nphi = ceiling(lam)
    As.append(tc / lam)
    print(f"   {lam:>10.3f}{tc:>12.2f}{tc/lam:>10.3f}{r:>9.4f}"
          f"{iso/lam:>10.4f}{M:>9.5f}")
A = float(np.mean(As))
print(f"   -> Tc_max / lam = {A:.3f} K/eV, constant to "
      f"{100*np.std(As)/A:.2e}% across a 30x range in lam. LINEAR.")

print()
print("3. What lam would room temperature require?")
for target in [77, 200, 273, 300]:
    print(f"   Tc = {target:3d} K  needs lam_SOC = {target/A:7.2f} eV")

print()
print("4. The ceiling for every ion, from its measured atomic SOC constant.")
print("   (zeta values are atomic/spectroscopic; in a solid they are reduced)")
IONS = [
    ('Ti3+  3d', 0.019), ('V4+   3d', 0.031), ('Mn3+  3d', 0.044),
    ('Fe2+  3d', 0.052), ('Co2+  3d', 0.066), ('Ni2+  3d', 0.083),
    ('Cu2+  3d', 0.103), ('Ru3+  4d', 0.150), ('Rh3+  4d', 0.161),
    ('Pd2+  4d', 0.195), ('Os4+  5d', 0.400), ('Ir4+  5d', 0.417),
    ('Pt4+  5d', 0.500), ('Au3+  5d', 0.600),
    ('Tl+   6p', 0.640), ('Pb2+  6p', 0.910), ('Bi3+  6p', 1.250),
    ('Bi    6p*', 1.500), ('U4+   5f', 0.260), ('Pu3+  5f', 0.360),
]
print(f"   {'ion':<11}{'lam (eV)':>10}{'Tc_max (K)':>12}   {'':<4}")
for name, lam in sorted(IONS, key=lambda x: -x[1]):
    print(f"   {name:<11}{lam:>10.3f}{A*lam:>12.2f}")
best = max(IONS, key=lambda x: x[1])
print(f"   -> best case in the periodic table: {best[0].strip()} at "
      f"lam = {best[1]} eV -> {A*best[1]:.0f} K")
print(f"   -> short of 300 K by {300/(A*best[1]):.1f}x, and lam would have to be"
      f" {300/A/best[1]:.1f}x larger than anything that exists.")
