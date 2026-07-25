"""
How large can A be?

The scaling theorem is settled: Tc_max = A * lam_SOC, A = 19.43 K/eV for the
cubic orbital-selective t2g model. A is a pure number fixed by the HOPPING
STRUCTURE, so the remaining question is how large A can be made.

That question is not open-ended. This project already proved a ceiling that no
geometry can exceed (ceiling.py, C3 + C6):

    Tc <= min( U/(4 n_phi) , coef * U * M_trg )      two caps
    U  <= Delta_iso / 2                              projection validity
    M_trg <= <v^2> / Delta_iso^2                     EXACT sum rule (C3)

  => Tc <= min( Delta_iso/(8 n_phi) , coef*<v^2>/(2 Delta_iso) )

The first rises with Delta_iso, the second falls. They cross at
Delta* = 2 sqrt(coef * n_phi * <v^2>), giving the geometry-independent ceiling

    Tc_ceiling = sqrt(coef * <v^2> / n_phi) / 4

For a cosine band <v^2> = 2 t^2, so Tc_ceiling = sqrt(coef/n_phi) * t / (2*sqrt(2)).

So A is bounded, and the gap between what the cubic t2g model ACHIEVES and what
the ceiling ALLOWS is the entire remaining question. This file measures it.

Caveats that apply to both sides equally, so the RATIO is more robust than
either absolute number: coef = 0.67 is mean-field and uncalibrated; the ceiling
assumes the sum rule is saturated AND Delta_iso can be tuned to Delta*
independently, which no specific lattice is guaranteed to allow.
"""
import numpy as np
import jeff

KB = 0.086173e-3
COEF3 = 0.55 * 0.67
NPHI = 3.0                       # j=1/2 doublet, orbitally isotropic


def achieved(lam=0.417, nk=32, npts=61):
    """Optimise the stiffness cap over t at fixed lam, in the cubic model."""
    best = None
    for r in np.linspace(0.05, 0.30, npts):
        t = r * lam
        out, _ = jeff.analyse(t=float(t), lam=lam, tp=0.0, nk=nk, bands=(4, 5))
        d = out[0]
        if d['iso'] <= 0:
            continue
        U = d['iso'] / 2
        tc = COEF3 * U * d['M_trg'] / KB
        if best is None or tc > best[0]:
            best = (tc, t, r, d['iso'], d['M_trg'], d['nphi'])
    return best


print(__doc__)
print("=" * 78)
tc, t, r, iso, M, nphi = achieved()
lam = 0.417
print(f"1. What the cubic t2g model achieves, at lam = {lam} eV")
print(f"   t* = {t:.4f} eV   t*/lam = {r:.4f}   Delta_iso = {iso:.4f} eV")
print(f"   M_trg = {M:.5f}   n_phi = {nphi:.3f}   Tc = {tc:.2f} K")
print(f"   A = Tc/lam = {tc/lam:.2f} K/eV")
print(f"   in units of the hopping:  Tc/t = {tc*KB/t:.5f}")

print()
print("2. What the ceiling theorem allows, at the SAME t")
vsq = 2 * t ** 2                              # cosine band
tc_ceil = np.sqrt(COEF3 * vsq / NPHI) / 4
dstar = 2 * np.sqrt(COEF3 * NPHI * vsq)
print(f"   <v^2> = 2t^2 = {vsq:.6f} eV^2")
print(f"   Tc_ceiling = sqrt(coef*<v^2>/n_phi)/4 = {tc_ceil/KB:.1f} K"
      f"   (Tc/t = {tc_ceil/t:.5f})")
print(f"   optimal gap Delta* = {dstar:.4f} eV = {dstar/t:.3f} t")
print(f"   the model's actual Delta_iso = {iso:.4f} eV = {iso/t:.3f} t")

print()
print("3. THE HEADROOM")
gap = (tc_ceil / KB) / tc
print(f"   achieved Tc/t = {tc*KB/t:.5f}   ceiling Tc/t = {tc_ceil/t:.5f}")
print(f"   the cubic t2g structure reaches {100/gap:.1f}% of the ceiling.")
print(f"   FACTOR AVAILABLE: {gap:.2f}x")
print(f"   => A_max <= {gap * tc/lam:.0f} K/eV  (vs {tc/lam:.1f} achieved)")

print()
print("4. Why it falls short: the two conditions are not met at the same point")
print(f"   sum rule would allow M_trg <= <v^2>/Delta_iso^2 = {vsq/iso**2:.4f}")
print(f"   actual M_trg = {M:.5f}  ->  saturates {100*M/(vsq/iso**2):.1f}% of it")
print(f"   and Delta_iso/Delta* = {iso/dstar:.3f}  ->  the gap is "
      f"{'too small' if iso < dstar else 'too large'} by {dstar/iso:.2f}x")
print("   SOC fixes Delta_iso at 3*lam/2 in the atomic limit; the optimum wants")
print("   Delta* = 2.97 t. Those two agree only at one t, and at that t the")
print("   metric has not grown yet. THAT is the whole deficit.")

print()
print("5. What lam room temperature needs, at each level of A")
A_ach = tc / lam
A_max = gap * A_ach
print(f"   {'A (K/eV)':>10}{'lam for 300 K':>16}   {'reachable by':<28}")
for A, tag in [(A_ach, 'nothing (max is ~2-3 eV, 7p)'),
               (A_max, 'Bi/Pb 6p at 1.25-1.5 eV?')]:
    print(f"   {A:>10.1f}{300/A:>16.2f}   {tag:<28}")

print()
print("6. Ceiling Tc for real ions IF A could be pushed to its own bound")
print(f"   {'ion':<12}{'lam (eV)':>10}{'Tc at A=' + f'{A_ach:.0f}':>14}"
      f"{'Tc at A=' + f'{A_max:.0f}':>14}")
for name, lm in [('Ir4+  5d', 0.417), ('Pt4+  5d', 0.50), ('Au3+  5d', 0.60),
                 ('Tl+   6p', 0.64), ('Pb2+  6p', 0.91), ('Bi3+  6p', 1.25),
                 ('Bi    6p', 1.50), ('7p (superheavy)', 2.50)]:
    print(f"   {name:<12}{lm:>10.2f}{A_ach*lm:>14.0f}{A_max*lm:>14.0f}")
print()
print("   The right-hand column is NOT a prediction. It is what the project's")
print("   own ceiling theorem permits if a hopping structure could be found")
print("   that saturates the sum rule at the SOC-limited gap. No such structure")
print("   is known, and the cubic one reaches 11% of it.")
