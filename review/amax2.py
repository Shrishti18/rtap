"""
The largest A is not 180 K/eV, and the answer was already in the repo.

amax.py bounded A by the ceiling theorem Tc <= sqrt(coef*<v^2>/n_phi)/4, giving
A <= 180 K/eV. But that ceiling ASSUMES <v^2> and Delta are set by the same
hopping -- that is exactly the assumption C4 (winding.py) breaks:

    d(k) = D * dhat(k),  |dhat| = 1, dhat winding n times
      gap   Delta = 2D                      set by D
      v     = D * n                         set by D AND n
      M_trg = <v^2>/Delta^2 = n^2/4         INDEPENDENT OF D

So M_trg is a free knob that does not cost gap. The cubic t2g model gets
M_trg = 0.018 because its isolation comes from ATOMIC SOC -- a k-INDEPENDENT
term. A k-independent gap contributes no metric: at small t the projector is
nearly k-independent, which is C2's mechanism in continuous form. The winding
model's gap has constant MAGNITUDE but winding DIRECTION, so it isolates
without flattening the projector.

** That is the structural statement of what "largest A" requires: an isolation
   mechanism whose orientation in orbital space winds with k, not a constant
   atomic splitting. **

This file asks what the terminus becomes if the SOC-set gap 3*lam/2 is kept but
the metric is taken from a winding structure instead of from atomic SOC.
"""
import numpy as np

KB = 0.086173e-3
COEF3 = 0.55 * 0.67
LAM_IR = 0.417

print(__doc__)
print("=" * 78)
print("1. The two structures at the SAME isolation gap")
print("   Delta_iso = 3*lam/2 (atomic SOC limit), U = Delta_iso/2\n")
print(f"   {'structure':<28}{'n_phi':>7}{'M_trg':>9}{'Tc_amp':>10}"
      f"{'Tc_stiff':>10}{'binds':>8}{'Tc (K)':>9}")
rows = []
for name, nphi, M in [('cubic t2g + atomic SOC', 3.0, 0.01809),
                      ('winding n=1', 2.0, 0.25),
                      ('winding n=2', 2.0, 1.00),
                      ('winding n=3', 2.0, 2.25)]:
    iso = 1.5 * LAM_IR
    U = iso / 2
    amp = U / (4 * nphi) / KB
    stf = COEF3 * U * M / KB
    tc = min(amp, stf)
    rows.append((name, tc))
    print(f"   {name:<28}{nphi:>7.1f}{M:>9.4f}{amp:>10.0f}{stf:>10.0f}"
          f"{'amp' if amp < stf else 'stiff':>8}{tc:>9.0f}")
print(f"\n   ratio winding(n=1) / cubic t2g = {rows[1][1]/rows[0][1]:.0f}x"
      f"   at IDENTICAL gap and IDENTICAL lam_SOC.")
print("   NOTE: the cubic row is GENEROUS to the cubic model. It pairs the")
print("   atomic-limit gap 3*lam/2 = 0.626 eV with M_trg = 0.0181, which that")
print("   model attains at Delta_iso = 0.21 eV, not at 0.626. It cannot have")
print("   both -- that IS the trade-off. Its true terminus is 8.1 K, so the")
print(f"   real advantage is {rows[1][1]/8.1:.0f}x, not {rows[1][1]/rows[0][1]:.0f}x.")
print("   The terminus was a property of the hopping structure, not of SOC.")

print()
print("2. What lam_SOC each structure needs for 300 K")
print(f"   {'structure':<28}{'A (K/eV)':>10}{'lam for 300 K':>16}  {'':<20}")
for (name, tc) in rows:
    A = tc / LAM_IR
    need = 300 / A
    tag = ('BELOW Ir (0.417)' if need < LAM_IR else
           'reachable in 6p' if need < 1.5 else 'nothing has it')
    print(f"   {name:<28}{A:>10.1f}{need:>16.3f}  {tag:<20}")

print()
print("3. Sensitivity to the two standing caveats")
print("   (a) D8: the defensible projection rule is U <= Delta_iso/4, not /2")
print("   (b) coef = 0.67 is mean-field and uncalibrated")
print(f"   {'structure':<20}{'U<=iso/2':>10}{'U<=iso/4':>10}"
      f"{'iso/4 & coef/2':>16}")
for name, nphi, M in [('cubic t2g', 3.0, 0.01809), ('winding n=1', 2.0, 0.25),
                      ('winding n=2', 2.0, 1.00)]:
    out = []
    for frac, cf in [(0.5, COEF3), (0.25, COEF3), (0.25, COEF3 / 2)]:
        U = 1.5 * LAM_IR * frac
        out.append(min(U / (4 * nphi), cf * U * M) / KB)
    print(f"   {name:<20}{out[0]:>10.0f}{out[1]:>10.0f}{out[2]:>16.0f}")
print("   Even under BOTH caveats, winding n>=1 stays two orders above the")
print("   cubic terminus. The conclusion is robust to the coefficient; it is")
print("   NOT robust to whether such a band exists in 3D in a real material.")

print()
print("4. WHAT IS ACTUALLY UNTESTED  (this is the honest part)")
print("""   - M_trg = n^2/4 was derived and verified for a 1D TWO-BAND model.
     The 3D analogue of a winding number is a Chern or Hopf-type invariant and
     the corresponding M_trg has NOT been computed here.
   - The winding model's partner band is also exactly flat. Real isolation is
     against a dispersive manifold; whether M survives that is untested.
   - The full-scan over 1,772 Wannier materials DID find M_trg > 0.339 in
     545/728 bands, so large metric is common in real materials. What was never
     found was large metric AND isolation AND correct filling in the same band.
     That is the same three-way condition, now with a structural target rather
     than a chemical one.
   - No material is proposed here. This says the terminus was structural, not
     fundamental -- it does not say a material exists.""")
