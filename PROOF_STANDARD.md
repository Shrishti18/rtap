# Problem statement: conditions under which an RTAP superconductor is PROVEN

**No theoretical condition set can establish this.** C1′–C5 are conditions on a
*model*. Even if a real material satisfied all five exactly, the T_c formula
producing "300 K" is mean-field, calibrated only against BdG — i.e. against
itself — and the one intended mean-field-vs-exact check did not work. Theory can
select a candidate and predict a number. It cannot certify a superconductor.

What follows is the standard that can. It is deliberately the one that recent
failed claims (LK-99, several hydride reports) would not have passed.

---

## Stage 0 — Pre-registration (what makes it a test rather than a claim)

Before measuring, record publicly:

- **P0.1** Composition, space group, lattice parameters, synthesis route.
- **P0.2** The predicted T_c **with its error bar**, from the theory, and the
  quantities it depends on (λ, M_min, U, Δ_iso) with their computed values.
- **P0.3** The predicted **H_c2(0)**, coherence length ξ, and gap ratio 2Δ/k_BT_c.
  A flat-band superconductor is a specific prediction, not a generic one — for
  T_c = 300 K the Pauli limit alone gives μ₀H_c2 ≳ 560 T and ξ ≲ 0.8 nm. If the
  measured ξ is tens of nm, the mechanism is not this one even if it superconducts.
- **P0.4** What result would **falsify** the mechanism, as distinct from
  falsifying superconductivity.

Without P0.2–P0.4 a positive result is a discovery of *something*, not
confirmation of *this*.

---

## Stage 1 — Transport: necessary, nowhere near sufficient

- **T1.1 Four-point resistivity, R = 0** within instrument noise, with the noise
  floor stated in Ω and the voltage resolution in nV. "Resistance drops sharply"
  is not R = 0.
- **T1.2 Current independence.** R = 0 must hold across ≥ 2 decades of measuring
  current, and J_c(T) must be measured. A percolating filament shows
  current-dependent apparent R.
- **T1.3 Both sweep directions, ≥ 3 thermal cycles**, no hysteresis beyond
  instrument lag. Excludes a structural or chemical transition masquerading as
  T_c — this is precisely what the LK-99 resistivity feature turned out to be.
- **T1.4 Multiple contact geometries and ≥ 3 independent samples.**
- **T1.5 Magnetic field suppression.** T_c must shift down monotonically with
  applied field, and the shift must be consistent with the H_c2 of P0.3. A
  feature that does not move in field is not superconductivity.

**Transport alone can never prove it.** A single connected filament of any
superconducting impurity phase gives R = 0 through a bulk sample.

---

## Stage 2 — Magnetism: the load-bearing evidence

- **M2.1 Meissner expulsion, field-cooled**, distinguished from zero-field-cooled
  shielding. Both curves reported. FC expulsion is the hard one to fake; ZFC
  shielding alone is consistent with a hollow shell of any origin.
- **M2.2 Volume fraction ≥ 50%**, from χ_v after a stated demagnetization
  correction with the sample geometry given. **This is the single most important
  number**, because it separates a bulk phase from a trace impurity. Report the
  raw moment, the mass, the density, and the correction — not just χ.
- **M2.3 Full hysteresis loops** M(H) at several T, showing type-II behaviour,
  H_c1 and H_c2, and flux pinning consistent with J_c from T1.2.
- **M2.4 Background and artifact control.** Empty-holder subtraction, and an
  explicit test against ferromagnetic-impurity and instrument-artifact
  explanations. A diamagnetic-looking signal at the 10⁻⁴ level of a sample with
  ppm magnetic impurities is not evidence.
- **M2.5 Levitation is not evidence.** Partial, orientation-dependent levitation
  is achieved by diamagnets and by ferromagnetic pinning. It has no place in the
  argument.

---

## Stage 3 — Thermodynamics: proof that it is a bulk phase transition

**This is the stage that separates a real discovery from every recent failed
claim, and it is the one most often skipped.**

- **H3.1 Specific-heat jump at T_c.** A bulk second-order transition *must* show
  one. Report ΔC/γT_c; BCS weak coupling gives 1.43, strong coupling 2–3. A
  flat-band condensate may deviate — **state the predicted value in P0.2** and
  compare.
- **H3.2 Entropy balance.** ∫(C_s − C_n)/T dT = 0 across the transition. This
  cannot be faked by an impurity phase at the few-percent level.
- **H3.3 The jump must scale with the superconducting volume fraction from M2.2.**
  If 60% of the sample is superconducting magnetically but the specific-heat jump
  corresponds to 2%, the two measurements are describing different things.

**A claim without a specific-heat jump is not established, regardless of how
clean the resistivity and magnetization look.**

---

## Stage 4 — Spectroscopy: that the gap exists and is what was predicted

- **S4.1 A gap in the density of states** by tunnelling/STS or ARPES, opening at
  T_c, with 2Δ/k_BT_c reported against the P0.3 prediction.
- **S4.2 Gap symmetry and structure** consistent with the proposed pairing.
- **S4.3 For this mechanism specifically:** the superfluid stiffness measured
  directly (μSR penetration depth, or optical spectral weight via the
  Ferrell–Glover–Tinkham sum rule) and compared against `D_s = c·U·M_min`. **This
  is the measurement that tests the mechanism rather than the phenomenon** — it
  is the only route to confirming that quantum geometry, not conventional
  electron–phonon coupling, is carrying the stiffness.

---

## Stage 5 — Independent reproduction

- **R5.1 Two independent groups**, independent synthesis, independent
  instruments, reproducing Stages 1–3.
- **R5.2 Same-sample verification** — the identical specimen measured in a second
  laboratory.
- **R5.3 Full data and synthesis release** sufficient for a third party to
  reproduce without contacting the authors.

---

## The standard, stated once

> **Proven** = R = 0 (Stage 1) **AND** ≥ 50% Meissner volume fraction (Stage 2)
> **AND** a specific-heat jump with entropy balance matching that volume fraction
> (Stage 3) **AND** a spectroscopic gap closing at T_c (Stage 4) **AND**
> independent reproduction (Stage 5) — all at ambient pressure, T ≥ 300 K.
>
> **Mechanism confirmed** additionally requires the measured superfluid stiffness
> to match `c·U·M_min` within its error bar (S4.3), and the measured ξ and H_c2
> to match P0.3.

Any subset is a candidate. Only the conjunction is a discovery.

---

## What C1′–C5 are for

They are **Stage 0**: they select what to synthesize and fix the numbers the
experiment must hit. Their value is that they make the claim *falsifiable in
advance* — P0.2 and P0.3 are predictions that can be wrong. That is the whole
contribution of the theory, and it is worth having.

Nothing before Stage 5 justifies the word "discovered."

---

## Standing caveats that survive into Stage 0

- `c = 0.67` is mean-field and **uncalibrated against exact methods**. The
  sawtooth ED validation did not work. Every predicted T_c is a ceiling on a
  ceiling. **Calibrating this against exact diagonalization or QMC is the
  cheapest thing that could still close the route, and it should be done before
  any synthesis.**
- `U ≤ Δ_iso/2` permits 50% relative corrections; `Δ_iso/4` is defensible and
  halves every T_c.
- C1′–C5 constrain band geometry and one interaction constant. They say nothing
  about **competing orders** — a flat band with strong attraction is equally
  unstable to charge density waves, phase separation, and magnetism. Nothing in
  the conditions shows superconductivity wins. That analysis is not done.
- No cRPA or spectroscopic negative-U value exists for Au²⁺ on an e_g doublet.
  The 158 meV requirement is compared against a **bismuthate analogy**, not a
  measurement.
