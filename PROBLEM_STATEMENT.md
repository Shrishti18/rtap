# Problem Statement — Room-Temperature Ambient-Pressure Superconductor

**Find, synthesize, and verify to the standard of `PROOF_STANDARD.md` a material
realizing a two-orbital flat band with non-removable quantum metric M_min ≥ 0.34
and an instantaneous on-site attraction U ≥ 158 meV on those same two orbitals,
doped strictly off half filling, stable at ambient pressure.**

Target: T_c = min(λ/4, 0.369·M_min)·U = **0.1644·U ≥ 25.9 meV (300 K)**.

---

## Epistemic discipline (syādvāda)

Every condition below carries a predication tag stating **in what respect** it
holds. No claim is made unconditionally. The respects used:

| tag | meaning here |
|---|---|
| **syād-asti** | established in the stated respect (verified computation or published measurement) |
| **syān-nāsti** | refuted in the stated respect |
| **syād-asti-nāsti** | holds in one stated respect, fails in another |
| **syād-avaktavyam** | presently undecidable — no measurement or valid computation exists |

A condition is **not satisfied** until it is *syād-asti in the respect of the
actual candidate material, measured* — model-level asti does not transfer.

---

## Part I — Band-structure conditions (solved at model level)

### C1′ — Pairing-weight optimum: maximize min(λ/4, 0.369·M_min)

Two orbitals per cell; **no constraint on their weight balance.**
*(The original "n_φ = 2, ⟨ρ⟩ = ½ each" is __syān-nāsti__ — disproved: λ = 1/n_φ
holds only for uniform weights, and the reference model achieves λ = 0.657 at
weights 0.35/0.65, better than the 0.5 that balance gives.)*

**Status: syād-asti** — in the respect of the reference model, λ = 0.657,
verified. **Syād-avaktavyam** in the respect of any real compound.

### C2 — Non-removable metric: M_min ≥ 0.34 (3D) / 0.84 (2D)

The metric **after minimization over orbital-position assignments** — M_naive is
inadmissible evidence. (A pure-winding model gave M_naive = n²/4 with M_min = 0
and D_s = 0 exactly: the whole metric was gauge, and the lattice was disconnected
dimers.)

**Status: syād-asti** in the respect of the models — e_g compass: M_min =
0.70–0.87, already minimized; higher-Chern reference: M_min = M_naive to 4
digits at C ≥ 2. **Syād-asti-nāsti** in the respect of real Wannier bands —
measured over 696: median retention M_min/M_naive = 0.957, 69.8% above 0.34, but
8/696 are C4-like collapses and 102/696 lose half. **Syād-avaktavyam** for the
target compound.

### C2′ — Connectivity (no folding)

Real-space hopping must occupy consecutive shells; a single n-th-shell bond
decouples the lattice (Read's obstruction forbids exactly-flat + finite-range +
C ≠ 0; exponential tails are the price and the protection).

**Status: syād-asti** in the model respect (shells 1–15, exponential decay).

### C3 — Isolation: Δ_iso ≥ 2U — **not a constraint**

Gap scales with |d|, metric depends only on d̂: verified M_trg = 0.28543
unchanged at gaps 0.21 → 100 eV. Any C2-passing model can be scaled to any gap.

**Status: syād-asti** (theorem-level in two-orbital models). Do not spend search
effort on gap size.

---

## Part II — Material conditions (open; these are the problem)

### C4 — Filling: **f ∉ {0, ½, 1}**, sharpened

E_F inside the band (727/728 real bands fail: **syān-nāsti** in the respect of
undoped known materials) — **and strictly off half filling**, because at f = ½
the attractive-U model has an exact pseudospin SU(2) degeneracy between
superconductivity and charge order; doping lifts it **toward SC**
([PRB 55, 1185](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.1185)).
The naive "ideally f = ½" of earlier drafts is **syān-nāsti** — half filling is
the one filling that must be avoided.

Geometry and filling are statistically independent (measured: joint 2 vs 0.46
expected), so this is engineered separately, and only within ~0.1–0.3 eV.

**Status: syād-avaktavyam** for any specific candidate until its DFT is done.

### C5 — Attraction: U ≥ 158 meV, instantaneous (ω₀ ≳ 2U), on the two orbitals

**Status: syād-asti-nāsti-avaktavyam**, and the respects matter:

- *asti*: electronic valence-skipping attraction is real and large in the Bi 6s
  case (~1.9 eV, BaBiO₃) — but there n_φ = 1 and tr g ≡ 0: right magnitude,
  wrong orbital, dead by theorem.
- *asti* (new, and the reason for this route): **the ambient-pressure charge
  order of Cs₂Au₂X₆ is itself direct evidence of negative U on Au d** — a
  disproportionated Au⁺/Au³⁺ lattice *is* the pair crystal of d⁹ → d¹⁰ + d⁸.
- *avaktavyam*: **no cRPA, spectroscopic, or transport determination of U for
  Au²⁺ on an e_g doublet exists.** The 158 meV requirement is compared against
  an analogy, not a measurement. This number is the project's F1-class gap:
  computable by one cRPA study.

### C6 — Competing orders: superconductivity must **win**, not merely be possible

The same U that pairs also crystallizes pairs. Requirements, each grounded:

1. **f ≠ ½ exactly** (degeneracy argument above — merged into C4).
2. **The known Au²⁺-halide family sits IN the competing phase at ambient
   pressure.** Cs₂Au₂X₆ (X = Cl, Br, I) is charge-ordered Au⁺/Au³⁺ at ambient;
   the homogeneous Au²⁺ state appears only at **~12.5 GPa** (valence + structural
   transition; metallization ~6 GPa)
   ([JACS 117, 1994](https://pubs.acs.org/doi/10.1021/ja00104a016),
   [J. Chem. Phys. 110, 9174](https://pubs.aip.org/aip/jcp/article-abstract/110/18/9174/528248/)).
   So the route runs *through* the CDW: the condition is a **doping or chemical
   substitution that suppresses charge order at P = 0** — precisely the
   Ba→K move that turned charge-ordered BaBiO₃ into superconducting BKBO.
3. Phase separation and magnetism checked at the same level (mean-field minimum;
   ED preferred).

**Status: syān-nāsti** at ambient in the known family (it *is* the competitor);
**syād-asti** under pressure (homogeneous valence exists); **syād-avaktavyam**
whether doping substitutes for pressure here — the BKBO precedent says it can.

### C7 — Disorder tolerance

Flat bands have no kinetic protection. Two graces, one grounded: disorder
destroys CDW **before** SC in the attractive model
([PRB 55, 1185](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.1185))
— disorder is partly on our side against C6; and geometric stiffness is
interaction-generated, not velocity-generated. Requirement: predicted D_s must
exceed the disorder broadening of the doped crystal.
**Status: syād-avaktavyam** quantitatively.

### C8 — Synthesizability at ambient

The family **exists** — Cs₂Au₂X₆, Rb analogues, (MA)₂Au₂I₆ — which is more than
any prior candidate in this project had. The condition: a member (or doping
x) whose **ambient-pressure** ground state is the homogeneous (or dynamically
fluctuating) valence state with the e_g doublet at E_F, thermodynamically stable
or metastable at 300 K.
**Status: syād-asti** that the chemistry exists; **syān-nāsti** that the needed
phase exists at ambient today; **syād-avaktavyam** for the doped series (not in
the literature we found — this is an open experimental question, stated as
such).

---

## Part III — Meta-conditions (the ones this project's history demands)

### C9 — Theory calibration **precedes synthesis**

`c = 0.67` is mean-field, calibrated against BdG — against itself. The one
mean-field-vs-exact check (sawtooth ED) **failed and was stood down**
(**syān-nāsti** in the respect of exact methods). Requirement: calibrate
`c` by ED/DMRG/QMC on a flat-band model **before** any material work; if the
true stiffness is k× below mean field, the required U is k× above 158 meV, and
the route may close on paper for the cost of a day of computing.
Additionally: every per-band descriptor uses a projector whose rank is selected
by the self-consistency rule (`screen/selector.py`), never by convention —
eight errors of one class stand behind this line.

### C10 — Survival prior and staged falsification (the honesty condition)

**The prior that any candidate passing C1′–C8 on paper survives to a verified
discovery is well under 10%.** This is a condition, not a mood: the search plan
must price it in.

- *Base rate, favourable respect (syād-asti)*: theory-led superconductor
  discovery has succeeded — H₃S and LaH₁₀ were **predicted, then measured** at
  203 K / ~250 K. Prediction-first is not fantasy.
- *Unfavourable respect (syān-nāsti)*: those successes used Eliashberg theory,
  quantitatively trusted for decades; ours is uncalibrated mean field (C9). And
  the ambient-pressure room-T claims of recent years (LK-99, retracted
  hydride/carbonaceous claims) failed exactly at the stages our
  `PROOF_STANDARD.md` gates.
- *Operational content*: pre-register kill criteria per stage **before** the
  stage runs (C9 kill: calibration collapse; DFT kill: no isolated doublet at
  E_F under any x; cRPA kill: U < 158 meV; synthesis kill: charge order robust
  to all attempted dopings; measurement kills: Stages 1–5). Proceed only while
  the expected information per unit cost remains positive. **Report every kill
  as a result** — this project's negatives (E_F gate, SOC ceiling, winding
  retraction) were its most reliable products.

---

## The decision chain, priced

| step | cost | kills on | tag today |
|---|---|---|---|
| 1. ED/QMC calibration of c (C9) | days, compute only | c ≪ 0.67 | syān-nāsti → must run first |
| 2. DFT of Cs₂Au₂X₆ doped series (C4, C8) | weeks, compute | no doublet at E_F at P=0 | syād-avaktavyam |
| 3. cRPA of U on Au²⁺ e_g (C5) | weeks, compute | U < 158 meV | syād-avaktavyam — also fills a genuine literature hole |
| 4. Competing-order ED on the fitted model (C6) | weeks, compute | CDW wins at all accessible f | syād-avaktavyam |
| 5. Synthesis of the doped series (C8) | months, lab | no ambient homogeneous phase | syān-nāsti today |
| 6. `PROOF_STANDARD.md` Stages 0–5 | months–years, lab | any gate | — |

Steps 1–4 are pure computation and can kill the route for under a month of
effort. **No synthesis before all four pass.**

---

## One-sentence form

> Dope a Cs₂Au₂X₆-class Au²⁺ halide off its charge order at ambient pressure so
> that the valence-skipping e_g doublet — whose model band structure already
> satisfies every geometric condition (λ = 0.657, M_min = 0.70–0.87, exactly
> flat, gap-free-of-metric-cost, connected) — sits at a filling strictly inside
> (0,1) and off ½, verify U ≥ 158 meV by cRPA and the stiffness constant by
> exact diagonalization first, and accept nothing as discovered short of
> `PROOF_STANDARD.md` Stages 0–5.

Every "syād-avaktavyam" above is an open task. There are six. That is the
problem.

---

## Sources

- SC–CDW degeneracy at half filling and its lifting by doping; disorder kills
  CDW before SC: [Phys. Rev. B 55, 1185 (1997)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.1185),
  [arXiv:cond-mat/9606038](https://arxiv.org/pdf/cond-mat/9606038)
- Cs₂Au₂X₆ P–T phase diagram, Au⁺/Au³⁺ → Au²⁺ valence transition ~12.5 GPa,
  metallization: [JACS 117 (1994)](https://pubs.acs.org/doi/10.1021/ja00104a016),
  [J. Chem. Phys. 110, 9174 (1999)](https://pubs.aip.org/aip/jcp/article-abstract/110/18/9174/528248/),
  [Physica B+C 1980 (electronic state under pressure)](https://www.sciencedirect.com/science/article/abs/pii/0375960180905101),
  [high-pressure XRD to 18 GPa](https://www.researchgate.net/publication/255023734)
- Hybrid mixed-valence gold iodide (MA)₂Au₂I₆:
  [Inorg. Chim. Acta (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0020169317303924)
- Internal: `PROOF_STANDARD.md`, `review/eg_triangular.py`, `review/higherC.py`,
  `review/winding_stiffness.py`, `screen/selector.py`, `screen/minmetric.csv`,
  `screen/FINDINGS.md`
