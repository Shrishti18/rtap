# Brief: RTAP superconductor search

It's time to do the RTAP superconductor breakthrough. I believe there is a full
**family** of RTAP superconductors out there that can actually be built, and
that survives being built.

Everything below is what we have already established, how we established it, and
how much to trust each item. **Check us — do not trust us.** Several items here
are corrections to things we ourselves got wrong, and the corrections are the
most valuable part.

---

## 0. The frame — read this before anything else

The target formula is

    Tc = min( U·λ/4 , 0.369·U·M_min )

- `λ = λ_max(A)`, `A_ab = ⟨ρ_a(k) ρ_b(k)⟩_BZ`, `ρ_a` = orbital weight of the band.
- `M_min` = quantum metric **minimized over orbital-position assignments**.
- For 300 K you need `Tc ≥ 25.9 meV`.

**Both caps are mean-field.** Every number is a ceiling on a ceiling. A ceiling
above 300 K licenses `Tc ≤ ceiling`, **not** `Tc ≥ 300 K`. This inversion is the
single most common error in this problem — do not make it.

---

## 1. What we derived, with method and confidence

### 1.1 `λ = 1/n_φ` holds ONLY for uniform orbital weights — do not force balance
**How:** computed `λ_max(⟨ρ_aρ_b⟩)` directly on candidate models.
**Result:** a model with weights 0.35/0.65 gave **λ = 0.657**, *better* than the
0.5 that perfectly balanced weights give.
**Why believe it:** `λ_max ≥ 1/n_φ` always, by the Rayleigh quotient; equality
only at uniform ρ. Verified numerically.
**Consequence:** the condition is *maximize* `min(λ/4, 0.369 M_min)`, not
"n_φ = 2 with ⟨ρ⟩ = ½". We spent effort optimizing the wrong proxy and lost 55%.
**Confidence: high (theorem + computation).**

### 1.2 Isolation gap and quantum metric are INDEPENDENTLY scalable
**How:** for any two-orbital band, `H = d₀𝟙 + d·σ` exactly, so `tr g = ¼|∂d̂|²`
(direction only) while `gap = 2|d|` (magnitude only).
**Result:** `M_trg = 0.28543` **unchanged** at gaps of 0.21, 1.0, 5.0, 50.2 and
100.0 eV, and at exactly zero bandwidth.
**Why believe it:** algebraic identity, confirmed to 5 digits across a 500× range.
**Consequence:** **do not hunt large-gap materials.** The famous "large metric
needs a small gap" trade-off is NOT a constraint for two-orbital systems. It is
real only when the gap comes from a k-*independent* term (see 1.4).
**Confidence: high.**

### 1.3 `M_naive` is worthless evidence. Only `M_min` counts.
**How:** we built a model with `d̂ = (cos nk, sin nk, 0)` giving `M_naive = n²/4`
— apparently a huge metric, gap-independent. Then three checks:
1. Fourier transform to real space: `H₁₂(R) = D·δ_{R,n}` — **one bond**, so the
   lattice decomposes into **decoupled dimers**.
2. Flux threading an N-site ring: energy spread **2×10⁻¹⁶**, i.e. `D_s = 0`
   exactly (each dimer's phase is removable by a local gauge).
3. Orbital shift to the true Wannier centre `d_B = −n`: `M_min = 3×10⁻²⁷`.
**Result:** the entire "metric" was the Wannier spread of a molecular orbital
across a disconnected dimer — pure gauge, zero current.
**Why believe it:** three independent methods agree.
**Consequence:** **any claimed metric must be reported after minimization over
orbital positions, and the lattice must be checked for connectivity in real
space.** This retracted our own headline result. Do not skip it.
**Confidence: high.**

### 1.4 A k-INDEPENDENT gap isolates without contributing metric
**How:** cubic orbital-selective t₂g model + atomic SOC, scanned over t/λ.
**Result:** in the isolated regime `M_trg ∝ t²` **exactly** (0.0007 / 0.0027 /
0.0060 at t = 0.01 / 0.02 / 0.03) and never exceeds 0.04. The model reaches only
**10.9%** of its own theoretical ceiling.
**Why:** a constant on-site splitting leaves the projector nearly k-independent.
This is the continuous version of "n_φ = 1 ⟹ tr g ≡ 0".
**Confidence: high.**

### 1.5 The atomic-SOC route has a hard ceiling: `Tc ≤ 19.4 K/eV × λ_SOC`
**How:** `H = t·H_hop + λ·H_soc` with both terms dimensionless implies
`H/λ = (t/λ)H_hop + H_soc`, so **eigenvectors depend only on t/λ** (hence M and
n_φ are functions of that ratio alone) while **eigenvalues scale linearly**
(hence `Δ_iso = λ·f(t/λ)`). The optimum therefore sits at a fixed ratio and
`Tc_max ∝ λ`.
**Result:** `Tc_max/λ = 19.43 K/eV`, constant to **five digits** over a 13× range
in λ. Best case in the periodic table is Bi 6p at ζ ≈ 1.5 eV → **29 K**. Room
temperature would need λ ≈ 15 eV; nothing has it.
**Why believe it:** exact scaling argument, numerically confirmed.
**Consequence:** **the entire atomic-SOC family is closed, for every ion.** Do
not propose j_eff = 1/2 iridates or any relative. This is a solved negative.
**Confidence: high.**

### 1.6 Lattice distortion buys isolation by destroying orbital multiplicity
**How:** added a tetragonal crystal field to the t₂g model and scanned.
**Result:** isolation rose, but `n_φ` fell **3.00 → 1.05**, `λ_pair` rose
1/3 → 0.95, and `M_trg` fell **170×**.
**Consequence:** distortion is not a route to isolation. It trades away exactly
the thing being protected.
**Confidence: high.**

### 1.7 The Fermi-level condition kills essentially everything, and is INDEPENDENT of geometry
**How:** screened all 1,772 JARVIS-WTB Wannier Hamiltonians; 728 candidate bands
in 179 materials survived the geometric filters.
**Result:** **727 of 728 fail the E_F condition.** Correlation of geometry with
"E_F inside" is +0.124 (gap uniformity), −0.087 (M_trg), −0.019 (1/n_φ) — i.e.
**statistically independent**; joint occurrence 2 observed vs 0.46 expected.
**Consequence:** search geometry first, then engineer E_F separately. But doping
moves E_F only ~0.1–0.3 eV; a 1 eV shift is a different compound, not doping.
**Confidence: high (large sample).**

### 1.8 Real Wannier bands KEEP their metric under gauge minimization
**How:** computed `M_min` in the physical Wannier-centre gauge (from
`wannier90.wout`, converted to fractional with the POSCAR lattice) for 696 bands.
**Result:** median `M_min/M_naive = 0.957`; **only 8 of 696** collapse below 0.10.
**Consequence:** the dimer pathology of 1.3 is a *constructed* pathology, not the
generic case. Real materials' metric is largely physical. Good news.
**Confidence: high.**

### 1.9 Large metric is common; it is NOT the scarce ingredient
**How:** same screen, using a rank-selection rule that compares the intra-group
splitting against the pairing scale (see 4.1).
**Result:** **43.5%** of candidate bands exceed the 3D threshold 0.34 (median M
= 0.279). *(This is a correction: we previously reported 77% and median 0.97,
using single-band descriptors paired with a group isolation gap — an
inconsistency that inflated it 3.5×.)*
**Consequence:** the scarce ingredients are **n_φ, filling, and attraction** —
not metric.
**Confidence: medium-high (method-dependent; sensitivity band 43–49%).**

### 1.10 No published cRPA U exists for any Ir non-oxide
**How:** literature search. Every ab-initio screened-interaction value for
Ir⁴⁺ 5d is an oxide; non-oxide work fits U empirically to a gap.
**Our estimate:** U(Ir⁴⁺, telluride) ≈ **0.95 eV**, from an oxide anchor of
1.9 eV times a heavy-ligand factor of 0.5 measured independently in three
families (Nb/Ta MX₂ across S/Se/Te; RuX₃ across Cl/Br/I; Fe pnictide vs
chalcogenide).
**Confidence: medium (extrapolation, explicitly not a calculation).** This is a
genuine hole in the literature and a publishable single calculation on its own.

---

## 2. What the literature establishes (cite, don't re-derive)

- **At exactly half filling, superconductivity and charge order are EXACTLY
  degenerate** in the attractive-U model, by a hidden SU(2) pseudospin symmetry.
  Doping lifts the degeneracy **toward superconductivity**.
  → *f = ½ is the one filling that must be AVOIDED.*
  [PRB 55, 1185](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.1185)
- **Disorder destroys the charge order before it destroys superconductivity** —
  disorder is partly an ally against the competing phase. *(Same reference.)*
- **Flat-band `D_s ∝ |U|` is confirmed beyond mean field** by QMC, DMRG and DMFT;
  BdG-vs-DMRG agreement for pairing and superfluid weight is good even in 1D
  where fluctuations are strongest. *(This partially closes our biggest caveat —
  see 3.1.)*
- **Read's obstruction:** exactly flat + strictly finite-range + C ≠ 0 is
  impossible. Something must give; giving up finite range costs only exponential
  tails, and those tails are what keep the lattice connected (cf. 1.3).
- **Cs₂Au₂X₆ (X = Cl, Br, I)** is charge-ordered Au⁺/Au³⁺ at ambient pressure and
  reaches the homogeneous Au²⁺ state only at **~12.5 GPa**, metallizing ~6 GPa.
  [JACS 117 (1994)](https://pubs.acs.org/doi/10.1021/ja00104a016),
  [JCP 110, 9174](https://pubs.aip.org/aip/jcp/article-abstract/110/18/9174/528248/)
- **Cs₄Au₃Cl₁₂** stabilizes Au²⁺ at ambient but as a **lattice-polaron crystal** —
  the ions relax. A "clamped-ion" attraction is fictitious if they do.
  [Nature Chem (2023)](https://www.nature.com/articles/s41557-023-01305-y),
  [arXiv:2602.11572](https://arxiv.org/abs/2602.11572)
- **Ambient-pressure record** is ~133–138 K (Hg-1223). A 2026 pressure-quenched
  151 K claim exists; **we have not verified it** — read it before citing.

---

## 3. What is open — this is the actual work

### 3.1 Is `c = 0.67` right? *(partially closed, finish it)*
The stiffness prefactor is mean-field, calibrated against BdG — i.e. against
itself. Our one intended mean-field-vs-exact check failed and was abandoned.
Literature says `D_s ∝ |U|` **is** confirmed by QMC/DMRG/DMFT, so the *scaling*
is safe; the **coefficient** is not pinned. **Do this first: it costs days of
compute and can close the whole route.** If the true stiffness is k× below mean
field, every required U rises k×.

### 3.2 Is there a real attraction, on the right orbitals?
Never verified in any material in this project. Fullerides give 17–20 meV (too
weak by ~10×). Bismuthates give ~1.9 eV but on Bi 6s where n_φ = 1 and
`tr g ≡ 0` **identically** — right magnitude, wrong orbital, dead by theorem.
**The open question:** a valence-skipping ion whose active manifold is a
*doublet*, so the metric isn't forced to zero. Au²⁺ (d⁹ → d¹⁰ + d⁸) is the
obvious candidate and **no cRPA or spectroscopic U exists for it.**

### 3.3 Does superconductivity beat charge order?
Unaddressed. A flat band with strong attraction is maximally unstable to CDW,
phase separation and magnetism too. Requires: f ≠ ½ (2.1), plus an actual
competing-order calculation on the fitted model.

### 3.4 Does any real crystal realize any of this?
No DFT has been run on any candidate. Every model in this project is a model.

---

## 4. Traps — we hit these eight times; each cost real work

### 4.1 Rank mismatch (**the big one — eight occurrences**)
A quantity must be computed on a projector whose **rank matches the degeneracy
structure** of what it projects onto, and *which rank is decided by Tc, not by
convention*:
- **manifold treatment** valid iff intra-group splitting `δ ≲ k_B·Tc`;
- **single-band** valid iff `δ ≫ k_B·Tc`, using the band's **own** nearest gap.
Mixing one treatment's `U` with the other's `M` inflated one candidate from 9 K
to 186 K.
**And the deep reason it kept recurring:** the single-band validity test is
**vacuous** — its own U-cap enforces `k_BTc ≤ U/8 ≤ d_near/8 ≤ δ/8` identically.
Any rule that asks "which test passes" inherits a tautology as a bias. *Check
whether your validity condition is entailed by your own definitions before you
use it to select.*

### 4.2 Non-periodic gauges
Check `max|H(k) − H(k+2π)| ≈ 0` before finite-differencing. A `2t·cos(k/2)`
convention is **not** 2π-periodic and injects a false discontinuity at the BZ
boundary. Tell-tale: the metric peaks where the gap is **largest**.

### 4.3 Rank-1 projectors at a degeneracy
For a Kramers pair the rank-1 projector is gauge-arbitrary: `⟨tr g⟩` diverges as
nk² and shifts under random U(2) rotation. Use the degenerate manifold — *unless*
a good quantum number (e.g. conserved `s_z`) labels the members, in which case
rank-1 is well defined.

### 4.4 `M` vs `M_trg`
All thresholds are on `M_trg = ⟨tr g⟩`. A separately-normalized `M` scaled by
`(2π)^(d−1)` differs by **39.5× in 3D**.

### 4.5 A ceiling is not a prediction
See §0.

---

## 5. Standard of proof

Nothing counts as discovered short of: **R = 0** (four-point, current-independent,
field-suppressed) **AND ≥50% Meissner volume fraction** (field-cooled, with
demagnetization correction) **AND a specific-heat jump with entropy balance
matching that volume fraction** **AND a spectroscopic gap closing at Tc** **AND
independent reproduction** — all at ambient pressure, T ≥ 300 K.

The specific-heat jump is the stage that separates real discoveries from the
recent failed claims, and it is the one most often skipped. Levitation is not
evidence.

For the **mechanism** additionally: measured superfluid stiffness (μSR or
Ferrell–Glover–Tinkham sum rule) matching `c·U·M_min`, and measured ξ and H_c2
matching the prediction. Note a 300 K superconductor implies
μ₀H_c2 ≳ 560 T and ξ ≲ 0.8 nm from the Pauli limit alone — if the measured
coherence length is tens of nm, the mechanism is wrong even if it superconducts.

---

## 6. Honest prior

**P(a candidate that passes every paper condition survives to a verified
discovery) is well under 10%.** Price that in: pre-register kill criteria per
stage, run the cheap computational kills first (§3.1 → §3.4), and **report every
kill as a result.** Our negatives — the E_F gate, the SOC ceiling, the dimer
retraction — were this project's most reliable products.

Prediction-first discovery is real: H₃S and LaH₁₀ were predicted, then measured.
But those used Eliashberg theory trusted quantitatively for decades. Ours is
uncalibrated mean field until §3.1 is done.

---

## 7. What I want

A **family**, not a single compound — materials that can actually be synthesized
and survive synthesis. So prefer:

- **robustness over optimality** — a design tolerating ~15% parameter error beats
  one needing exact cancellations;
- **known chemistry over inverse design** — a family that exists and needs doping
  beats a hypothetical structure with ideal hoppings;
- **isotropic stiffness** — check `g_xx`, `g_yy`, `g_zz` separately. A model with
  `g_yy = 0` is quasi-2D and its 3D Tc is set by the weak direction;
- **anything that survives §3.1–§3.4.**

Do whatever is needed. Tell me what dies and why — that is worth as much as what
lives.
