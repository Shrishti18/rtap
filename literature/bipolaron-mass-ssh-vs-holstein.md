# Bipolaron effective mass: bond/Peierls/SSH vs. Holstein coupling

Literature pull compiled 2026-07-25. Numbers are quoted or read directly from the cited sources
(figure/table/equation locations given where relevant); anything not explicitly stated in a source
is flagged as such rather than estimated. Free PDFs are mirrored in `pdfs/` (arXiv preprints /
Gold-OA journal copies only — see per-entry OA status).

---

## Headline answer: polynomial or exponential at U_eff ~ 2-4 ω0?

**Exponential — never polynomial — but exponential in √λ rather than in λ.**

- The anchor paper (Zhang et al., PRX 2023) never fits a power law anywhere. It describes growth
  *below* an "optimal" coupling λ_op as merely "weak-to-moderate," and states explicitly that
  above λ_op bond-Peierls bipolarons "enter a regime of exponential mass enhancement" (Sec. III).
  Their own numbers at U/t=8, t/Ω=10/3 already show it: m\*/m₀ ≈ 2 (λ≈0.3) → 5.4–6.3 (λ=0.5,
  two independent fits) → ≈15–20 (λ≈0.6) — a ~2–3× jump per Δλ≈0.1–0.15, i.e. exponential-looking
  behavior even inside the nominally "light" window, not mild power-law growth.
- Sous's 2026 perspective (arXiv:[2605.16625](https://arxiv.org/abs/2605.16625)) makes this
  quantitative with an instanton-derived asymptotic scaling:
  **Tc ~ exp[−c(t/Ω)√λ]** (bond-Peierls) vs. **Tc ~ exp[−c′(t/Ω)λ]** (Holstein).
  Since Tc ∝ 1/m\* in their BEC treatment, this is **m\*_BP ~ exp(c(t/Ω)√λ)** — a *stretched
  exponential*: asymptotically much slower-growing than Holstein's plain exp(λ), but categorically
  not a power law.
- **Unit-conversion caveat:** the paper's λ = α²/(2Ωt) is bandwidth-referenced, not
  phonon-frequency-referenced like a U_eff/ω0 ratio — I could not find a clean, literature-stated
  conversion between the two conventions, so I'm not forcing one. Physically, though, U_eff ~ 2–4ω0
  is unambiguously moderate-to-strong coupling, which places it at or past every λ_op the paper
  reports (λ_op ranges ~0.3–0.9 depending on t/Ω, see §1.2) — i.e. expect to already be inside the
  stretched-exponential regime, not sitting comfortably in the flat weak-coupling shoulder.
- **The result is also adiabaticity-gated.** Marijanović, Qu & Demler (arXiv:2502.13005, 1D SSH,
  independent group) find the light regime only in the *anti-adiabatic* limit (ω0≳t); in the
  *adiabatic* limit (ω0<t — the regime most of the PRX paper's own headline numbers actually use),
  mass grows exponentially in the ordinary Holstein-like sense (exp(λ), not exp(√λ)) once coupling
  is significant. So which exponent you get depends on which adiabaticity regime your ω0/t sits in.

---

## 1. Core paper: Zhang, Sous, Reichman, Berciu, Millis, Prokof'ev, Svistunov, PRX 13, 011010 (2023)

**Citation correction:** first author is **C. Zhang**, not Y. Zhang. Full cite: C. Zhang, J. Sous,
D. R. Reichman, M. Berciu, A. J. Millis, N. V. Prokof'ev, B. V. Svistunov, "Bipolaronic
high-temperature superconductivity," *Phys. Rev. X* **13**, 011010 (2023).
**arXiv:[2203.07380](https://arxiv.org/abs/2203.07380)** (not 2207.05087). DOI:
[10.1103/PhysRevX.13.011010](https://doi.org/10.1103/PhysRevX.13.011010). PRX is Gold OA (CC BY) —
confirmed via Semantic Scholar (`isOpenAccess: true`, `status: GOLD`). PDF: [`pdfs/2203.07380.pdf`](pdfs/2203.07380.pdf).

Model: bond-Peierls (bond-SSH) electron-phonon coupling, 2D square lattice, lattice constant a≡1,
on-site Hubbard U. Mass unit m₀ = 2mₑ = 1/t. Coupling definitions:
- λ (bond-Peierls) = α²/(2Ωt)
- λ (Holstein, rescaled for their Fig. 1 comparison) = α²_H/(8Ωt)
- λ_ME (Migdal–Eliashberg) via m\*/m|_FS = 1+λ_ME (Eq. A3)
- Adiabaticity variable used throughout: **t/Ω** (reciprocal of ω0/t)

### 1.1 m\*/m₀ vs. coupling and adiabaticity

| t/Ω | λ | U/t | m\*/m₀ | Basis |
|---|---|---|---|---|
| 10/3 | 0.5 | 8 | **6.33 ± 1.0** | Fig. 6 caption (dispersion fit) — exact, stated |
| 10/3 | 0.5 | 8 | **5.4 ± 0.3** | Fig. 7 caption (R²(τ) fit, independent method) — exact, stated |
| 10/3 | ≈0.3 | 8 | ≈2 | Fig. 3c / Fig. 8 — visual read |
| 10/3 | ≈0.6 | 8 | ≈15–20 | Fig. 3c / Fig. 8 — visual read |
| 10/3 | ≈0.65 | 0 | ≈90–100 | Fig. 3c — visual read |
| 2 | ≈0.4 | 8 | ≈2 | Fig. 9(b2) — visual read |
| 2 | ≈1.0 | 8 | ≈20–30 | Fig. 9(b2) — visual read |
| 1 | ≈0.5 | 8 | ≈1–1.5 | Fig. 9(b3) — visual read |
| 1 | ≈1.0 | 8 | ≈4–5 | Fig. 9(b3) — visual read |
| 2 (Holstein) | 0.25→0.55 | 8 | rises steeply (log axis) | Fig. 12a — low confidence |

"Exact" = numeric value printed as text in a caption. "Visual read" = estimated off the rendered
figure (log-scale axes; the paper does not publish a data table) — treat as approximate only.

### 1.2 The "light" coupling window

No single universal boundary: the Tc/Ω-maximizing λ_op is **t/Ω-dependent**, shifting to larger λ
as t/Ω decreases (Appendix C). Only one explicit number is given: *"for t/Ω = 10/3, the optimal λ
is ∼0.5"* (App. C, p.10). Fig. 2's "intermediate coupling" points: λ=0.64 (t/Ω=2), λ=0.5 (t/Ω=10/3),
λ=0.5 (t/Ω=5). Beyond λ_op the text states bipolarons *"enter a regime of exponential mass
enhancement that becomes prominent for λ>λ_op and larger t/Ω"* (Sec. III, p.5) — the light regime is
explicitly bounded, not asymptotic.

Tc/Ω dome peaks (U/t=8, visual reads except where noted):

| t/Ω | Tc/Ω peak | near λ |
|---|---|---|
| 1 | ≈0.17–0.18 | ≈0.8–0.9 |
| 2 | ≈0.12–0.13 | ≈0.6–0.7 |
| 10/3 | ≈0.08–0.09 | ≈0.5 |
| 5 | ≈0.05–0.06 | ≈0.3–0.4 |
| Holstein (t/Ω=2) | **≈0.05Ω** (stated explicitly in text) | ≈0.25 |

### 1.3 Holstein vs. Peierls contrast, as the paper actually states it

- *"This behavior is completely absent in the standard Holstein model in which bipolarons rapidly
  become heavy in a manner that depends exponentially on the electron-phonon coupling strength
  [12, 13]"* (Sec. III, p.5).
- *"it appears that in all models the polaronic mass enhancement grows exponentially in λ"*
  (Intro, p.2).
- **Not found:** no fitted exponential formula or numeric exponent anywhere in text/captions/
  appendices — the exponential claim is qualitative, backed only by citations [12,13] (Bonča–
  Katrašnik–Trugman 2000; Macridin–Sawatzky–Jarrell 2004), not re-derived in this paper.
- **Not found:** the words "polynomial" and "algebraic" never appear in this paper — that
  dichotomy is external framing, not the authors'. Their own contrast is "weak-to-moderate
  enhancement" (bond-Peierls) vs. "increases very rapidly" / exponential (Holstein); no power-law
  fit is given for either.
- "Lang-Firsov" is never mentioned in this paper.

### 1.4 Tc estimates and their assumptions

- **The only absolute-temperature number in the paper: Tc ∼ 70 K** (Sec. IV Outlook, p.6–7), from
  Tc/Ω≈0.2 (peak at t/Ω∼1, Fig. 1) × Ω≈0.03 eV ("typical value of phonon frequency"), explicitly
  conditioned on *"if and only if the unusual limit of t∼Ω can be achieved."* Not tied to a specific
  measured material's phonon frequency (their own cited FeSe/FeTe values are 22/17 meV, App. H —
  the gap between that and the assumed 30 meV is not addressed in the text).
- **Formula** (Eq. 3, Sec. III p.4): 2D dilute/hard-core interacting Bose-gas BKT condensation
  temperature, Tc ≈ 1.84 ρ_BP/m\*_BP (constant from Pilati–Giorgini–Prokof'ev, PRL 100, 140405
  (2008)), reducing to Tc ≈ 0.5/(m\*_BP R²_BP) if R²_BP≥1, else 0.5/m\*_BP.
- **Density** n=ρ_BP is not a physical carrier density — it's set to the theoretical maximum
  non-overlap value ρ_BP = min{1/(πR²_BP), 1/π} ("largest Tc from this mechanism").
- **Lattice constant a≡1; dimensionality: 2D square lattice.** Bipolaron mass plugged into the 70 K
  estimate is not restated numerically at that point in the text.
- **Material target:** none computed quantitatively. Qualitative candidate = iron
  pnictides/chalcogenides (Sec. IV; App. H; Fig. 5): Fe–Fe hopping ≈50 meV, transverse phonon
  ≈22 meV (FeSe) / ≈17 meV (FeTe) ⇒ t/Ω≈2–3; λ≈0.5 "in one member of this family" (their Ref. 60).
  Explicit caveat in-text: *"the model studied here is not directly applicable to the pnictides."*
- Other Tc/Ω bounds quoted (not Kelvin): McMillan bound Tc/Ω∼0.05 at λ=1, μ\*=0.12 (Eq. F1);
  Migdal-Eliashberg theory ≲0.05 (Fig. 1 empty symbols).

---

## 2. Follow-ups and rebuttals, 2023–2026

**No 2D or 3D rebuttal of the light-bipolaron claim was found.** Every 2D/3D extension located
(all by the original author group) supports rather than refutes it. The only genuine qualification
found is 1D and adiabaticity-regime-dependent (Marijanović et al., below).

### 2.1 Supporting extensions (same author group)

| Result | Citation | Free PDF |
|---|---|---|
| **3D, explicit.** 3D cubic lattice (L=140), U=8t + long-range Coulomb V=U/10. Small bipolarons (R²<3): "weak mass enhancement m\*_BP<10" (units 2mₑ). Tc/Ω ≈0.07 (t/Ω=2) up to ≈0.27 (t/Ω=1), exceeding McMillan bound ≈0.05. | Sous, Zhang, Berciu, Reichman, Svistunov, Prokof'ev, Millis, "Bipolaronic superconductivity out of a Coulomb gas," *Phys. Rev. B* **108**, L220502 (2023). arXiv:[2210.14236](https://arxiv.org/abs/2210.14236) | [`pdfs/2210.14236.pdf`](pdfs/2210.14236.pdf) |
| Realistic phonon dispersion added; Tc "remains relatively high ... and even keeps increasing in the deep adiabatic regime both in two and three dimensions." | C. Zhang, N. Prokof'ev, B. Svistunov, "Effects of phonon dispersion on the bond-bipolaron superconductivity," *Phys. Rev. B* **111**, 184513 (2025). arXiv:[2409.14132](https://arxiv.org/abs/2409.14132) | [`pdfs/2409.14132.pdf`](pdfs/2409.14132.pdf) |
| **Synthesis/perspective**, explicitly 2D+3D. Gives the asymptotic instanton scaling Tc ~ exp[−#(t/Ω)√λ] (bond-Peierls) vs. exp[−#(t/Ω)λ] (Holstein) — the number behind the headline answer above. | J. Sous, "Bipolaronic High-Temperature Superconductivity from Phonon-Modulated Hopping: A Perspective," arXiv:[2605.16625](https://arxiv.org/abs/2605.16625) (2026, solo) | [`pdfs/2605.16625.pdf`](pdfs/2605.16625.pdf) |
| Semiclassical/instanton mechanism: strong-coupling bond-Peierls bipolaron slides within a degenerate manifold whose barrier vanishes as coupling grows (vs. Holstein's deepening self-trapping well) — the analytic reason for sub-exponential (√λ) mass growth. | K.-S. Kim, Z. Han, J. Sous, "Semiclassical theory of bipolaronic superconductivity in a bond-modulated electron-phonon model," *Phys. Rev. B* **109**, L220502 (2024). arXiv:[2308.01961](https://arxiv.org/abs/2308.01961) | [`pdfs/2308.01961.pdf`](pdfs/2308.01961.pdf) |

Further same-group 2D extensions (C. Zhang, solo; arXiv only, not mirrored locally — tangential to
the mass question): triangular lattice outperforms square, arXiv:[2507.07662](https://arxiv.org/abs/2507.07662)
(2025); long-range Coulomb, Tc stays high even at ω/t=0.5, arXiv:[2407.10444](https://arxiv.org/abs/2407.10444)
(2024); Holstein+SSH mixing leaves mass "nearly unchanged" while compressing pair size,
arXiv:[2511.06350](https://arxiv.org/abs/2511.06350) (2025).

### 2.2 Independent confirmation (non-group), with a mild tempering

- A. T. Ly, B. Cohen-Stead, S. Malkaruge Costa, S. Johnston, "Comparative study of the
  superconductivity in the Holstein and optical SSH models," *Phys. Rev. B* **108**, 184501 (2023).
  arXiv:[2307.10809](https://arxiv.org/abs/2307.10809) — [`pdfs/2307.10809.pdf`](pdfs/2307.10809.pdf).
  2D square lattice, **determinant QMC at finite density** (independent method from the anchor
  paper's diagrammatic MC). Confirms SSH bipolarons stay light "out to relatively large values of
  λ" vs. Holstein, but cautions "resulting Tc gains are small in the weak coupling limit."

### 2.3 Closest thing to a rebuttal — 1D, adiabaticity-regime-dependent

- F. Marijanović, Y.-F. Qu, E. Demler, "Bipolaron dynamics in the one-dimensional SSH model,"
  arXiv:[2502.13005](https://arxiv.org/abs/2502.13005) (2025, ETH Zürich, independent group) —
  [`pdfs/2502.13005.pdf`](pdfs/2502.13005.pdf). Light mass holds only in the **anti-adiabatic**
  regime (ω0≳t); in the **adiabatic** regime (slow phonons — arguably the more physically relevant
  one), "increasing electron-phonon interactions results in an exponential increase of the
  bipolaron mass," qualitatively like Holstein. Directly qualifies the light-mass claim as
  regime-dependent, even in 1D. This is a genuine, citable qualification — not a rebuttal of the
  2D/3D BEC-Tc results, which use different (adiabatic, t/Ω~1–10/3) parameters throughout.

### 2.4 Tangential (different mechanism, worth flagging but off-axis)

- Z. Zhang, Kuklov, Prokof'ev, Svistunov, "Superconductivity of bipolarons from quadratic
  electron-phonon interaction," *Phys. Rev. B* **111**, 134504 (2024/25), arXiv:[2408.03266](https://arxiv.org/abs/2408.03266)
  — quadratic coupling, not SSH/Peierls; by 2 of the 7 original authors as an alternate route.
- M. Grundner, T. Blatz, J. Sous, U. Schollwöck, S. Paeckel, "Cooper-Paired Bipolaronic
  Superconductors," arXiv:[2308.13427](https://arxiv.org/abs/2308.13427) (2023) — 1D DMRG many-body
  check finds a "fragmented condensate" rather than a simple BEC, qualifying the dilute-limit Tc
  estimates.
- G. D. Adebanjo, J. P. Hague, P. E. Kornilovitch, *Annalen der Physik* (2026),
  arXiv:[2507.17398](https://arxiv.org/abs/2507.17398) — 3D FCC/BCC "superlight" bipolarons, but
  Holstein/extended-Holstein (density) coupling, not SSH/Peierls.

---

## 3. Bipolarons in flat bands: does Lang-Firsov exponential scaling survive?

**This exact intersection is a thin, largely open field.** No paper found directly computes a
Holstein/SSH bipolaron mass on a lattice flat band and reports a replacement scaling law for
exp(g²). What exists: (a) flat-band Holstein models studied for other observables, (b) a solid body
of work proving quantum metric sets a *finite* composite-particle mass for non-phonon pairing
(Cooper pairs, excitons) in flat bands, (c) evidence from ordinary dispersive-band SSH/Peierls
models that off-diagonal coupling itself evades exponential Lang-Firsov scaling. Nobody has yet
visibly stitched (b) and (c) together into an explicit flat-band bipolaron mass formula.

**My synthesis of why quantum geometry should matter (reasoning, not a quoted result):** Lang-Firsov
band narrowing is multiplicative on the *bare* hopping, t_eff = t·exp(−S). On a lattice flat band
built from destructive interference (compact localized states), t=0 by construction — so the LF
exponential has nothing to multiply, and mass renormalization cannot come from single-particle
hopping suppression at all. The pair/composite mass that does emerge (per the quantum-metric papers
below) comes from second-order, geometry-controlled virtual processes instead, with 1/m\* set by
the quantum metric rather than by any bare t. This is a plausible bridge, not a literature-confirmed
result — flagging it as such.

| Finding | Citation | Free PDF |
|---|---|---|
| Holstein coupling on Lieb lattice, QMC. Derives U_eff = −λ²/ω0², finds CDW order at 1/3, 2/3 filling. **Does not compute polaron/bipolaron mass or address Lang-Firsov scaling** — a clear open gap. | Feng, Scalettar, "Interplay of flat electronic bands with Holstein phonons," *Phys. Rev. B* **102**, 235152 (2020). arXiv:[2009.05595](https://arxiv.org/abs/2009.05595) | [`pdfs/2009.05595.pdf`](pdfs/2009.05595.pdf) |
| Not flat-band, but mechanistically central (see §2.3): anti-adiabatic SSH bipolaron stays light because momentum-dependent coupling enables phonon-assisted hopping, breaking naive exp(g²) scaling. | Marijanović, Qu, Demler, arXiv:[2502.13005](https://arxiv.org/abs/2502.13005) (2025) | [`pdfs/2502.13005.pdf`](pdfs/2502.13005.pdf) |
| Pure SSH/off-diagonal coupling: mass stays near the bare value for essentially all coupling, diverging only at a critical line — qualitatively unlike Holstein's smooth exponential growth. | Marchand, Stamp, Berciu, "The polaron paradigm: a dual coupling effective band model," (2016). arXiv:[1609.03096](https://arxiv.org/abs/1609.03096) | [`pdfs/1609.03096.pdf`](pdfs/1609.03096.pdf) |
| **Foundational quantum-geometry template** (non-phonon): [1/m\*]ᵢⱼ = −(λ/N_cN_orb)Σₖ gᵢⱼ(k) — inverse two-body bound-state mass set purely by the quantum metric g, finite despite a flat band. Sawtooth ladder, Lieb, Harper model. | Törmä, Liang, Peotta, "Quantum metric and effective mass of a two-body bound state in a flat band," *Phys. Rev. B* **98**, 220511(R) (2018). arXiv:[1810.09870](https://arxiv.org/abs/1810.09870) | [`pdfs/1810.09870.pdf`](pdfs/1810.09870.pdf) |
| Flat-band excitons on Lieb lattice: 1/Mμν ∝ gμν/U. Electrons/holes individually infinitely heavy; bound composite is light, mass set by quantum metric. | Ying, Law, "Flat band excitons and quantum metric," arXiv:[2407.00325](https://arxiv.org/abs/2407.00325) (2024) | [`pdfs/2407.00325.pdf`](pdfs/2407.00325.pdf) |

**Adjacent, not directly on-topic (linked, not mirrored):**
- Pimenov, *Phys. Rev. B* **109**, 195153 (2024), arXiv:[2401.15155](https://arxiv.org/abs/2401.15155)
  — **caution:** this is a Fermi-*polaron* (mobile impurity in a fermion bath), not electron-phonon;
  finds band geometry reduces effective interaction (correction ∝ E_F·Σ quantum metrics).
- Yu, Ciccarino, Bianco, Errea, Narang, Bernevig, *Nature Physics* (2024), arXiv:[2305.02340](https://arxiv.org/abs/2305.02340)
  — quantum (Fubini-Study) metric contributes directly to the dimensionless e-ph coupling λ itself
  (~50% in graphene, ~90% in MgB₂); general bands, not flat-band-specific, but suggests g² in
  Lang-Firsov could itself inherit geometric structure.

**Unverified — flag explicitly, do not cite as fact:** a search snippet surfaced "Single-Polaron
Physics in Flat Bands: Exact Solution and Spectral Properties" (Kienesberger, Herzog-Arbeitman, Yu
et al., tentatively APS 2026), which would claim an exact single-polaron solution with e-ph coupling
in a uniform-quantum-geometry flat band — exactly the sought result, if real. **Could not confirm**
via arXiv listings, Semantic Scholar, or direct search; may be a very recent/unindexed abstract
(e.g. APS March Meeting 2026) or a search artifact. Worth checking arXiv new-submissions directly.

**Terminology false-friends to avoid when searching further:** "bipolaron flat bands" in the
Alexandrov school (e.g. *Phys. Rev. B* **53**, 2863) means an *emergent, narrow bipolaron band
arising from strong coupling*, not bipolarons on a pre-existing lattice-geometric flat band.
"Spin polarons in flat-band ferromagnets" (arXiv:2510.26798) and kagome RVB polarons
(arXiv:2606.12204) are magnetic/doped-Mott polarons, unrelated to electron-phonon coupling.

---

## 4. Empirical anchors: U_eff in A3C60 and Ba1-xKxBiO3

### 4.1 A3C60 fullerides — Jahn-Teller-mediated U_eff

Consistent pattern across independent ab initio and model-fit approaches: **bare** intramolecular
Jahn-Teller coupling J_JT ≈ 60–150 meV is largely cancelled by Hund's coupling, leaving a **net**
attractive J_eff of only **≈17–20 meV**.

| Value | Material/method | Citation | Free PDF |
|---|---|---|---|
| U_ph(0) ≈ −0.10 to −0.15 eV (density channel), J_ph(0) ≈ −51 meV (exchange channel); bare Hund's J≈31–36 meV, Hubbard U=0.82–1.07 eV; **net J_eff ≈ −17 meV**. Max JT phonon ≈0.2 eV, W≈0.5 eV. | K3C60→Cs3C60, cRPA+cDFPT ab initio | Nomura, Sakai, Capone, Arita, *J. Phys.: Condens. Matter* **28**, 153001 (2016), arXiv:[1512.05755](https://arxiv.org/abs/1512.05755); underlying values in Nomura, Sakai, Capone, Arita, *Sci. Adv.* **1**, e1500568 (2015), arXiv:1505.05849 | [`pdfs/1512.05755.pdf`](pdfs/1512.05755.pdf) |
| J_H ≈ −0.02 eV (−20 meV), λ≈0.13, inferred from K3C60/K4C60 spin-gap fit (0.07–0.1 eV); vibron scale ~0.1 eV. | Model fit | Capone, Fabrizio, Castellani, Tosatti, *Science* **296**, 2364 (2002), arXiv:[cond-mat/0207058](https://arxiv.org/abs/cond-mat/0207058) | [`pdfs/cond-mat_0207058.pdf`](pdfs/cond-mat_0207058.pdf) |
| Bare J_H≈0.03–0.1 eV, J_JT≈0.06–0.12 eV; **net J=J_JT−J_H≈0.02 eV**, λ_eff≈0.16–0.2; Hg phonon ħω≈90 meV; U≈1 eV, W≈0.6 eV. | Review | Capone, Fabrizio, Castellani, Tosatti, *Rev. Mod. Phys.* **81**, 943 (2009), arXiv:[0809.0910](https://arxiv.org/abs/0809.0910) | [`pdfs/0809.0910.pdf`](pdfs/0809.0910.pdf) |
| Same J_H/J_JT range independently cited | Cs3C60 context | Klupp, Matus, Kamarás, Ganin, McLennan, Rosseinsky, Takabayashi, McDonald, Prassides, *Nat. Commun.* **3**, 912 (2012) (OA) | [`pdfs/klupp_natcommun3_912_2012.pdf`](pdfs/klupp_natcommun3_912_2012.pdf) |
| J=0.07 eV attractive JT term, fit to reproduce the Cs3C60 insulating gap. | Cs3C60, DFT+DMFT | Baldassarre et al. (incl. Capone), *Sci. Rep.* **5**, 15240 (2015) (OA) | [`pdfs/baldassarre_screp5_15240_2015.pdf`](pdfs/baldassarre_screp5_15240_2015.pdf) |
| λ~0.5–1, ω_ph~0.2 eV, bare U~1–1.5 eV, JT energy ~tenths of eV, μ\*~0.3–0.4. | Review | Gunnarsson, *Rev. Mod. Phys.* **69**, 575 (1997), arXiv:[cond-mat/9611150](https://arxiv.org/abs/cond-mat/9611150) | [`pdfs/cond-mat_9611150.pdf`](pdfs/cond-mat_9611150.pdf) |
| λ~0.5–1, U/W~1.5–2.5, ω_ph/W~0.1–0.25. | Model | Han, Gunnarsson, Crespi, *Phys. Rev. Lett.* **90**, 167006 (2003), arXiv:[cond-mat/0208454](https://arxiv.org/abs/cond-mat/0208454) | [`pdfs/cond-mat_0208454.pdf`](pdfs/cond-mat_0208454.pdf) |

### 4.2 Ba1-xKxBiO3 — negative-U from Bi charge disproportionation

Context: x≈0.4, Tc≈29.8–30.5 K onset (Cava, Batlogg, Krajewski et al., *Nature* **332**, 814 (1988)
— no OA found, pre-arXiv).

| Value | Method | Citation | Free PDF |
|---|---|---|---|
| Foundational breathing-mode real-space-pairing (negative-U) model; no numeric U retrievable from primary text. | Theory | Rice, Sneddon, *Phys. Rev. Lett.* **47**, 689 (1981). DOI:[10.1103/PhysRevLett.47.689](https://doi.org/10.1103/PhysRevLett.47.689) — **no free PDF found** (paywalled, pre-arXiv) | — |
| Electronic (non-phonon) negative-U via nonlinear screening; no single headline number confirmed. | Theory | Varma, *Phys. Rev. Lett.* **61**, 2713 (1988). DOI:[10.1103/PhysRevLett.61.2713](https://doi.org/10.1103/PhysRevLett.61.2713) — **no free PDF found** | — |
| **\|U_eff\| ≈ 1.9–2.0 eV**, intersite 2zV≈0.13–0.24 eV — fit to BaBiO3/K-doped optical gap (~2 eV) and transport gap (~0.24 eV); W~1.6 eV from Batlogg DOS. | Review/fit | Taraphder, Pandit, Krishnamurthy, Ramakrishnan (1996), arXiv:[cond-mat/9601068](https://arxiv.org/abs/cond-mat/9601068) | [`pdfs/cond-mat_9601068.pdf`](pdfs/cond-mat_9601068.pdf) |
| λ=0.34 (harmonic+anharmonic e-ph), called **"too small"** for conventional Tc alone. | LDA, Ba0.6K0.4BiO3 | Meregalli, Savrasov, *Phys. Rev. B* **57**, 14453 (1998), arXiv:[cond-mat/9801251](https://arxiv.org/abs/cond-mat/9801251) | [`pdfs/cond-mat_9801251.pdf`](pdfs/cond-mat_9801251.pdf) |
| U=1.21±0.12 eV, ω_ph=70 meV (breathing mode); **U\*=U−2g\*²/ω_ph=(1.21−3.06) eV=−1.85 eV**, λ\*≈0.46. | Ab initio | Lukyanov et al., *Phys. Rev. B* **105**, 045131 (2022), arXiv:[2110.00084](https://arxiv.org/abs/2110.00084) | [`pdfs/2110.00084.pdf`](pdfs/2110.00084.pdf) |

**Disagreement flagged:** BKBO estimates span more than an order of magnitude — λ=0.34 from
harmonic-phonon LDA linear-response theory (explicitly called insufficient on its own) vs.
\|U_eff\| ~1.85–2.0 eV independently recovered by two very different methods 26 years apart
(1996 optical/transport-gap fit; 2022 ab initio nonlinear coupling) — suggesting ~2 eV is a fairly
robust cross-check, while pure harmonic-phonon linear-response theory alone underdelivers.

---

## Appendix: coupling-constant notation across sources (read before comparing numbers)

Different papers in this pull use genuinely different dimensionless coupling conventions — do not
equate them without checking definitions:
- Zhang et al. PRX 2023: λ = α²/(2Ωt) (bond-Peierls), bandwidth(t)-referenced.
- Sous 2026 perspective: same λ, used in the asymptotic exp[−c(t/Ω)√λ] scaling.
- A3C60 papers (Capone et al., Gunnarsson): λ = (coupling)/(bandwidth W or Hubbard U), phonon
  energy ω_ph quoted separately in eV/meV.
- BKBO papers: mix of λ (LDA linear-response, Meregalli-Savrasov) and U/U\* in eV (Taraphder;
  Lukyanov) — these are not the same quantity and the disagreement noted in §4.2 is partly a
  disagreement of *convention*, not only of magnitude.
- User's U_eff/ω0: phonon-frequency-referenced, closest in spirit to a Lang-Firsov g²=E_p/ω0 style
  coupling. No paper in this pull states a direct λ ↔ U_eff/ω0 conversion; treat any such
  conversion as an assumption you supply, not a literature-sourced number.

## PDFs not available free

Rice & Sneddon, PRL 47, 689 (1981); Varma, PRL 61, 2713 (1988); Cava et al., Nature 332, 814 (1988).
All three predate arXiv and sit behind APS/Nature paywalls with no OA copy located.
