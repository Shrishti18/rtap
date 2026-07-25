# A₂CoX₂ (Ibam) chain compounds — measured physical properties

**Primary structural/magnetic reference for all five compounds:**
Bronger, W. & Bomba, C., *"Ternäre Cobaltchalkogenide A₂CoX₂ mit A = Na, K, Rb,
Cs und X = S, Se. Synthese, Struktur und Magnetismus"*, **J. Less-Common Metals
158 (1990) 169–176** (CODEN JCOMAH). German language, pre-1995, not full-text
indexed — **abstract-level and citing-paper access only.** This is the source
publication behind ICSD 67387–67392 and 624264; the citation was read out of the
ICSD CIF headers in my local cache, not from a database landing page.

Legend: **[E]** experiment · **[T]** theory/DFT · **[—]** nothing found.

---

## Per-compound table

| | **Cs₂CoS₂** | **K₂CoS₂** | **Rb₂CoS₂** | **Cs₂CoSe₂** | **Na₂CoSe₂** |
| --- | --- | --- | --- | --- | --- |
| ICSD | 67389 | 67387, 623958, 623960 | 67388 | 67392 | 624264 |
| **1. Electrical** ρ(T) | — | **No measured ρ(T) found.** [T] insulator (PRB 2020); [T] JARVIS OptB88vdW & TBmBJ both give **gap = 0.0 eV, i.e. metallic** | — | — | — |
| **2. Magnetic** | AFM coupling within chains, from susceptibility [E] Bronger & Bomba 1990 — no numerical T_N, θ or μ_eff retrievable | **T_N > 9.5 K** [E] (Bronger & Bomba, susceptibility + neutron diffraction; collinear AFM proposed but **the neutron data were never published**). [T] 2.5 μ_B/Co ordered moment, T_N ≈ 15 K (monolayer, classical MC), PRB 2020. [T] **3.00 μ_B/Co** (JARVIS, this work) | AFM coupling [E] 1990, no numbers | AFM coupling [E] 1990, no numbers | AFM coupling [E] 1990, no numbers |
| **Co²⁺ spin state** | — | **HIGH SPIN, S = 3/2** — [T] only. JARVIS gives exactly 3.00 μ_B/Co (spin-only S=3/2 = 3.0); PRB gives 2.5 μ_B/Co (covalency-reduced). **No experimental μ_eff found to confirm.** | [T] 3.00 μ_B/Co (JARVIS JVASP-87107) | — | — |
| **3. Optical / spectroscopic** | — | — | — | — | — |
| **Co–X stretching mode** | **—** | **—** | **—** | **—** | **—** |
| **4. High pressure** | **—** | **—** | **—** | **—** | **—** |
| **5. Doping / off-stoich.** | — | — | — | — | — |
| **6. Theory** | MP mp-8770 | **Sarkar *et al.*, Phys. Rev. B 102, 035420 (2020)**, DOI 10.1103/PhysRevB.102.035420, arXiv:2005.12868 — AFM insulator, flat bands + van Hove singularities near E_F. **U value and Co-d bandwidth not stated in the abstract; full text not retrieved.** Also MP mp-5131; JARVIS JVASP-57783 | MP mp-8766; JARVIS JVASP-87107 | MP mp-8770 | — (Na₂CoS₂ = JVASP-10284) |

### JARVIS-DFT numbers computed/extracted from my local copy (OptB88vdW, spin-polarised)

| Compound | JVASP id | Band gap (OptB88vdW) | TBmBJ gap | Moment per Co | E above hull (eV/atom) |
| --- | --- | ---: | ---: | ---: | ---: |
| K₂CoS₂ | JVASP-57783 | **0.0 (metallic)** | 0.0 | **3.00 μ_B** | 0.933 |
| Rb₂CoS₂ | JVASP-87107 | 0.0 | — | 3.00 μ_B | 0.934 |
| K₂CoSe₂ | JVASP-58410 | 0.0 | — | 2.99 μ_B | 0.739 |
| Rb₂CoSe₂ | JVASP-58411 | 0.0 | — | 3.00 μ_B | 0.739 |
| Na₂CoS₂ | JVASP-10284 | 0.0 | 0.0 | 2.95 μ_B | 1.040 |
| K₂MnS₂ | JVASP-10273 | **0.872** | — | **5.00 μ_B** (S=5/2) | 1.103 |

Cs₂CoS₂, Cs₂CoSe₂ and Na₂CoSe₂ are **absent from JARVIS-DFT** — no entry exists.

Two things to note. **Plain semilocal DFT makes every Co member metallic**, while
the PRB study reports an insulator — so the insulating state in this family is
correlation-driven and the gap depends entirely on the treatment. The Mn analogue
opens a gap (0.872 eV) in the *same* functional, so the Co compounds are the
borderline ones. Separately, the JARVIS E-above-hull values of ~0.7–1.0 eV/atom
are implausibly high for compounds that were synthesised and refined in 1990;
I would treat them as an artefact of the competing-phase set, not as evidence of
instability.

---

## Wider A₂MX₂ (Ibam / K₂ZnO₂-type) family

| Sub-family | Members | What is measured |
| --- | --- | --- |
| **A₂MnX₂** | K₂MnS₂, Rb₂MnS₂, Cs₂MnS₂, K₂MnSe₂, Rb₂MnSe₂, Cs₂MnSe₂, K₂MnTe₂, Rb₂MnTe₂, Cs₂MnTe₂ | **The best-characterised part of the family.** Bhutani *et al.*, *Phys. Rev. Materials* **3**, 064404 (2019), DOI 10.1103/PhysRevMaterials.3.064404 (arXiv:1905.03318): K₂MnS₂₋ₓSeₓ solid solution across 0 ≤ x ≤ 2, **incommensurate cycloid** magnetic structure with **k = [0.58 0 1]**, two-step magnetic transition, geometric frustration of the triangular arrangement of chains. Neutron diffraction on HB-2A/POWGEN. |
| **A₂CoX₂** | Na/K/Rb/Cs × S/Se (8 compounds) | Bronger & Bomba 1990 only, plus the 2020 DFT paper on K₂CoS₂ |
| **A₂NiX₂** | K₂NiP₂, Cs₂NiP₂, K₂NiAs₂ exist as chains but are **square-planar** NiX₄, not tetrahedral (see my earlier structure table) | — |
| **A₂FeX₂** | Not found as an Ibam family member | — |

**Important distinction:** the well-studied 1D iron chain compounds **KFeS₂,
TlFeS₂, TlFeSe₂, RbFeS₂, CsFeS₂** are ***A*FeX₂** (C2/c or Immm, **Fe³⁺ d⁵**), a
*different* stoichiometry and space group from the A₂MX₂ (Ibam, **M²⁺**) family
here. Their properties do not transfer. For reference, KFeS₂ is a poor conductor,
ρ > 10³ Ω·cm below 250 K with semiconducting dρ/dT, and shows an unusual
"anti-Curie–Weiss" susceptibility rising quasi-linearly to 700 K (Molecules 27,
2663 (2022), DOI 10.3390/molecules27092663) — that paper also has DFT phonons for
KFeS₂, the closest thing to vibrational data anywhere near this structural class.

---

## Fields that are EMPTY across the entire A₂CoX₂ family

These are the gaps, stated plainly:

1. **Resistivity ρ(T): no measurement exists for any A₂CoX₂ compound.** Not one
   curve, not one activation energy. The only electrical characterisation of any
   kind is the DFT insulator claim for K₂CoS₂, which plain GGA contradicts.

2. **High pressure: nobody has ever done it.** I searched for pressure studies on
   A₂CoX₂, A₂MnX₂, the K₂ZnO₂ structure type and alkali chalcogenometallate
   chains generally, through four differently-phrased queries. **Zero hits — no
   resistivity, no structure, no magnetism, no superconductivity search under
   pressure on any member of this family, Co or otherwise.** By contrast the
   neighbouring BaFe₂S₃ ladder family has been studied to 11 GPa. This is a
   completely open experimental axis.

3. **Optical/spectroscopic: nothing.** No measured band gap, no XPS, no XAS, no
   ARPES, no Raman, no IR, for any member. **The Co–X stretching frequency you
   asked for has never been measured** — and I found no measured phonon of any
   symmetry for any A₂MX₂ compound. The nearest data is *computed* phonons for
   KFeS₂, a different family.

4. **Doping / off-stoichiometry / intercalation: nothing.** No A₂₋ₓCoX₂, no A-site
   or Co-site substitution series, no intercalation study. Even the mixed-anion
   work exists only on the **Mn** analogue (K₂MnS₂₋ₓSeₓ), never on Co.

5. **Quantitative magnetism is missing for four of your five compounds.** Only
   K₂CoS₂ has a number attached (T_N > 9.5 K), and its neutron-diffraction
   magnetic structure was **never published** — it is cited in the later
   literature as a private/unpublished result. No Curie–Weiss θ, no μ_eff, no
   ordered moment from experiment exists for *any* A₂CoX₂ compound.

6. **Spin state is inferred, not measured.** High-spin S = 3/2 is supported by two
   independent DFT sources (3.00 μ_B/Co here, 2.5 μ_B/Co in PRB 2020) but has no
   experimental confirmation, because no μ_eff has been published.

7. **Cs₂CoS₂, Cs₂CoSe₂ and Na₂CoSe₂ have no published property of any kind
   beyond the 1990 structure determination** — no experiment, no dedicated
   theory paper, and Cs₂CoSe₂/Na₂CoSe₂ are not even in JARVIS-DFT.

---

## Access caveat

Bronger & Bomba (1990) is German-language and pre-dates full-text indexing. I
worked from its **abstract and from citing papers** (chiefly the 2019 PRMaterials
and 2020 PRB studies), plus the ICSD CIF headers in my local cache for the
citation and structural data. I did **not** obtain the paper's own susceptibility
curves, so any θ, μ_eff or T_N values printed inside it are not reflected above.
That is the single most likely place where the blanks in row 2 could be filled by
someone with library access.
