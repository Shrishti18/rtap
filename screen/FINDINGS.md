# RTAP flat-band screen — findings

Screen of **1,772 JARVIS Wannier tight-binding Hamiltonians** (the complete
JARVIS-WTB set) for the derived specification: an isolated, narrow, partially
filled two-band manifold whose quantum geometry could support a high
superfluid stiffness.

Pipeline: `manifold_bands` locates isolated 2-band groups → `descriptors` per
single band → `tc_max` → Fermi-level gate.

---

## Finding 1 — Geometry is abundant, not scarce

The quantum metric was expected to be the bottleneck. It is not.

| quantity | median | max |
| --- | ---: | ---: |
| `M_trg` (per band, 728 bands) | **0.95** | 11.18 |
| `lam` (per band, = 1/n_phi) | 0.117 | 0.993 |
| `lam · M_trg` | 0.087 | 1.682 |

`M_trg` exceeds the 3D R=1.00 threshold of 0.339 in **545 of 728 bands (75%)**.
The two 2-orbital reference models predicted 0.19 and 0.37; real 44–72 orbital
Hamiltonians give the projector far more room to rotate, and the distribution
peaks well above the threshold rather than below it.

**The predicted anti-correlation is confirmed on real materials:**
corr(lam, M_trg) = **−0.29** Pearson, **−0.37** Spearman, **−0.45** log–log,
against the ~−0.4 predicted from random band structures. The 2D histogram shows
it directly: of 49 bands with lam > 0.8, 33 have M_trg < 0.1; the M_trg > 4
column is almost entirely lam < 0.2.

**Dilution, not stiffness, is the binding geometric constraint.** `tc_max`
reports `binds = 'amp'` for **75–80% of bands at every R** (586/728 at R=1.00,
539/728 at R=0.55). The amplitude cap `U·lam/4` is below the stiffness cap
`coef·R·U·M_trg` in three quarters of cases.

**The geometry was never the limiting factor.** Restricted to bands within
0.5 eV of E_F — the only ones with any physical claim — the maximum
`lam · M_trg` is **0.473**. Through
`Tc ≤ U·√(0.369·lam·M/4)`, that gives **Tc ≤ 0.208·U**, so even at the most
pessimistic interaction strength in the screen (U = 0.18 eV, the R=0.55
antiadiabatic cap) the geometric bound permits **Tc > 400 K**. Nothing in the
band geometry of these materials forbids room-temperature pairing.

One caveat on that same number: 0.473 sits just above the random-structure
maximum of 0.40 and below the ~0.5 "structurally anomalous" line. The near-E_F
geometry is *sufficient*, but it is not *exceptional*.

### What actually limits Tc_max

| R | max Tc_max | caps binding U | binds |
| --- | ---: | --- | --- |
| 1.00 | 2469 K | iso 658, phys 70, **adiab 0** | amp 586, stiff 142 |
| 0.85 | 1172 K | iso 546, **adiab 182** | amp 580, stiff 148 |
| 0.55 | 421 K | iso 374, **adiab 354** | amp 539, stiff 189 |

The entire R-dependence is the antiadiabatic cap on U. At R=1.00 it never
binds; at R=0.55 it binds half the set. The fall from 2469 K to 421 K is
mediator speed, not geometry.

---

## Finding 2 — The joint condition has ONE occupant in 1,772 materials

Isolated **and** narrow **and** partially filled is satisfied by a single band
in the entire dataset.

| stage | survivors |
| --- | ---: |
| materials scanned | 1,772 |
| with ≥1 isolated 2-band group (`W_band ≤ 0.5 eV`, `iso ≥ 0.05 eV`) | 179 (10.1%) |
| candidate bands | 728 |
| **after the Fermi-level gate** | **1** |

Breakdown: 716 bands fail because E_F is not inside the group, 11 more fail the
filling cut, 1 passes.

**The failure is structural, not statistical.** Isolation is *created by a gap*,
and E_F sits either inside that gap or in the dense part of the spectrum — the
two conditions are close to mutually exclusive in an undoped ground state:

* E_F misses the candidate group by a **median of 5.15 eV**, which is
  **21× the median group width** (0.241 eV). Only 12 of 728 groups contain
  E_F at all.
* Filling is almost perfectly **bimodal**: **200 bands at f < 0.05**,
  **527 at f > 0.95**, and **one** band in between. Nothing sits near
  half-filling.

The high `lam·M` values belong entirely to unfillable bands. Sr3BiSb's 1.33 —
more than 3× the random-structure maximum — is a property of an empty
conduction manifold. Restricting to bands near E_F drops the maximum from
1.682 to 0.473. The A3BX antiperovskite cluster (Ba3AsN, Sr3BiN, Ca3Bi2,
Sr3BiSb, Ba3PbO, …) that dominated the pre-gate ranking fails entirely on
`ef_outside`; its large `d_iso` was the signature of a band gap.

**The correct reading is "not an undoped ground state", not "impossible".**
Every cuprate parent compound would fail this same gate — they are
antiferromagnetic insulators whose superconductivity appears only on doping.
A screen over undoped stoichiometric ground states cannot see a material whose
useful filling is reached chemically. What this result rules out is finding the
target *ready-made*; it says nothing about doped or gated versions of the same
structures.

### Validation of the gate

E_F is read from each material's own `OUTCAR`, so it shares the absolute energy
scale of the `_hr.dat` eigenvalues. On a random sample of 12 materials the
Wannier spectrum brackets E_F **12/12** — including the two headline losers,
Ca3Bi2 (E_F = 3.30, spectrum [−8.16, 8.41]) and ZnPd5Se (E_F = 7.79, spectrum
[−6.97, 19.13]). The Wannier models do span the Fermi level; the isolated
narrow groups simply are not there.

---

## Finding 3 — The sole survivor is unremarkable

**HfFeCl6** — JVASP-6766, 2D, SG 149, group 41-42, band 41.

| | value |
| --- | ---: |
| filling f | 0.917 |
| lam | 0.206 |
| n_phi | 4.86 |
| M_trg | 0.497 |
| **lam · M_trg** | **0.102** |
| W_band | 0.091 eV |
| d_iso | 0.151 eV |
| Tc_max (R=1.00 / 0.55) | 45.0 K / 35.6 K |

`lam·M = 0.102` is **exactly the random-structure median of 0.10**. This is not
a discovery — it is one ordinary band that happened to clear a gate everything
else failed. Its n_phi of 4.86 is nowhere near the doublet target of 2, its
filling of 0.917 is at the edge of the accepted window, and its Tc_max is below
the boiling point of nitrogen at every mediator speed. It should not be read as
a candidate.

---

## Caveats

**Grid.** The sweep used nk=24 in 3D (nk=40 in 2D). `M_trg` is therefore a
**lower bound** — about 2% low on a smooth band, worse near a touching. `d_iso`
on a coarse grid is an **upper bound**, the unsafe direction, since a coarse
grid can miss a touching point entirely. The planned nk=64 confirmation was
started but abandoned once the E_F gate eliminated every band it was checking.

**M_naive / M_min: not measured (affects eliminated bands only).** The
gauge-invariant minimal metric was not obtained. Every `M_trg` here is the naive
value and therefore an upper bound on the true gauge-invariant metric — the one
prior datapoint (CrSe2, ratio 2.27) suggests the overestimate can be a factor
of 2–3. This does not affect any conclusion: it can only reorder bands the E_F
gate has already removed, and Finding 1's argument is that geometry is
*abundant*, which a downward correction of 2–3× does not overturn (0.473/2.27 =
0.21, still permitting Tc ≤ 0.14·U). Where a bounded Powell refinement was
attempted, its budget cap means M_min would itself have been an upper bound on
the true minimum, so the ratio would have been a lower bound.

**Scope of the Wannierisation.** JARVIS-WTB Wannierises a chosen energy window
per material. A qualifying group falling outside that window is invisible to
this screen. The 10.1% group-discovery rate is a property of the dataset's
disentanglement choices as much as of the materials.

**Search parameters.** Only `nb = 2` groups were sought, with
`W_band ≤ 0.5 eV` and `iso ≥ 0.05 eV`. Wider manifolds (n_phi = 3 t2g-derived
groups, for instance) and looser windows were not searched.

**No attractive U was verified anywhere.** This is the most important
limitation. The screen *grants* U up to `min(d_iso/2, U_cap, 2·w0)` and asks
only whether the band geometry could exploit it if it existed. It never tested
whether an attractive pairing interaction is present in any of these materials.
**Every Tc_max in this report is conditional on a pairing interaction that was
never checked.** The geometric screen is now saturated on this dataset; pairing
chemistry is untouched and is the binding open question.

---

## Data

| file | contents |
| --- | --- |
| `full_bands.csv` | 728 candidate bands, all Phase-A descriptors |
| `full_bands_ef.csv` | the same, plus E_F, group extrema, filling, gate status |
| `stage1_sites.csv` | 37,001 site rows from the CIF symmetry pre-filter |
| `full_scan.py`, `ef_gate.py`, `analyze_full.py` | the sweep and its analysis |
| `fastops.py` | rank-r reformulations, validated against `harvest` to ~1e-14 |
| `validate_fastops.py` | that validation |

---

## Terminus

The project closes on `review/LIGAND_U.md`. The final chain:

> n_φ = 2 forced → M ≤ ⟨v²⟩/Δ² → winding breaks the tie → isolation vs.
> metallicity → SOC is the only non-hopping scale → isolating it quenches n_φ
> → g → 0

The last link is the first one: C1/C2 (n_φ → 1 ⇒ tr g ≡ 0) is what closes the
tetragonal-crystal-field loophole in the SOC route. The cRPA U question that was
expected to decide it does invert — U(Ir⁴⁺, telluride) ≈ 0.95 eV < 3λ = 1.25 eV —
and is irrelevant: U is capped by Δ_iso/2 = 0.106 eV at the metric optimum, 9×
below the telluride estimate. Ceiling **4–16 K** depending on the projection
condition (8.1 K under the project's own U ≤ Δ_iso/2), three orders of magnitude
short.

**No viable route to room-temperature superconductivity through flat-band
quantum geometry in d-electron systems.** Every link derived.

Two limitations carried on the record: there is **no published cRPA U for any Ir
chalcogenide or halide** (the 0.95 eV is an extrapolation from Ta/Nb/Ru/Fe ligand
factors, and is a genuine hole in the literature), and `c = 0.67` is mean-field
and uncalibrated, so the Kelvin figures are ceilings on a ceiling.

### The three root assumptions

The framework assumed singlet pairing from an instantaneous attraction in an
isolated narrow band. All three were tested.

| assumption | status |
| --- | --- |
| (i) singlet pairing | **branch.** The pairing channel moves only the amplitude cap, which sits 12.6× above the binding one. `review/triplet.py`: the on-site orbital-antisymmetric channel is exactly empty without SOC (0.000e+00), is repulsive in every real d-electron system (Kanamori U − 3J > 0 for measured J/U = 0.10–0.25), and at λ = 1/3, 2/3 or 1 the terminus stays 8.1 K. |
| (ii) instantaneous attraction | **costs a factor, does not move the ceiling.** |
| (iii) isolated narrow band | **root, and closed both ways.** Gapped: M ∝ t² exactly, ceiling 8.1 K. Ungapped: `<tr g>` has no limit — `review/lieb.py` gives `<tr g> = ln(nk)/2π + const`, increment per octave 0.11033 against ln2/2π = 0.11032, with the whole divergence localised at the touching (excision restores convergence to 0.999). And D_s at the touching is *smaller* than gapped anyway: 0.100 vs 0.125 at U = 0.5. |

Isolation caps the metric; removing isolation removes the metric's *meaning*,
not merely its usefulness. There is no third position.

**The terminus and the project's worst error are the same statement.** A rank-1
projector at a degeneracy is not a property of the Hamiltonian — that is D1, and
it is also why "the metric of the flat band" at a touching does not exist. The
first theorem closes the last door.

### Retracted

- **`c_eff = 0.047` for the touching flat band, and `M_Lieb = 4.3016`.** Both
  were computed in the non-periodic `2t·cos(k/2)` gauge, where `H(0) = −H(2π)`,
  so differencing across the BZ boundary injected a false discontinuity (the
  same term `harvest._trg_shifted` was fixed for; it shows up as a metric peak
  at Γ, where the gap is *largest*). Corrected, `<tr g>` is log-divergent, so a
  conversion coefficient at a touching is not a constant — it falls like
  1/ln(nk) without bound. The M-free comparison D_s(touching) < D_s(gapped)
  replaces it and carries the same conclusion.
- **Pairing cap "202 K, 25×".** Correct values: 102.3 K at U ≤ Δ_iso/2, ratio
  12.6× and constant across projection conditions.

---

## RETRACTION: "geometry is abundant" degrades by 3.5× under self-consistency

The claim — *median `M_trg` = 0.97, 542/708 bands above the 0.339 threshold* —
is **single-band `M` throughout**, paired with the group's isolation gap. That
mismatch is the D1 error in the primary data product (see
`screen/CONSTD_FINDINGS.md`).

### The standing rule that replaces it

A quantity must be computed on a projector whose rank matches the degeneracy
structure of what it projects onto — and **which rank that is is decided by Tc,
not by convention**:

- **rank-r manifold** treatment is valid iff the intra-group splitting δ ≲ k_B·Tc
  (the bands then pair coherently): `U ≤ d_iso/2`, `M = M_manifold`, `lam = lam_manifold`.
- **single-band** treatment is valid iff δ ≫ k_B·Tc (the partner is spectrally
  inert): `U ≤ δ_nearest/2` — the band's **own** nearest gap — `M = M_trg`, `lam = lam`.

Using one treatment's `U` with the other's `M` is what inflated LiPbAu₂ from 9 K
to 186 K.

### Result (`selfconsist.py`, 708 bands, 179 materials, 599 s)

Group ranks: 688 of 708 bands sit in rank-2 groups, 20 in rank-1.
Self-consistency admits **single** for 222, **either** for 486, and **neither**
for none. The two resolutions of "either" agree:

| | median `M` | above 0.339 |
|---|---|---|
| stored (single-band throughout) | **0.9657** | **542/708 (76.6%)** |
| applied, either → min | 0.2797 | 318/708 (44.9%) |
| applied, either → manifold | 0.2895 | 322/708 (45.5%) |

**The median falls below the threshold it was being compared against.**

### The degradation is uniform, not a tail effect

| percentile | stored | applied |
|---|---|---|
| 10th | 0.0829 | 0.0182 |
| 25th | 0.3668 | 0.0704 |
| 50th | 0.9657 | 0.2797 |
| 75th | 2.2008 | 0.5662 |
| 90th | 3.4893 | 1.0919 |
| 99th | 7.6737 | 2.5280 |

A roughly constant ~3.5× across the whole distribution. This is **unlike** the Tc
ranking, where the median was near-correct (1.06×) and only the tail was
corrupted (up to 374×). Predicting "survives in the median, degrades at the top"
was wrong in both directions: here it degrades in the median and the top alike.

Cause: among the 486 either-cases, `M_manifold < M_single` in **444**, median
0.281 against 1.483 — a 5.3× inflation. A rank-2 projector removes the
intra-pair channel, and for a near-degenerate pair that channel is most of the
single-band metric. In 44/708 bands `M_manifold < 0.01` outright — the
"both bands span the same orbitals" case.

### What survives

**Geometry is common but not abundant**: ~45% of candidate bands exceed the 3D
threshold, not 77%, and the typical band sits just below it. The qualitative
conclusion that large quantum metric is *not* the scarce ingredient still holds
— `n_φ`, isolation and filling remain scarcer — but the margin is 3.5× smaller
than every strategic discussion has assumed, including the geometry-first
reframe.

### Not affected

The 727/728 E_F gate (band energies and E_F only, no metric anywhere), the four
theorems (analytic, independently verified), and F1 (chemistry, no band structure).

*(Review note: the self-consistency criterion uses bare inequalities and the
BZ-minimum intra-group splitting, both generous to manifold validity. With a
factor-3 safety margin the fraction above threshold is 48.6% and the median
0.329 — so the abundance figure is **45–49%** under any reasonable margin, and
the retraction stands. Stored manifold columns were spot-checked against raw
hr.dat on three materials including AgBr: 6-digit agreement.)*
