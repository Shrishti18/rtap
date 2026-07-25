# cRPA U for Ir⁴⁺ in a telluride or iodide host

The question posed: *if U comes back below ~1.25 eV, C8 inverts and the SOC route
is open; if it comes back at 1.8–2.0 like the oxides, the terminus is final.*

**Answer: U ≈ 0.7–1.0 eV, so C8's inequality does invert — and the route stays
closed anyway, because U was never the binding constraint.** The ceiling on the
j_eff = 1/2 route is **8.1 K**, and it is independent of U.

---

## 1. The literature number

### 1a. There is no published cRPA U for any Ir chalcogenide or halide

Searched for cRPA/constrained-RPA U on IrTe₂, Ir sulfides/selenides, IrX₃,
A₂IrX₆. **Nothing.** Every published cRPA U for Ir⁴⁺ 5d t₂g is an oxide. This is
a real gap, not a search failure — the non-oxide Ir literature uses empirical
DFT+U with U chosen to reproduce a gap, not ab-initio screening.

So the number has to be built from (oxide anchor) × (measured ligand factor).

### 1b. Oxide anchor

| system | U (eV) | J (eV) | window | source |
|---|---|---|---|---|
| Sr₂IrO₄ | 1.82 | 0.22 | Ir 5d t₂g | RP-iridate cRPA |
| Sr₂IrO₄ | 1.93 | 0.16 | Ir 5d t₂g | Arita, Kuneš, Kozhevnikov, Eguiluz, Imada, PRL 108, 086403 |
| Ir 5d t₂g (range in the review) | 1.96–2.54 | — | — | D6 |

Anchor: **U(Ir⁴⁺, oxide) = 1.9 eV.**

### 1c. Ligand factor, measured at fixed metal, fixed subspace, fixed structure

The cleanest data is the MX₂ cRPA sweep (Nb, Ta × S, Se, Te; monolayer,
`arXiv:2402.01273`), which holds everything but the chalcogen fixed:

| metal | phase | subspace | U(S) | U(Se) | U(Te) | Te/S |
|---|---|---|---|---|---|---|
| Ta (5d) | 1T | full d | 1.89 | 1.52 | 0.94 | **0.50** |
| Ta (5d) | 1H | full d | 2.40 | 2.45 | 1.14 | **0.475** |
| Ta (5d) | 1T | d_z² | 1.15 | 1.14 | 0.78 | 0.68 |
| Nb (4d) | 1H | full d | 2.30 | 2.25 | 0.94 | 0.41 |
| Nb (4d) | 1T | full d | 1.77 | 1.29 | 0.93 | 0.53 |

The bare V is lowest for Te in both Ta phases (1T full-d: 7.87 / 8.85 / **7.62**;
1H full-d: 8.34 / 8.96 / **7.88**) — not monotonic in the ligand, but the Te end
is the least bare-repulsive — which the authors attribute to the Wannier
functions spilling onto the chalcogen: *"the coupling to
neighbouring chalcogen p states gets stronger in Te-based TMDs, which makes the
Wannier functions more extended."* Same mechanism, halide series, Ru 4d
(`npj Quantum Materials 7, 75`, `arXiv:2203.01626`): U_eff **2.7 (Cl) → 2.1 (Br)
→ 1.4 (I)**, ratio **0.52**, explicitly following their cRPA U_avg trend. And
Fe 3d: U is 24% larger in FeSe than FeTe (Se→Te alone = 0.81).

Three independent families, two different metals rows, halide and chalcogen
series: **the heavy-ligand factor is 0.5, reproducibly.**

### 1d. The number

U(Ir⁴⁺, telluride) ≈ 1.9 × 0.5 ≈ **0.95 eV**, bracket **0.7–1.0 eV**.

Two reasons this is an over-estimate, i.e. the bound is conservative in the
direction that helps the inversion scenario:
- the factor 0.5 is S→Te and Cl→I. **O→Te is a bigger step** and is not included.
- the MX₂ values are **monolayer**, which screens badly; bulk U is lower still.

**C8's inequality U_rep < 3λ = 1.25 eV is therefore satisfied by a telluride, and
marginally by an iodide.** The user's stated inversion condition is met.

---

## 2. It does not help, and the reason is D2

### 2a. The chemistry: Ir⁴⁺ is not attainable in a telluride

IrTe₂ is a **negative charge-transfer** system. Te 5p sits *above* Ir t₂g; the
measured Ir 5d → Te 5p charge-transfer energy is ≈ **0.1 eV**; the nominal
valence is **Ir³⁺ (t₂g⁶)** with Te ligand holes, not Ir⁴⁺
(`Nature Communications 6, 7342`). The j_eff = 1/2 construction requires the hole
to sit on the metal in a t₂g state of essentially pure d character. In a
telluride it does not — it goes to the ligand. **The same p–d hybridization that
lowers U is what removes the d⁵ configuration.** These are not two effects to be
traded off; they are one effect.

Iodides are the same story one step milder: [IrI₆]²⁻ is reported as a solution
species (IrI₄ + KI), but the solid-state A₂IrX₆ vacancy-ordered perovskites are
established only for **X = Cl, Br**. K₂IrCl₆ and K₂IrBr₆ are the closest realized
"Ir⁴⁺ in a polarizable halide", and they reach near-ideal j = 1/2 precisely
because the IrX₆ octahedra are **isolated** — no inter-octahedral dispersion.
That is the t → 0 limit, where the projector is k-independent and tr g → 0. The
one realized member of the family sits at the dead end of the metric axis.

### 2b. The model: U is not the binding constraint

`ligand_scan.py` / `ligand_scan2.py`, λ = 0.417 eV (RIXS, fixed by the ion and
unaffected by the ligand), rank-2 projector, nk = 32.

Large-t (covalent, telluride-like) branch — U inverts, isolation dies:

| ligand | t (eV) | λ/t | W=8t | iso (eV) | M_trg | U_cRPA | U<3λ? |
|---|---|---|---|---|---|---|---|
| oxide | 0.26 | 1.60 | 2.08 | **−1.444** | 0.332 | 1.90 | no |
| sulfide | 0.36 | 1.16 | 2.88 | **−2.241** | 0.524 | 1.35 | no |
| selenide | 0.44 | 0.95 | 3.52 | **−2.878** | 0.678 | 1.15 | **yes** |
| telluride | 0.52 | 0.80 | 4.16 | **−3.514** | 0.829 | 0.95 | **yes** |

`iso < 0` throughout: the j = 1/2 doublet **overlaps** the j = 3/2 quartet. There
is no isolated manifold to project onto, so `U ≤ Δ_iso/2` is not merely violated,
it is unsatisfiable. M_trg is large (0.83, far above the 0.339 threshold) — and
meaningless, because it is the metric of a manifold that isn't separated.

Small-t (ionic, molecular-halide-like) branch — isolation returns, metric dies:

| t (eV) | λ/t | iso (eV) | M_trg | U ≤ iso/2 | Tc_geo (K) |
|---|---|---|---|---|---|
| 0.02 | 20.9 | 0.466 | 0.0027 | 0.233 | 3 |
| 0.03 | 13.9 | 0.387 | 0.0060 | 0.193 | 5 |
| **0.052** | **8.0** | **0.212** | **0.0180** | **0.106** | **8.1** |
| 0.068 | 6.1 | 0.084 | 0.0305 | 0.042 | 5.5 |
| 0.078 | 5.3 | 0.005 | 0.0399 | 0.002 | 0 |

M_trg ∝ t² exactly in this regime (0.0007 / 0.0027 / 0.0060 at t = 0.01 / 0.02 /
0.03), as it must: the metric measures inter-site hybridization, and isolation is
the statement that there isn't any.

**Ceiling under the project's own projection rule: Tc = 8.1 K at t = 52 meV.**
Relaxing to U ≤ Δ_iso gives 16 K; tightening to D8's U ≤ Δ_iso/4 gives 4 K.

At the optimum, U is capped at **0.106 eV** — **9× below** the telluride cRPA
estimate of 0.95 eV, and 18× below the oxide value. **Whatever cRPA returns for a
telluride, it is not the limiter.** The limiter is Δ_iso, which is set by λ, which
is set by the ion.

### 2c. The tetragonal-crystal-field loophole is closed by C1

The obvious objection to §2b: the cubic model isolates the doublet by SOC alone,
yet Sr₂IrO₄ *is* isolated at t = 0.26 eV where the model says iso = −1.44 eV.
Real j = 1/2 materials get isolation from a tetragonal crystal field too. If
Δ_tet can hold the manifold apart at large t, then large M and isolation coexist
and the 8 K ceiling is a model artefact.

`ligand_scan3.py` adds Δ_tet on the xy orbital and scans (t, Δ_tet):

| Δ_tet (eV) | t (eV) | iso | M_trg | **n_phi** | lam_pair | Tc_geo (K) |
|---|---|---|---|---|---|---|
| 0.00 | 0.05 | 0.227 | 0.0167 | **3.000** | 0.333 | **8.1** |
| 0.40 | 0.05 | 0.222 | 0.0114 | 2.078 | 0.483 | 5.4 |
| 0.70 | 0.05 | 0.370 | 0.0047 | 1.467 | 0.683 | 3.7 |
| 1.00 | 0.05 | 0.589 | 0.0017 | 1.230 | 0.813 | 2.2 |
| 1.50 | 0.05 | 1.022 | 0.0004 | **1.097** | 0.912 | 0.9 |
| 2.00 | 0.05 | 1.488 | 0.0001 | **1.052** | 0.950 | 0.4 |

Δ_tet does buy isolation — up to 1.5 eV of it — and it buys it by **quenching the
orbital multiplicity**. n_phi runs 3.000 → 1.05 and lam_pair runs 1/3 → 0.95.
That is the n_phi → 1 limit the project already knows is **dead by theorem**: a
single-orbital band has a k-independent projector and tr g ≡ 0 (C2). M_trg falls
by 170× across the table.

**The cubic point is the optimum.** Every deviation from it trades the one thing
that makes the j = 1/2 doublet interesting — n_phi = 3 exactly, protected by cubic
symmetry — for isolation. The 8.1 K ceiling stands.

---

## 3. What this does to the claims

| claim | status |
|---|---|
| **C8** as stated (window open iff λ > U_rep/3) | **inverts** for a telluride/iodide host: U ≈ 0.95 eV < 1.25 eV. The inequality is satisfied. |
| **C8's conclusion** (the SOC route is closed) | **survives, by a different and stronger argument.** The inequality was never the binding constraint. D5 already showed C8's inequality is 4× looser than the project's own rule; §2b shows the project's own rule caps U at 0.106 eV, and §2c shows nothing can be done about it. |
| **D2** (isolation and stiffness trade off) | **confirmed and quantified**: a hard 8.1 K ceiling, with the tetragonal loophole closed. |
| **C1/C2** (n_phi → 1 is dead) | **load-bearing in a new place**: it is what closes the Δ_tet loophole. |

**The terminus is final.** Not because U came back at 1.8–2.0 — it came back at
~0.95 — but because the SOC route's ceiling does not depend on U at all.

---

## 4. Caveats

- **No cRPA U for any Ir non-oxide exists.** The 0.95 eV is an extrapolation from
  Nb/Ta/Ru/Fe ligand series, not a calculation on an Ir compound. It is a
  *bound* in the useful direction (§1d), not a measurement.
- The 8.1 K ceiling is computed in `jeff.hk`, the orbital-selective cubic t₂g
  model, at **tp = 0**, which makes each t₂g orbital strictly 2D (Attack 2). The
  full 6×6 is 3D through SOC, but a real 3D hopping structure has not been tested.
- The Tc's use **c = 0.67, mean-field, uncalibrated against exact methods** (the
  standing caveat from the stood-down fix #3). They are upper bounds within mean
  field. 8.1 K is a ceiling on a ceiling.
- λ_SOC is taken as ligand-independent. Heavy ligands add their own SOC (the RuX₃
  paper finds ligand SOC becomes important for Br and I and that the two sources
  *do not* add constructively). Including it would move λ by tens of meV, not the
  factor of ~7 that §2b would require.
- `iso` here is the gap to the j = 3/2 quartet in a 6-state model. Real materials
  have e_g, ligand p, and further bands that can only reduce it.

---

## 5. On the record

### F1. There is no published cRPA U for any Ir chalcogenide or halide.

Stated as a finding, not a caveat. Every ab-initio screened-interaction value for
Ir⁴⁺ 5d t₂g in the literature is an oxide. The non-oxide Ir literature — IrTe₂,
IrX₃, A₂IrX₆ — parameterises U empirically, by fitting a gap. The 0.95 eV in §1d
is an **extrapolation** from Ta/Nb (MX₂), Ru (RuX₃) and Fe (pnictide/chalcogenide)
ligand factors applied to an Ir oxide anchor. It is a bound in the useful
direction, not a calculation on an Ir compound. Anyone who needs this number for
real has to compute it.

### F2. The 8 K ceiling is a ceiling on a ceiling.

`c = 0.67` is mean-field and uncalibrated against exact methods (the standing
caveat from the stood-down ED validation). Every Tc here is an upper bound
*within* mean field, and the true stiffness is expected to be lower.

### F3. The Tc bracket, by projection condition

| condition | U at t = 52 meV | Tc |
|---|---|---|
| U ≤ Δ_iso (no safety factor) | 0.212 eV | **16.3 K** |
| U ≤ Δ_iso/2 (the project's own rule) | 0.106 eV | **8.1 K** |
| U ≤ Δ_iso/4 (D8's defensible rule) | 0.053 eV | **4.1 K** |

The 8.1 K quoted in §2b is the **Δ_iso/2** number — the project's own rule, not
D8's tightening. The factor 2 between it and 16.3 K is the safety factor on the
projection, not an amplitude-vs-stiffness branch: stiffness binds throughout
(at U = 0.106 eV the pairing cap U·λ_pair/4 is 202 K, 25× above the geometric
cap). Applying D8 would give 4.1 K.

Three orders of magnitude short under every one of them.

## Sources

- MX₂ cRPA sweep: <https://arxiv.org/pdf/2402.01273>
- RuX₃ (X = Cl, Br, I): <https://www.nature.com/articles/s41535-022-00481-3>, <https://arxiv.org/pdf/2203.01626>
- Fe pnictide/chalcogenide cRPA: <https://arxiv.org/pdf/0911.3705>
- Sr₂IrO₄ cRPA and j = 1/2 review: <https://arxiv.org/pdf/2110.12877>
- IrTe₂ negative charge transfer / Ir³⁺: <https://www.nature.com/articles/ncomms8342>
- K₂IrCl₆ isolated octahedra, fcc j = 1/2: <https://arxiv.org/pdf/1903.01660>
