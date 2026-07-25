# Work order: can a winding gap come from symmetry instead of fine-tuning?

## State of the problem

The 8.1 K terminus has been **retracted for the winding class**. The mechanism,
now confirmed in 1D and 3D:

- An isolation gap that is **k-independent** (atomic SOC — a constant on-site
  splitting) isolates a band *without contributing quantum metric*. The
  projector stays nearly k-independent, so `tr g` stays small. This is C2's
  mechanism in continuous form, and it is why cubic t2g + atomic SOC reached
  only 11% of its own ceiling.
- A **winding** gap does not have that defect. Write a two-band manifold as
  `H(k) = d0(k)·I + d(k)·σ`. If `|d(k)|` is constant (= D) but the direction
  `d̂(k)` winds `n` times, the gap is `2D` everywhere, both bands are exactly
  flat, and `M_trg = |n|²/4` — **independent of D**. Verified: 1D gives
  0.25000 / 1.00000 / 2.24998 / 3.99995 for n = 1..4 at D = 1 and D = 5 alike;
  3D adds across axes, `n = (1,1,1)` → `M_trg = 0.74`; and it survives
  hybridization with a dispersive manifold (`M_trg` 0.7375 → 0.7436 as
  hybridization grows to the winding scale itself).
- At Ir's own λ = 0.417 eV: cubic t2g gives **8.1 K**, winding n=1 gives
  **334 K**. λ needed for 300 K falls to **0.374 eV, below Ir**.

**The one thing blocking this from being a route**: the winding model needs the
real-space hopping `H(R)` nonzero **only** at `R = ±n` — an n-th neighbour hop
with all shorter-range hops absent. In 3D with `n = (1,1,1)`, a ⟨111⟩ hop and no
⟨100⟩ hops. No mechanism has been shown. One attempt (`target3d.py`) failed on Γ
degeneracy. This is exactly the gap where Cs₂CoS₂ died: "this structure works"
is not "this structure exists".

## The question

**Can `|d(k)| ≈ const` with winding `d̂(k)` arise from something other than
fine-tuned long-range hopping?**

### Lead 1 — non-symmorphic symmetry (highest priority, never tried here)

Screw axes and glide planes involve a **fractional translation**, so their
representation carries a k-dependent phase and sublattice/orbital components are
*forced* to interchange as k moves. That is structurally a winding d-vector, and
it would be **enforced by the space group rather than fine-tuned**. Investigate:

- In non-symmorphic space groups, what is the generic form of a two-band `H(k)`
  near the zone-boundary planes where non-symmorphic degeneracies are enforced?
  Does `|d(k)|` tend to constant magnitude with winding direction, and over how
  much of the BZ — a whole plane, or only along lines?
- Known flat-band or narrow-band systems in non-symmorphic groups. Any published
  quantum-metric or superfluid-weight calculation for one?
- Non-symmorphic band sticking, hourglass fermions, Dirac nodal lines — does any
  produce a manifold with near-uniform gap to everything else?
- Which non-symmorphic space groups are most promising, and are there real
  materials in them with narrow bands near E_F?

**The sharp version**: non-symmorphic symmetry guarantees *degeneracy* at the
zone boundary. What is needed here is the opposite — a *uniform gap* with a
winding direction. Do those coexist, or does the enforced degeneracy at the
boundary destroy the uniform gap? Settle this first; it may kill the lead
immediately.

### Lead 2 — other mechanisms for constant-|d| winding

- **SOC in a non-centrosymmetric crystal.** Rashba/Dresselhaus d-vectors wind by
  construction. How close to constant magnitude can `|d|` get, and over what
  region? This is the closest thing to an already-realized winding gap.
- **Chern insulators / flux lattices**, where Berry curvature is spread uniformly
  — is `|d|` near-constant in any standard model (Haldane, Qi-Wu-Zhang)? QWZ has
  an explicit `d(k)`; compute its `|d|` uniformity and `M_trg` directly.
- **Moiré systems** — the isolation scale is independent of intra-layer hopping.
- Any published lattice model with an **exactly flat band and a uniform gap**.

## Calibration you must use — a constant gap alone is NOT the signature

This nearly caused a false positive here. Measured `unif_gap = min_k Δ / max_k Δ`:

| structure | unif_gap | M_trg |
|---|---|---|
| winding, n = 1, 2, 3 | **1.000000** | 0.25 / 1.00 / 2.25 |
| cubic t2g + atomic SOC, at its optimum | **0.911** | **0.018** |
| same, t → 0 | **0.980** | → 0 |

Atomic SOC produces a nearly constant gap too — of course it does, a
k-independent term gives a k-independent gap. **The discriminant is a constant
gap TOGETHER WITH a large metric.** Report both for every candidate, and do not
report a uniform gap as a hit on its own.

## What to return

For each mechanism: (i) does it give `|d|` near-constant with winding `d̂`
*without* fine-tuning; (ii) what enforces it — symmetry, topology, or accident;
(iii) `M_trg` and gap uniformity, computed where you can rather than argued;
(iv) any real material or published model that realizes it, with numbers.

Cite real papers. Say plainly where the literature is silent. The prior is that
symmetry does not hand this over for free — but if it does, it is the whole
question.

## Running in parallel — do not duplicate

`screen/constd.py` is screening all 728 candidate bands from the 1,772-material
JARVIS-WTB sweep for the constant-gap + large-metric signature in real Wannier
data. That covers "does it already exist in a known material". Yours is "can
symmetry enforce it in principle".
