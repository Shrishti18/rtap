# Does the winding structure exist in real materials?

Screen of all 605 unique candidate bands (179 materials) from the 1,772-material
JARVIS-WTB sweep, for the signature that turns the 8.1 K terminus into 334 K.

## What the signature actually is

A constant gap **alone is not it**. Measured `unif_gap = min_k Δ / max_k Δ`:

| structure | unif_gap | M_trg |
|---|---|---|
| winding, n = 1, 2, 3 | **1.000000** | 0.25 / 1.00 / 2.25 |
| cubic t2g + atomic SOC, at its optimum | **0.911** | **0.018** |
| same, t → 0 | **0.980** | → 0 |

Atomic SOC gives a nearly constant gap too — a k-independent term gives a
k-independent gap. **The discriminant is a constant gap TOGETHER WITH a large
metric.** Screening on uniformity alone would have returned the whole
atomic-SOC family as false positives.

## Defect found and fixed mid-screen

`constd.py` takes the gap on whichever side has the smaller mean, which for a
near-degenerate pair is the splitting **inside** the doublet, not the isolation
gap. LiPbAu₂ 24/25 has `gap_mean = 0.024` against `d_iso = 0.390`. The corrected
isolation-gap uniformity (gap from the whole group to the nearest band outside
it, pointwise minimum over both sides) is in the table below as `unif_iso`. It
moves several candidates substantially: FePSe₃ 0.94 → 0.45, CoCl₂O₈ 0.82 → 0.32,
RbW₃Cl₉ 0.83 → 0.34, K₄IrO₄ 0.93 → 0.54; and up for LiPbAu₂ 2/3, 0.73 → 0.92.
**`constd.csv`'s `unif_gap` column is the uncorrected quantity.**

## Result 1 — the trade-off is measurable in real data

```
corr(unif_gap, M_trg) = -0.42        over 605 bands
```

Constant gap and large metric are **anti-correlated** in real materials. This is
the isolation-vs-metricity trade-off appearing empirically for the first time in
this project, rather than as a model result.

Counts: `unif > 0.5`: 95/605. `M_trg > 0.339`: 472/605. **Both: 23/605.**

## Result 2 — the structure exists, and n_φ is what is scarce

Of the 23, almost all have `n_φ` between 10 and 26 — the pair is diluted over
10–26 orbitals, so the amplitude cap `Tc = U/(4 n_φ)` binds and kills them. The
winding ideal is `n_φ = 2`.

**Exactly two bands in 605 have large metric and `n_φ < 4`:**

| | LiPbAu₂ band 24 | band 25 |
|---|---|---|
| space group | 71 (Immm) | 71 |
| `n_φ` | **3.05** | **2.97** |
| `M_trg` | 0.547 | 0.540 |
| `d_iso` | 0.390 eV | 0.390 eV |
| `unif_iso` | 0.608 | 0.608 |
| `W_band` | 0.240 eV | 0.242 eV |
| Tc (amplitude cap, binds) | **185 K** | **190 K** |
| Tc (geometric cap) | 456 K | 450 K |

That is the first real material in this project to clear 100 K on its own
descriptors.

## Result 3 — and it fails on filling, like everything else

`E_F = −3.569 eV`; the band runs `[−4.793, −4.552]`. **E_F sits 0.98 eV above the
band top — 4.1 band widths.** `filling = 1.0`.

All 23 candidates fail the same way, missing E_F by **4 to 180 band widths**.
Every one is completely full or completely empty.

## Result 4 — geometry and filling are statistically independent

| quantity | corr with `ef_inside` |
|---|---|
| `unif_gap` | +0.124 |
| `M_trg` | −0.087 |
| `1/n_φ` | −0.019 |

```
P(geometry)  = 23/605 = 0.038
P(at E_F)    = 12/605 = 0.020
P(both) if independent = 0.46 bands       observed = 2
```

Only **1 band of 605** has filling in (0.05, 0.95): HfFeCl₆ JVASP-6766 band 41,
`f = 0.917`, `M_trg = 0.497`, `n_φ = 4.86`, `d_iso = 0.151` → **45 K** (amplitude
cap binds). That is the lone survivor of the original E_F gate.

## What this changes

The original screen reported the E_F gate as a death: 727/728 eliminated. That
framing treats a **tunable** condition as a fatal one. The independence above
says the correct reading is different:

- **Geometry is not tunable** — it is set by the crystal structure and orbital
  content. Finding `n_φ ≈ 3`, `M_trg > 0.5`, `d_iso = 0.39 eV` is rare (2/605).
- **Filling is tunable** — by doping, substitution, gating, or pressure.

So the search should be **geometry-first, then engineer E_F onto it**, not
screen for both at once and report the joint failure. LiPbAu₂ is the concrete
instance: the geometry that gives 185 K exists in a real, catalogued compound,
1.0 eV from where it would need to be.

## What is NOT established

- **The winding number of LiPbAu₂ 24/25 has not been computed.** `M_trg = 0.547`
  is consistent with winding (n = 1 gives 0.25 in 1D, n = (1,1,1) gives 0.74 in
  3D) but the large metric may have another origin. This is the next check and
  it is cheap.
- `unif_iso = 0.608` is not 1.0. This is not a clean winding structure, only a
  partially uniform one.
- Moving E_F by 0.98 eV in an Au-based metal means heavy doping, and LiPbAu₂ has
  many other bands at E_F that would carry the current instead. No claim is made
  that this compound superconducts.
- All standing caveats apply: `c = 0.67` is mean-field and uncalibrated; `U` is
  granted, never verified; `U ≤ d_iso/2` is the loose projection rule (D8 argues
  for `/4`, which halves every Tc above).

## Files

- `constd.py` — the screen (eigenvalues only, 605 bands, 179 materials, 623 s)
- `constd.csv` — per-band `unif_gap` (uncorrected), `M_trg`, `n_φ`, `d_iso`, gap

---

# RETRACTION: LiPbAu₂ is not 185 K, and is not winding

Both claims in the section above fail on direct computation (`lipbau2.py`).

## The 185 K was a rank-1/rank-2 mismatch — D1's fifth appearance

The stored Tc pairs `M_trg` and `lam` from the **single-band** descriptors of
band 24 with `d_iso = 0.390 eV`, the isolation gap of the two-band **group**.
Different manifolds. Band 24's own nearest neighbour is band 25 at 0.024 eV.

| treatment | U (eV) | M | lam | amp | geo | **Tc** |
|---|---|---|---|---|---|---|
| (a) rank-2 manifold, U ≤ iso/2 | 0.1952 | 1.0531 | 0.1931 | 109 | 879 | **109 K** |
| (b) single band, U ≤ pair_gap/2 | 0.0097 | 0.5472 | 0.3285 | 9 | 23 | **9 K** |
| STORED — (a)'s U with (b)'s M, lam | 0.1952 | 0.5472 | 0.3285 | 186 | 457 | ~~186 K~~ |

**Self-consistency picks (b).** Treatment (a) assumes both bands pair coherently,
which requires the intra-pair splitting ≲ Tc: 0.024 eV = 274 K against Tc = 109 K
fails by 2.5×. Treatment (b) assumes band 25 is negligible, requiring splitting
≫ Tc: 274 K against 9 K holds. **LiPbAu₂ gives 9 K.** The claim that it is the
first material here to clear 100 K is withdrawn entirely.

## The error is systematic across the whole stored scan

Every stored Tc uses single-band `M_trg`/`lam` with the group's `d_iso`.
Recomputing all 708 unique bands in the consistent manifold treatment:

```
stored   median 41.0 K   max 7615 K
manifold median 33.2 K   max 3044 K
ratio stored/manifold: median 1.063, 90th pct 4.67, max 374
inflated by >1.5x : 239/708      deflated (<0.67) : 44/708
```

Median impact is small (1.06), **but the tail is where every candidate lives** —
AgBr band 0 is inflated 26×, ZnPd₅Se 7.5×, Ba₃AsN 6.9×. The top of the stored
ranking is the part most corrupted by this.

## And it is not a winding structure

For any two-band `H = d₀I + d·σ` the rank-1 metric is exactly
`tr g = ¼Σ|∂d̂|²`. So if band 24's metric came from winding inside the 24/25
pair, an effective two-band model built from that pair would reproduce it:

```
|d(k)|  min 0.0097  mean 0.0118  max 0.0140   uniformity 0.694  (winding: 1.000)
M_trg, effective 2-band model : 0.1676
M_trg, full Hamiltonian       : 0.5472        ratio 0.306
```

**70% of the metric comes from bands outside the group**, ≥ 0.39 eV away —
consistent with the manifold cross term being only −0.034 (3%), since a rank-2
projector removes exactly the intra-pair channel. `|d|` is not constant either.
LiPbAu₂ is not an empirical instance of integer winding.

## Where the winding class now stands

- Integer winding **cannot be symmetry-enforced**: a non-symmorphic operation
  with fractional translation 1/m gives winding n = 1/m, and half-integer
  winding *is* the band-sticking degeneracy. Symmetry protection and isolation
  are mutually exclusive.
- So integer winding requires **fine-tuned hopping** — suppressed short-range,
  dominant long-range — which is why it is rare.
- The best of 605 real bands is not it.

**The winding class is unrealized in any known material.** It is not closed by
theorem — the 334 K for the model stands — but there is no material and no
symmetry route to one.

## Correction to the strategic reframe

"Filling is dopable" was stated too broadly. It holds for shifts of ~0.1–0.3 eV.
LiPbAu₂ needs E_F moved 0.98 eV — 4 band widths — in a metal with other bands
already at E_F. That is not doping, it is a different compound.
