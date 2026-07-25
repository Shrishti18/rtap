# Adversarial review — defects ranked by whether they change a conclusion

Setup reproduced: `vh.py` gives M = 5.0679 / 1.9073, lam = 0.5000 / 0.6125 exactly.
`VALIDATE_ALL.py` = 17 PASS / 0 FAIL. All findings below are in the gaps that
suite does not cover.

---

## TIER 1 — changes a conclusion

### D1. `jeff.analyse` measures a gauge-dependent quantity of a degenerate pair. C7's M_trg is not unconverged — it is **divergent**.

**What is wrong.** Bands 4 and 5 are an *exactly* degenerate Kramers pair:
`max|E5 - E4| = 2.7e-15` over the whole BZ (also true of 0-1 and 2-3). For a
degenerate subspace the individual eigenvector returned by LAPACK is arbitrary,
so the rank-1 projector `P = u u†` built in `analyse` is not a property of the
Hamiltonian. Only the rank-2 projector onto the pair is.

Two independent symptoms confirm it:

| nk | 16 | 20 | 24 | 32 | 40 | 48 | 64 |
|---|---|---|---|---|---|---|---|
| M_trg single band | 2.496 | 3.824 | 5.346 | 9.452 | 14.69 | 21.12 | **37.25** |
| M_trg rank-2 manifold | 0.4965 | 0.5184 | 0.5312 | 0.5447 | 0.5512 | 0.5549 | 0.5585 |

Log-log slope of the single-band series is **+1.95**, i.e. M_trg ∝ nk², the
signature of an eigenvector that re-randomises between adjacent k-points
(dP/dk ~ 1/dk, tr g ~ nk²). It does not converge and has no limit.

It is also gauge-dependent at *fixed* nk. Applying a random U(2) rotation inside
the degenerate pair at nk=20: M_trg = 3.750, 3.813, 3.828, 3.750, 3.807 —
while `lam_pair` and `nphi` stay at **exactly** 0.333333 and 3.000000.

**Converged value.** The rank-2 manifold metric does converge, at order q≈2
(Richardson stable to 4 digits): **M_trg = 0.563**.

**Affects C7.** Split verdict:
- `lam = 1/3` and `n_phi = 3` — **survive**, exactly gauge-invariant (the j=1/2
  doublet is orbitally isotropic, 1/3 on each t2g orbital, so any mixing
  preserves them). VALIDATE §8 tests these and is right to pass them.
- `M_trg ~ 4-5` — **does not survive**. It is an artefact of the grid. The
  correct number is 0.563, **7-9x smaller** than quoted.

**Does the conclusion survive?** Partially. 0.563 still clears the 3D threshold
0.339, so the sign of the verdict is unchanged — but the margin collapses from
~15x to **1.66x**, which is inside the uncertainty of `c = 0.67` itself. Any
claim resting on M_trg being comfortably above threshold is unsupported.

### D2. The j=1/2 manifold is **not isolated** at the physical Ir⁴⁺ parameters. C7's premise fails before its metric matters.

**What is wrong.** `analyse` computes the metric of bands (4,5) as though they
were a separated manifold. They are not. At t = 0.4 eV, λ = 0.45 eV:

```
iso = E4.min() - E3.max() = -2.50 eV      (NEGATIVE: j=1/2 and j=3/2 overlap)
```

The cause is a scale mismatch the model never checks: the atomic SOC splitting is
3λ/2 = 0.675 eV, while the orbital bandwidth is W = 8t = 3.2 eV. Confirmed
numerically (single-orbital W = 3.1901 at nk=40 → 3.2). **C7 quotes the atomic
splitting 3λ_SOC/2 as if it were the band gap.** In the solid it is not.

Scanning for a regime where the manifold *is* isolated shows the tradeoff is
fatal, not merely tight:

| t | λ | 3λ/2 | iso | M_trg (rank-2) |
|---|---|---|---|---|
| 0.4 | 0.45 | 0.675 | **-2.498** | 0.531 |
| 0.4 | 1.0 | 1.500 | -1.673 | 0.161 |
| 0.4 | 2.0 | 3.000 | -0.173 | 0.045 |
| 0.4 | 3.0 | 4.500 | **+1.327** | **0.020** |
| 0.05 | 0.45 | 0.675 | **+0.278** | **0.014** |

Isolation needs roughly 3λ/2 > 8t, i.e. λ ≳ 5.3t. Everywhere isolation exists,
M_trg has collapsed to 0.014–0.045, **8–24x below** the 0.339 threshold. There is
no point in this model's parameter space with both an isolated j=1/2 manifold and
M_trg above threshold.

**Affects C7 and the premise of C8.** The docstring of `jeff.py` asks exactly the
right question — "is SOC buying isolation at the cost of stiffness?" — and the
answer the code contains, but does not report, is **yes, catastrophically**.

**Does the conclusion survive?** No. C7 as stated (an isolated j=1/2 manifold with
M_trg ~ 4-5) describes no accessible parameter regime.

### D3. `cs2cos2.hk_chain` omits the ligand p-p overlap *and* the entire π channel. Its flat bands are an artefact.

**What is wrong.** Two distinct omissions in the inter-site hopping
`tAB = -J Σ_S outer(wsig(S-A), wsig(S-B))`:

1. **No ⟨p_A|p_B⟩ factor.** A two-step d→p→d hop must carry the overlap between
   the ligand p orbital pointing at A and the one pointing at B. That factor is
   ~cos(θ_bridge). The computed bridge angle here is **83.88°**, cos = 0.107.
2. **π channel absent.** `wpi()` is defined at lines 59-62 and **never called**.
   π enters only the on-site crystal field via `epi`, never the hopping.

Quantified against a proper sum over all three ligand p orbitals (‖t‖_F, J = 1):

| construction | ‖t_AB‖_F | vs code |
|---|---|---|
| code (σ, no p-overlap) | 1.0000 | — |
| proper σ only | 0.1066 | 0.107 |
| proper σ + π (t_π/t_σ = 0.5) | 0.7542 | 0.754 |

The two errors **partially cancel** in magnitude — net overestimate only ~1.33x
rather than 9.4x — which is presumably why this survived. But they do not cancel
in *structure*:

| construction | rank of t_AB | singular values |
|---|---|---|
| code | **2** | 1.483, 0.517, 0, 0, 0 |
| proper σ only | 2 | 0.158, 0.055, 0, 0, 0 |
| proper σ + π | **5** | 0.711, 0.577, 0.512, 0.335, 0.322 |

**Three of the five d orbitals get exactly zero inter-site dispersion in the code
and are therefore flat by construction.** With π restored the hopping is full
rank and no orbital is artificially flat. At 84° — squarely Goodenough-Kanamori
near-90° territory — π is the *dominant* channel (0.754 vs 0.107), which is the
worst possible place to omit it.

**Affects** any flat band, metric or n_phi extracted from the Cs₂CoS₂ model.
**Does it survive?** No. Flatness in this model must be re-derived with the π
channel before it means anything.

---

## TIER 2 — weakens support without flipping a verdict

### D4. The sawtooth ED-vs-BdG agreement does not reproduce. C5's `c = 0.67` has never been checked against an exact method.

The docstring predicts `S_ED = D_raw / N`. Running the shipped code (U = 2.0,
mu = E_flat = -2.0, D_raw = 0.36832):

| N | S_ED | S_ED/(D_raw/N) | with f(1-f), f=1/N | S_ED·N²/D_raw |
|---|---|---|---|---|
| 4 | 0.050201 | 0.545 | 2.908 | 2.181 |
| 5 | 0.028923 | 0.393 | 2.454 | 1.963 |
| 6 | 0.020750 | 0.338 | 2.434 | 2.028 |
| 8 | 0.011599 | 0.252 | 2.303 | 2.015 |

Neither normalisation gives ≈1, and **both drift monotonically with N** — so the
reported "ratio 0.98–1.03 at N=4" is not reproducible from these files, and even
if it were at one N it would be coincidence, since the ratio is not N-stable.

What the data *does* show is a clean power law: log-log slope **-2.09**, and
`S_ED·N²/D_raw = 2.02 ± 0.06` across N = 4..8. So the true relation is
**S_ED = 2·D_raw/N²**, not D_raw/N. That is the expected scaling for a *single
pair* on an N-cell ring (the pair carries charge 2e and sees flux 2θ, giving
E ~ (2θ/N)²/2M*), whereas `bdg_draw` is a thermodynamic half-filled quantity.
The comparison is between two different fillings and different limits.

**Affects C5.** The `c = 0.67` calibration itself still passes VALIDATE §5
(0.6659 in 2D, 0.6750 in 3D) — but that check is BdG-against-geometry, i.e.
entirely *within* mean field. The sawtooth ED was the only intended
mean-field-vs-exact validation, and it does not work. **C5 survives as a
mean-field statement; it is unsupported as a statement about the true stiffness.**

### D5. C8's inequality is not the one the rest of the project uses.

With Δ_iso = 3λ/2 (atomic), the project's own projection condition U ≤ Δ_iso/2
gives **λ ≥ 4U/3**. C8 instead uses **λ > U/3** — looser by a factor of **4**.
The stated inequality is not derivable from the stated projection rule, and no
derivation for it appears in the shipped files (`ceiling.py` does not contain it;
VALIDATE §10 is the bare arithmetic `2.0/3 = 0.667`).

**Does the conclusion survive?** Yes, and it is strengthened. Under the project's
own rule, Ir⁴⁺ would need λ ≥ 2.67 eV against a measured ~0.42 eV — short by
6.4x rather than 1.6x. The verdict "window closed" holds under either inequality.

### D6. C8's input numbers check out; the reviewer's inversion scenario does not.

Literature check on both sides of the comparison:

- **λ_SOC(Ir⁴⁺)**: RIXS gives **417 meV** (Sr₂IrO₄) and **400 meV** (Sr₃Ir₂O₇);
  the general quoted range is 0.4–0.5 eV. The project uses 0.45 eV — about
  **8% high** relative to the directly measured values. Using 0.417 makes the
  window slightly *more* closed.
- **U_rep(Ir 5d t₂g)**: cRPA gives **2.54 eV** and **1.96 eV** in different
  studies. The project's 2.0 eV sits at the bottom of that range, i.e. it is the
  conservative choice.

Inverting C8 requires U_rep < 3λ = 3 × 0.417 = **1.25 eV**, which is **36% below
the lowest published cRPA value**. The "if U_rep is really ~1.2 eV" scenario is
not supported by the literature. **C8's conclusion is robust** — more robust than
the author claims.

### D7. `Tc = U·lam/4` is correct, but it is the half-filling result and nothing enforces that.

Re-derived independently from the linearised gap equation. For a flat band at μ,
the projected pair interaction is V_kk' = -U Σ_a ρ_a(k) ρ_a(k'), and
tanh(ξ/2T)/2ξ → 1/4T as ξ → 0. Writing Δ(k) = Σ_a c_a ρ_a(k):

```
c_a = (U/4Tc) Σ_b A_ab c_b ,  A_ab = <rho_a rho_b>   =>   Tc = U·lam_max(A)/4
```

which is exactly `harvest.spec`. And for uniform k-independent weights over n_phi
orbitals, λ_max = n_phi·(1/n_phi²) = 1/n_phi, confirming **C1**. Both are sound.

The gap: the step `tanh(ξ/2T)/2ξ → 1/4T` requires ξ = 0, i.e. the flat band sits
*exactly* at μ (half filled). Away from half filling the pair susceptibility is
cut off and Tc falls. **No filling factor appears anywhere in `spec` or
`tc_max`.** Every Tc in the project is therefore an at-half-filling upper bound.
This is the same physics D4's f(1-f) was groping at.

### D8. `U ≤ Δ_iso/2` is too loose to justify the projection.

The projected model's leading correction from virtual excitation out of the
manifold is O(U²/Δ_iso), i.e. **O(U/Δ_iso) = 50%** relative at U = Δ_iso/2. That
is not a controlled expansion. A defensible criterion is U ≲ Δ_iso/4 to Δ_iso/10.
Since Tc ∝ U linearly, tightening to Δ_iso/4 **halves every isolation-limited
Tc**. In the full-sweep data the `iso` cap binds most candidates, so this is not
a corner case.

---

## TIER 3 — no defect found (reported as requested)

- **Attack 1 — `jeff.hk` spin-orbital ordering: CLEAN.** `np.kron(Ho, I2)` and
  `np.einsum('...ij,ab->...iajb').reshape(6,6)` agree to **0.0** exactly. Both are
  orbital-major/spin-minor, matching `h_soc`'s `np.kron(L,S)`. The transposition
  bug hypothesised does not exist. Atomic splitting = 0.675 = 3λ/2 exactly.
- **Attack 2 — W = 8t: CONFIRMED** (3.1901 at nk=40 → 3.2). But note the default
  `tp=0.0` makes **each t2g orbital strictly 2D** (d_xy has exactly zero
  kz-dependence, verified 0.0). The docstring says tp "is what makes the lattice
  genuinely 3D for every orbital" — yet `analyse` defaults to tp=0. The full 6×6
  is still 3D because SOC mixes the orbitals, so this is not itself an error, but
  the default contradicts the stated intent.
- **Attack 3 — `_trg_shifted`: CLEAN.** Derivation verified independently: V and D
  are both diagonal so they commute, giving V(dP + i[D,P])V† = V dP V† + i[D,P']
  = dP'. Numerically A_i is Hermitian to **0.0**, [D,P] anti-Hermitian to **0.0**,
  tr g real to **0.0**. The commutator form converges O(h²)
  (1.2889 → 1.3779 over nk = 32..1024); the explicit phase-then-roll form differs
  by an amount that **grows ∝ nk and is localised at the BZ boundary** (interior
  discrepancy falls 0.168 → 0.0009 as nk goes 32 → 512), exactly the defect the
  commutator form was introduced to remove. The README's claim is correct.
- **Attack 4 — twist convention: CLEAN.** Factor chain re-derived: pair momentum
  2q ⇒ ∇θ = 2q ⇒ F = F₀ + 2ρ_s q² ⇒ d²E/dq² = 4ρ_s, and `Draw` returns
  `2·polyfit(...)[0]` = d²E/dq², so **ρ_s = D_raw/4** ✓. Lattice-XY constants
  0.893 (2D BKT) and 2.202 (3D) are the standard values, and using the bare
  stiffness with 0.893 is consistent with the universal jump. Numerically E(q) is
  even to **1e-15**; the author's one-sided fit differs from a symmetric fit by
  **8.3e-4** relative; c = 0.6696. **Every Tc is correctly scaled.**
- **Attack 7 — I/O: essentially clean.** `write_hr`/`read_hr` round-trip: R and
  deg **exact**, H to **6.7e-9** — limited by the `%14.8f` write format, not a
  logic error. The (a,b) transpose convention is self-consistent (file order is
  b-outer/a-inner, `reshape` then `transpose(0,2,1)` inverts it correctly).
- **C3 sum rule: derivation correct.** M = Σ_m |v_nm|²/(E_n-E_m)² ≤ (1/Δ²)Σ_m
  |v_nm|² = ⟨v²⟩/Δ² follows immediately from |E_n - E_m| ≥ Δ.
- **C2: correct.** n_phi = 1 forces ρ_a = δ_{a,a₀} for a.e. k, so P is a constant
  projector and dP ≡ 0. (VALIDATE §2 tests this by *constructing* a constant
  projector, which is true by construction rather than a test of a real band —
  weak, but the claim itself is sound.)

**Weak tests worth noting:** VALIDATE §9 (`0.5/(4·0.55·0.67) = 0.339`) and §10
(`2.0/3 = 0.667`) are arithmetic on definitions, not independent checks. §10 in
particular validates neither of the two numbers C8 actually rests on.

---

## Summary table

| Claim | Verdict |
|---|---|
| **C1** lam = 1/n_phi | **Survives.** Re-derived from the linearised gap equation. |
| **C2** n_phi=1 ⇒ tr g ≡ 0 | **Survives.** Analytically sound; test is vacuous but the claim is right. |
| **C3** M ≤ ⟨v²⟩/Δ² | **Survives.** |
| **C4** M = n²/4 | **Survives.** |
| **C5** D_raw = 0.67·U·M | **Survives as mean-field.** Unsupported as exact — the ED cross-check (D4) does not reproduce. |
| **C6** thresholds 0.339 / 0.839 | **Survives** as arithmetic given C5; inherits C5's mean-field-only status. |
| **C7** j_eff=1/2: lam=1/3, M_trg~4-5, gap=3λ/2 | **lam=1/3 survives. M_trg~4-5 FAILS** (divergent; true value 0.563). **gap=3λ/2 FAILS** — it is the atomic splitting, and the real band gap is **-2.5 eV**. |
| **C8** window closed for Ir⁴⁺ | **Conclusion survives** and is robust to the literature values. **The stated inequality λ > U/3 is unsupported** and inconsistent with the project's own U ≤ Δ_iso/2 (which is 4x stricter and closes the window harder). |

**The single most consequential finding is D2**: at the physical Ir⁴⁺ point the
j=1/2 manifold overlaps j=3/2 by 2.5 eV, and every parameter choice that opens a
real gap drives M_trg 8–24x below threshold. C7 describes no reachable regime.
