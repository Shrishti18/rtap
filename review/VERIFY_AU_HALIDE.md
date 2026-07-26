# Verification of the Au-halide survivor-programme document

## 1. Their §1.1 correction to me is CORRECT, and sharper than they state

**Claim:** for two orbitals with `rho_1 = 1/2 + x`, `rho_2 = 1/2 - x`, `<x> = 0`,
`lam_max = 1/2` for **every mean-balanced texture**, not only pointwise-uniform
ones.

**Verified numerically** across five textures including extreme ones:

| texture | `<rho_1>` | `<x^2>` | `lam_max` |
|---|---|---|---|
| pointwise uniform | 0.5000 | 0.0000 | **0.500000** |
| 0.3 cos k | 0.5000 | 0.0450 | **0.500000** |
| 0.49 sign(cos k) | 0.5000 | 0.2401 | **0.500000** |
| uniform random ±0.5 | 0.5000 | 0.0828 | **0.500000** |
| 0.45 sin 3k | 0.5000 | 0.1012 | **0.500000** |

My brief said `lam = 1/n_phi` holds "only for uniform orbital weights". The true
condition is only that the **means** are balanced — strictly weaker, and it
matters for design.

**Extension (derived here, not in their document).** For mean imbalance `w` and
k-variance `v = <x^2>`, exactly:

    det A = v ,   tr A = w^2 + (1-w)^2 + 2v ,   lam_max = [T + sqrt(T^2 - 4v)]/2

| w | v = 0 | v = 0.05 | v = 0.15 |
|---|---|---|---|
| 0.50 | 0.5000 | 0.5000 | 0.5000 |
| 0.40 | 0.5200 | 0.5247 | 0.5445 |
| 0.35 | 0.5450 | 0.5549 | 0.5913 |
| 0.25 | 0.6250 | 0.6478 | 0.7153 |

So **k-variance is worthless on its own and only helps once the means are
imbalanced.** Design rule: imbalance first, texture second.

Their §1.2 (algebraic vs chemical separability), §1.3 (the SOC theorem closes
only two-scale Hamiltonians `t*H_hop + lam*H_soc`, not systems with independent
ligand or molecular scales) and §1.4 (declare rank before quoting anything) are
all fair and I accept them.

## 2. Their Gate 0 is the right test — and it has already been run

They define `eta` = exact/mean-field stiffness coefficient and preregister:

    eta >= 0.7  : programme plausible
    0.4 - 0.7   : need the 2x headroom envelope
    eta <  0.4  : "extremely unlikely"

**Measured this session, against two independent materials:**

| system | predicted | measured | factor | eta_eff |
|---|---|---|---|---|
| MATBG (BM model built here) | 16.3 K | **1.7 K** | 9.6x | **0.104** |
| BKBO (`ceiling.py`) | 235 K | **30 K** | 7.7x | **0.130** |

Both bind on the amplitude cap `U*lam/4`, so this calibrates the framework
slightly more broadly than their `eta`, which is defined on the geometric cap.

**`eta_eff = 0.10-0.13`, below their own 0.4 kill line.**

## 3. Their design vector fails under its own criterion

`U = 0.40 eV, lam = 0.60, M_min = 0.75, Delta_iso = 0.80 eV`:

| eta | amp (meV) | geo (meV) | Tc | |
|---|---|---|---|---|
| 1.000 raw formula | 60.0 | 110.7 | 696 K | |
| 0.500 their assumption | 30.0 | 55.4 | 348 K | pass |
| **0.130 BKBO-calibrated** | 7.8 | 14.4 | **91 K** | **fail** |
| **0.104 MATBG-calibrated** | 6.2 | 11.5 | **72 K** | **fail** |

Their 2x headroom requirement was the right instinct and is **5x too small**.

## 4. The consequence points away from their chosen family

Required U at `lam = 0.6`, calibrated:

    eta = 0.130  ->  U >= 1.33 eV
    eta = 0.104  ->  U >= 1.66 eV

That is **bismuthate scale**, and 5-9x above the 0.24-0.40 eV an Au d9 valence
skipper can plausibly supply. The calibration does not close the programme — it
collapses it onto the single mechanism with a ~1.9 eV attraction, the **s2
valence skippers**, which the sublattice-degeneracy argument reopened (one s
orbital on two symmetry-related sites gives n_phi = 2 without needing orbital
degeneracy).

So: right framework, right gates, right first calculation — **wrong family.**
A2Au2X6 should be demoted to the control experiment it is well suited to be
(Au(II) is chemically accessible; Cs4Au3Cl12 is a clean polaron negative
control), and the lead should move to an s2 skipper on two equivalent sites.

## 5. What I could not verify

- The `P4_2/mnm` Cs2Au2Cl6 claim (18.5 meV/f.u. above ground state,
  dynamically stable, 3D-connected Au-halide conduction states). Requires
  reading arXiv:2404.08465. **If true it is the strongest structural fact in
  their document** and survives the family demotion, because the same
  phase-retention argument would apply to an s2 analogue.
- The ambient-retention report for a quenched cubic phase (Bull. Chem. Soc.
  Jpn. 73, 1445).
