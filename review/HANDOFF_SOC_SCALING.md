# Work order: is the ceiling exactly linear in λ_SOC?

Self-contained. Needs numpy only. Nothing from this repo is required.

## Why this matters

The project terminated at **8.1 K for Ir⁴⁺**. That was written up as if it
settled the class. It does not — it settles one ion. If the claim below holds,
the same computation settles **every ion in the periodic table** and states
exactly what spin-orbit constant room temperature would require. If it fails,
the terminus is ion-specific and the search space reopens.

## The claim

The model Hamiltonian is

```
H(k) = t * H_hop(k)  +  λ * H_soc
```

with `H_hop` and `H_soc` both dimensionless and fixed. Two consequences that
should be exact:

**(a)** `H/λ = (t/λ)·H_hop + H_soc`, so the *eigenvectors* depend only on the
ratio `t/λ`. Every quantity built from eigenvectors alone — `M_trg = <tr g>`,
`n_phi`, `lam_pair` — is a function of `t/λ` **only**, invariant under
`(t, λ) → (s·t, s·λ)`.

**(b)** Eigenvalues scale linearly, so `Δ_iso = λ · f(t/λ)` for a fixed `f`.

The binding constraint is the stiffness cap `Tc_geo = coef · U · M_trg` with
`coef = 0.55 × 0.67` (3D) and `U ≤ Δ_iso/2`. So

```
Tc_geo(t, λ) = coef · (λ/2) · f(t/λ) · M_trg(t/λ)
```

Maximising over `t` at fixed `λ`, the optimum sits at a **fixed ratio**
`(t/λ)*`, and

```
Tc_max = A · λ_SOC       A = (coef/2) · f((t/λ)*) · M_trg((t/λ)*)
```

**The ceiling is exactly linear in the spin-orbit constant**, with a single
universal coefficient `A` in K/eV.

## The model, written out

t2g basis order `xy, yz, zx`, each × spin up/down → 6 states. Orbital-major,
spin-minor.

```python
import numpy as np

LX = np.array([[0,0,0],[0,0,-1j],[0,1j,0]], complex)
LY = np.array([[0,0,1j],[0,0,0],[-1j,0,0]], complex)
LZ = np.array([[0,-1j,0],[1j,0,0],[0,0,0]], complex)
SX = np.array([[0,1],[1,0]], complex)
SY = np.array([[0,-1j],[1j,0]], complex)
SZ = np.array([[1,0],[0,-1]], complex)

def h_soc(lam):                      # -(lam/2) l_eff . sigma on t2g x spin
    return -0.5*lam*(np.kron(LX,SX) + np.kron(LY,SY) + np.kron(LZ,SZ))

def hop_t2g(k, t, tp=0.0):           # orbital-selective, cubic, a=1
    kx, ky, kz = k[...,0], k[...,1], k[...,2]
    H = np.zeros(k.shape[:-1]+(3,3), complex)
    H[...,0,0] = -2*t*(np.cos(kx)+np.cos(ky)) - 4*tp*np.cos(kx)*np.cos(ky)
    H[...,1,1] = -2*t*(np.cos(ky)+np.cos(kz)) - 4*tp*np.cos(ky)*np.cos(kz)
    H[...,2,2] = -2*t*(np.cos(kz)+np.cos(kx)) - 4*tp*np.cos(kz)*np.cos(kx)
    return H

def hk(k, t, lam, tp=0.0):
    Ho = hop_t2g(k, t, tp)
    H6 = np.einsum('...ij,ab->...iajb', Ho, np.eye(2)).reshape(Ho.shape[:-2]+(6,6))
    return H6 + h_soc(lam)
```

Descriptors for the j=1/2 doublet, **bands (4,5), rank-2 projector** — this
matters, see the trap list:

```python
def analyse(t, lam, tp=0.0, nk=32, bands=(4,5)):
    kk = 2*np.pi*(np.arange(nk)+0.5)/nk
    K = np.stack(np.meshgrid(kk,kk,kk,indexing='ij'), axis=-1)
    dk = 2*np.pi/nk
    E, V = np.linalg.eigh(hk(K, t, lam, tp))
    bl = list(bands)
    Um = V[..., :, bl]
    P = np.einsum('...ia,...ja->...ij', Um, Um.conj())
    trg = np.zeros(P.shape[:-2])
    for i in range(3):
        dP = (np.roll(P,-1,i) - np.roll(P,1,i)) / (2*dk)
        trg += 0.5*np.real(np.einsum('...ij,...ji->...', dP, dP))
    r = len(bl)
    rho = np.real(np.einsum('...ii->...i', P)).reshape(-1,6)
    orb = rho.reshape(-1,3,2).sum(-1)          # trace over spin -> xy,yz,zx
    Am = (orb.T @ orb)/orb.shape[0]
    w = orb.mean(0)
    gl = float(E[...,bl].min() - E[...,bl[0]-1].max())
    gu = float(E[...,bl[-1]+1].min() - E[...,bl].max()) if bl[-1] < 5 else np.inf
    return dict(M_trg=float(trg.mean()), iso=float(min(gl,gu)),
                nphi=float(r**2/np.sum(w**2)),
                lam_pair=float(np.linalg.eigvalsh(Am)[-1])/r**2)
```

## What to compute

1. **Invariance.** Confirm `M_trg`, `n_phi`, `iso/λ` are identical at fixed
   `t/λ` across at least a 20× range in λ. Report the residual.
2. **The coefficient.** For each λ, maximise `Tc_geo = coef·(iso/2)·M_trg/k_B`
   over `t` (scan `t/λ ∈ [0.05, 0.30]`). Report `Tc_max/λ` for λ spanning
   0.1 → 3.0 eV. Is it constant? To what precision? **That constancy is the
   whole claim.**
3. **The number.** Report `A` in K/eV, the optimal ratio `(t/λ)*`, and
   `iso/λ` and `M_trg` at the optimum.
4. **Room temperature.** What λ_SOC does Tc = 77 / 200 / 300 K require?
5. **The periodic table.** Apply `A·λ` to real atomic spin-orbit constants.
   **Include 6p — Tl, Pb, Bi.** The `t2g` manifold maps to an effective `l = 1`,
   which is exactly the structure of an atomic `p` shell under SOC, so the same
   model applies with λ → ζ_6p. Bi 6p is believed to be ~1.25–1.5 eV, roughly
   3× Ir. **Scoping this project to d-electrons was never justified, and this is
   the check that was skipped.** Source every ζ; mark estimates as estimates.

## What would falsify the claim

- `Tc_max/λ` not constant → the reduction to `t/λ` is wrong somewhere, and the
  terminus is ion-specific. **Report this loudly.**
- The optimal ratio `(t/λ)*` drifting with λ → same conclusion.
- `tp ≠ 0` (second-neighbour hopping) breaking the two-parameter structure.
  **Please test `tp/t = 0.1, 0.3` explicitly** — with `tp ≠ 0` there are three
  scales, and the invariance only survives if `tp` scales with `t`. Say what
  happens when it does not.

## Traps this project has already fallen into — do not repeat them

1. **Rank-1 projectors at a degeneracy are meaningless.** Bands 4 and 5 are a
   Kramers pair, degenerate to ~3e-15. A rank-1 projector built from one member
   is not a property of the Hamiltonian: `<tr g>` diverges as `nk²` and shifts
   under a random U(2) rotation inside the pair. Use the **rank-2** projector as
   coded above. An earlier `M_trg = 5.3` was this artefact.
2. **Non-periodic gauges.** If you build any other lattice for comparison, check
   `max|H(k) − H(k+2π)| ≈ 0` before differencing. A `2t·cos(k/2)` convention is
   *not* 2π-periodic and injects a false discontinuity at the BZ boundary — it
   shows up as a metric peak where the gap is *largest*, which is the tell.
3. **`M` vs `M_trg`.** All thresholds are on `M_trg = <tr g>`. A separate `M`
   scaled by `(2π)^(d−1)` exists elsewhere and mixing them is a 39.5× error in 3D.
4. **Check `iso > 0` before using it.** At large `t/λ` the j=1/2 doublet overlaps
   the j=3/2 quartet and `iso` goes negative; `U ≤ iso/2` is then unsatisfiable,
   not merely violated. Skip those points, don't clip them.

## Reference points to reproduce

At `λ = 0.417 eV` (RIXS value for Ir⁴⁺) the optimum should land near
`t = 0.052 eV`, giving `iso = 0.2115`, `M_trg = 0.0180`, `n_phi = 3.000`,
`lam_pair = 0.3333`, and `Tc_geo = 8.1 K`. If you do not reproduce those, stop
and report the discrepancy before going further.

## Not in scope for you

A separate agent is attacking the *premise* — whether atomic SOC is really the
only isolation scale that does not track hopping (moiré/superlattice, molecular
and CDW gaps, Kondo hybridization, p-block). Do not duplicate that. Your job is
the scaling law **given** the premise.
