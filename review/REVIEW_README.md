# Review package: flat-band superconductivity derivation

Run `python3 VALIDATE_ALL.py` first. Requires numpy + scipy only.
17 checks, all currently passing, each against an ANALYTIC reference.

## Load-bearing claims and where they live

| | claim | implemented in |
|---|---|---|
| C1 | `lam = 1/n_phi` -> forces n_phi=2, Tc <= U/8 | `harvest.descriptors`, `geom.py`, `nphi.py` |
| C2 | `n_phi=1 => tr g == 0` identically | `gauge.py`, VALIDATE §2 |
| C3 | `M <= <v^2>/Delta^2` (exact sum rule) | `sumrule.py` |
| C4 | `M = n^2/4` from winding number | `winding.py` |
| C5 | `D_raw = c*U*M`, c=0.67, dimension-independent | `cdmodel.py`, `ed2.py`, `sawtooth.py` |
| C6 | thresholds `M>=0.339` (3D) / `0.839` (2D) | `harvest.spec`, `harvest.tc_max` |
| C7 | j_eff=1/2: lam=1/3, M_trg ~4-5, gap=3*lam_SOC/2 | `jeff.py` |
| C8 | window open iff `lam_SOC > U_rep/3`; Ir4+ fails | `ceiling.py`, arithmetic only |

## KNOWN OPEN DEFECT
`jeff.analyse` M_trg = 5.3 at nk=24 but 3.75 at nk=20. UNCONVERGED.
This is the only quantity in the project with no analytic reference and it
was quoted to three figures. Converge it.

## Files by role
- **core library**: `harvest.py` (descriptors, minimal metric, spec/tc_max, I/O),
  `geom.py`, `cdmodel.py`, `bdg2.py`
- **claims**: `sumrule.py` C3, `winding.py` C4, `jeff.py` C7, `ceiling.py` C8
- **BdG vs exact-diagonalisation cross-checks**: `ed.py`, `ed2.py`, `sawtooth.py`, `audit*.py`
- **materials models**: `cs2cos2.py` (real chain from CIF), `target3d.py`
  (constructed lattice), `eg3d.py`, `fbmodel.py`
- **chemistry / crystal field**: `cfield.py`, `aom.py`, `chem.py`, `shortlist.py`,
  `cobaltate.py`, `spinstate.py`, `bridgeU.py`
- **database screen** (largely superseded): `cifscan.py`, `ptgrp.py`, `mkcif.py`
- **superseded / historical, do not review unless tracing an error**:
  `bdg.py` (pre-TRS-fix), `mech.py` (superseded by `lf.py`), `main.py`,
  `final.py`, `opt.py`, `robust.py`
- **test fixtures**: `*_hr.dat`, `mktest.py`, `vh.py`, `demo.py`

## Conventions that have bitten before
- `harvest.descriptors` returns BOTH `M` (scaled by (2pi)^(d-1)) and `M_trg`
  (raw <tr g>). **All thresholds are on `M_trg`.** Mixing them is a 39.5x
  error in 3D and caused a real bug earlier.
- Orbital shift convention: `u_a -> exp(+i k.d_a) u_a`, so seeding from
  Wannier centres r_a needs `d_a = -r_a`.
- `_trg_shifted` uses the commutator form, NOT phase-then-roll. The old form
  had a BZ-boundary term and depended on the k-grid offset.
- BdG twist: pair momentum is `2q`, hence `J = D_raw/4`.
