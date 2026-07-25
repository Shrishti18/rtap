# RTAP flat-band spec harvester

Screens Wannier90 `_hr.dat` against the derived specification.

## Spec (derived, not assumed)
    Tc = min( U*lam/4 , stiffness bound )
    lam = lam_max(A),  A_ab = <rho_a rho_b>_k,  rho_a(k) = P_aa(k)

Exact result: under uniform pairing lam = 1/2 identically, so Tc_MF = U/8
regardless of band geometry. Breaking uniform pairing raises lam toward 1
(ceiling Tc = U/4) but destroys the metric integral faster -> Tc always falls.
Optimum is AT uniform pairing.

Targets for Tc = 300 K:
    U      >= 0.21 eV        (dimension-independent; excludes moire)
    d_iso  >= 2U ~ 0.42 eV   (projection antecedent)
    M      >= ~4.7  in 2D    (M = int tr g d^dk / 2pi, in units of the Chern floor)
                              ~5x looser in 3D (coefficient UNCALIBRATED)
    unif   ~ 0               (orbitals symmetry-equivalent)

Topology guarantees only M >= |C|, i.e. 1. That is ~5x too weak.
Target large-metric, symmetry-equivalent-orbital, NON-topological
(obstructed / RSI) flat bands -- not Chern bands.

## Usage
    import harvest
    R,H,deg = harvest.read_hr("wannier90_hr.dat")
    Hk,dks  = harvest.hk_grid(R,H,deg, nk=60, dim=3)
    for b,W,iso in harvest.auto_bands(Hk, wmax=0.5, iso_min=0.2):
        d  = harvest.descriptors(Hk,dks,b)
        sp = harvest.spec(d, Tc_K=300, dim=3)

Gauge-safe: P(k) is gauge invariant, so g_ij = (1/2)Tr[d_iP d_jP] by finite
difference needs no Wannier gauge fixing, no branch tracking, no smoothing.

## Validation (vh.py, demo.py)
Synthetic checkerboard `_hr.dat` reproduces the direct BdG run exactly:
    M   5.0679 vs 5.068   |  lam 0.5000 vs 0.5000
    M   1.9073 vs 1.907   |  lam 0.6125 vs 0.6125
Spec back-derives U_req = 0.207 eV from lam=0.5, matching the hand derivation.

## Scaling / caveats
- Memory ~ nk^dim * norb^2 * 16 B. 60^3 x 20 orb ~ 1.4 GB. Slice k3 if larger.
- C2D = 0.068 (D_s ~ C2D*U*M) is calibrated in-class in 2D. The 3D coefficient
  is NOT calibrated -- spec() flags this via `calibrated=False`. Run a 3D BdG
  calibration before trusting 3D M_req.
- No network here: point read_hr at local JARVIS/Wannier90 files.

## Pipeline (cheap -> expensive)
1. `cifscan.scan(cif)`   -- symmetry only, ~ms/structure, no DFT.
     PASS  = one site per PRIMITIVE cell carrying a 2-fold orbital doublet
     PASS* = same, and the site is a 4d/5d transition metal
     Excludes ALL triclinic/monoclinic/orthorhombic sites outright.
2. `harvest.descriptors/spec` -- needs a Wannier _hr.dat. Verifies on the
     BANDS what symmetry only permitted: is the manifold at E_F actually the
     doublet, is it narrow, is d_iso >= 2U, is M large enough.
3. `bdg2.Ds` -- full BdG, only for survivors.

## Known limits of step 1
- Symmetry cannot see the Fermi level. Sr in SrTiO3 and Na in Ba2NaOsO6 both
  PASS structurally and are closed-shell: step 2 kills them.
- Electron count must place E_F inside the doublet. For 4d/5d low spin that
  means d7-d9 (t2g filled, e_g partly filled). Ba2NaOsO6 is Os(7+) = 5d1,
  i.e. t2g -> n_phi = 3 -> diluted, despite PASS* at step 1.
- Needs explicit _symmetry_equiv_pos_as_xyz in the CIF (most have it); there
  is no space-group-number lookup table here.
