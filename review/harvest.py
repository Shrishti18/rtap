"""
Spec harvester for RTAP flat-band candidates.

Screens Wannier90 _hr.dat against the derived specification:
    Tc = min( U*lam/4 ,  BKT/3D-XY stiffness bound )
    lam = lam_max(A),  A_ab = <rho_a rho_b>_k,  rho_a(k) = P_aa(k)
    lam <= 1/2 at uniform pairing (exact); breaking it always loses.

Quantities per band:
    M        = int tr g dk /(2pi)^d      (gauge-invariant, = Omega_I/V_cell)
    lam      = lam_max(A)                (amplitude channel; Tc_MF = U*lam/4)
    W        = bandwidth
    d_iso    = min(gap above, gap below)
    unif     = max_a |w_a - 1/nphi|      (uniform-pairing deviation)

g_ij = (1/2) Tr[d_i P  d_j P]. P(k) is gauge-invariant, so finite-differencing
P avoids all gauge/branch problems. No Wannierization gauge fixing required.
"""
import re
import numpy as np

KB_MEV = 0.086173  # meV/K
C_STIFF = 0.67     # D_raw = C_STIFF * U * <tr g>; ED-calibrated, 2D and 3D alike.


# ---------------- I/O ----------------

def read_hr(path):
    """Parse Wannier90 _hr.dat -> (R (nr,3) int, H (nr,n,n) complex, deg (nr,))."""
    with open(path) as f:
        lines = [l for l in f.read().splitlines()]
    n = int(lines[1].split()[0])
    nr = int(lines[2].split()[0])
    deg, i = [], 3
    while len(deg) < nr:
        deg += [int(x) for x in lines[i].split()]
        i += 1
    deg = np.array(deg[:nr], float)
    body = np.array([[float(x) for x in l.split()] for l in lines[i:] if l.strip()])
    R = body[:, :3].astype(int).reshape(nr, n * n, 3)[:, 0, :]
    hh = (body[:, 5] + 1j * body[:, 6]).reshape(nr, n, n)
    H = np.transpose(hh, (0, 2, 1))          # file is column-major in (m,n)
    return R, H, deg


def write_hr(path, R, H, deg=None, comment="synthetic"):
    nr, n, _ = H.shape
    deg = np.ones(nr) if deg is None else deg
    with open(path, "w") as f:
        f.write(comment + "\n%d\n%d\n" % (n, nr))
        for j in range(0, nr, 15):
            f.write(" ".join("%d" % d for d in deg[j:j + 15]) + "\n")
        for r in range(nr):
            for b in range(n):
                for a in range(n):
                    z = H[r, a, b]
                    f.write("%5d%5d%5d%5d%5d%14.8f%14.8f\n"
                            % (R[r, 0], R[r, 1], R[r, 2], a + 1, b + 1, z.real, z.imag))


# ---------------- Bloch Hamiltonian ----------------

def hk_grid(R, H, deg, nk, dim=2):
    """H(k) on a uniform grid. Returns (shape..., n, n)."""
    ks = [2 * np.pi * (np.arange(nk) + 0.5) / nk for _ in range(dim)]
    mesh = np.meshgrid(*ks, indexing="ij")
    phase = np.zeros(mesh[0].shape, complex)
    out = np.zeros(mesh[0].shape + H.shape[1:], complex)
    for r in range(R.shape[0]):
        phase = np.exp(1j * sum(mesh[d] * R[r, d] for d in range(dim)))
        out += phase[..., None, None] * (H[r] / deg[r])
    return out, [2 * np.pi / nk] * dim


# ---------------- Geometry ----------------

def projector(Hk, bands):
    """P(k) for the band manifold `bands` (list of indices)."""
    E, V = np.linalg.eigh(Hk)
    U = V[..., :, bands]                       # (...,n,nb)
    return np.einsum('...ia,...ja->...ij', U, U.conj()), E


def metric_trace(P, dks):
    """tr g(k) with g_ij = (1/2) Tr[d_i P d_j P], central differences."""
    d = len(dks)
    trg = np.zeros(P.shape[:-2])
    for i in range(d):
        dP = (np.roll(P, -1, i) - np.roll(P, 1, i)) / (2 * dks[i])
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', dP, dP))
    return trg


def descriptors(Hk, dks, band):
    P, E = projector(Hk, [band])
    trg = metric_trace(P, dks)
    d = len(dks)
    # int tr g d^dk / (2pi)  -> reads directly as multiples of the Chern floor |C|
    M = trg.mean() * (2 * np.pi) ** (d - 1)
    rho = np.real(np.einsum('...ii->...i', P))              # (...,norb)
    rf = rho.reshape(-1, rho.shape[-1])
    A = (rf.T @ rf) / rf.shape[0]
    lam = float(np.linalg.eigvalsh(A)[-1])
    w = rf.mean(0)
    nphi = 1.0 / np.sum(w ** 2)
    unif = float(np.max(np.abs(w - 1.0 / len(w))))
    Eb = E[..., band]
    W = float(Eb.max() - Eb.min())
    nb = E.shape[-1]
    gl = float(Eb.min() - E[..., band - 1].max()) if band > 0 else np.inf
    gu = float(E[..., band + 1].min() - Eb.max()) if band < nb - 1 else np.inf
    return dict(M=float(M), M_trg=float(trg.mean()), lam=lam, W=W,
                d_iso=float(min(gl, gu)), unif=unif, nphi=float(nphi), w=w,
                Emid=float(Eb.mean()))


# ---------------- Spec ----------------

def spec(dsc, Tc_K=300.0, dim=2, R=1.0):
    """
    Required U and metric threshold for a target Tc. Energies in eV.

    ALL thresholds are on M_trg = <tr g>_BZ (raw BZ mean), NOT on the
    (2pi)^(d-1)-scaled `M`. The two differ by 6.28x in 2D and 39.5x in 3D.

    Calibrated (flat-band BdG, constant-|d| model, both 2D and 3D):
        D_raw = c * U * M_trg ,  c = 0.67  (dimension independent)
    XY mapping (pair momentum 2q => grad theta = 2q => J = D_raw/4):
        2D lattice XY: Tc = 0.89 J = 0.2225 D_raw
        3D lattice XY: Tc = 2.20 J = 0.55   D_raw
    Amplitude cap: Tc = U * lam / 4, with lam = 1/n_phi (= 1/2 at n_phi = 2).

    R in (0,1] is the mediator recovery factor: R = 1 for an instantaneous
    (electronic) attraction; for a boson of frequency w0, R is set by
    lam_LF = U/(2 w0). Measured by ED: R = 0.58 / 0.81 / 0.93 / 0.98
    at lam_LF = 1.0 / 0.5 / 0.25 / 0.125.
    """
    Tc = Tc_K * KB_MEV / 1000.0
    lam = dsc['lam']
    U_req = 4 * Tc / lam
    coef = 0.2225 * C_STIFF if dim == 2 else 0.55 * C_STIFF
    M_req = lam / (4 * coef * R)
    Mt = dsc['M_trg'] if 'M_trg' in dsc else dsc['M'] / (2*np.pi)**(dim-1)
    return dict(Tc_K=Tc_K, U_req=U_req, iso_req=2 * U_req,
                M_req_trg=M_req, M_trg=Mt, R=R,
                pass_iso=dsc['d_iso'] >= 2 * U_req,
                pass_M=Mt >= M_req,
                calibrated=True,
                **tc_max(dsc, dim, R))


def tc_max(dsc, dim=3, R=1.0, U_cap=None, w0=None):
    """
    Best achievable Tc for a band or manifold.

        Tc = min( U*lam/4 , coef*R*U*M_trg )
        U  <= min( d_iso/2 ,  U_cap ,  2*w0 )

    THREE caps on U, not one. My earlier version used d_iso/2 alone, which
    lets U run to 2 eV whenever the manifold is well isolated and produces
    nonsense (thousands of K). The other two:

      U_cap : the physically available screened interaction. cRPA values run
              ~2-4 eV bare for 3d, ~1-2 eV for 4d/5d. Default 1.5 eV.
      2*w0  : antiadiabaticity. A mediator of frequency w0 cannot supply an
              effectively instantaneous U much beyond ~2*w0 -- and this is
              the SAME physics R encodes, so leaving it out double-counted
              the optimism. Defaults are tied to R:
                  R=1.00 -> electronic, w0 ~ 1.5 eV  (cap inactive)
                  R=0.85 -> hydrogen phonon, 250 meV -> U <= 0.50 eV
                  R=0.55 -> oxide phonon,     90 meV -> U <= 0.18 eV
    """
    coef = 0.2225 * C_STIFF if dim == 2 else 0.55 * C_STIFF
    Mt = dsc['M_trg'] if 'M_trg' in dsc else dsc['M'] / (2*np.pi)**(dim-1)
    if w0 is None:
        w0 = {1.00: 1.5, 0.85: 0.250, 0.55: 0.090}.get(round(R, 2), 1.5)
    if U_cap is None:
        U_cap = 1.5
    caps = dict(iso=dsc['d_iso'] / 2.0, phys=U_cap, adiab=2.0 * w0)
    which = min(caps, key=caps.get)
    U_use = caps[which]
    amp = dsc['lam'] / 4.0
    stf = coef * R * Mt
    Tc_eV = U_use * min(amp, stf)
    return dict(Tc_max_K=Tc_eV / (KB_MEV / 1000.0),
                binds='amp' if amp < stf else 'stiff',
                U_used=U_use, U_capped_by=which)


def auto_bands(Hk, wmax=0.5, iso_min=0.2):
    """Candidate isolated narrow SINGLE bands. See manifold_bands for why this
    is usually the wrong object for an n_phi=2 doublet."""
    E = np.linalg.eigvalsh(Hk)
    nb = E.shape[-1]
    out = []
    for b in range(nb):
        Eb = E[..., b]
        W = Eb.max() - Eb.min()
        gl = Eb.min() - E[..., b - 1].max() if b > 0 else np.inf
        gu = E[..., b + 1].min() - Eb.max() if b < nb - 1 else np.inf
        iso = min(gl, gu)
        if W <= wmax and iso >= iso_min:
            out.append((b, float(W), float(iso)))
    return out


# ---------------- minimal quantum metric (Huhtinen et al.) ----------------

def _trg_shifted(P, dks, kgrids, shifts):
    """tr g under orbital shifts d_a, commutator form (no phase applied to P).
    dP' = V(dP + i[D,P])V^dag  =>  g_ij = (1/2)Tr[A_i A_j], A_i = d_iP + i[D_i,P].
    A_i is periodic because P is: no zone-boundary wrap, no grid-offset
    dependence. `kgrids` is unused and kept for signature compatibility."""
    d = len(dks)
    A = []
    for i in range(d):
        dP = (np.roll(P, -1, i) - np.roll(P, 1, i)) / (2 * dks[i])
        Di = np.diag(shifts[:, i]).astype(complex)
        comm = (np.einsum('ab,...bc->...ac', Di, P)
                - np.einsum('...ab,bc->...ac', P, Di))
        A.append(dP + 1j * comm)
    trg = np.zeros(P.shape[:-2])
    for i in range(d):
        trg += 0.5 * np.real(np.einsum('...ij,...ji->...', A[i], A[i]))
    return trg.mean()


def minimal_metric(Hk, dks, band, nk, dim, centres=None, restarts=6, seed=0):
    """M_min over orbital-position assignments. Naive value (shifts=0) is only
    an UPPER BOUND. Convention: u_a -> e^{+i k.d_a} u_a, so seeding from
    Wannier centres r_a requires d_a = -r_a (convention I -> II)."""
    from scipy.optimize import minimize as _min
    bands = band if isinstance(band, (list, tuple)) else [band]
    P, _ = projector(Hk, list(bands))
    norb = P.shape[-1]
    naive = float(_trg_shifted(P, dks, None, np.zeros((norb, dim))))

    def obj(x):
        s = np.zeros((norb, dim))
        s[1:] = x.reshape(norb - 1, dim)
        return float(_trg_shifted(P, dks, None, s))

    best, bestx = naive, np.zeros((norb - 1) * dim)
    starts = [np.zeros((norb - 1) * dim)]
    if centres is not None:
        c = np.asarray(centres, float)
        rel = (c[1:] - c[0]).ravel()
        starts += [-rel, rel]
    rng = np.random.default_rng(seed)
    starts += [rng.uniform(-1, 1, (norb - 1) * dim) for _ in range(restarts)]
    for x0 in starts:
        r = _min(obj, x0, method='Powell',
                 options=dict(xtol=1e-6, ftol=1e-8, maxiter=20000))
        if r.fun < best:
            best, bestx = float(r.fun), r.x
    return dict(M_min=best, M_naive=naive, shifts=bestx.reshape(norb - 1, dim))


def read_wout_centres(path):
    """Final Wannier centres (Cartesian, Angstrom) from wannier90.wout.
    Convert to FRACTIONAL and negate before passing to minimal_metric."""
    txt = open(path, errors='ignore').read()
    blocks = re.findall(r'Final State(.*?)(?:Sum of centres|\Z)', txt, re.S)
    if not blocks:
        return None
    rows = re.findall(r'WF centre and spread\s+\d+\s*\(\s*([-\d.]+)\s*,\s*'
                      r'([-\d.]+)\s*,\s*([-\d.]+)\s*\)', blocks[-1])
    return np.array([[float(a) for a in r] for r in rows]) if rows else None


# ---------------- multi-band manifolds (the actual target) ----------------

def manifold_bands(Hk, nb=2, wmax=0.5, iso_min=0.05, ef=None, window=2.0):
    """
    Contiguous nb-band manifolds isolated from the REST of the spectrum.

    `wmax` gates the WIDEST INDIVIDUAL BAND in the manifold, not the manifold's
    total span. This matters and my first version had it wrong: the target is a
    SPLIT doublet -- two narrow bands separated by a gap -- because an unsplit
    manifold spanning the same orbitals has P = subspace identity and tr g = 0
    identically. Gating on the span (hi-lo) therefore rejects exactly the
    structures we want. The constant-|d| reference has two exactly flat bands
    split by 2D: span 2D (large), max band width 0 (ideal).
    """
    E = np.linalg.eigvalsh(Hk)
    nbnd = E.shape[-1]
    out = []
    for b in range(nbnd - nb + 1):
        sub = E[..., b:b + nb]
        lo, hi = sub.min(), sub.max()
        if ef is not None and not (lo - window <= ef <= hi + window):
            continue
        gl = lo - E[..., b - 1].max() if b > 0 else np.inf
        gu = E[..., b + nb].min() - hi if b + nb < nbnd else np.inf
        iso = min(gl, gu)
        wb = max(float(E[..., j].max() - E[..., j].min())
                 for j in range(b, b + nb))
        if wb <= wmax and iso >= iso_min:
            out.append((tuple(range(b, b + nb)), wb, float(iso),
                        float(hi - lo)))
    return out


def descriptors_manifold(Hk, dks, bands):
    """descriptors() for a rank-len(bands) projector."""
    E, V = np.linalg.eigh(Hk)
    U = V[..., :, list(bands)]
    P = np.einsum('...ia,...ja->...ij', U, U.conj())
    trg = metric_trace(P, dks)
    rho = np.real(np.einsum('...ii->...i', P))
    rf = rho.reshape(-1, rho.shape[-1])
    A = (rf.T @ rf) / rf.shape[0]
    w = rf.mean(0)
    Es = E[..., list(bands)]
    nbnd = E.shape[-1]
    b0, b1 = bands[0], bands[-1]
    gl = float(Es.min() - E[..., b0 - 1].max()) if b0 > 0 else np.inf
    gu = float(E[..., b1 + 1].min() - Es.max()) if b1 + 1 < nbnd else np.inf
    r = len(bands)
    lam_raw = float(np.linalg.eigvalsh(A)[-1])
    # rho_a sums to r, so A ~ r^2: lam_raw is inflated by r^2 vs the
    # single-band case. lam_eff = lam_raw/r^2 reduces to 1/n_phi at r=1.
    return dict(M=float(trg.mean() * (2 * np.pi) ** (len(dks) - 1)),
                M_trg=float(trg.mean()), M_trg_per_band=float(trg.mean() / r),
                lam=lam_raw / r ** 2, lam_raw=lam_raw, rank=r,
                nphi=float(r ** 2 / np.sum(w ** 2)),
                W=max(float(E[..., j].max() - E[..., j].min()) for j in bands),
                W_span=float(Es.max() - Es.min()), d_iso=float(min(gl, gu)),
                unif=float(np.max(np.abs(w - r / len(w)))),
                w=w, Emid=float(Es.mean()))
