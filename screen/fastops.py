"""Exact drop-in replacements for the harvest routines that dominate cost.

Validated against the shipped implementations in validate_fastops.py to ~1e-14;
these change the arithmetic, not the result.

hk_grid_fast
    harvest.hk_grid accumulates `out += phase[...,None,None] * H[r]` over every
    R-vector, allocating a full (nk^dim, norb, norb) temporary each time
    (133 x 1.9 GB of traffic for a 60-orbital material at nk=32). The same sum
    is one matmul: phases (nk^dim, nR) @ H (nR, norb*norb).

trg_from_u / minimal_metric_fast
    For a SINGLE band the projector is rank one, P = u u^dag, so tr g never
    needs the (nk^dim, norb, norb) projector -- only the (nk^dim, norb)
    eigenvector. Each Powell evaluation becomes a reduction over that array,
    milliseconds rather than tens of seconds, which is what makes optimising
    ~150 shift parameters feasible at all.

    Follows harvest v4's commutator form: A_i = d_iP + i[D_i, P]. No phase is
    applied to P, so there is no zone-boundary wrap term and no dependence on
    the k-grid offset. (The v3 phase-then-roll form needed an explicit
    2h - 2pi correction on the wrapped slices; that whole issue is gone.)
"""
import numpy as np


def hk_grid_fast(R, H, deg, nk, dim=2):
    nr, n, _ = H.shape
    ks = [2 * np.pi * (np.arange(nk) + 0.5) / nk for _ in range(dim)]
    mesh = np.meshgrid(*ks, indexing="ij")
    kflat = np.stack([m.ravel() for m in mesh], 1)            # (nk^dim, dim)
    ph = np.exp(1j * (kflat @ R[:, :dim].T.astype(float)))    # (nk^dim, nr)
    Hf = (H / deg[:, None, None]).reshape(nr, n * n)
    out = (ph @ Hf).reshape(mesh[0].shape + (n, n))
    return out, [2 * np.pi / nk] * dim


def band_vectors(Hk, band):
    """Normalised eigenvector of `band` on the whole grid -> (..., norb)."""
    _, V = np.linalg.eigh(Hk)
    return V[..., :, band]


def trg_from_u(u, dks, shifts=None):
    """<tr g> for a single band from the rank-1 identity, commutator form.

    With A_i = d_iP + i[D_i, P], D_i = diag(shifts[:, i]), P = u u^dag, and
        ov  = <u(k+e_i) | u(k-e_i)>
        s_p = <u | u(k+e_i)>,   t_p = <u(k+e_i)| D_i |u>   (likewise for -)
        var = <d^2> - <d>^2  over the weights |u_a|^2
    the three pieces of Tr[A_i A_i] are

        Tr[(d_iP)^2]       = (2 - 2|ov|^2) / (2h)^2
        2i Tr[d_iP [D,P]]  = -4 (Im(s_p t_p) - Im(s_m t_m)) / (2h)
        -Tr[[D,P]^2]       = 2 var

    and tr g = sum_i (1/2) Tr[A_i A_i]. Nothing larger than (nk^dim, norb) is
    ever formed.
    """
    tot = 0.0
    for i, h in enumerate(dks):
        up = np.roll(u, -1, i)
        um = np.roll(u, 1, i)
        ov = np.einsum('...a,...a->...', np.conj(up), um)
        acc = (2.0 - 2.0 * np.abs(ov) ** 2) / (2 * h) ** 2
        if shifts is not None:
            d = np.asarray(shifts, float)[:, i]
            du = d * u
            sp = np.einsum('...a,...a->...', np.conj(u), up)
            tp = np.einsum('...a,...a->...', np.conj(up), du)
            sm = np.einsum('...a,...a->...', np.conj(u), um)
            tm = np.einsum('...a,...a->...', np.conj(um), du)
            acc = acc - 4.0 * (np.imag(sp * tp) - np.imag(sm * tm)) / (2 * h)
            w = np.abs(u) ** 2
            m1 = np.einsum('...a,a->...', w, d)
            m2 = np.einsum('...a,a->...', w, d * d)
            acc = acc + 2.0 * (m2 - m1 ** 2)
        tot += float(np.mean(0.5 * acc))
    return tot


def minimal_metric_fast(Hk, dks, band, dim, centres=None, restarts=2, seed=0,
                        maxiter=20000, u=None, maxfev=None):
    """M_min / M_naive in raw <tr g> units. Same contract as
    harvest.minimal_metric; `u` avoids repeating eigh, `maxfev` bounds Powell."""
    from scipy.optimize import minimize as _min
    if u is None:
        u = band_vectors(Hk, band)
    norb = u.shape[-1]
    naive = trg_from_u(u, dks)

    def obj(x):
        s = np.zeros((norb, dim))
        s[1:] = x.reshape(norb - 1, dim)
        return trg_from_u(u, dks, s)

    best, bestx = naive, np.zeros((norb - 1) * dim)
    starts = [np.zeros((norb - 1) * dim)]
    if centres is not None:
        # v4's D = diag(d_a) is the sign-flip of gauge.py's `shift`: on the
        # decoupled dimer the exact zero sits at d_B = -1, not +1. Seed BOTH
        # signs so the physical start cannot land on the wrong side.
        c = np.asarray(centres, float)
        rel = (c[1:] - c[0]).ravel()
        starts.append(-rel)
        starts.append(rel)
    rng = np.random.default_rng(seed)
    starts += [rng.uniform(-1, 1, (norb - 1) * dim) for _ in range(restarts)]
    for x0 in starts:
        opts = dict(xtol=1e-6, ftol=1e-8, maxiter=maxiter)
        if maxfev:
            opts["maxfev"] = maxfev
        r = _min(obj, x0, method="Powell", options=opts)
        if r.fun < best:
            best, bestx = float(r.fun), r.x
    return dict(M_min=best, M_naive=naive, shifts=bestx.reshape(norb - 1, dim))


def descriptors_fast(E, V, dks, band):
    """Same contract as harvest.descriptors, from a precomputed eigh."""
    u = V[..., :, band]
    d = len(dks)
    trg = trg_from_u(u, dks)          # raw <tr g>_BZ -- the screening quantity
    M = trg * (2 * np.pi) ** (d - 1)  # legacy scaled value
    rho = np.abs(u) ** 2
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
    return dict(M=float(M), M_trg=float(trg), lam=lam, W=W,
                d_iso=float(min(gl, gu)), unif=unif, nphi=float(nphi), w=w,
                Emid=float(Eb.mean()))
