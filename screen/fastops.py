"""Exact drop-in replacements for the two harvest routines that dominate cost.

Both are validated against the shipped implementations in `validate_fastops.py`
to ~1e-12; they change the arithmetic, not the result.

hk_grid_fast
    harvest.hk_grid accumulates `out += phase[...,None,None] * H[r]` over every
    R-vector, allocating a full (nk^dim, norb, norb) temporary each time
    (133 x 1.9 GB of traffic for a 60-orbital material at nk=32). The same sum
    is one matmul: phases (nk^dim, nR) @ H (nR, norb*norb).

trg_shifted_fast / minimal_metric_fast
    For a SINGLE band the projector is rank one, P = u u^dag, so

        Tr[dP_i dP_i] = 2 - 2 |<u(k+e_i)|u(k-e_i)>|^2

    and tr g never needs the (nk^dim, norb, norb) projector at all -- only the
    (nk^dim, norb) eigenvector. Under an orbital shift u_a -> e^{i k.d_a} u_a
    the overlap picks up a phase set by the finite-difference step:

        <u+|u-> -> sum_a conj(u+_a) u-_a e^{-i Dk_i d_{a,i}}

    with Dk_i = k(idx+1) - k(idx-1) = 2h_i away from the zone boundary and
    2h_i - 2pi on the two wrapped slices (see _delta_k -- dropping that term
    does NOT reproduce the shipped result). Each Powell evaluation is then a
    reduction over a precomputed (nk^dim, norb) array -- milliseconds rather
    than tens of seconds.
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


def _overlap_products(u, dks):
    """conj(u(k+e_i)) * u(k-e_i) elementwise, per direction -> list of (...,norb)."""
    return [np.conj(np.roll(u, -1, i)) * np.roll(u, 1, i) for i in range(len(dks))]


def _delta_k(shape, dks):
    """Actual k(idx+1) - k(idx-1) per direction on the rolled grid.

    This is 2h everywhere except the two boundary slices, where the periodic
    roll wraps and the true difference is 2h - 2pi. The shipped _trg_shifted
    evaluates exp(i k.d) at the stored k values and then rolls, so it carries
    exactly this boundary term; ignoring it changes the answer.
    """
    out = []
    for i, h in enumerate(dks):
        nk = shape[i]
        ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk
        d = np.roll(ks, -1) - np.roll(ks, 1)      # (nk,)
        sh = [1] * len(dks)
        sh[i] = nk
        out.append(d.reshape(sh))
    return out


def trg_from_u(u, dks, shifts=None):
    """<tr g> for a single band, from the rank-1 identity. `shifts` fractional."""
    prods = _overlap_products(u, dks)
    dk = _delta_k(u.shape[:-1], dks)
    tot = 0.0
    for i, h in enumerate(dks):
        c = prods[i]
        if shifts is not None:
            c = c * np.exp(-1j * dk[i][..., None] * shifts[:, i])
        ov = c.sum(-1)
        tot += np.mean(1.0 - np.abs(ov) ** 2) / (2 * h) ** 2
    return float(tot)


def minimal_metric_fast(Hk, dks, band, dim, centres=None, restarts=2, seed=0,
                        maxiter=20000, u=None, maxfev=None):
    """M_min / M_naive, same contract as harvest.minimal_metric.

    `u` lets a caller pass the already-computed band eigenvector so eigh is not
    repeated; `maxfev` bounds the Powell budget for large orbital counts."""
    from scipy.optimize import minimize as _min
    if u is None:
        u = band_vectors(Hk, band)
    norb = u.shape[-1]
    prods = _overlap_products(u, dks)
    dk = _delta_k(u.shape[:-1], dks)
    inv = [1.0 / (2 * h) ** 2 for h in dks]

    def obj_from(shifts):
        tot = 0.0
        for i in range(len(dks)):
            ph = np.exp(-1j * dk[i][..., None] * shifts[:, i])
            ov = (prods[i] * ph).sum(-1)
            tot += np.mean(1.0 - np.abs(ov) ** 2) * inv[i]
        return float(tot)

    naive = obj_from(np.zeros((norb, dim)))

    def obj(x):
        s = np.zeros((norb, dim))
        s[1:] = x.reshape(norb - 1, dim)
        return obj_from(s)

    best, bestx = naive, np.zeros((norb - 1) * dim)
    starts = [np.zeros((norb - 1) * dim)]
    if centres is not None:
        c = np.asarray(centres, float)
        starts.append((c[1:] - c[0]).ravel())
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
    """Same contract as harvest.descriptors, but from a precomputed eigh.

    Avoids re-running eigh once per candidate band, and never forms the
    (nk^dim, norb, norb) projector: rho_a = |u_a|^2 and tr g comes from the
    rank-1 identity.
    """
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
