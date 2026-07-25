import numpy as np

# ---------- Bloch-vector models (two-band) ----------

def qwz(kx, ky, m):
    """Qi-Wu-Zhang. C=+1 for 0<m<2."""
    hx = np.sin(kx)
    hy = np.sin(ky)
    hz = m + np.cos(kx) + np.cos(ky)
    return np.stack([hx, hy, hz], axis=-1)

def haldane(kx, ky, t1=1.0, t2=0.15, phi=np.pi/2, M=0.0):
    a1 = np.array([1.0, 0.0])
    a2 = np.array([-0.5, np.sqrt(3)/2])
    a3 = np.array([-0.5, -np.sqrt(3)/2])
    b1 = a2 - a3; b2 = a3 - a1; b3 = a1 - a2
    hx = t1*sum(np.cos(kx*a[0] + ky*a[1]) for a in (a1, a2, a3))
    hy = t1*sum(np.sin(kx*a[0] + ky*a[1]) for a in (a1, a2, a3))
    hz = M - 2*t2*np.sin(phi)*sum(np.sin(kx*b[0] + ky*b[1]) for b in (b1, b2, b3))
    return np.stack([hx, hy, hz], axis=-1)

def checkerboard(kx, ky, t=1.0, tp=0.5, M=0.0):
    hx = 2*t*np.cos((kx+ky)/2)*np.cos((kx-ky)/2)
    hy = 2*t*np.sin((kx+ky)/2)*np.sin((kx-ky)/2)
    hz = M + 2*tp*(np.cos(kx) - np.cos(ky))
    return np.stack([hx, hy, hz], axis=-1)

# ---------- Mobius boost: conformal on S^2, preserves C and int tr g ----------

def mobius_boost(n, b):
    """Stereographic w -> b*w. Conformal => int tr g and C invariant."""
    n = n/np.linalg.norm(n, axis=-1, keepdims=True)
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    denom = 1.0 - nz
    small = denom < 1e-12
    denom = np.where(small, 1e-12, denom)
    w = (nx + 1j*ny)/denom
    w = b*w
    aw2 = np.abs(w)**2
    nz2 = (aw2 - 1)/(aw2 + 1)
    t = 2*w/(aw2 + 1)
    nx2, ny2 = np.real(t), np.imag(t)
    nz2 = np.where(small, 1.0, nz2)
    nx2 = np.where(small, 0.0, nx2)
    ny2 = np.where(small, 0.0, ny2)
    out = np.stack([nx2, ny2, nz2], axis=-1)
    return out/np.linalg.norm(out, axis=-1, keepdims=True)

# ---------- Geometry from n-hat ----------

def geometry(n, dk):
    """tr g = (1/4)|grad n|^2 ; Omega = (1/2) n.(dx n x dy n)."""
    dnx = (np.roll(n, -1, 0) - np.roll(n, 1, 0))/(2*dk)
    dny = (np.roll(n, -1, 1) - np.roll(n, 1, 1))/(2*dk)
    trg = 0.25*((dnx**2).sum(-1) + (dny**2).sum(-1))
    cross = np.cross(dnx, dny)
    omega = 0.5*(n*cross).sum(-1)
    return trg, omega

def rho(n):
    """Orbital weights rho_1,2 = (1 +/- n_z)/2 for |u> = (cos(th/2), sin(th/2)e^{i ph})."""
    nz = n[...,2]
    return np.stack([(1+nz)/2, (1-nz)/2], axis=-1)

def descriptors(n, dk, Nk):
    trg, om = geometry(n, dk)
    cell = dk*dk
    intg = trg.sum()*cell/(2*np.pi)          # = |C| lower bound
    C = om.sum()*cell/(2*np.pi)
    kappa = (trg**2).mean()/(trg.mean()**2)  # Cauchy-Schwarz >= 1
    r = rho(n)
    ipr = (r**2).sum(-1).mean()
    # pairing matrix M_ab = <rho_a rho_b>
    M = np.einsum('xya,xyb->ab', r, r)/(Nk*Nk)
    lam = np.linalg.eigvalsh(M)[-1]          # in units of U
    return dict(intg=intg, C=C, kappa=kappa, ipr=ipr, lam=lam, meanrho=r.mean((0,1)))
