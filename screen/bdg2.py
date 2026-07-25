import numpy as np
from scipy.optimize import minimize
import geom

PAULI = [np.array([[0,1],[1,0]], complex),
         np.array([[0,-1j],[1j,0]], complex),
         np.array([[1,0],[0,-1]], complex)]

def hflat(nvec, Egap):
    nd = sum(nvec[...,i,None,None]*PAULI[i] for i in range(3))
    return Egap*(np.eye(2, dtype=complex) - nd)/2.0

def bdg_eigs(hp, hm, Delta):
    """TRS singlet BdG: M = [[h(k+q), D],[D, -h(k-q)]]. Spin-down fixed by TRS."""
    sh = hp.shape[:-2]
    M = np.zeros(sh+(4,4), complex)
    M[...,:2,:2] = hp
    M[...,2:,2:] = -hm
    D = np.diag(Delta).astype(complex)
    M[...,:2,2:] = D
    M[...,2:,:2] = D.conj().T
    return np.linalg.eigvalsh(M)

def make_h(nfun, KX, KY, q, Egap, mu):
    I = np.eye(2)
    hp = hflat(nfun(KX+q, KY), Egap) - mu*I
    hm = hflat(nfun(KX-q, KY), Egap) - mu*I
    return hp, hm

def free_energy(Delta, hp, hm, U):
    E = bdg_eigs(hp, hm, Delta)
    return -0.5*np.abs(E).sum(-1).mean() + (Delta**2).sum()/U

def solve(nfun, KX, KY, q, Egap, U, mu=0.0, x0=None):
    hp, hm = make_h(nfun, KX, KY, q, Egap, mu)
    x0 = np.array([0.3, 0.3]) if x0 is None else x0
    r = minimize(free_energy, x0, args=(hp, hm, U), method='Nelder-Mead',
                 options=dict(xatol=1e-10, fatol=1e-13, maxiter=20000, maxfev=20000))
    return r.x, r.fun

def Ds(nfun, KX, KY, Egap, U, mu=0.0, qs=(0.0, 0.015, 0.03, 0.045)):
    es, x = [], None
    for q in qs:
        x, e = solve(nfun, KX, KY, q, Egap, U, mu, x0=x)
        es.append(e)
    c = np.polyfit(np.array(qs), np.array(es), 2)
    return 2*c[0], x
