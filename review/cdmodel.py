import numpy as np
from scipy.optimize import minimize
# |d(k)| = const  =>  H = d.sigma has H^2 = D^2 I  =>  TWO EXACTLY FLAT BANDS at +-D,
# gap 2D, no touching anywhere, d-hat free to wind => nonzero metric.
# d = D(cos k1, sin k1 cos k2, sin k1 sin k2): finite Fourier, |d|=D exactly.
SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex); SZ=np.array([[1,0],[0,-1]],complex)
def dvec(K):
    k1,k2=K[...,0],K[...,1]
    return np.stack([np.cos(k1), np.sin(k1)*np.cos(k2), np.sin(k1)*np.sin(k2)],-1)
def H(K,D=1.0,shift=True):
    d=D*dvec(K)
    h=d[...,0,None,None]*SX+d[...,1,None,None]*SY+d[...,2,None,None]*SZ
    return h+ (D*np.eye(2) if shift else 0)   # lower band at 0, upper at 2D
def geom(K,dk):
    E,V=np.linalg.eigh(H(K))
    u=V[...,:,0]; P=np.einsum('...i,...j->...ij',u,u.conj())
    trg=np.zeros(P.shape[:-2])
    for i in range(K.shape[-1]):
        dP=(np.roll(P,-1,i)-np.roll(P,1,i))/(2*dk)
        trg+=0.5*np.real(np.einsum('...ij,...ji->...',dP,dP))
    rho=np.real(np.einsum('...ii->...i',P)).reshape(-1,2)
    A=(rho.T@rho)/rho.shape[0]
    return dict(M=float(trg.mean()),lam=float(np.linalg.eigvalsh(A)[-1]),w=rho.mean(0),
                W=float(E[...,0].max()-E[...,0].min()),gap=float(np.min(E[...,1]-E[...,0])))
def E0(Delta,Hp,Hm,U):
    sh=Hp.shape[:-2]; M=np.zeros(sh+(4,4),complex)
    M[...,:2,:2]=Hp; M[...,2:,2:]=-Hm
    Dm=np.diag(Delta).astype(complex); M[...,:2,2:]=Dm; M[...,2:,:2]=Dm
    return -0.5*np.abs(np.linalg.eigvalsh(M)).sum(-1).mean()+(Delta**2).sum()/U
def Draw(K,U,D=12.5,qs=(0.,0.02,0.04,0.06)):
    es=[]; x=np.array([.3,.3]); e=np.zeros(K.shape[-1]); e[0]=1.
    for q in qs:
        Hp=H(K+q*e,D); Hm=H(K-q*e,D)
        r=minimize(E0,x,args=(Hp,Hm,U),method='Nelder-Mead',
                   options=dict(xatol=1e-10,fatol=1e-13,maxiter=8000,maxfev=8000))
        x=r.x; es.append(r.fun)
    return 2*np.polyfit(np.array(qs),np.array(es),2)[0], x
