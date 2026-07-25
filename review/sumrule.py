import numpy as np
# EXACT: for an isolated band n,  g_ij = sum_{m!=n} Re[v^i_nm v^j_mn]/(E_n-E_m)^2
# with v_nm = <u_n|dH/dk|u_m>. So the metric is NOT free -- it is velocity
# matrix elements divided by gaps squared. Hence the bound
#     M = <tr g>  <=  <(dH/dk)^2>_n / Delta_iso^2
# LARGE METRIC REQUIRES SMALL GAP. And the projection needs U <= Delta/2.
# The two demands are in direct opposition. That is the pattern.
def check(hfun, dhfun, nk=2000, band=0):
    k=2*np.pi*(np.arange(nk)+0.5)/nk
    H=hfun(k); dH=dhfun(k)
    E,V=np.linalg.eigh(H)
    n=band; nb=E.shape[-1]
    u=V[...,:,n]
    # sum rule
    trg_sr=np.zeros(nk); vsq=np.zeros(nk)
    for m in range(nb):
        if m==n: continue
        um=V[...,:,m]
        v=np.einsum('ki,kij,kj->k',np.conj(u),dH,um)
        trg_sr+=np.abs(v)**2/(E[...,n]-E[...,m])**2
        vsq+=np.abs(v)**2
    # direct finite-difference metric
    P=np.einsum('ki,kj->kij',u,u.conj())
    dP=(np.roll(P,-1,0)-np.roll(P,1,0))/(2*(2*np.pi/nk))
    trg_fd=0.5*np.real(np.einsum('kij,kji->k',dP,dP))
    gap=np.min([np.abs(E[...,m]-E[...,n]).min() for m in range(nb) if m!=n])
    return trg_sr.mean(), trg_fd.mean(), vsq.mean(), gap
# --- sawtooth (t'=sqrt2 t), flat band ---
T=1.0; TP=np.sqrt(2)
def h_saw(k):
    H=np.zeros(k.shape+(2,2),complex)
    H[...,0,0]=2*T*np.cos(k); H[...,0,1]=TP*(1+np.exp(-1j*k))
    H[...,1,0]=np.conj(H[...,0,1]); return H
def dh_saw(k):
    H=np.zeros(k.shape+(2,2),complex)
    H[...,0,0]=-2*T*np.sin(k); H[...,0,1]=TP*(-1j)*np.exp(-1j*k)
    H[...,1,0]=np.conj(H[...,0,1]); return H
# --- constant-|d| 1D (two flat bands, TRS) ---
D=5.0
SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex)
def h_cd(k): return D*(np.cos(k)[:,None,None]*SX+np.sin(k)[:,None,None]*SY)+D*np.eye(2)
def dh_cd(k): return D*(-np.sin(k)[:,None,None]*SX+np.cos(k)[:,None,None]*SY)
print(f"{'model':>14} {'M(sumrule)':>11} {'M(findiff)':>11} {'<v^2>':>9} {'gap':>7} {'<v^2>/gap^2':>12}")
for nm,hf,df,b in [('sawtooth',h_saw,dh_saw,0),('const-|d| 1D',h_cd,dh_cd,0)]:
    a,bb,v,g=check(hf,df,band=b)
    print(f"{nm:>14} {a:11.6f} {bb:11.6f} {v:9.3f} {g:7.3f} {v/g**2:12.6f}")
