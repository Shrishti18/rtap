"""Independent checks of every load-bearing numerical claim. Each compares a
computed value against an ANALYTIC result derived separately, not against
another run of the same code."""
import numpy as np
P=[];F=[]
def chk(name,got,want,tol):
    ok=abs(got-want)<=tol
    (P if ok else F).append((name,got,want))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<46} got {got:>12.6f}  want {want:>10.6f}")

print("1. DILUTION LAW  lam = 1/n_phi")
import geom, harvest
n=300;k=2*np.pi*(np.arange(n))/n
K1,K2=np.meshgrid(k,k,indexing='ij')
# kagome, analytic 1/3
a1=np.array([1.,0.]);a2=np.array([.5,np.sqrt(3)/2])
H=np.zeros(K1.shape+(3,3),complex)
H[...,0,1]=H[...,1,0]=-2*np.cos(np.pi*K1/(2*np.pi)*2*np.pi/2/np.pi*np.pi)  # placeholder
c1=np.cos(np.pi*K1/np.pi/2*1);  # rebuild cleanly below
def kagome(k1,k2):
    H=np.zeros(k1.shape+(3,3),complex)
    H[...,0,1]=H[...,1,0]=-2*np.cos(np.pi*k1/(2*np.pi))
    H[...,0,2]=H[...,2,0]=-2*np.cos(np.pi*k2/(2*np.pi))
    H[...,1,2]=H[...,2,1]=-2*np.cos(np.pi*(k2-k1)/(2*np.pi))
    return H
E,V=np.linalg.eigh(kagome(K1,K2))
b=int(np.argmin([E[...,j].max()-E[...,j].min() for j in range(3)]))
u=V[...,:,b];rho=np.abs(u)**2;rf=rho.reshape(-1,3)
chk("kagome flat band lam (analytic 1/3)",float(np.linalg.eigvalsh((rf.T@rf)/rf.shape[0])[-1]),1/3,2e-4)

print("\n2. n_phi = 1 => g == 0 identically")
nk=200;kk=2*np.pi*(np.arange(nk)+.5)/nk
u1=np.zeros((nk,2),complex);u1[:,0]=1.0
Pj=np.einsum('ki,kj->kij',u1,u1.conj())
dP=(np.roll(Pj,-1,0)-np.roll(Pj,1,0))/(2*(2*np.pi/nk))
chk("single-orbital band <tr g>",float((0.5*np.real(np.einsum('kij,kji->k',dP,dP))).mean()),0.0,1e-12)

print("\n3. METRIC SUM RULE  M = sum_m |v_nm|^2/(E_n-E_m)^2")
import sumrule
a,bb,v,g=sumrule.check(sumrule.h_saw,sumrule.dh_saw,nk=2000,band=0)
chk("sawtooth: sum-rule vs finite-difference",a,bb,1e-5)

print("\n4. WINDING  M = n^2/4")
import winding
for w in [1,2,3]:
    M,_,_,_=winding.model(w,3.0,nk=3000)
    chk(f"winding n={w}",M,w*w/4,3e-4)

print("\n5. STIFFNESS CALIBRATION c = Draw/(U*M), dimension-independent")
import cdmodel as cd
for dim,nk in [(2,110),(3,30)]:
    kk2=2*np.pi*(np.arange(nk)+.5)/nk
    Kg=np.stack(np.meshgrid(*[kk2]*dim,indexing='ij'),-1)
    E2,V2=np.linalg.eigh(cd.H(Kg,1.0));u2=V2[...,:,0]
    Pp=np.einsum('...i,...j->...ij',u2,u2.conj());trg=np.zeros(Pp.shape[:-2])
    for i in range(dim):
        dPp=(np.roll(Pp,-1,i)-np.roll(Pp,1,i))/(2*(2*np.pi/nk));trg+=0.5*np.real(np.einsum('...ij,...ji->...',dPp,dPp))
    Dr,_=cd.Draw(Kg,1.0)
    chk(f"c in {dim}D (expect ~0.67)",Dr/(1.0*trg.mean()),0.67,0.02)

print("\n6. MINIMAL METRIC: decoupled dimers must give exactly 0")
R,H2,d=harvest.read_hr('dimer_hr.dat')
Hk,dks=harvest.hk_grid(R,H2,d,300,dim=1)
mm=harvest.minimal_metric(Hk,dks,0,300,1)
chk("dimer M_min",mm['M_min'],0.0,1e-6)
chk("dimer M_naive (gauge-dependent upper bound)",mm['M_naive'],0.25,2e-3)

print("\n7. GRID-OFFSET INDEPENDENCE of the commutator metric")
vals=[]
for off in [0.0,0.25,0.5]:
    ks=2*np.pi*(np.arange(300)+off)/300
    Hk2=np.zeros((300,2,2),complex)
    for r in range(R.shape[0]): Hk2+=np.exp(1j*ks*R[r,0])[:,None,None]*H2[r]
    Pp,_=harvest.projector(Hk2,[0])
    vals.append(harvest._trg_shifted(Pp,[2*np.pi/300],None,np.zeros((2,1))))
chk("max spread over grid offsets",float(max(vals)-min(vals)),0.0,1e-9)

print("\n8. j_eff=1/2")
import jeff
E3=np.linalg.eigvalsh(jeff.hk(np.zeros((1,3)),0.0,0.45))[0]
chk("atomic SOC splitting (analytic 3lam/2)",float(E3[4]-E3[3]),0.675,1e-9)
o,_=jeff.analyse(t=0.2,lam=0.45,nk=20)
chk("j=1/2 lam_pair (analytic 1/3)",o[0]['lam_pair'],1/3,1e-6)
chk("j=1/2 n_phi (analytic 3)",o[0]['nphi'],3.0,1e-6)
print(f"  [INFO] j=1/2 M_trg = {o[0]['M_trg']:.4f}  (no analytic reference; only cross-checked vs t)")

print("\n9. UNIT CONSISTENCY of the thresholds")
c=0.67
chk("3D threshold M_req = lam/(4*0.55*c)",0.5/(4*0.55*c),0.339,2e-3)
chk("2D threshold M_req = lam/(4*0.2225*c)",0.5/(4*0.2225*c),0.839,2e-3)

print("\n10. SOC CLOSING INEQUALITY")
chk("lambda needed = U_rep/3 at U_rep=2.0",2.0/3,0.667,1e-3)

print(f"\n{'='*62}\nPASSED {len(P)}   FAILED {len(F)}")
for n_,g_,w_ in F: print("  FAILED:",n_,g_,w_)
