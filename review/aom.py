import numpy as np
# Angular Overlap Model, MX4 tetrahedron vs flattening, WITH pi bonding.
# E_i = sum_L [ e_sig |S_sig(i,L)|^2 + e_pi (|S_pix|^2+|S_piy|^2) ]
# Overlap factors from Wigner rotations of d-orbitals; standard AOM closed forms.
NAMES=['z2','x2-y2','xy','xz','yz']
def aom_energies(ligs, epi):
    # ligand at (theta,phi): standard AOM angular factors (Schaffer-Gliemann)
    E=np.zeros(5)
    for (th,ph) in ligs:
        c,s=np.cos(th),np.sin(th); C,S=np.cos(ph),np.sin(ph)
        c2p,s2p=np.cos(2*ph),np.sin(2*ph)
        # sigma overlaps
        Fs=np.array([ (3*c*c-1)/2,
                      (np.sqrt(3)/2)*s*s*c2p,
                      (np.sqrt(3)/2)*s*s*s2p,
                      np.sqrt(3)*s*c*C,
                      np.sqrt(3)*s*c*S ])
        # pi overlaps (two ligand pi orbitals)
        Fp1=np.array([ -np.sqrt(3)*s*c,
                       s*c*c2p,
                       s*c*s2p,
                       c*(2*c*c-1)/1*C - 0*S,   # standard: (2c^2-1) -> cos(2th)? use exact below
                       0 ])
        # Use exact Schaffer tables instead (safer): build from rotation matrices
        E+=0  # placeholder
    return E
# The hand-tabulated factors are error-prone; do it EXACTLY via rotation matrices.
def dmat(axis_angles):
    # real-d rotation matrix for ligand direction: rotate z-axis to (th,ph)
    th,ph=axis_angles
    # rotation R = Rz(ph) Ry(th); act on quadratic forms
    def Ry(t): return np.array([[np.cos(t),0,np.sin(t)],[0,1,0],[-np.sin(t),0,np.cos(t)]])
    def Rz(t): return np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
    R=Rz(ph)@Ry(th)
    # d-basis as symmetric traceless matrices
    s2,s6=np.sqrt(2),np.sqrt(6)
    B=[np.diag([-1,-1,2])/s6,                       # z2
       np.diag([1,-1,0])/s2,                        # x2-y2
       np.array([[0,1,0],[1,0,0],[0,0,0]])/s2,      # xy
       np.array([[0,0,1],[0,0,0],[1,0,0]])/s2,      # xz
       np.array([[0,0,0],[0,0,1],[0,1,0]])/s2]      # yz
    D=np.zeros((5,5))
    for j,Q in enumerate(B):
        Qr=R@Q@R.T
        for i,P in enumerate(B): D[i,j]=np.sum(P*Qr)
    return D
def energies(ligs,epi,esig=1.0):
    H=np.zeros((5,5))
    # in the LIGAND frame, sigma couples only to z2 (index 0), pi to xz,yz (3,4)
    diag=np.diag([esig,0,0,epi,epi])
    for L in ligs:
        D=dmat(L)
        H+=D@diag@D.T
    return H
def tet(beta):
    # 4 ligands at polar angle beta from z, alternating phi; beta_reg=arccos(-1/3)/... 
    # standard: ligands at (beta,45),(beta,225) up-pair and (pi-beta,135),(pi-beta,315)
    b=np.radians(beta)
    return [(b,np.radians(45)),(b,np.radians(225)),
            (np.pi-b,np.radians(135)),(np.pi-b,np.radians(315))]
breg=np.degrees(np.arccos(1/np.sqrt(3)))  # 54.7356: regular tetrahedron half-angle
print("MX4 vs flattening. beta = ligand polar angle; regular = %.2f deg."%breg)
print("beta < reg: FLATTENED toward xy-plane... (check signs below)  epi in units of esig\n")
for epi in [0.0,0.15,0.25,0.4]:
    print("e_pi/e_sig = %.2f"%epi)
    print(f"  {'beta':>6} {'ordering (low->high)':<46} {'{xz,yz} isolated?':>18} {'d-window':>9}")
    for beta in [45,50,54.74,60,65,70]:
        H=energies(tet(beta),epi)
        w,V=np.linalg.eigh(H)
        # group degenerate levels, label by dominant orbital
        order=np.argsort(w); labs=[]
        i=0; groups=[]
        while i<5:
            j=i
            while j+1<5 and abs(w[order[j+1]]-w[order[i]])<1e-6: j+=1
            mem=[NAMES[int(np.argmax(V[:,order[k]]**2))] for k in range(i,j+1)]
            groups.append(sorted(set(mem),key=mem.index)); i=j+1
        lab=' < '.join('{'+','.join(g)+'}' if len(g)>1 else g[0] for g in groups)
        # is {xz,yz} an isolated 2-fold group, and how many electrons sit below it?
        dbl=None; below=0
        for g in groups:
            if set(g)=={'xz','yz'}: dbl=g; break
            below+=2*len(g)
        win=("d%d-d%d"%(below+1,below+3)) if dbl else '-'
        print(f"  {beta:>6.2f} {lab:<46} {('YES' if dbl else 'no'):>18} {win:>9}")
    print()
