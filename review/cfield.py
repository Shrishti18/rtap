import numpy as np
# Point-charge / sigma-antibonding ordering: orbital i is destabilised in
# proportion to sum_L |Y_i(Omega_L)|^2 over ligand directions. Completeness
# gives sum_i |Y_i|^2 = 1 per ligand, so the trace is fixed and only the
# SPLITTING pattern matters -- which is all we need.
def Ys(v):
    x,y,z=v/np.linalg.norm(v)
    s3=np.sqrt(3)
    return np.array([(3*z*z-1)/2, (s3/2)*(x*x-y*y), s3*x*y, s3*x*z, s3*y*z])
NAMES=['z2','x2-y2','xy','xz','yz']
def levels(ligands,weights=None):
    w=np.ones(len(ligands)) if weights is None else np.asarray(weights,float)
    E=np.zeros(5)
    for L,wl in zip(ligands,w): E+=wl*Ys(np.asarray(L,float))**2
    return E
GEOM={
 'linear (2, +-z)':            ([[0,0,1],[0,0,-1]],None),
 'octahedral':                 ([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],None),
 'oct. ELONGATED along z':     ([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],[1,1,1,1,.5,.5]),
 'oct. COMPRESSED along z':    ([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],[1,1,1,1,1.6,1.6]),
 'square planar (xy)':         ([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]],None),
 'tetrahedral':                ([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],None),
 'tet. FLATTENED (z-squash)':  ([[1,1,.6],[1,-1,-.6],[-1,1,-.6],[-1,-1,.6]],None),
 'tet. ELONGATED (z-stretch)': ([[1,1,1.7],[1,-1,-1.7],[-1,1,-1.7],[-1,-1,1.7]],None),
}
print("Which orbitals form the 2-fold level, and what d-count puts it PARTIAL?\n")
print(f"{'geometry':>28}  {'ordering (low -> high)':<44} {'doublet':>9} {'d-window':>9}")
for g,(lig,w) in GEOM.items():
    E=levels(lig,w); o=np.argsort(E)
    lab=[]; i=0
    groups=[]
    while i<5:
        j=i
        while j+1<5 and abs(E[o[j+1]]-E[o[i]])<1e-6: j+=1
        groups.append([NAMES[o[k]] for k in range(i,j+1)]); i=j+1
    lab=' < '.join('{'+','.join(g2)+'}' if len(g2)>1 else g2[0] for g2 in groups)
    # find the first 2-fold group and the electrons below it
    dbl=None; below=0
    for g2 in groups:
        if len(g2)==2 and set(g2)=={'xz','yz'}: dbl=g2; break
        below+=2*len(g2)
    if dbl is None:
        dbl2=[g2 for g2 in groups if len(g2)==2]
        if dbl2:
            below=0
            for g2 in groups:
                if len(g2)==2: dbl=g2; break
                below+=2*len(g2)
    win=f"d{below+1}-d{below+3}" if dbl else "none"
    print(f"{g:>28}  {lab:<44} {str(dbl) if dbl else '-':>9} {win:>9}")
print()
print("=> The window is NOT universal. My d5-d7 assumed LINEAR coordination.")
print("   Real chain/ladder pnictides are TETRAHEDRAL-derived, which moves it.")
