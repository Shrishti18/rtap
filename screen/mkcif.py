import numpy as np
def rot(ax,a):
    v=np.array(ax,float); v/=np.linalg.norm(v)
    K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    return np.rint(np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K).astype(int)
def close(gens):
    G=[np.eye(3,dtype=int)]
    for _ in range(10):
        new=[]
        for g in G:
            for h in gens:
                p=g@h
                if not any((p==x).all() for x in G+new): new.append(p)
        if not new: break
        G+=new
    return G
def s(W,w):
    out=[]
    for i in range(3):
        t=''
        for j,c in enumerate('xyz'):
            if W[i,j]==1: t+='+'+c
            elif W[i,j]==-1: t+='-'+c
        f=w[i]%1.0
        if abs(f)>1e-9:
            from fractions import Fraction
            fr=Fraction(f).limit_denominator(12); t+='+%d/%d'%(fr.numerator,fr.denominator)
        out.append(t.lstrip('+') if t[0]=='+' else t)
    return ','.join(out)
def cif(fn,cell,gens,cent,sites,title):
    G=close(gens); ops=[(W,np.array(c,float)) for W in G for c in cent]
    L=["data_"+title,"_cell_length_a %.5f"%cell[0],"_cell_length_b %.5f"%cell[1],"_cell_length_c %.5f"%cell[2],
       "_cell_angle_alpha %.3f"%cell[3],"_cell_angle_beta %.3f"%cell[4],"_cell_angle_gamma %.3f"%cell[5],
       "loop_","_symmetry_equiv_pos_as_xyz"]
    L+=["  '%s'"%s(W,w) for W,w in ops]
    L+=["loop_","_atom_site_label","_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z"]
    L+=["  %s %.5f %.5f %.5f"%(n,x,y,z) for n,(x,y,z) in sites]
    open(fn,'w').write('\n'.join(L)+'\n')
C4z=rot([0,0,1],np.pi/2); C3d=rot([1,1,1],2*np.pi/3); INV=-np.eye(3,dtype=int)
C2z=rot([0,0,1],np.pi); C2y=rot([0,1,0],np.pi)
# hexagonal setting: integral only in the FRACTIONAL basis
C3h=np.array([[0,-1,0],[1,-1,0],[0,0,1]]); C2h_=np.array([[1,-1,0],[0,-1,0],[0,0,-1]])
P=[(0,0,0)]; F=[(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]
cif('t_perov.cif',(3.905,3.905,3.905,90,90,90),[C4z,C3d,INV],P,
    [('Ti',(0,0,0)),('Sr',(.5,.5,.5)),('O1',(.5,0,0)),('O2',(0,.5,0)),('O3',(0,0,.5))],'SrTiO3_Pm-3m')
cif('t_dperov.cif',(8.28,8.28,8.28,90,90,90),[C4z,C3d,INV],F,
    [('Os',(0,0,0)),('Na',(.5,.5,.5)),('Ba',(.25,.25,.25)),('O',(.235,0,0))],'Ba2NaOsO6_Fm-3m')
cif('t_ortho.cif',(5.4,5.9,7.7,90,90,90),[C2z,C2y,INV],P,
    [('M',(0,0,0)),('X',(.3,.25,.1))],'orthorhombic_Pmmm')
cif('t_trig.cif',(3.4,3.4,5.8,90,90,120),[C3h,C2h_,INV],P,
    [('M',(0,0,0)),('X',(1/3.,2/3.,.25))],'trigonal_P-3m1')
print("wrote 4 test CIFs")
