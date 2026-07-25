import numpy as np, harvest
# H(k) = D[cos k sigma_x + sin k sigma_z].  |d|=D const -> two EXACTLY flat bands at +-D.
# theta(k) varies (d_z = D sin k) and shifts only move the AZIMUTHAL angle,
# so the metric should NOT be removable -> M_min > 0.  rho = (1 +- sin k)/2, <rho>=1/2.
SX=np.array([[0,1],[1,0]],complex); SZ=np.array([[1,0],[0,-1]],complex)
def hr(D=1.0):
    R=np.array([[0,0,0],[1,0,0],[-1,0,0]])
    H=np.zeros((3,2,2),complex)
    H[1]=(D/2)*(SX-1j*SZ); H[2]=(D/2)*(SX+1j*SZ)
    return R,H,np.ones(3)
if __name__=="__main__":
    R,H,d=hr(); harvest.write_hr('fb_hr.dat',R,H,d,'flat pair, theta-winding')
    R,H,d=harvest.read_hr('fb_hr.dat')
    Hk,dks=harvest.hk_grid(R,H,d,200,dim=1)
    E=np.linalg.eigvalsh(Hk)
    ds=harvest.descriptors(Hk,dks,0)
    mm=harvest.minimal_metric(Hk,dks,0,200,1)
    print("bands: lower [%.6f,%.6f]  upper [%.6f,%.6f]  gap=%.4f"%(
        E[...,0].min(),E[...,0].max(),E[...,1].min(),E[...,1].max(),
        float((E[...,1]-E[...,0]).min())))
    print("<rho> = %s   lam = %.6f"%(np.round(ds['w'],5),ds['lam']))
    print("M_naive = %.6f   M_MIN = %.6f   shift=%s"%(mm['M_naive'],mm['M_min'],np.round(mm['shifts'].ravel(),4)))
    print("  -> analytic tr g = 1/4 = 0.25 uniformly (kappa=1)")
