import numpy as np, harvest
print("\nSPEC CHECK @ Tc=300K (energies eV)")
hdr=f"{'case':>6} {'band':>4} {'M':>7} {'lam':>6} {'W':>6} {'d_iso':>6} {'unif':>6} | {'U_req':>6} {'iso_req':>7} {'M_req':>6} | {'iso':>4} {'M':>4} {'unif':>4}"
print(hdr); print('-'*len(hdr))
for tag in ['M0.0','M0.9']:
    R,H,d=harvest.read_hr(f"cb_{tag}_hr.dat")
    Hk,dks=harvest.hk_grid(R,H,d,240,dim=2)
    for b,W,iso in harvest.auto_bands(Hk,wmax=10.0,iso_min=0.0):
        ds=harvest.descriptors(Hk,dks,b); sp=harvest.spec(ds,300.0,dim=2)
        print(f"{tag:>6} {b:>4} {ds['M']:7.3f} {ds['lam']:6.3f} {ds['W']:6.2f} {ds['d_iso']:6.2f} {ds['unif']:6.3f} | "
              f"{sp['U_req']:6.3f} {sp['iso_req']:7.3f} {sp['M_req']:6.2f} | "
              f"{'OK' if sp['pass_iso'] else 'FAIL':>4} {'OK' if sp['pass_M'] else 'FAIL':>4} {'OK' if sp['pass_unif'] else 'FAIL':>4}")
