import numpy as np, harvest
print("VALIDATION: harvester vs direct calculation (final.py / opt.py)")
print(f"{'case':>8} {'M_harv':>8} {'M_ref':>7} {'lam_harv':>9} {'lam_ref':>8} {'unif':>7} {'nphi':>5}")
ref={'M0.0':(5.068,0.5000),'M0.9':(1.907,0.6125)}
for tag in ['M0.0','M0.9']:
    R,H,d=harvest.read_hr(f"cb_{tag}_hr.dat")
    Hk,dks=harvest.hk_grid(R,H,d,240,dim=2)
    ds=harvest.descriptors(Hk,dks,band=0)
    mr,lr=ref[tag]
    print(f"{tag:>8} {ds['M']:8.4f} {mr:7.3f} {ds['lam']:9.4f} {lr:8.4f} {ds['unif']:7.4f} {ds['nphi']:5.2f}")
