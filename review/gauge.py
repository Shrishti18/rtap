import numpy as np
# Fully dimerized chain: H_AB(k) = D e^{-ik}.  Real space: A in cell i binds B in cell i-1.
# The dimers are DECOUPLED -> a pair physically cannot move -> D_s MUST be 0.
# But the naive metric from this H(k) is nonzero. Which is right?
n=400; k=2*np.pi*np.arange(n)/n
def M_of(shift):
    # shift = position assigned to orbital B (in units of the lattice constant)
    hAB=np.exp(-1j*k)*np.exp(1j*k*shift)      # convention change: B at position `shift`
    th=np.pi/2*np.ones(n); ph=np.angle(hAB)
    u=np.stack([np.cos(th/2)*np.ones(n), np.sin(th/2)*np.exp(1j*ph)],-1)
    P=np.einsum('...i,...j->...ij',u,u.conj())
    dP=(np.roll(P,-1,0)-np.roll(P,1,0))/(2*(2*np.pi/n))
    trg=0.5*np.real(np.einsum('...ij,...ji->...',dP,dP))
    return float(trg.mean())
print("dimerized chain, metric vs assumed position of orbital B:")
for s in [0.0,0.25,0.5,0.75,1.0,1.25]:
    print(f"   B at {s:5.2f}a   <tr g> = {M_of(s):.6f}")
print("\n-> M is NOT gauge invariant. Minimum is 0 at s=1 (the physically correct")
print("   assignment: B really does sit one cell over). Decoupled dimers, D_s=0. Correct.")
print("   Wannier90 _hr.dat uses lattice-vector phases only (convention I),")
print("   so harvest.descriptors returns the NAIVE metric = an UPPER BOUND on M_min.")
print("   => pass_M can give FALSE POSITIVES. This is the Huhtinen minimal-metric issue.")
print("\nRobust to the convention (u_a -> e^{ik.d_a} u_a leaves |u_a| alone):")
print("   rho_a, lam, uniform pairing, n_phi, d_iso, W   -- all unaffected.")
print("   Only M is affected.")
