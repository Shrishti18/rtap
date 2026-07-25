import numpy as np
# Where does <v^2> ~ 2t^2 come from? It assumes the gap scale and the velocity
# scale are the SAME hopping. Decouple them:
#     d(k) = D * dhat(k)  with |dhat|=1 and dhat winding n times
#  => gap  Delta = 2D          (set by D)
#     v    = D |d dhat/dk| = D n   (set by D AND n)
#  => M = <v^2>/Delta^2 = (D n)^2/(2D)^2 = n^2/4      -- INDEPENDENT OF D
# So the winding number is a free knob that raises M without touching the gap.
SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex)
def model(n,D,nk=4000):
    k=2*np.pi*(np.arange(nk)+0.5)/nk
    H=D*(np.cos(n*k)[:,None,None]*SX+np.sin(n*k)[:,None,None]*SY)
    E,V=np.linalg.eigh(H); u=V[:,:,0]
    P=np.einsum('ki,kj->kij',u,u.conj())
    dP=(np.roll(P,-1,0)-np.roll(P,1,0))/(2*(2*np.pi/nk))
    trg=0.5*np.real(np.einsum('kij,kji->k',dP,dP))
    rho=np.abs(u)**2; A=(rho.T@rho)/nk
    return float(trg.mean()), float(np.linalg.eigvalsh(A)[-1]), float((E[:,1]-E[:,0]).min()), float(E[:,0].max()-E[:,0].min())
print("d(k) = D[cos(nk) sx + sin(nk) sy]: two EXACTLY flat bands, gap 2D, winding n")
print(f"{'n':>3} {'D':>5} {'M_trg':>9} {'n^2/4':>8} {'lam':>7} {'gap':>7} {'W':>9}")
for n in [1,2,3,4]:
    for D in [1.0,5.0]:
        M,lam,gap,W=model(n,D)
        print(f"{n:>3} {D:>5.1f} {M:9.5f} {n**2/4:8.4f} {lam:7.4f} {gap:7.3f} {W:9.1e}")
print()
print("=> M = n^2/4 EXACTLY, and INDEPENDENT of D (hence of the gap).")
print("   The 'universal ceiling' Tc <= 0.152 t assumed n=1, i.e. velocity and")
print("   gap set by the same hopping. Winding breaks that tie.")
print()
coef=0.55*0.67; KB=0.086173e-3
print("Consequence.  Tc = min(U/(4 n_phi), coef*U*M), U <= Delta/2 = D, n_phi=2:")
print(f"{'n':>3} {'M=n^2/4':>8} {'amp=U/8':>9} {'stiff':>9} {'binds':>6} {'Tc/D':>7}")
for n in [1,2,3,4]:
    M=n**2/4; amp=1/8.; stf=coef*M
    print(f"{n:>3} {M:8.3f} {amp:9.4f} {stf:9.4f} {'amp' if amp<stf else 'stiff':>6} {min(amp,stf):7.4f}")
ncrit=np.sqrt(1/(2*coef))
print(f"\n   amplitude takes over at n >= sqrt(1/(2*coef)) = {ncrit:.2f}")
print(f"   for n >= 2 the ceiling is Tc = U/8 with U bounded ONLY by Delta/2 and chemistry")
print(f"   300 K needs U = {8*300*KB*1000:.0f} meV, i.e. gap Delta = {16*300*KB*1000:.0f} meV. No t-ceiling.")
