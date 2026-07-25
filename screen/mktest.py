import numpy as np, harvest
# checkerboard 2-band model as real-space hoppings:
# H = t(cos kx+cos ky) sx + t(cos ky-cos kx) sy + [M+2tp(cos kx-cos ky)] sz
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex)
def build(t=1.0,tp=0.5,M=0.0):
    R=np.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]])
    Hx=(t/2)*sx-(t/2)*sy+tp*sz     # coeff of cos kx  -> split 1/2 to +-x
    Hy=(t/2)*sx+(t/2)*sy-tp*sz     # coeff of cos ky
    H=np.array([M*sz,Hx,Hx,Hy,Hy])
    return R,H,np.ones(5)
for tag,M in [("M0.0",0.0),("M0.9",0.9)]:
    R,H,d=build(M=M); harvest.write_hr(f"cb_{tag}_hr.dat",R,H,d,f"checkerboard tp=0.5 M={M}")
print("wrote cb_M0.0_hr.dat cb_M0.9_hr.dat")
