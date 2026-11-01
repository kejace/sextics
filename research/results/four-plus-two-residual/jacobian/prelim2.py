from sympy import *
X0,X1,X2,Y0,Y1,Y2,U,V,z = symbols('X0 X1 X2 Y0 Y1 Y2 U V z')
B = X2**6 + Rational(15,7)*X2**4*Y2**2 + Rational(10,7)*X2**2*Y2**4 + Y2**6
gt = 1 + Rational(15,7)*V + Rational(10,7)*V**2 + V**3      # gtilde(V) = p(1,V)
print("B(1,t) == gt(t^2):", expand(B.subs({X2:1,Y2:V}) - gt.subs(V,V**2)) == 0)
print("disc(7 gt):", factorint(discriminant(expand(7*gt),V)))
h5 = sum(z**k for k in range(5))
cH = Poly(expand(h5**4), z).all_coeffs()[::-1]
print("Hilb(R_H):", cH, "sum", sum(cH))
RW = Poly(expand(h5**4*(1+z**2)**2), z).all_coeffs()[::-1]
R = lambda d: RW[d] if 0 <= d < len(RW) else 0
print("Hilb(R_W):", RW, "sum", sum(RW))
print("W_H (d=6, sum w=8, n=4): h^{4-p,p}_prim = R_{6(p+1)-8}:", [(4-p,p,R(6*(p+1)-8)) for p in range(5)])
RT = Poly(expand(h5**4*(1+z**2)), z).all_coeffs()[::-1]
RTd = lambda d: RT[d] if 0 <= d < len(RT) else 0
print("T_H (d=6, sum w=6, n=3): h^{3-p,p} = R_{6(p+1)-6}:", [(3-p,p,RTd(6*(p+1)-6)) for p in range(4)])
print("  split by u-power: deg6 = R_H6 + R_H4 =", cH[6], "+", cH[4])
RB = Poly(expand(h5**2*(1+z**2)), z).all_coeffs()[::-1]
RBd = lambda d: RB[d] if 0 <= d < len(RB) else 0
print("T_B (d=6, sum w=4, n=1): h^{1,0}=R_2 =", RBd(2), " h^{0,1}=R_8 =", RBd(8), " (R_2 = 3 quadratic monomials [char zeta] + u [char zeta^2])")
# C_B: w^6 = B in P^2: R_B (x) C[w]/(w^5), degree 6m-3
RC = Poly(expand(h5**2*sum(z**k for k in range(5))), z).all_coeffs()[::-1]
RCd = lambda d: RC[d] if 0 <= d < len(RC) else 0
print("C_B genus-10 check: h^{1,0}=R_3 =", RCd(3), " by w-power:", [ (k, sum(1 for i in range(5) for j in range(5) if i+j==3-k)) for k in range(4)])
# sextic fourfold X: h^{3,1} = R_6 of F
RX = Poly(expand(h5**6), z).all_coeffs()[::-1]
print("X sextic fourfold: h^{4,0}=R_0=",RX[0]," h^{3,1}=R_6=",RX[6]," h^{2,2}prim=R_12=",RX[12])
