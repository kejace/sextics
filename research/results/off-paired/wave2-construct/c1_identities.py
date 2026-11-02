# Exact sanity checks of the structural identities of part 1 (sympy, light).
from sympy import *
X0,X1,X2,Y0,Y1,Y2 = symbols('X0 X1 X2 Y0 Y1 Y2')
F = X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6 + X0**3*X1*Y0*Y1/2 + X1**3*X2*Y1*Y2/2 + X0**2*X1*X2*Y0*Y2/3
H = X0**6+X1**6+Y0**6+Y1**6 + X0**3*X1*Y0*Y1/2
C = X1**3*Y1/2 + X0**2*X1*Y0/3
m = X2*Y2
out = {}
out['block_form'] = expand(F - (H + C*m + X2**6 + Y2**6)) == 0
# t1 = i(X2^3+Y2^3), t2 = i(X2^3-Y2^3): t1^2 = H + C m - 2 m^3 and t2^2 = H + C m + 2 m^3 modulo F
t1sq = -(X2**3+Y2**3)**2; t2sq = -(X2**3-Y2**3)**2
out['t1_identity_modF'] = simplify(expand(t1sq - (H + C*m - 2*m**3)) + F) == 0  # difference should be -F ... check
out['t1_diff'] = factor(expand(t1sq - (H + C*m - 2*m**3)))
out['t2_diff'] = factor(expand(t2sq - (H + C*m + 2*m**3)))
# v = X2^6 - Y2^6: v^2 = (H + C m)^2 - 4 m^6 mod F
v = X2**6 - Y2**6
out['v_diff'] = factor(expand(v**2 - ((H+C*m)**2 - 4*m**6)))
# genus-2 branch factorisation
u = symbols('u')
out['branch_factor'] = factor((H + C*u)**2 - 4*u**6)
# Weierstrass discriminant of t^2 = 2 m^3 - C m - H  (cubic in m)
disc = discriminant(2*u**3 - C*u - H, u)
out['weierstrass_disc'] = factor(disc)
# Sigma rationality: lambda = (1+w)(1-s^2)/(2 w^3 (1+s^2)), v' = 2 s (1+w)/(1+s^2) satisfy v'^2 = (1+w)^2 - 4 lambda^2 w^6
w,s = symbols('w s')
lam = (1+w)*(1-s**2)/(2*w**3*(1+s**2)); vp = 2*s*(1+w)/(1+s**2)
out['sigma_param'] = simplify(vp**2 - ((1+w)**2 - 4*lam**2*w**6)) == 0
# g-action: (lambda,w,s) -> (-lambda, w, 1/s) preserves the curve (1+w)(1-s^2) = 2 lambda w^3 (1+s^2)
Lm = symbols('L')
E = (1+w)*(1-s**2) - 2*Lm*w**3*(1+s**2)
Eg = E.subs({Lm:-Lm, s:1/s}, simultaneous=True)
out['g_preserves_E'] = simplify(Eg*s**2 + E) == 0 or simplify(Eg*s**2 - E) == 0, factor(Eg*s**2)
for k,v_ in out.items(): print(k, '=', v_)
