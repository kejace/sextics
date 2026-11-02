import sys; sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
X0, X1, X2, Y0, Y1, Y2 = VARS
Fs, As, q = sympy_F()
L = X0 + 2*X1 + 3*X2; M = Y0 + 3*Y1 + 5*Y2; c = sp.Rational(5, 7**7)
eX = X0**6 + X1**6 + X2**6; eY = Y0**6 + Y1**6 + Y2**6
ft = eX + 2*c*L**6; gt = eY + 2*c*M**6
f1 = sp.expand(Fs - (ft + gt + 30*c*L**2*M**2*(L**2 + M**2))) == 0
S = L**2 + M**2; P = L*M
f1b = sp.expand(c*((L+M)**6 + (L-M)**6) - 2*c*S*(S**2 + 12*P**2)) == 0
t = sp.symbols('t')
f1c = sp.expand(2*c*(1 + 15*t**2 + 15*t**4 + t**6) - 2*c*(t**2+1)*(t**4+14*t**2+1)) == 0
print('fact1 TS form:', f1, ' B = 2cS(S^2+12P^2):', f1b, ' B(1,t) factorization:', f1c)
# fact 2: cubic in w = (M/L)^2 with eps_X = e(X)/L^6, eps_Y = e(Y)/M^6
epsX, epsY, w = sp.symbols('epsX epsY w')
# F/L^6 with e(X) = epsX L^6, e(Y) = epsY M^6, M^2 = w L^2
expr = sp.expand((epsX + epsY*w**3 + 2*c*(1 + 15*w + 15*w**2 + w**3)) - ((epsX+2*c) + 30*c*w + 30*c*w**2 + (epsY+2*c)*w**3))
print('fact2 cubic:', expr == 0, '; pure cube iff epsY = 8c (and then epsX enters as constant):', sp.solve([sp.Eq(3*(epsY+2*c), 30*c)], epsY))
# fact 5: e(X)M^6 - e(Y)L^6 = (L^2+M^2) D_res mod F
Dres = eX*(L**4 - L**2*M**2 + M**4) + 2*c*L**6*(L**4 + 14*L**2*M**2 + M**4)
Dres2 = eY*(M**4 - L**2*M**2 + L**4) + 2*c*M**6*(M**4 + 14*L**2*M**2 + L**4)
lhs = sp.expand(eX*M**6 - eY*L**6 - (L**2+M**2)*Dres)
# should be a multiple of F: lhs = -(L^6)*(e(Y) + e(X) + B) ... check divisibility
qq, rr = sp.div(sp.Poly(lhs, *VARS), sp.Poly(Fs, *VARS))
print('fact5: e(X)M^6 - e(Y)L^6 - (L^2+M^2) D_res divisible by F:', rr.is_zero, ' quotient =', sp.factor(qq.as_expr()))
print('fact5b: D_res + D_res_mirror == (L^4-L^2M^2+M^4) F:', sp.expand(Dres + Dres2 - (L**4 - L**2*M**2 + M**4)*Fs) == 0)
print('deg D_res =', sp.Poly(Dres, *VARS).total_degree(), ' Y-parity even:', all((m[3]+m[4]+m[5]) % 2 == 0 for m in sp.Poly(Dres, *VARS).monoms()))
logjson(claim='construct facts 1,2,5 (TS form, cubic in w, D_res factorization) symbolic identities', verdict='CONFIRMED' if (f1 and f1b and f1c and expr == 0 and rr.is_zero) else 'REFUTED')
