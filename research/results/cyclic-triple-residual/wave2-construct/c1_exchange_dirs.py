# C1: directions v in which F has the section-2.7 shape  F = ell^6 + ell*Lambda + f  (D_v^k F = 0 for k = 2..5)
# and the cubic shape F = ell^6 + ell^3*Gamma + Phi (checked for coordinate directions only).  Exact, sympy.
import sympy as sp, json, time
X0,X1,X2,Y0,Y1,Y2 = V = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
t = sp.Rational(5,7)
F = (X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6
     + t*(X0**3*X1*Y0*Y1 + X0**3*X2**3 + X0**2*Y0*Y1*Y2**2 + X0*X1*Y1**3*Y2 + X0*X2**3*Y1*Y2 + X1*X2**3*Y0*Y1))
v = sp.symbols('v0:6')
def Dv(P):
    return sum(v[i]*sp.diff(P, V[i]) for i in range(6))
t0=time.time()
D2 = sp.expand(Dv(Dv(F)))
eqs = sp.Poly(D2, *V).coeffs()          # quadrics in v
G = sp.groebner(eqs, *v, order='grevlex')
print('number of quadrics in v:', len(eqs))
print('Groebner basis of {D_v^2 F = 0}:', G.exprs)
# the variety: expected  v0*v1=..., all products zero and v0=v2=v4=v5=0  -> v in span(e_X1) or span(e_Y0)
# coordinate check of exponent patterns
pats = {}
for i,x in enumerate(V):
    degs = sorted(set(sp.Poly(F, x).monoms()[k][0] for k in range(len(sp.Poly(F,x).monoms()))))
    pats[str(x)] = degs
print('exponent sets of F in each coordinate:', pats)
json.dump({'quadrics': len(eqs), 'groebner': [str(g) for g in G.exprs], 'exponent_sets': pats, 'seconds': round(time.time()-t0,1)},
          open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c1_exchange_dirs.json','w'), indent=1)
