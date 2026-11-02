import sys; sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, P
from fractions import Fraction
import itertools, sympy
F = qcp.F
# 1. diagonal symmetry group: exponent lattice differences
sup = list(F.keys())
e0 = sup[0]
diffs = [tuple(e[i]-e0[i] for i in range(6)) for e in sup[1:]]
M = sympy.Matrix(diffs)
# Smith normal form to get the finite diagonal symmetry group (mod scalars handled by including degree row)
from sympy.matrices.normalforms import smith_normal_form
S = smith_normal_form(M, domain=sympy.ZZ)
print("Smith diag of exponent-difference lattice:", [S[i,i] for i in range(min(S.shape))])
# 2. bidegree parts
def bideg(e): return (e[0]+e[1]+e[2], e[3]+e[4]+e[5])
parts = {}
for e,c in F.items(): parts.setdefault(bideg(e), {})[e]=c
for k in sorted(parts): print(k, len(parts[k]), 'terms')
F60, F42, F24, F06 = parts[(6,0)], parts[(4,2)], parts[(2,4)], parts[(0,6)]
# 3. two-plane feasibility at degree 8: Y-even quadric a = a20(X)+a02(Y) with (a*F)_{(4,4)} = a20*F24 + a02*F42 = 0
monX2 = [e for e in core.monomials(2,0) if bideg(e)==(2,0)]
monY2 = [e for e in core.monomials(2,0) if bideg(e)==(0,2)]
cols = []
for m in monX2: cols.append(core.shift(F24, m))
for m in monY2: cols.append(core.shift(F42, m))
rows = sorted(set(e for c in cols for e in c))
A = sympy.Matrix([[c.get(e,0) for c in cols] for e in rows])
print("two-plane (4,4)-cancellation system: unknowns", len(cols), "rank", A.rank(), "-> nullity", len(cols)-A.rank())
# also (X=0 plane with m1=6, Y=0 with m2=4): need (aF) to have no bidegree (k,l) with k<=5 and l<=3: only (4,4)? (5,3) parity-odd -> same. 
# 4. Is F_{60} smooth as a plane sextic? and F_{06}?
X0,X1,X2,Y0,Y1,Y2 = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
def topoly(Pd): return sum(sympy.Rational(c.numerator,c.denominator)*X0**e[0]*X1**e[1]*X2**e[2]*Y0**e[3]*Y1**e[4]*Y2**e[5] for e,c in Pd.items())
f60 = topoly(F60); f06 = topoly(F06)
for nm,f,vs in (('F60',f60,(X0,X1,X2)),('F06',f06,(Y0,Y1,Y2))):
    G = sympy.groebner([f]+[sympy.diff(f,v) for v in vs], *vs, order='grevlex')
    print(nm, 'singular-locus groebner basis has', len(G.exprs), 'elements; is (1)?', G.exprs==[1] or any(g.is_number for g in G.exprs))
# 5. F42 as quadratic form in Y with quartic-in-X coefficients: rank along the curve? print its structure
f42 = topoly(F42)
H = sympy.Matrix(3,3, lambda i,j: sympy.diff(f42,(Y0,Y1,Y2)[i],(Y0,Y1,Y2)[j]))
print('F42 Y-Hessian det (degree 12 in X), nonzero?', sympy.expand(H.det())!=0)
print('F42 =', f42)
