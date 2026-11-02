"""Equivariant Hilbert series of S/J under the Klein group for a smooth Klein-symmetric sextic (independent)."""
import sys, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
import sympy as sp
t = sp.symbols('t')
# group elements as sign vectors on (X0,X1,X2,Y0,Y1,Y2)
G = {'1': (1,1,1,1,1,1), 'r': (1,1,1,1,1,-1), 'rho': (1,1,1,-1,-1,1), 'tau': (1,1,1,-1,-1,-1)}
def trace_series(s):
    expr = 1
    for si in s:
        expr *= (1 - si*t**5) / (1 - si*t)
    return sp.series(expr, t, 0, 25).removeO()
tr = {g: sp.Poly(trace_series(s), t) for g, s in G.items()}
# characters of the Klein group indexed by numerator character c=(a,b): a = Y2-parity, b = (Y0+Y1)-parity
# a monomial with character (a,b) has eigenvalue (-1)^a under r, (-1)^b under rho, (-1)^(a+b) under tau
def dim_sector(d, c):
    a, b = c
    vals = {'1': 1, 'r': (-1)**a, 'rho': (-1)**b, 'tau': (-1)**(a+b)}
    return sp.Rational(sum(vals[g] * tr[g].coeff_monomial(t**d) for g in G), 4)
rows = {}
for d in (6, 8, 9, 10, 12, 14, 15, 16, 18):
    rows[d] = {str(c): int(dim_sector(d, c)) for c in [(0,0),(0,1),(1,0),(1,1)]}
    tot = sum(rows[d].values())
    print(d, rows[d], 'total', tot, ' Fermat total', fermat_count(d),
          ' Fermat by char', {str(c): fermat_count(d, c) for c in [(0,0),(0,1),(1,0),(1,1)]})
# monomial counts
print('Y-odd Y2-even sextic monomials:', len(monos_chi(6,(0,1))), ' Y-even Y2-even sextic monomials:', len(monos_chi(6,(0,0))))
# J_d^chi = S_d^chi - (S/J)_d^chi
for d in (8,):
    for c in [(0,0),(0,1),(1,0),(1,1)]:
        print(f'J_{d}^{c} = {len(monos_chi(d,c)) - int(dim_sector(d,c))}')
logjson(VLOG, {'check': 'sector-dims', 'dims': rows})
