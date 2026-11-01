"""Experiment 0: sanity checks, symmetry, D=1, D8* and its factors."""
import sys, time
from common import *
import sympy as sp

print('F terms', len(F), 'A', core.to_str(A), 'q', q)
print('chi(F) sectors:', {chi(e) for e in F})
print('chi(A):', CHI_A)
# reflections: swap X0<->Y0, swap X1<->Y1 preserve F?
def swap(P, i, j):
    out = {}
    for e, c in P.items():
        e2 = list(e); e2[i], e2[j] = e2[j], e2[i]
        out[tuple(e2)] = c
    return out
print('F swap01 invariant:', swap(F, 0, 3) == F, ' F swap14 invariant:', swap(F, 1, 4) == F, ' F swap25:', swap(F, 2, 5) == F)
print('A swap01 invariant:', swap(A, 0, 3) == A)

# structure: F = S^3 + 2 P^2 S + 5 S Q1 + 5/7 P Q2 + R
S = add(mul(X0, X0), mul(Y0, Y0)); Pm = mul(X0, Y0)
Q1 = add(add(mul(mul(X1, X1), mul(Y1, Y1)), mul(mul(X1, X2), mul(Y1, Y2))), mul(mul(X2, X2), mul(Y2, Y2)))
Q2 = mul(mul(X1, Y1), add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
R = add(add(power(X1, 6), power(X2, 6)), add(power(Y1, 6), power(Y2, 6)))
Fcheck = add(add(power(S, 3), scale(mul(mul(Pm, Pm), S), 2)), add(scale(mul(S, Q1), 5), add(scale(mul(Pm, Q2), Fraction(5, 7)), R)))
print('F == S^3 + 2P^2S + 5SQ1 + 5/7 PQ2 + R :', Fcheck == F)

# G_S, G_P, D8* = 4 G_S^2 - G_P^2 ; check F_x = 2x G_S + y G_P
GS = add(add(scale(mul(S, S), 3), scale(mul(Pm, Pm), 2)), scale(Q1, 5))
GP = add(scale(mul(Pm, S), 4), scale(Q2, Fraction(5, 7)))
Fx = pderiv(F, 0); Fy = pderiv(F, 3)
print('F_x == 2x G_S + y G_P:', Fx == add(scale(mul(X0, GS), 2), mul(Y0, GP)))
print('F_y == 2y G_S + x G_P:', Fy == add(scale(mul(Y0, GS), 2), mul(X0, GP)))
D8 = sub(scale(mul(GS, GS), 4), mul(GP, GP))
delta = sub(mul(X0, X0), mul(Y0, Y0))
print('(F_x^2 - F_y^2) == delta * D8*:', sub(mul(Fx, Fx), mul(Fy, Fy)) == mul(delta, D8))
# factors: D8a = (F_x - F_y)/(x - y), D8b = (F_x + F_y)/(x + y)
D8a = exact_div(sub(Fx, Fy), sub(X0, Y0))
D8b = exact_div(add(Fx, Fy), add(X0, Y0))
print('D8* == D8a * D8b:', D8 == mul(D8a, D8b), ' deg', core.degree(D8a), core.degree(D8b))
print('chi(D8*)=', chi_of(D8), 'chi(D8a)=', chi_of(D8a), 'chi(D8b)=', chi_of(D8b))
# 2 G_S - G_P and 2 G_S + G_P vs D8a, D8b
print('D8a == 2G_S - G_P:', D8a == sub(scale(GS, 2), GP), ' D8b == 2G_S + G_P:', D8b == add(scale(GS, 2), GP))

# E1 identity of fqlx: 2 x y D8* = (2y G_S - x G_P) F_x + (2x G_S - y G_P) F_y
lhs = scale(mul(mul(X0, Y0), D8), 2)
rhs = add(mul(sub(scale(mul(Y0, GS), 2), mul(X0, GP)), Fx), mul(sub(scale(mul(X0, GS), 2), mul(Y0, GP)), Fy))
print('fqlx E1 identity:', lhs == rhs)

# ---- tests ----
t0 = time.time()
r = test_D(ONE, 'constant', note='D=1: is A in J(F)?')
print('D=1 ->', {k: v for k, v in r.items() if k not in ('system', 'res')})
for name, D in [('D8*', D8), ('D8a', D8a), ('D8b', D8b), ('G_S', GS), ('G_P', GP)]:
    r = test_D(D, 'fqlx-degree8', note=name)
    print(name, '->', {k: v for k, v in r.items() if k not in ('system', 'res')})
# cross-check the generic engine on D8* (should also be inconsistent)
t0 = time.time()
out = fixed_d.solve_fixed_D_fast(F, A, q, D8, verbose=True, want_exact=False)
print('generic engine on D8*:', out, f'({time.time()-t0:.1f}s)')
log(core.to_str(D8), 8, 'fqlx-degree8', True, out['consistent'], out['rank'], out['nullity'], time.time() - t0, 'D8* generic engine cross-check', rows=out['rows'], cols=out['cols'])
