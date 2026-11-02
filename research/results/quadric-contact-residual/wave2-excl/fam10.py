"""Structured candidate denominators at degrees 10-12 (sector-restricted full fixed-D solves at P2).
Families: D8* x (quadratic monomials), D8* x (natural invariant quadrics), the 13 Gamma-invariant non-Jacobian E1 directions at
degree 10 (wave-1 closure_family_d10.json), D8* x (natural quartics) at degree 12, (2GS-GP)^a (2GS+GP)^b products.
Usage: python3 fam10.py <part...>   parts: q10 inv10 nat10 q12 pow12 cub11"""
import sys, time, json
from w2lib import *
parts = sys.argv[1:] or ['q10', 'inv10', 'nat10', 'q12', 'pow12']
done = set()
try:
    for line in open(LOG):
        r = json.loads(line)
        if r.get('exp') == 'fixedD' and r.get('prime') == p: done.add(r['D'])
except FileNotFoundError:
    pass
def run(D, name, family, note=''):
    if name in done:
        print('skip (done)', name, file=sys.stderr); return
    return test_D(D, name, family, note=note)
Zn = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
if 'q10' in parts:
    for m in monomials(2):
        nm = '*'.join(f'{Zn[i]}^{m[i]}' if m[i] > 1 else Zn[i] for i in range(NV) if m[i])
        run(mul(D8, {m: Fraction(1)}), f'D8*{nm}', 'D8*quadratic-monomial', note='degree 10, E1 automatic')
if 'nat10' in parts:
    nat = [('S0', S0), ('P0', P0), ('delta', DELTA), ('S1', S1), ('P1', P1), ('S2', S2), ('P2', P2),
           ('X1^2-Y1^2', sub(mul(X1, X1), mul(Y1, Y1))), ('X2^2-Y2^2', sub(mul(X2, X2), mul(Y2, Y2))), ('X2^2+2Y2^2', add(mul(X2, X2), scale(mul(Y2, Y2), 2))),
           ('S0+S1', add(S0, S1)), ('S0+S2', add(S0, S2)), ('S1+S2', add(S1, S2)), ('S0+S1+S2', add(add(S0, S1), S2)), ('S0-S1', sub(S0, S1)),
           ('P1+P2', add(P1, P2)), ('P1-P2', sub(P1, P2)), ('P0+P1', add(P0, P1)), ('P0+P1+P2', add(add(P0, P1), P2)),
           ('X0X1+Y0Y1', add(mul(X0, X1), mul(Y0, Y1))), ('X0Y1+Y0X1', add(mul(X0, Y1), mul(Y0, X1))), ('X1X2+Y1Y2', add(mul(X1, X2), mul(Y1, Y2))),
           ('X0X2+Y0Y2', add(mul(X0, X2), mul(Y0, Y2))), ('X0^2+Y1^2', add(mul(X0, X0), mul(Y1, Y1)))]
    for nm, Qd in nat:
        run(mul(D8, Qd), f'D8*({nm})', 'D8*natural-quadric', note='degree 10')
if 'inv10' in parts:
    fam = json.load(open(LANE + '/closure_family_d10.json'))
    Ds = fam['D']; ne = fam['res']['n_extra']; nj = fam['res']['n_jac']
    extras = Ds[len(Ds) - ne:]
    for k, s in enumerate(extras):
        D = core.parse(s)
        run(D, f'inv10-extra-{k}', 'Gamma-invariant-E1-extra-d10', note='one of the 13 Gamma-invariant non-Jacobian E1 directions at degree 10 (wave-1 basis, mod 2^31-1 coefficients)')
if 'q12' in parts:
    quart = [('GS', GS), ('GP', GP), ('2GS-GP', D8A), ('2GS+GP', D8B), ('delta^2', mul(DELTA, DELTA)), ('S0^2', mul(S0, S0)), ('P0^2', mul(P0, P0)), ('S0*P0', mul(S0, P0)),
             ('Q1', Q1), ('Q2', Q2), ('S0*S1', mul(S0, S1)), ('S1^2', mul(S1, S1)), ('P1^2', mul(P1, P1)), ('S0^2+2P0^2', add(mul(S0, S0), scale(mul(P0, P0), 2))),
             ('X0^4+Y0^4', add(power(X0, 4), power(Y0, 4))), ('X1^4+Y1^4', add(power(X1, 4), power(Y1, 4))), ('S2^2', mul(S2, S2)), ('P2^2', mul(P2, P2)), ('S0*S2', mul(S0, S2))]
    for nm, Qd in quart:
        run(mul(D8, Qd), f'D8*({nm})', 'D8*natural-quartic', note='degree 12')
if 'pow12' in parts:
    for a, b in ((3, 0), (2, 1), (1, 2), (0, 3)):
        D = mul(power(D8A, a), power(D8B, b))
        run(D, f'(2GS-GP)^{a}(2GS+GP)^{b}', 'contour-powers', note='degree 12 products of the two polar quartics')
    run(mul(D8A, D8A), '(2GS-GP)^2', 'contour-powers', note='degree 8'); run(mul(D8B, D8B), '(2GS+GP)^2', 'contour-powers', note='degree 8')
    run(D8A, '2GS-GP', 'contour-powers', note='degree 4'); run(D8B, '2GS+GP', 'contour-powers', note='degree 4')
if 'cub11' in parts:
    for m in monomials(3):
        nm = '*'.join(f'{Zn[i]}^{m[i]}' if m[i] > 1 else Zn[i] for i in range(NV) if m[i])
        run(mul(D8, {m: Fraction(1)}), f'D8*{nm}', 'D8*cubic-monomial', note='degree 11, E1 automatic')
print('DONE', file=sys.stderr)
