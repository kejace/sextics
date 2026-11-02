"""Cone / hyperplane-union families in the reduced model (core P^3, surface S_D = V(D) cap Y_H)."""
import sys, random, itertools
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
from threefold import *
random.seed(5)
p = core.PRIMES[0]
P = core.parse
Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
f = reduce_mod(P('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6'), p)
g = reduce_mod(P('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6'), p)
# roots of f(t,1), g(s,1) mod p
def roots_binary(poly, i, j):
    from flint import nmod_poly
    coeffs = [0] * 7
    for e, c in poly.items(): coeffs[e[i]] = c % p
    pol = nmod_poly(coeffs, p)
    return [int(r) for r, m in pol.roots()] if hasattr(pol, 'roots') else []
tf = roots_binary(f, 0, 1); sg = roots_binary(g, 3, 4)
print('roots f(t,1):', tf, ' roots g(s,1):', sg, flush=True)
def lin(c):  # dict index->coeff
    return {tuple(1 if j == i else 0 for j in range(6)): c % p for i, c in c.items() if c % p}
def pw(Q, k): return core.power(Q, k, p)
def addp(*Qs):
    out = {}
    for Q in Qs: out = add(out, Q, p)
    return out
def sc(Q, c): return scale(Q, c % p, p)
X0, X1, Y0, Y1 = [reduce_mod(P(v), p) for v in ['X0', 'X1', 'Y0', 'Y1']]
# (a) is a*L^6 in J(H)?  (E1 for D = L^6)  via the kernel: check rank
Hi = [pderiv(Hp, i, p) for i in CORE]
def e1_core(D):
    Jcols = []
    k = core.degree(D)
    for idx, i in enumerate(CORE):
        for m2 in cmon(k + 1):
            Jcols.append(shift(Hi[idx], m2, 1, p))
    M, rows = build(Jcols, p); r1 = M.rank()
    M2, _ = build(Jcols + [mul(ap, D, p)], p, dict(rows)); r2 = M2.rank()
    return r1 == r2
cands = []
for t in tf: cands.append((lin({0: 1, 1: -t}), f'L=X0-{t}X1 (f(t,1)=0)'))
for s in sg: cands.append((lin({3: 1, 4: -s}), f'L=Y0-{s}Y1 (g(s,1)=0)'))
cands += [(X0, 'X0'), (Y0, 'Y0'), (lin({0: 1, 3: 1}), 'X0+Y0'), (lin({0: 1, 1: 2, 3: 3, 4: 5}), 'X0+2X1+3Y0+5Y1')]
for L, name in cands:
    ok = e1_core(pw(L, 6))
    print(f'E1 for D=L^6, {name}: {ok}', flush=True)
    log({'D': f'L^6, {name}', 'd': 6, 'family': 'threefold-L6-E1', 'e1': ok, 'consistent': None, 'rank': None, 'nullity': None, 'seconds': 0, 'note': 'is a L^6 in J(H)? (reduced model E1)', 'engine': 'threefold-v1', 'p': p})
# (b) cones with vertex on S_H
tests = []
if tf:
    t = tf[0]; l = lin({0: 1, 1: -t})
    tests += [(addp(pw(l, 6), pw(Y0, 6), pw(Y1, 6)), f'cone l^6+Y0^6+Y1^6, l=X0-{t}X1, vertex (t:1:0:0) on S_H'),
              (addp(pw(l, 6), g), 'cone l^6+g(Y0,Y1), vertex on S_H'), (addp(pw(l, 6), sc(g, -1)), 'cone l^6-g, vertex on S_H'),
              (addp(pw(l, 4), pw(Y0, 4), pw(Y1, 4)), 'cone quartic l^4+Y0^4+Y1^4, vertex on S_H'),
              (addp(pw(l, 3), pw(Y0, 3), pw(Y1, 3)), 'cone cubic l^3+Y0^3+Y1^3, vertex on S_H'),
              (addp(pw(l, 2), pw(Y0, 2), pw(Y1, 2)), 'quadric cone l^2+Y0^2+Y1^2, vertex on S_H'),
              (addp(mul(pw(l, 2), mul(Y0, Y1, p), p), pw(Y0, 4), pw(Y1, 4)), 'cone quartic l^2 Y0Y1 + Y0^4+Y1^4, vertex on S_H')]
    for k in range(2):
        D = {}
        for a_ in range(7):
            for b in range(7 - a_):
                c = 6 - a_ - b
                if (b + c) % 2: continue
                D = add(D, sc(mul(mul(pw(l, a_), pw(Y0, b), p), pw(Y1, c), p), random.randint(1, 50)), p)
        tests.append((D, f'random even sextic cone in (l,Y0,Y1), vertex on S_H #{k}'))
    lp = lin({0: 1, 1: -(t + 1)})
    tests.append((addp(pw(lp, 6), pw(Y0, 6), pw(Y1, 6)), 'control: cone l\'^6+Y0^6+Y1^6, vertex off S_H'))
if sg:
    s = sg[0]; m = lin({3: 1, 4: -s})
    tests += [(addp(pw(X0, 6), pw(X1, 6), pw(m, 6)), f'cone X0^6+X1^6+m^6, m=Y0-{s}Y1, vertex (0:0:s:1) on S_H'),
              (addp(f, pw(m, 6)), 'cone f(X0,X1)+m^6, vertex on S_H'), (addp(f, sc(pw(m, 6), -1)), 'cone f-m^6, vertex on S_H'),
              (addp(pw(X0, 4), pw(X1, 4), pw(m, 4)), 'cone quartic X0^4+X1^4+m^4, vertex on S_H'),
              (addp(pw(X0, 3), pw(X1, 3), pw(m, 3)), 'cone cubic X0^3+X1^3+m^3 (Y-odd), vertex on S_H')]
# (c) cones with vertex at a general point of S_H (mod p): find v with H(v)=0, all coords nonzero
v = None
for _ in range(2000):
    cand = [random.randrange(1, p) for _ in range(3)]
    # solve for Y1 coordinate: H(x0,x1,y0,Y1)=0 univariate
    from flint import nmod_poly
    coeffs = [0] * 7
    for e, c in Hp.items():
        coeffs[e[4]] = (coeffs[e[4]] + c * pow(cand[0], e[0], p) * pow(cand[1], e[1], p) * pow(cand[2], e[3], p)) % p
    rts = nmod_poly(coeffs, p).roots()
    if rts:
        v = (cand[0], cand[1], 0, cand[2], int(rts[0][0]), 0); break
if v:
    # three linear forms vanishing at v: X0 - (v0/v4) Y1, X1 - (v1/v4) Y1, Y0 - (v3/v4) Y1
    iv4 = pow(v[4], -1, p)
    L1 = lin({0: 1, 4: -v[0] * iv4}); L2 = lin({1: 1, 4: -v[1] * iv4}); L3 = lin({3: 1, 4: -v[3] * iv4})
    assert core.evaluate(Hp, v, p) == 0
    tests += [(addp(pw(L1, 6), pw(L2, 6), pw(L3, 6)), 'Fermat sextic cone, vertex general point of S_H'),
              (addp(pw(L1, 4), pw(L2, 4), pw(L3, 4)), 'Fermat quartic cone, vertex general point of S_H'),
              (addp(pw(L1, 3), pw(L2, 3), pw(L3, 3)), 'Fermat cubic cone, vertex general point of S_H')]
# (d) six hyperplanes through a plane sextic: D = H + lam L^6, and binary forms in two linear forms
for L, name in cands[:2] + cands[6:8] + cands[12:14]:
    for lam in [1, -1, 7]:
        tests.append((addp(Hp, sc(pw(L, 6), lam)), f'H + ({lam}) L^6, {name}'))
tests += [(addp(pw(X0, 6), pw(Y0, 6)), 'X0^6+Y0^6 (6 planes through X0=Y0=0)'),
          (addp(pw(X0, 6), sc(pw(Y0, 6), -1)), 'X0^6-Y0^6'),
          (f if False else addp(pw(lin({0: 1, 3: 1}), 6), pw(lin({0: 1, 3: -1}), 6)), '(X0+Y0)^6+(X0-Y0)^6'),
          (addp(pw(X0, 2), pw(Y0, 2)), 'X0^2+Y0^2'), (mul(f, pw(Y0, 0), p), 'f(X0,X1) again (6 planes through X0=X1=0)')]
for D, note in tests:
    run3(D, note, 'threefold-cones', p)
