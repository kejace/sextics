#!/usr/bin/env python3
"""Rigorous 'serves no primitive target' test by E1 alone: for a denominator D (degree d, pure parity),
   T(D) := {A' in S_6^- : A' D in J_{d+6}}  (every served target satisfies E1, so served(D) ⊆ T(D)).
J_6^- (dim 18) is always contained.  If dim T(D) mod p = 18 then T(D) = J_6^- over Q (saturated-lattice argument:
the rational solution lattice reduces injectively, and dim J_6^- = 18 is the exact Koszul count), hence D serves no
primitive target over C — the strong form of exclusion needed for the sum principle served(D1 D2) = served(D1)+served(D2).
Usage: python3 wave2_excl_e1dual.py prime workers [names...]"""
import sys, time, json, random, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
from multiprocessing import Pool
p = int(sys.argv[1]); workers = int(sys.argv[2]); only = sys.argv[3:]
F, A, q = load()
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
Z = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
def blocks(F):
    out = {0: {}, 2: {}, 4: {}, 6: {}}
    for e, c in F.items(): out[e[3] + e[4] + e[5]][e] = c
    return out
B = blocks(F)
X = lambda s: core.parse(s)
def EX(P_): return add(add(mul(X('X0'), pderiv(P_, 0)), mul(X('X1'), pderiv(P_, 1))), mul(X('X2'), pderiv(P_, 2)))
def EY(P_): return add(add(mul(X('Y0'), pderiv(P_, 3)), mul(X('Y1'), pderiv(P_, 4))), mul(X('Y2'), pderiv(P_, 5)))
cands = {}
for i in range(NV): cands[f'F_{Z[i]}'] = Fi[i]
cands['F_X0+F_X1+F_X2'] = add(add(Fi[0], Fi[1]), Fi[2]); cands['F_X0+2F_X1+3F_X2'] = add(add(Fi[0], scale(Fi[1], 2)), scale(Fi[2], 3))
cands['F_Y0+2F_Y1+3F_Y2'] = add(add(Fi[3], scale(Fi[4], 2)), scale(Fi[5], 3))
cands['E_X(F)'] = EX(F); cands['E_Y(F)'] = EY(F); cands['E_X(F)-E_Y(F)'] = sub(EX(F), EY(F))
cands['X0F_X0+Y0F_Y0'] = add(mul(X('X0'), Fi[0]), mul(X('Y0'), Fi[3])); cands['X0F_X0'] = mul(X('X0'), Fi[0]); cands['Y0F_Y0'] = mul(X('Y0'), Fi[3])
cands['X1F_X0'] = mul(X('X1'), Fi[0]); cands['Y1F_Y0'] = mul(X('Y1'), Fi[3]); cands['Y0F_X0'] = mul(X('Y0'), Fi[0]); cands['X0F_Y0'] = mul(X('X0'), Fi[3])
cands['a=F60'] = B[0]; cands['d=F06'] = B[6]; cands['b=F42'] = B[2]; cands['c=F24'] = B[4]
cands['a+d'] = add(B[0], B[6]); cands['a-d'] = sub(B[0], B[6]); cands['b+c'] = add(B[2], B[4]); cands['b+2c+3d'] = add(add(B[2], scale(B[4], 2)), scale(B[6], 3))
cands['3a+2b+c'] = add(add(scale(B[0], 3), scale(B[2], 2)), B[4])
for I in itertools.combinations(range(6), 3):
    cands['fermat-' + ''.join(Z[i] for i in I)] = {tuple(6 if k == i else 0 for k in range(6)): Fraction(1) for i in I}
cands['X0^6+Y0^6'] = X('X0**6+Y0**6'); cands['sum-sixth-powers'] = X('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
cands['(X0^2+X1^2+X2^2)^3'] = power(X('X0**2+X1**2+X2**2'), 3); cands['(X0^2+Y0^2)^3'] = power(X('X0**2+Y0**2'), 3)
random.seed(7)
for t in range(2): cands[f'random-even-sextic-{t}'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(6, 0)}
cands['random-odd-sextic'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(6, 1)}
cands['random-even-quintic'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(5, 0)}
cands['random-even-quartic'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(4, 0)}
cands['random-even-cubic'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(3, 0)}
# degree 7 / 8 structured J-elements
cands['X0^2F_X0'] = mul(X('X0**2'), Fi[0]); cands['X0E_X(F)'] = mul(X('X0'), EX(F)); cands['Y0^2F_X0'] = mul(X('Y0**2'), Fi[0]); cands['X0Y0F_Y0'] = mul(X('X0*Y0'), Fi[3])
cands['X0^3F_X0'] = mul(X('X0**3'), Fi[0]); cands['X0^3F_X0+Y0^3F_Y0'] = add(mul(X('X0**3'), Fi[0]), mul(X('Y0**3'), Fi[3]))
cands['Hess-minor-X0X0.Y0Y0'] = sub(mul(Fij[0][0], Fij[3][3]), mul(Fij[0][3], Fij[0][3])); cands['Hess-minor-X0X0.X1X1'] = sub(mul(Fij[0][0], Fij[1][1]), mul(Fij[0][1], Fij[0][1]))
cands['E_X(F)(X0^2+Y0^2)'] = mul(EX(F), X('X0**2+Y0**2')); cands['a(X)(Y0^2+Y1^2+Y2^2)'] = mul(B[0], X('Y0**2+Y1**2+Y2**2'))
def one(nm):
    D = cands[nm]; d = core.degree(D); par = core.parity_of(D)
    t0 = time.time()
    Dp = core.reduce_mod(D, p); Fp = core.reduce_mod(F, p); Fip = [pderiv(Fp, i, p) for i in range(NV)]
    tpar = (par + 1) % 2                       # A' is Y-odd; A' D has parity par+1
    Amons = monomials(6, 1)
    cols = [mul(Dp, {m: 1}, p) for m in Amons]
    nA = len(cols)
    for i in range(NV):
        for m in monomials(d + 1, (tpar - par_i(i)) % 2): cols.append(shift(Fip[i], m, 1, p))
    rows = {e: j for j, e in enumerate(monomials(d + 6, tpar))}
    M = nmod_mat(len(rows), len(cols), p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = v
    Xn, nullity = M.nullspace()
    proj = nmod_mat(nullity, nA, [int(Xn[j, k]) for k in range(nullity) for j in range(nA)], p)
    dimT = proj.rank()
    rec = {'exp': 'e1-dual-served-targets', 'D': nm, 'd': d, 'parD': par, 'prime': p, 'rows': len(rows), 'cols': len(cols),
           'dim_T_E1': dimT, 'dim_T_mod_J6': dimT - 18, 'seconds': round(time.time() - t0, 1),
           'result': 'serves NO primitive target (rigorous over C: E1 alone)' if dimT == 18 else f'E1 allows a {dimT-18}-dim space of primitive targets (E2 not tested here)'}
    log2(**rec); print(json.dumps(rec), flush=True)
    return rec
if __name__ == '__main__':
    names = only if only else [n for n in cands if core.parity_of(cands[n]) is not None]
    with Pool(workers) as pool: res = pool.map(one, names)
    with open(f'{OUT}/e1dual_p{p}.json', 'w') as f: json.dump(res, f, indent=1)
