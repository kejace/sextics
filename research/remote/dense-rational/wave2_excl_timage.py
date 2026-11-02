#!/usr/bin/env python3
"""Served-target spaces T(D) = {A' in S_6^- : some (N, Q) certifies A' with denominator D} (mod p) for named
degree-5/6 denominators, via fixed_d.target_image.  dim T(D) = 18 means D serves only J_6^- (nothing primitive).
Reading rule: dim T_p(D) is modular evidence; it equals dim T_Q(D) unless p divides the index of the projected
solution lattice, which we do not certify here (we report two primes when available).  For the sum principle
served(D1 D2) = served(D1) + served(D2) (mod J) the useful strong form is "T(D) = J_6^-".
Usage: python3 wave2_excl_timage.py prime workers [names...]"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
from multiprocessing import Pool
import random
p = int(sys.argv[1]); workers = int(sys.argv[2]); only = sys.argv[3:]
F, A, q = load()
Fi = [pderiv(F, i) for i in range(NV)]
Z = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
def blocks(F):
    out = {0: {}, 2: {}, 4: {}, 6: {}}
    for e, c in F.items(): out[e[3] + e[4] + e[5]][e] = c
    return out
B = blocks(F)
def EX(P_): return add(add(mul(core.parse('X0'), pderiv(P_, 0)), mul(core.parse('X1'), pderiv(P_, 1))), mul(core.parse('X2'), pderiv(P_, 2)))
def EY(P_): return add(add(mul(core.parse('Y0'), pderiv(P_, 3)), mul(core.parse('Y1'), pderiv(P_, 4))), mul(core.parse('Y2'), pderiv(P_, 5)))
cands = {}
for i in range(3): cands[f'F_{Z[i]}'] = Fi[i]
cands['F_X0+F_X1+F_X2'] = add(add(Fi[0], Fi[1]), Fi[2])
cands['F_X0+2F_X1+3F_X2'] = add(add(Fi[0], scale(Fi[1], 2)), scale(Fi[2], 3))
cands['E_X(F)'] = EX(F); cands['E_Y(F)'] = EY(F); cands['E_X(F)-E_Y(F)'] = sub(EX(F), EY(F))
cands['X0F_X0+Y0F_Y0'] = add(mul(core.parse('X0'), Fi[0]), mul(core.parse('Y0'), Fi[3]))
cands['X0F_X0'] = mul(core.parse('X0'), Fi[0]); cands['Y0F_Y0'] = mul(core.parse('Y0'), Fi[3])
cands['X1F_X0'] = mul(core.parse('X1'), Fi[0]); cands['Y1F_Y0'] = mul(core.parse('Y1'), Fi[3])
cands['a=F60'] = B[0]; cands['d=F06'] = B[6]; cands['b=F42'] = B[2]; cands['c=F24'] = B[4]
cands['a+d'] = add(B[0], B[6]); cands['a-d'] = sub(B[0], B[6]); cands['b+c'] = add(B[2], B[4]); cands['b+2c+3d'] = add(add(B[2], scale(B[4], 2)), scale(B[6], 3))
cands['fermat-X0X1X2'] = core.parse('X0**6+X1**6+X2**6'); cands['fermat-X0X1Y1'] = core.parse('X0**6+X1**6+Y1**6'); cands['fermat-Y0Y1Y2'] = core.parse('Y0**6+Y1**6+Y2**6')
cands['X0^6+Y0^6'] = core.parse('X0**6+Y0**6'); cands['sum-sixth-powers'] = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
random.seed(7)
for t in range(2):
    cands[f'random-even-sextic-{t}'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(6, 0)}
random.seed(8)
cands['random-even-quintic'] = {m: Fraction(random.randrange(-9, 10)) for m in monomials(5, 0)}
def one(nm):
    D = cands[nm]
    t0 = time.time()
    try:
        dim, basis, S = fixed_d.target_image(F, q, D, p=p)
        rec = {'exp': 'target-image', 'D': nm, 'd': core.degree(D), 'prime': p, 'target_image_dim': dim,
               'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 1),
               'result': 'serves only J_6^- (nothing primitive) mod p' if dim == 18 else f'serves a {dim-18}-dim space of primitive targets mod p'}
    except Exception as ex:
        rec = {'exp': 'target-image', 'D': nm, 'd': core.degree(D), 'prime': p, 'error': repr(ex)[:300], 'seconds': round(time.time() - t0, 1)}
    log2(**rec); print(json.dumps(rec), flush=True)
    return rec
if __name__ == '__main__':
    names = only if only else list(cands)
    with Pool(workers) as pool:
        res = pool.map(one, names)
    with open(f'{OUT}/timage_p{p}.json', 'w') as f: json.dump(res, f, indent=1)
