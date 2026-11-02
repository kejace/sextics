import sys, time, random
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
F, A, q = load(); p = P0
ctx = Ctx(F, A, p)
random.seed(77)
def rnd(): return random.randint(1, p - 1)
def combo(polys, coeffs):
    acc = {}
    for P, c in zip(polys, coeffs):
        if c: acc = add(acc, scale(P, c, p), p)
    return acc
Fi = [core.pderiv(F, i) for i in range(6)]
X = [core.parse(v) for v in core.VARS]
# Y-odd polars
for i in range(3, 6):
    e = full_test(F, A, Fi[i], f'F_{i}', 'polar', note='Y-odd partial derivative', verbose=True)
    print(f'F_{i}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]}', flush=True)
# Z_j F_i, all 36 (both parities)
for i in range(6):
    for j in range(6):
        D = core.mul(X[j], Fi[i])
        e = full_test(F, A, D, f'{core.VARS[j]}*F_{i}', 'polar-times-coordinate')
        print(f'{core.VARS[j]}*F_{i}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
# random elements of J_d per character and parity, d = 6,7,8
for d in [6, 7, 8]:
    for parD in [0, 1]:
        for k in range(7):
            jk, jp = jac_block_cols(ctx, d, parD, k)
            if not jp: continue
            for trial in range(3):
                D = combo(jp, [rnd() for _ in jp])
                e = full_test(F, A, D, f'random J_{d} parD={parD} k={k} trial {trial}', 'J_d-random',
                              note=f'random F_p point of J(F)_d block ({len(jp)} generators)')
                print(f'J_{d} parD={parD} k={k} trial {trial}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} rows={e.get("rows")} cols={e.get("cols")} {e["seconds"]}s', flush=True)
# products of two polars (degree 10)
for i in range(6):
    for j in range(i, 6):
        D = core.mul(Fi[i], Fi[j])
        if core.parity_of(D) is None: continue
        e = full_test(F, A, D, f'F_{i}*F_{j}', 'polar-product')
        print(f'F_{i}*F_{j}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
