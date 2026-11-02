"""Served-target scan: which structured eigenvector denominators carry ANY nontrivial G-invariant (3,1) class?"""
import oplib, time, sys, itertools, random, json
from hodge import core
from hodge.core import mul, pderiv, add, monomials, shift
F, A, q = oplib.load()
p = oplib.P0
Fi = [pderiv(F, i) for i in range(6)]
Fij = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
D0 = {}
for e, c in Fi[0].items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = c / 6
V = core.VARS
def mstr(e): return '*'.join(f'{V[i]}**{k}' if k > 1 else V[i] for i, k in enumerate(e) if k) or '1'
cands = []
for (i, j) in itertools.combinations(range(6), 2):
    for (k, l) in itertools.combinations_with_replacement(range(6), 2):
        Dm = core.sub(mul(Fij[i][k], Fij[j][l]), mul(Fij[i][l], Fij[j][k]))
        if Dm: cands.append((f'Hess[{i}{j}|{k}{l}]', 8, Dm))
for r in (1, 2, 3):
    for I in itertools.combinations(range(6), r):
        EI = {}
        for i in I: EI = add(EI, shift(Fi[i], tuple(1 if j == i else 0 for j in range(6))))
        for m in monomials(1) + monomials(2):
            cands.append((f'E_{"".join(V[i] for i in I)}(F)*{mstr(m)}', 6 + sum(m), shift(EI, m)))
for i in range(6):
    cands.append((f'D0*F_{V[i]}', 9, mul(D0, Fi[i])))
for (i, j) in itertools.combinations_with_replacement(range(6), 2):
    cands.append((f'F_{V[i]}*F_{V[j]}', 10, mul(Fi[i], Fi[j])))
# random colon-space elements per character at d=7,8
rng = random.Random(7)
for d in (7, 8):
    for c in [(a, b) for a in range(6) for b in range(6)]:
        dim, basis, info = oplib.colon_space(F, A, q, d, c, p)
        for t in range(min(2, dim)):
            cands.append((f'random colon element d={d} char={c} #{t}', d, oplib.rand_combo(basis, p, seed=rng.randrange(10**9))))
print('candidates', len(cands), file=sys.stderr)
seen = set()
out = open('served_scan.jsonl', 'a')
for name, d, D in cands:
    key = tuple(sorted(core.reduce_mod(D, p).items()))
    if key in seen or not D: continue
    seen.add(key)
    if oplib.char_of(core.reduce_mod(D, p)) is None: continue
    t0 = time.time()
    r = oplib.served_targets(F, q, D, p)
    dt = time.time() - t0
    rec = {'name': name, 'd': d, 'char': oplib.char_of(core.reduce_mod(D, p)), 'dim_served': r['dim_served'], 'dim_nontrivial': r['dim_nontrivial'], 'seconds': round(dt, 1)}
    out.write(json.dumps(rec) + '\n'); out.flush()
    print(rec, file=sys.stderr)
    if r['dim_nontrivial'] > 0:
        print('!!! NONTRIVIAL SERVED CLASSES:', name, file=sys.stderr)
        for b in r['basis']:
            print('   ', core.to_str(oplib.lift_poly_small(b, p) or {}), file=sys.stderr)
print('done', file=sys.stderr)
