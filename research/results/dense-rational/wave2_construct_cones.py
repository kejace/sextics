#!/usr/bin/env python3
"""E1 test for ternary cone denominators D = g(X) (sector X) or D = g(Y) (sector Y), g of degree k:
E1 <=> A*g in J(F)_{k+6}.  Computes rank_p J_{k+6}^{pi} (pi = parity of A*g), compares with the Hilbert-series cap
(rigorous rank over Q when attained), then the rank of [J | A*x^beta for all ternary beta].  If the augmented rank
rises by the number of cone monomials, no nonzero cone of degree k satisfies E1 (rigorous over Q when the cap is attained).
Usage: python3 wave2_construct_cones.py k sector [prime]"""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from flint import nmod_mat
k, sector = int(sys.argv[1]), sys.argv[2]
p = int(sys.argv[3]) if len(sys.argv) > 3 else core.PRIMES[0]
LOG = '/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl'
t0 = time.time()
F, A, q, rec = core.load_instance('dense-rational', mod=p)
Fi = [core.pderiv(F, i, p) for i in range(6)]
m = k + 6
# parity bookkeeping: Y-parity of A = 1; of g(X) = 0; of g(Y) = k mod 2
par_g = 0 if sector == 'X' else (k % 2)
par = (1 + par_g) % 2
cols = core.monomials(m, par)              # target monomials of degree m, parity par
cidx = {e: j for j, e in enumerate(cols)}
# Hilbert series caps with parity: S and R = C[z]/(regular sequence of six quintics), parity grading s (s^2=1)
def series(deg):
    # returns dict (d, parity) -> dim for S and R up to degree deg
    import collections
    S = collections.Counter(); R = collections.Counter()
    for e in itertools.product(range(deg+1), repeat=6):
        d = sum(e)
        if d > deg: continue
        pe = (e[3]+e[4]+e[5]) % 2
        S[(d, pe)] += 1
        if all(x <= 4 for x in e): R[(d, pe)] += 1
    return S, R
S, R = series(m)
cap = S[(m, par)] - R[(m, par)]
assert S[(m, par)] == len(cols)
rows = []
for i in range(6):
    pf = 1 if i >= 3 else 0
    for mono in core.monomials(m - 5, (par - pf) % 2):
        rows.append(core.shift(Fi[i], mono, 1, p))
nJ = len(rows)
cones = []
cone_vars = (0, 1, 2) if sector == 'X' else (3, 4, 5)
for e3 in itertools.product(range(k+1), repeat=3):
    if sum(e3) != k: continue
    e = [0]*6
    for t, v in enumerate(cone_vars): e[v] = e3[t]
    cones.append(core.shift(A, tuple(e), 1, p))
M = nmod_mat(nJ + len(cones), len(cols), p)
for r, poly in enumerate(rows + cones):
    for e, c in poly.items():
        M[r, cidx[e]] = c
rkJ = nmod_mat(nJ, len(cols), p)
# rank of J part: use a sub-matrix copy
MJ = nmod_mat(nJ, len(cols), p)
for r, poly in enumerate(rows):
    for e, c in poly.items():
        MJ[r, cidx[e]] = c
rJ = MJ.rank(); t1 = time.time()
rA = M.rank(); t2 = time.time()
ncone = len(cones)
res = {'exp': 'cone-E1-test', 'sector': sector, 'k': k, 'm': m, 'parity_target': par, 'prime': p, 'cols': len(cols),
       'gens_J': nJ, 'rank_J': rJ, 'cap_J': cap, 'cap_attained': rJ == cap, 'n_cones': ncone, 'rank_aug': rA,
       'dim_cones_in_JA': ncone - (rA - rJ), 'rigorous': (rJ == cap),
       'result': ('NO cone of degree %d in sector %s satisfies E1' % (k, sector)) if rA - rJ == ncone else ('%d-dim space of E1-admissible cones' % (ncone - (rA - rJ))),
       'seconds': [round(t1-t0,1), round(t2-t1,1)], 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
print(json.dumps(res), flush=True)
with open(LOG, 'a') as f: f.write(json.dumps(res) + '\n')
