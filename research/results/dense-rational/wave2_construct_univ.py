#!/usr/bin/env python3
"""E1 test for the universal degree-12 denominators built from the Y-degree pieces a=F60, b=F42, c=F24, d=F06:
   the residual-tangency divisor R' = closure{p in X : the line through [X(p)],[Y(p)] is tangent to X elsewhere},
   D_R' = (c + d)^2 + 4 a d  (discriminant of the residual quadratic of a + b t + c t^2 + d t^3 after removing t = 1),
   and the whole 10-dim space spanned by the products of two pieces.  E1 <=> A*D in J_18^- (Y-odd part).
Reports rank_p J_18^- against the Hilbert cap (rigorous over Q when attained) and the E1-admissible dimension.
Usage: python3 wave2_construct_univ.py [prime]"""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else core.PRIMES[0]
LOG = '/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl'
t0 = time.time()
F, A, q, rec = core.load_instance('dense-rational', mod=p)
Fi = [core.pderiv(F, i, p) for i in range(6)]
pieces = {}
for e, cf in F.items():
    yd = e[3] + e[4] + e[5]
    pieces.setdefault(yd, {})[e] = cf
a, b, c, d = pieces[0], pieces[2], pieces[4], pieces[6]
names = ['a', 'b', 'c', 'd']; P4 = [a, b, c, d]
prods, pnames = [], []
for i, j in itertools.combinations_with_replacement(range(4), 2):
    prods.append(core.mul(P4[i], P4[j], p)); pnames.append(names[i] + names[j])
cd = core.add(c, d, p)
DR = core.add(core.mul(cd, cd, p), core.scale(core.mul(a, d, p), 4, p), p)
m = 18; par = 1
cols = core.monomials(m, par); cidx = {e: j for j, e in enumerate(cols)}
import collections
S = collections.Counter(); R = collections.Counter()
for e in itertools.product(range(m+1), repeat=6):
    dd = sum(e)
    if dd > m: continue
    pe = (e[3]+e[4]+e[5]) % 2
    S[(dd, pe)] += 1
    if all(x <= 4 for x in e): R[(dd, pe)] += 1
cap = S[(m, par)] - R[(m, par)]
rows = []
for i in range(6):
    pf = 1 if i >= 3 else 0
    for mono in core.monomials(m - 5, (par - pf) % 2):
        rows.append(core.shift(Fi[i], mono, 1, p))
nJ = len(rows)
extra = [core.mul(A, DR, p)] + [core.mul(A, pr, p) for pr in prods]
def rank_of(polys):
    M = nmod_mat(len(polys), len(cols), p)
    for r, poly in enumerate(polys):
        for e, cf in poly.items(): M[r, cidx[e]] = cf
    return M.rank()
rJ = rank_of(rows); t1 = time.time()
rR = rank_of(rows + extra[:1]); t2 = time.time()
rAll = rank_of(rows + extra[1:]); t3 = time.time()
res = {'exp': 'universal-deg12-E1-test', 'prime': p, 'm': m, 'parity': par, 'cols': len(cols), 'gens_J': nJ, 'rank_J': rJ,
       'cap_J': cap, 'cap_attained': rJ == cap, 'rank_with_A_DRprime': rR, 'DRprime_E1': (rR == rJ),
       'n_products': len(prods), 'rank_with_A_products': rAll, 'dim_E1_admissible_products': len(prods) - (rAll - rJ),
       'DRprime_formula': '(c+d)^2 + 4 a d, a=F60 b=F42 c=F24 d=F06', 'seconds': [round(t1-t0,1), round(t2-t1,1), round(t3-t2,1)],
       'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
print(json.dumps(res), flush=True)
with open(LOG, 'a') as f: f.write(json.dumps(res) + '\n')
