"""Served-target spaces mod p for a list of denominators; union span; membership of A; decomposition.

usage: python3 serves.py <listfile> [tag]
listfile lines:  name | <D expression>   (same env as run_list2.py)
Writes served_<tag>.json with, per D, the served basis (mod p) and the union analysis.
"""
import sys, time, json, pickle, os
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
from flint import nmod_mat

p = core.PRIMES[0]
WD = qcp.WD
listfile = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(listfile).replace('.txt', '')
F_ = qcp.F
A_ = qcp.A
delta = P(open(f'{WD}/delta.txt').read())


def Fpart(xdeg):
    return {e: c for e, c in F_.items() if e[0] + e[1] + e[2] == xdeg}


def EI(I, G=None):
    G = F_ if G is None else G
    out = {}
    for i in I:
        out = core.add(out, core.mul(qcp.var(i), core.pderiv(G, i)))
    return out


env = {'P': P, 'Q': qcp.Q, 'core': core, 'delta': delta, 'F': F_, 'Fpart': Fpart, 'EI': EI}
odd6 = core.monomials(6, 1)
col = {e: j for j, e in enumerate(odd6)}
Ap = core.reduce_mod(A_, p)
Avec = [0] * len(odd6)
for e, v in Ap.items():
    Avec[col[e]] = v

results = []
allrows = []
owners = []
items = [l.strip() for l in open(listfile) if l.strip() and not l.startswith('#')]
for line in items:
    parts = line.split('|')
    name, Ds = parts[0].strip(), parts[1].strip()
    D = eval(Ds, env)
    if isinstance(D, str):
        D = P(D)
    t0 = time.time()
    dim, basis, S_ = fixed_d.target_image(F_, 1, D, p=p)
    dt = time.time() - t0
    rows = []
    for b in basis:
        r = [0] * len(odd6)
        for e, v in b.items():
            r[col[e]] = v % p
        rows.append(r)
    # is A in this D's served space?
    M = nmod_mat(len(rows), len(odd6), [x for r in rows for x in r], p) if rows else None
    inA = False
    if rows:
        M2 = nmod_mat(len(rows) + 1, len(odd6), [x for r in rows for x in r] + Avec, p)
        inA = (M2.rank() == M.rank())
    rec_ = {'name': name, 'D': Ds, 'd': core.degree(D), 'imgdim': dim, 'A_served': inA, 'seconds': round(dt, 1)}
    print(json.dumps(rec_), flush=True)
    results.append(rec_)
    qcp.log(D, name, None, None, -1, -1, dt, f'served-space: imgdim {dim}, A in image: {inA}', extra={'imgdim': dim, 'A_served': inA})
    allrows.extend(rows)
    owners.extend([name] * len(rows))
    # union so far
    if allrows:
        U = nmod_mat(len(allrows), len(odd6), [x for r in allrows for x in r], p)
        ru = U.rank()
        U2 = nmod_mat(len(allrows) + 1, len(odd6), [x for r in allrows for x in r] + Avec, p)
        inU = (U2.rank() == ru)
        print(f'   union rank so far = {ru}   A in union: {inU}', flush=True)
        results[-1]['union_rank'] = ru
        results[-1]['A_in_union'] = inU
with open(f'{WD}/served_{tag}.json', 'w') as f:
    json.dump({'p': p, 'results': results}, f, indent=1)
with open(f'{WD}/served_rows_{tag}.pkl', 'wb') as f:
    pickle.dump({'p': p, 'odd6': odd6, 'rows': allrows, 'owners': owners, 'items': items}, f)
print('done', flush=True)
