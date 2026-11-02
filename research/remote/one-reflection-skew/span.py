"""Span analysis: accumulate served spaces (mod p, mod J_f) of all tested D and test A membership."""
import sys, pickle, os, glob
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from flint import nmod_mat
p = P1
tm = mons5(6, 1); idx = {m: k for k, m in enumerate(tm)}; nT = len(tm)
def vec(P):
    v = [0] * nT
    for e, c in P.items(): v[idx[e]] = int(c) % p
    return v
def sigma(P):  # (X1,Y1) -> -(X1,Y1)
    return {e: (c if (e[1] + e[4]) % 2 == 0 else (-c) % p if isinstance(c, int) else -c) for e, c in P.items()}
Jrows = [vec(P) for P in jacobian6_odd(p)]
rJ = nmod_mat(len(Jrows), nT, [x for r in Jrows for x in r], p).rank()
recs = []
for fn in sorted(glob.glob('served/*.pkl')):
    with open(fn, 'rb') as fh: recs.append(pickle.load(fh))
active = [r for r in recs if r['served_mod_J'] > 0]
print(f'{len(recs)} tested, {len(active)} active (served_mod_J>0); dim J6odd={rJ}, nT={nT}, (S/J)_6^odd = {nT - rJ}')
for r in sorted(active, key=lambda r: (-r['served_mod_J'], r['d'])):
    print(f"  {r['name'][:60]:60s} d={r['d']} served_mod_J={r['served_mod_J']} A={r['serves_A']} A+={r['serves_A+']} A-={r['serves_A-']}")
rows = list(Jrows)
for r in active:
    for b in r['served_basis']:
        rows.append(vec(b)); rows.append(vec(sigma(b)))
M = nmod_mat(len(rows), nT, [x for r in rows for x in r], p)
rk = M.rank()
print('rank of J + sum of all served spaces (with sigma images):', rk, ' => span mod J has dim', rk - rJ)
for name, Aq in (('A', A), ('A+', A_plus), ('A-', A_minus)):
    M2 = nmod_mat(len(rows) + 1, nT, [x for r in rows for x in r] + vec(core.reduce_mod(Aq, p)), p)
    print(f'  {name} in span: {M2.rank() == rk}')
# sigma-even part of the target sector: project
