import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
F, A, q = load(); p = P0
G = GROUPS['triv']()
mons = [core.power(var(0), a) for a in range(1)]  # placeholder
# rebuild the ternX d=7 family exactly as w2_struct does
def monomials_in(coords, d):
    n = len(coords); out = []
    def rec(i, left, cur):
        if i == n - 1: out.append(cur + [left]); return
        for k in range(left + 1): rec(i + 1, left - k, cur + [k])
    rec(0, d, [])
    polys = []
    for ex in out:
        P = {(0,) * 6: Fraction(1)}
        for c, k in zip(coords, ex):
            if k: P = core.mul(P, core.power(c, k))
        polys.append(P)
    return polys
mons = [P for P in monomials_in([var(0), var(1), var(2)], 7) if core.parity_of(P) == 0]
basis = restricted_family(mons, G, p)
print('basis', len(basis), 'invariant?', [G.is_semi_invariant(B) for B in basis])
dimK, kvecs, _, _ = e1_family_modp(F, A, q, basis, G, 0, p)
print('dimK', dimK, 'kv', kvecs)
D = {}
for j, Bj in enumerate(basis):
    if kvecs[0][j]:
        D = core.add(D, core.scale(core.reduce_mod(Bj, p), kvecs[0][j] % p, p), p)
print('D terms', len(D), 'semi-inv?', G.is_semi_invariant(D))
# inspect: which perm fails
from charsolve import act
for perm, s in G.perms[:6]:
    Dg = act(D, perm)
    diff = [(e, D.get(e), Dg.get(e)) for e in set(D) | set(Dg) if D.get(e) != Dg.get(e)]
    print(perm, 'mismatches', len(diff), diff[:2])
