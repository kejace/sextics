#!/usr/bin/env python3
"""referee_membership.py -- replay of the Jacobian-membership table (dossier 4.2/4.6, excl ledger) at the unused
prime 1000000009, with the regular-sequence rank check that makes NON-membership rigorous over Q
(rank_Q[J|g] >= rank_p[J|g] = rank_p J + 1 = rank_Q J + 1).  Membership mod p is only evidence.
Shares the polynomial helpers of referee_replay.py (own arithmetic, own parser)."""
import sys, json, time
from fractions import Fraction
sys.argv = [sys.argv[0]]  # referee_replay runs nothing when imported
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
import referee_replay as R
from flint import nmod_mat

p = 1000000009
P = R.parse
f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
L = P('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2')
M = P('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
L2M = R.mul(R.mul(L, L), M)
assert R.sub(R.F, R.add(R.add(f, g), R.scale(L2M, Fraction(5, 7)))) == {}, 'F != f + g + (5/7) L^2 M'
EX = {}; EY = {}
for i in range(3): EX = R.add(EX, R.mul(R.mono(tuple(1 if j == i else 0 for j in range(6))), R.Fi[i]))
for i in range(3, 6): EY = R.add(EY, R.mul(R.mono(tuple(1 if j == i else 0 for j in range(6))), R.Fi[i]))
assert R.sub(R.scale(R.sub(EX, R.scale(EY, 2)), Fraction(1, 6)), R.sub(f, R.scale(g, 2))) == {}, 'f - 2g != (E_X F - 2 E_Y F)/6'
print('[exact] F = f + g + (5/7) L^2 M; f - 2g = (E_X F - 2 E_Y F)/6 in J: True (exact identities)')

tests = {'L2M': L2M, 'f': f, 'g': g, 'f-g': R.sub(f, g), 'A': R.A,
         'A*L': R.mul(R.A, L), 'A*M': R.mul(R.A, M), 'A*L^2': R.mul(R.A, R.mul(L, L)), 'A*M^2': R.mul(R.A, R.mul(M, M)), 'A*L*M': R.mul(R.A, R.mul(L, M)),
         'f*g': R.mul(f, g), 'f*L2M': R.mul(f, L2M), 'g*L2M': R.mul(g, L2M), '(L2M)^2': R.mul(L2M, L2M), 'f^2': R.mul(f, f), 'g^2': R.mul(g, g),
         'E_X F': EX, 'E_Y F': EY, 'F': R.F}
out = {}
cache = {}
for name, G in tests.items():
    k = R.degree(G)
    if k not in cache:
        rows = R.monomials(k); ridx = {e: i for i, e in enumerate(rows)}
        gens = [R.mul(R.mono(m), R.Fi[i]) for i in range(6) for m in R.monomials(k - 5)] if k >= 5 else []
        cache[k] = (rows, ridx, gens)
    rows, ridx, gens = cache[k]
    def build(colist):
        Rn, C = len(rows), len(colist); flat = [0] * (Rn * C)
        for j, Q in enumerate(colist):
            for e, v in Q.items(): flat[ridx[e] * C + j] = R.modp(v, p)
        return nmod_mat(Rn, C, flat, p)
    rJ = build(gens).rank() if gens else 0
    expJ = len(rows) - R.fermat_hf(k)
    rA = build(gens + [G]).rank()
    inJ = (rA == rJ)
    out[name] = dict(degree=k, rank_J=rJ, rank_J_expected=expJ, rankJ_ok=(rJ == expJ), in_J_modp=inJ,
                     verdict=('in J mod p (evidence only)' if inJ else ('NOT in J: rigorous over Q' if rJ == expJ else 'NOT in J mod p (J rank off!)')))
    print(f'[membership p={p}] {name:8s} deg {k:2d}: rank J = {rJ} (expected {expJ}), in J mod p = {inJ} -> {out[name]["verdict"]}', flush=True)
json.dump(out, open('/tmp/hodge/results/rank-one-biquartic/wave2-excl/referee_membership.json', 'w'), indent=1)
