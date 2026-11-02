"""Rigorous E2-mod-J test for D = D7^+ alone: with an EXACT rational N0 (N0(F) = A*D7), any certificate with
denominator D7 forces B = D7*div N0 - N0(D7) in J_14^- (derivation in the verify note).  Non-membership mod p with
rank_p J_14^- = Hilbert value (5052) is a proof over Q."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs); A = poly_from_sympy(As)
D7 = poly_from_sympy(sp.sympify(json.load(open('/tmp/hodge/results/short-power-presentation/e1vec_d7_p0.json'))[0]['D'],
                     locals=dict(zip([str(v) for v in VARS], VARS))))
raw = json.load(open('d7_N_exact.json'))['N']
N0 = [{tuple(int(x) for x in k.strip('()').split(',')): Fraction(v) for k, v in Ni.items()} for Ni in raw]
# re-verify exactly
Fi = [pderiv(F, i) for i in range(6)]
tot = {}
for i in range(6): tot = add(tot, mul(N0[i], Fi[i]))
assert add(tot, scale(mul(A, D7), -1)) == {}, 'N0 does not solve N0(F) = A D7'
print('N0 exact: OK', flush=True)
divN = {}
for i in range(6): divN = add(divN, pderiv(N0[i], i))
ND = {}
for i in range(6): ND = add(ND, mul(N0[i], pderiv(D7, i)))
B = add(mul(D7, divN), scale(ND, -1))
print('deg B =', degree(B), 'parity', parity_of(B), 'terms', len(B), flush=True)
H = hilb()
for p in (1000000007, 998244353):
    t0 = time.time()
    Fp = reduce_mod(F, p); Bp = reduce_mod(B, p)
    J = jac_gens(Fp, 14, 1, p)
    rows = {}
    for P in J + [Bp]:
        for e in P: rows.setdefault(e, len(rows))
    rJ = nmod_from_polys(J, p, rows)[0].rank()
    rJB = nmod_from_polys(J + [Bp], p, rows)[0].rank()
    exp = H[14]['J_minus']
    verdict = 'B NOT in J_14^- (rigorous over Q since rank J matches Hilbert)' if (rJB == rJ + 1 and rJ == exp) else ('B in J mod p' if rJB == rJ else 'rank mismatch')
    print(f'p={p}: rank J_14^- = {rJ} (expected {exp}), rank [J|B] = {rJB} -> {verdict} ({time.time()-t0:.0f}s)', flush=True)
    logjson(claim='D7^+ alone is not a denominator: B_{D7}(N0) in J_14^- fails (exact N0)', verdict='CONFIRMED' if rJB == rJ + 1 and rJ == exp else 'UNVERIFIED', p=p, rankJ=rJ, expected=exp, rankJB=rJB)
