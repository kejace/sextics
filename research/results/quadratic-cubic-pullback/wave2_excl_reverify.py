"""Light local re-verification of wave-1 rigorous structural claims at a second prime (and exact where cheap)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, Qs, monomials, mul, add, sub, scale, power, pderiv, logline, rank_of, hilb_R, hilb_S, dim_J
from hodge import jacobian, symmetry
from flint import nmod_mat
LOG = w2x.LANE + '/wave2-excl.log.jsonl'
p = core.PRIMES[1]
t0 = time.time()
# (1) Hilbert function sanity: totals vs known 1,6,21,56,126,246,426,666,951,1246,1506
tot = [hilb_R(k) for k in range(0, 13)]
print('hilb_R totals', tot)
print('hilb_R parity 6:', hilb_R(6,0), hilb_R(6,1), ' 10:', hilb_R(10,0), hilb_R(10,1), ' 16:', hilb_R(16,0), hilb_R(16,1))
logline(LOG, dict(exp='hilbert', note='parity-split Hilbert function of S/J(F) from Koszul series', totals=tot,
                  split={k: [hilb_R(k,0), hilb_R(k,1)] for k in range(0, 25)}))
# (2) counts 226 / 18 / 208 and A, delta independent mod J_6 (prime p)
amons = monomials(6, 1)
Jb = w2x.jacobian_part(1, p)
rJ = rank_of(Jb, amons, p)
delta = w2x.delta_poly()
Ap, dp = core.reduce_mod(A, p), core.reduce_mod(delta, p)
rA = rank_of(Jb + [Ap], amons, p); rD = rank_of(Jb + [dp], amons, p); rAD = rank_of(Jb + [Ap, dp], amons, p)
rec = dict(exp='reverify-h31', prime=p, odd_monomials=len(amons), J6_odd_rank=rJ, hilbJ6_odd=dim_J(6,1), h31_plus=len(amons)-rJ,
           A_notin_J=(rA == rJ+1), delta_notin_J=(rD == rJ+1), A_delta_independent_modJ=(rAD == rJ+2),
           rigorous='ranks mod p are lower bounds; rank J6_odd = Hilbert value certifies; independence statements rigorous')
print(rec); logline(LOG, rec)
# (3) symmetry group via engine (independent of struct1's Smith form)
G = symmetry.diagonal_group(F)
rec = dict(exp='reverify-symmetry', group_order=G.order if hasattr(G,'order') else None, note=str(G)[:300])
print(rec); logline(LOG, rec)
# (4) F not a join: rank of 21 second partials at prime p
second = [pderiv(pderiv(core.reduce_mod(F,p), i, p), j, p) for i in range(6) for j in range(i, 6)]
r2 = rank_of(second, monomials(4), p)
rec = dict(exp='reverify-notjoin', prime=p, rank_second_partials=r2, note='21 = max; a join G(l)+H(l\') has rank <= 12; rigorous (lower bound attains max)')
print(rec); logline(LOG, rec)
# (5) cubics in Q meet J_6^even only in C.F
import itertools
cubQ = []
for c in itertools.combinations_with_replacement(range(6), 3):
    cubQ.append(core.reduce_mod(mul(mul(Qs[c[0]], Qs[c[1]]), Qs[c[2]]), p))
emons = monomials(6, 0)
Je = w2x.jacobian_part(0, p)
rJe = rank_of(Je, emons, p); rC = rank_of(cubQ, emons, p); rU = rank_of(Je + cubQ, emons, p)
rec = dict(exp='reverify-cubicsQ', prime=p, dim_cubicsQ=rC, J6_even_rank=rJe, hilbJ6_even=dim_J(6,0), dim_intersection_modp=rC + rJe - rU,
           note='intersection dim mod p >= rational (certified J rank => rational intersection <= modular); F itself is in both so rational dim >= 1')
print(rec); logline(LOG, rec)
# (6) D in J checks for the wave-1 degree-6 families at prime p (rigorous non-membership when J rank certified)
def bideg(e): return (e[0]+e[1]+e[2], e[3]+e[4]+e[5])
parts = {}
for e, c in F.items(): parts.setdefault(bideg(e), {})[e] = c
tests = {'delta': delta, 'join012': add(add(power(Qs[0],3), power(Qs[1],3)), power(Qs[2],3)),
         'join345': add(add(power(Qs[3],3), power(Qs[4],3)), power(Qs[5],3)), 'F60': parts[(6,0)], 'F06': parts[(0,6)],
         'fermatX': core.parse('X0**6+X1**6+X2**6'), 'Q0cube+Q3cube': add(power(Qs[0],3), power(Qs[3],3)), 'EX_F': w2x.euler([1,1,1,0,0,0]),
         'Q0sq': power(Qs[0], 2), 'F-Q0cube': sub(F, power(Qs[0],3))}
res = {}
for nm, D in tests.items():
    inJ, info = jacobian.in_jacobian(F, D, p)
    par = core.parity_of(D)
    res[nm] = dict(in_J=inJ, rankJ=info['rank'], hilbJ=dim_J(6, par), certified=(info['rank'] == dim_J(6, par)))
rec = dict(exp='reverify-J6-membership', prime=p, results=res)
print(json.dumps(rec)); logline(LOG, rec)
# (7) two-plane (4,4)-cancellation system nullity (exact over Q with Fractions via flint fmpq)
from flint import fmpq_mat, fmpq
F24, F42 = parts[(2,4)], parts[(4,2)]
monX2 = [e for e in monomials(2,0) if bideg(e)==(2,0)]; monY2 = [e for e in monomials(2,0) if bideg(e)==(0,2)]
cols = [core.shift(F24, m, 1) for m in monX2] + [core.shift(F42, m, 1) for m in monY2]
rows = sorted(set(e for c in cols for e in c))
M = fmpq_mat(len(rows), len(cols))
for j, c in enumerate(cols):
    for e, v in c.items():
        M[rows.index(e), j] = fmpq(v.numerator, v.denominator)
rk = M.rank()
rec = dict(exp='reverify-twoplane44', unknowns=len(cols), rank_exact=rk, nullity_exact=len(cols)-rk, note='exact over Q (fmpq)')
print(rec); logline(LOG, rec)
print('total', round(time.time()-t0,1), 's')
