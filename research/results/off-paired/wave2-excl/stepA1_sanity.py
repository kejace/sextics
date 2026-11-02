"""A1: independent re-verification of the basic rigorous facts (exact over Q where possible, two new primes otherwise)."""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import exlib as X
from exlib import *
from hodge import core, fixed_d
from fractions import Fraction

F, A, q = load()
assert q == 1
D0 = core.parse('X0**4 + X0*X1*Y0*Y1/4 + X1*X2*Y0*Y2/9')
out = {}

# 0. builder self-check against core.defect (sign convention)
rng = random.Random(5)
cols, colpolys, rhs = fixedD_columns(F, A, q, D0, p=None, sector=False)
N = [{} for _ in range(NV)]; Q = {}
x = {}
for j, (kind, i, m) in enumerate(cols):
    v = Fraction(rng.randrange(-5, 6))
    if v == 0: continue
    x[j] = v
    if kind == 'N': N[i][m] = N[i].get(m, 0) + v
    else: Q[m] = Q.get(m, 0) + v
E = core.defect(F, A, D0, q, N)
E = sub(E, mul(mul(F, F), Q))
S = mul(A, mul(D0, D0))
for j, v in x.items(): S = add(S, scale(colpolys[j], v))
assert sub(E, S) == {}, 'builder sign mismatch!'
out['builder_matches_core_defect'] = True
print('builder matches core.defect')

# 1. A not in J(F): exact over Q
t0 = time.time()
Fi = [pderiv(F, i) for i in range(NV)]
gens = [shift(Fi[i], m, 1) for i in range(NV) for m in monomials(1)]
mons = monomials(6); idx = {m: k for k, m in enumerate(mons)}
from flint import fmpq_mat, fmpq
M = fmpq_mat(len(mons), len(gens)+1)
for j, g in enumerate(gens):
    for e, v in g.items(): M[idx[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
for e, v in A.items(): M[idx[e], len(gens)] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
M0 = fmpq_mat(len(mons), len(gens))
for j, g in enumerate(gens):
    for e, v in g.items(): M0[idx[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
r0, r1 = M0.rank(), M.rank()
out['A_in_J'] = {'rank_J6_odd': r0, 'rank_aug': r1, 'A_in_J': r0 == r1, 'exact': True, 'seconds': round(time.time()-t0, 2)}
print(out['A_in_J'])

# 2. colon dimensions d=4..8 (sum over characters) for two new primes
for p in PRIMES[:2]:
    dims = {}
    for d in range(4, 9):
        tot = 0; totg = 0
        for a in range(6):
            for b in range(6):
                chosen, syz, info = colon_sector_and_lifts(F, A, q, d, (a, b), p=p)
                tot += info['colon_dim']; totg += info['gauge_dim']
        dims[d] = (tot, totg)
    out[f'colon_dims_p{p}'] = dims
    print('colon dims (colon, gauge) mod', p, dims)

# 3. D0: exact over Q ranks (sector-restricted and full) and rational certificate
t0 = time.time()
cols, colpolys, rhs = fixedD_columns(F, A, q, D0, p=None, sector=True)
r, rb, R, C = exact_ranks_Q(colpolys, rhs)
out['D0_sector_exact'] = {'rank': r, 'rank_aug': rb, 'rows': R, 'cols': C, 'seconds': round(time.time()-t0, 1)}
print(out['D0_sector_exact'])
cert = inconsistency_certificate(colpolys, rhs, verbose=True)
out['D0_sector_cert'] = {k: v for k, v in cert.items() if k != 'y'}
out['D0_sector_cert']['reverified'] = verify_certificate_exact(colpolys, rhs, cert['y']) if cert['ok'] else None
json.dump(cert, open('cert_D0_sector.json', 'w'))
print(out['D0_sector_cert'])
t0 = time.time()
cols, colpolys, rhs = fixedD_columns(F, A, q, D0, p=None, sector=False)
r, rb, R, C = exact_ranks_Q(colpolys, rhs)
out['D0_full_exact'] = {'rank': r, 'rank_aug': rb, 'rows': R, 'cols': C, 'seconds': round(time.time()-t0, 1)}
print(out['D0_full_exact'])

# 4. Fermat control: F0 = f6+g6, D = f6 -> consistent (exact)
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6'); f6 = core.parse('X0**6+X1**6+X2**6')
cols, colpolys, rhs = fixedD_columns(F0, A, q, f6, p=None, sector=False)
r, rb, R, C = exact_ranks_Q(colpolys, rhs)
out['fermat_control_exact'] = {'rank': r, 'rank_aug': rb, 'consistent': r == rb, 'rows': R, 'cols': C}
print(out['fermat_control_exact'])
# explicit witness check
Nw = [core.parse('X0**4*Y0**3/6'), core.parse('X0**3*X1*Y0**3/6'), core.parse('X0**3*X2*Y0**3/6'), {}, {}, {}]
ok, info = core.verify_witness(F0, A, q, f6, Nw, {})
out['fermat_witness_verified'] = ok
print('fermat witness', ok)
# and f6 on off-paired fails e1 (exact over Q): A*f6 in J(F)_12 ?
t0 = time.time()
gens = [shift(Fi[i], m, 1) for i in range(NV) for m in monomials(7)]
tgt = mul(A, f6)
mons = sorted({e for g in gens for e in g} | set(tgt)); idx = {m: k for k, m in enumerate(mons)}
for p in PRIMES[:2]:
    from flint import nmod_mat
    flat = [0]*(len(mons)*(len(gens)+1))
    for j, g in enumerate(gens):
        for e, v in g.items(): flat[idx[e]*(len(gens)+1)+j] = (Fraction(v).numerator * pow(Fraction(v).denominator, -1, p)) % p
    for e, v in tgt.items(): flat[idx[e]*(len(gens)+1)+len(gens)] = (Fraction(v).numerator * pow(Fraction(v).denominator, -1, p)) % p
    Mb = nmod_mat(len(mons), len(gens)+1, flat, p)
    M0 = nmod_mat(len(mons), len(gens), [int(Mb[i, j]) for i in range(len(mons)) for j in range(len(gens))], p)
    out[f'A_f6_in_J_p{p}'] = {'rank': M0.rank(), 'rank_aug': Mb.rank()}
    print('A*f6 in J mod', p, out[f'A_f6_in_J_p{p}'], round(time.time()-t0, 1), 's')

json.dump(out, open('stepA1_out.json', 'w'), indent=1)
log({'step': 'A1', 'family': 'sanity', 'summary': {k: v for k, v in out.items()}})
print('DONE')
