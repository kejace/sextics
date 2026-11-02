import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from flint import nmod_mat, fmpq, fmpq_mat
p2 = core.PRIMES[0]
res = hess_quadrics(6, 0, p=p2, verbose=True)
write_msolve(res, f'{WD}/hess6p2-quadrics-d6-p0.ms')
print('second prime', p2, 'rank', res['rank_quadrics'])
# rational reconstruction test of the annihilator functionals (RREF basis of null space of J_12^- generators)
Fp = reduce_mod(F, p2)
gens = jacobian.jacobian_generators(Fp, 12, p2, parity=1)
Rs = Space(12, 1)
ng = len(gens)
flat = [0] * (ng * Rs.n)
for r, (_, _, P) in enumerate(gens):
    for e, v in P.items():
        flat[r * Rs.n + Rs.idx[e]] = v
G = nmod_mat(ng, Rs.n, flat, p2)
Cn, dimR = G.nullspace()
# RREF of the annihilator (as rows): transpose and rref
Ct = nmod_mat(dimR, Rs.n, [int(Cn[i, j]) for j in range(dimR) for i in range(Rs.n)], p2)
Cr, rk = Ct.rref()
ok = 0; bad = 0; maxden = 0
t0 = time.time()
for i in range(rk):
    row = [core.ratrecon(int(Cr[i, j]), p2) for j in range(Rs.n)]
    if any(v is None for v in row):
        bad += 1; continue
    ok += 1
    maxden = max(maxden, max(v.denominator for v in row))
print(f'annihilator RREF rows: {rk}; rationally reconstructible rows: {ok}, failed: {bad}, max denominator {maxden} ({time.time()-t0:.1f}s)')
