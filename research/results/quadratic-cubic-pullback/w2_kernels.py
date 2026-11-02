#!/usr/bin/env python3
"""w2_kernels.py -- wave 2, item 1: exact E1-admissible denominator spaces for quadratic-cubic-pullback.

K_d^par = {D in S_d^par : A*D in J(F)}  (mod p = 2^31-1), computed with hodge.jacobian.mult_kernel.
For each (d, parity) requested: appends one JSON line to w2_kernels.jsonl (dims, sizes, timing) and writes
w2_kernel_<d>_<par>.json with the full kernel basis (vectors over the monomials of S_d^par, mod p) and the
basis of K_d/J_d in standard monomials.  One process per (d, parity).

usage: python3 w2_kernels.py --jobs 9:0,10:0,11:0,12:0,7:1,8:1,9:1,10:1,11:1,12:1 --procs 8
"""
import sys, os, json, time, argparse
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import multiprocessing as mp

LANE = 'quadratic-cubic-pullback'
OUT = 'w2_kernels.jsonl'


def run_one(d, par):
    from hodge import core, jacobian
    F, A, q, _ = core.load_instance(LANE)
    t0 = time.time()
    p = core.PRIMES[0]
    # patch: mult_kernel returns only the quotient basis; we also want the full kernel vectors -> re-implement the
    # projection part here (same code path) so the full basis is saved.
    from hodge.core import monomials, shift, reduce_mod
    from flint import nmod_mat
    Ap = reduce_mod(A, p)
    Dm = monomials(d, par)
    K = len(Dm)
    target_par = (core.parity_of(Ap) + par) % 2
    gens = jacobian.jacobian_generators(F, d + 6, p, target_par)
    cols = [shift(Ap, m, 1, p) for m in Dm] + [g[2] for g in gens]
    rows = {}
    for col in cols:
        for e in col:
            if e not in rows:
                rows[e] = len(rows)
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    del flat
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    rankM = C - nullity
    del M
    vecs = []
    if nullity:
        proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], p)
        del X
        Pr, rk = proj.rref()
        for i in range(rk):
            vecs.append({Dm[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])})
    kernel_dim = len(vecs)
    red = jacobian._Reducer(F, d, p, par)
    reduced = [red.normal_form(v) for v in vecs]
    reduced = [v for v in reduced if v]
    basis = []
    if reduced:
        sm = red.standard
        sidx = {m: i for i, m in enumerate(sm)}
        flat = [0] * (len(reduced) * len(sm))
        for r, v in enumerate(reduced):
            for e, c in v.items():
                flat[r * len(sm) + sidx[e]] = c
        Q = nmod_mat(len(reduced), len(sm), flat, p)
        Qr, rk = Q.rref()
        for i in range(rk):
            basis.append({sm[j]: int(Qr[i, j]) for j in range(len(sm)) if int(Qr[i, j])})
    t3 = time.time()
    rec = {'lane': LANE, 'item': 'kernel', 'd': d, 'parity': par, 'p': p, 'rows': R, 'cols': C, 'nD': K,
           'rank_stack': rankM, 'nullity_stack': nullity, 'kernel_dim': kernel_dim, 'jacobian_dim': red.jacobian_dim,
           'ring_dim': red.ring_dim, 'dim_K_mod_J': len(basis), 'ring_dim_target': None,
           't_build': round(t1 - t0, 1), 't_null': round(t2 - t1, 1), 't_reduce': round(t3 - t2, 1),
           'seconds': round(t3 - t0, 1), 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    # ring dim in the target degree = rows - rank(J block) ; rank(J block) = rankM - (kernel-free A columns)
    rec['rank_J_target'] = rankM - (K - kernel_dim)
    rec['ring_dim_target'] = R - rec['rank_J_target']
    with open(OUT, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    def ser(v):
        return [[list(e), c] for e, c in v.items()]
    with open(f'w2_kernel_{d}_{par}.json', 'w') as f:
        json.dump({'d': d, 'parity': par, 'p': p, 'monomials': [list(m) for m in Dm],
                   'kernel_basis': [ser(v) for v in vecs], 'quotient_basis': [ser(v) for v in basis],
                   'standard_monomials': [list(m) for m in red.standard], 'record': rec}, f)
    print(json.dumps(rec), flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', default='9:0,10:0,7:1,8:1,9:1,10:1,11:0,12:0,11:1,12:1')
    ap.add_argument('--procs', type=int, default=8)
    a = ap.parse_args()
    jobs = []
    for tok in a.jobs.split(','):
        d, par = tok.split(':'); jobs.append((int(d), int(par)))
    ctx = mp.get_context('fork')
    with ctx.Pool(a.procs) as pool:
        for rec in pool.starmap(run_one, jobs, chunksize=1):
            pass


if __name__ == '__main__':
    main()
