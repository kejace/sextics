"""Block Macaulay closure for a G-graded homogeneous quadric ideal (quads JSON with per-variable characters):
I_k = S_k  (=> V(I) = {0})  is tested block by block (character psi of degree-k monomials): rank of the matrix whose rows
are m*q (deg m = k-2) and columns the degree-k monomials of character psi.  Exact mod p (python-flint), one process per block.
Usage: macaulay.py quads.json k [workers] [subset_json]   (subset_json: optional list of variable indices to keep)"""
import sys, json, time, os, itertools
from multiprocessing import Pool
from flint import nmod_mat

def chi_add(a, b):
    return tuple((x + y) % m for x, y, m in zip(a, b, (2, 2, 6, 2)))

def monos(vars_, k):
    return list(itertools.combinations_with_replacement(vars_, k))

def block_rank(args):
    psi, rows, cols, p = args
    t0 = time.time()
    cidx = {c: i for i, c in enumerate(cols)}
    nc = len(cols)
    flat = [0] * (len(rows) * nc)
    for r, row in enumerate(rows):
        for key, v in row.items():
            flat[r * nc + cidx[key]] = v
    M = nmod_mat(len(rows), nc, flat, p)
    rk = M.rank()
    return {'psi': psi, 'rows': len(rows), 'cols': nc, 'rank': rk, 'full': rk == nc, 'seconds': round(time.time() - t0, 1)}

def closure(fn, k, workers=8, subset=None, log=print):
    J = json.load(open(fn))
    p = J['p']; chi = [tuple(c) for c in J['chi']]
    quads = [{(a, b): v for a, b, v in cf} for cf in J['quads']]
    n = len(chi)
    vars_ = list(range(n)) if subset is None else sorted(subset)
    vs = set(vars_)
    # restrict quadrics to the subset
    Q = []
    for cf in quads:
        cf2 = {(a, b): v for (a, b), v in cf.items() if a in vs and b in vs}
        if cf2: Q.append(cf2)
    chi_q = [chi_add(chi[a], chi[b]) for cf in Q for (a, b) in [next(iter(cf))]]
    # columns per block
    cols = {}
    for mono in monos(vars_, k):
        c = (0, 0, 0, 0)
        for a in mono: c = chi_add(c, chi[a])
        cols.setdefault(c, []).append(mono)
    rows = {}
    for r, cf in enumerate(Q):
        for m in monos(vars_, k - 2):
            c = chi_q[r]
            for a in m: c = chi_add(c, chi[a])
            row = {}
            for (a, b), v in cf.items():
                key = tuple(sorted(m + (a, b)))
                row[key] = (row.get(key, 0) + v) % p
            rows.setdefault(c, []).append(row)
    tasks = [(str(psi), rows.get(psi, []), cols[psi], p) for psi in cols]
    tasks.sort(key=lambda t: -len(t[2]))
    log(f'macaulay k={k}: {len(vars_)} vars, {len(Q)} quadrics, {len(tasks)} blocks, max block {max(len(t[2]) for t in tasks)} cols x {max(len(t[1]) for t in tasks)} rows')
    out = []
    with Pool(workers) as pool:
        for res in pool.imap_unordered(block_rank, tasks):
            out.append(res); log(json.dumps(res))
    full = all(r['full'] for r in out)
    summary = {'k': k, 'n_vars': len(vars_), 'n_quadrics': len(Q), 'n_blocks': len(out), 'dim_S_k': sum(r['cols'] for r in out), 'dim_I_k': sum(r['rank'] for r in out),
               'dim_quotient_k': sum(r['cols'] - r['rank'] for r in out), 'I_k_equals_S_k': full, 'blocks': out}
    return summary

if __name__ == '__main__':
    fn, k = sys.argv[1], int(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    subset = json.load(open(sys.argv[4])) if len(sys.argv) > 4 else None
    s = closure(fn, k, workers, subset, log=lambda m: print(m, file=sys.stderr, flush=True))
    print(json.dumps({k: v for k, v in s.items() if k != 'blocks'}))
