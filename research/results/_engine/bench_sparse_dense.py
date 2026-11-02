"""sparse elimination on a dense-F instance (dense-rational d=8, 15M nnz): does it blow up?"""
import sys, time, json, random, os
os.environ['HODGE_QUIET'] = '1'
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, sge
random.seed(1)
for iid, d, nterms in [('cyclic-triple-residual', 8, 12), ('klein-coupled-character', 8, 12), ('dense-rational', 8, 12)]:
    F, A, q, rec = core.load_instance(iid)
    mons = core.monomials(d, 0)
    D = {m: core.Fraction(random.randint(-3, 3) or 1) for m in random.sample(mons, nterms)}
    p = core.PRIMES[0]
    t0 = time.time(); S = fixed_d.FixedDSystem(F, A, q, D, p); tb = time.time() - t0
    M = fixed_d._sparse_matrix(S)
    for max_mult in (3, 30):
        t0 = time.time()
        r = sge.solve_sparse_system(M, p, stop_cols=400, max_nnz=max_mult * M.nnz, need_rank=True)
        print(json.dumps({'instance': iid, 'd': d, 'rows': S.nrows, 'cols': S.ncols, 'nnz': int(M.nnz), 't_build': round(tb, 1), 'max_nnz_mult': max_mult,
                          **{k: r.get(k) for k in ('consistent', 'rank', 'nullity', 'status_sparse', 'npiv_sparse', 'active_rows', 'active_cols', 'active_nnz', 'fill_events', 't_sparse', 't_dense', 'k_dense', 't_total', 'verified')}}), flush=True)
        if r.get('status_sparse') in (0, 1):
            break
