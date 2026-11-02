import sys, time, json, random, os
os.environ['HODGE_QUIET'] = '1'
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, sge
random.seed(1)
iid, d = 'cyclic-triple-residual', 8
F, A, q, rec = core.load_instance(iid)
mons = core.monomials(d, 0)
D = {m: core.Fraction(random.randint(-3, 3) or 1) for m in random.sample(mons, 12)}
p = core.PRIMES[0]
S = fixed_d.FixedDSystem(F, A, q, D, p)
M = fixed_d._sparse_matrix(S)
for stop in (400, 100, 0):
    t0 = time.time()
    r = sge.solve_sparse_system(M, p, stop_cols=stop, need_rank=True)
    print(json.dumps({'instance': iid, 'd': d, 'stop_cols': stop, **{k: r.get(k) for k in ('consistent', 'rank', 'nullity', 'status_sparse', 'npiv_sparse', 'active_rows', 'active_cols', 'active_nnz', 'fill_events', 't_sparse', 't_dense', 'k_dense', 't_total')}}), flush=True)
# accepted witness at d=6 (consistent) with stop 0: exact lift from pivots
F, A, q, rec = core.load_instance('hesse-triangle'); d, D, N, Q = core.load_witness('hesse-triangle')
S = fixed_d.FixedDSystem(F, A, q, D, p); M = fixed_d._sparse_matrix(S)
for stop in (400, 0):
    t0 = time.time(); r = sge.solve_sparse_system(M, p, stop_cols=stop)
    sol = fixed_d.exact_solution_from_pivots(F, A, q, D, r['pivots'])
    print(json.dumps({'instance': 'hesse-triangle', 'stop_cols': stop, 'consistent': r['consistent'], 'rank': r['rank'], 'verified': r.get('verified'), 'exact': sol is not None, 't_total': r['t_total'], 'npiv_sparse': r['npiv_sparse'], 'active_nnz': r['active_nnz']}), flush=True)
