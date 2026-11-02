"""solve_sparse vs solve_compressed on accepted witnesses (consistent; exact lift) and perturbed D (inconsistent)."""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
p = core.PRIMES[0]
for iid in ['collision-h3', 'hesse-triangle', 'xtype-sparse-rational-1m142', 'contact46006', 'three-block-rational-d4-strict']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    S = fixed_d.FixedDSystem(F, A, q, D, p)
    t0 = time.time(); rs = fixed_d.solve_sparse(S); ts = time.time() - t0
    t0 = time.time(); rc = fixed_d.solve_compressed(S); tc = time.time() - t0
    ok = rs['consistent'] == rc['consistent'] and rs['rank'] == rc['rank'] and rs.get('verified')
    # exact lift from sparse pivots
    sol = fixed_d.exact_solution_from_pivots(F, A, q, D, rs['pivots'])
    print(json.dumps({'instance': iid, 'd': d, 'rows': S.nrows, 'cols': S.ncols, 'sparse': {k: rs[k] for k in ('consistent','rank','nullity','npiv_sparse','active_rows','active_cols','fill_events','t_sparse','t_dense','t_total')},
                      'compressed': {k: rc[k] for k in ('consistent','rank','nullity','t_compress','t_rref')}, 'agree': ok, 'exact_from_sparse_pivots': sol is not None,
                      't_sparse_total': round(ts, 2), 't_compressed_total': round(tc, 2)}), flush=True)
    D2 = dict(D); e = next(iter(D2)); D2[e] = D2[e] + 1
    S2 = fixed_d.FixedDSystem(F, A, q, D2, p)
    rs2 = fixed_d.solve_sparse(S2); rc2 = fixed_d.solve_compressed(S2)
    print(json.dumps({'instance': iid, 'perturbed': True, 'sparse_consistent': rs2['consistent'], 'compressed_consistent': rc2['consistent'],
                      'sparse_rank_lb': rs2.get('rank_lower_bound', rs2.get('rank')), 'compressed_rank': rc2['rank'], 't_sparse': rs2['t_total']}), flush=True)
