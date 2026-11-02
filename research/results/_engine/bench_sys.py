"""Structure + timing of real fixed-D systems on open instances (engine v2 dense/compressed path)."""
import sys, time, json, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
import numpy as np
random.seed(1)
cases = [('dense-rational', 8), ('eight-plane-configuration', 8), ('cycle-q2', 6)]
for iid, d in cases:
    F, A, q, rec = core.load_instance(iid)
    # random Y-even denominator of degree d, sparse (12 terms)
    mons = core.monomials(d, 0)
    D = {m: core.Fraction(random.randint(-3, 3) or 1) for m in random.sample(mons, 12)}
    p = core.PRIMES[0]
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, D, p)
    tb = time.time() - t0
    colw = np.array([len(c) for c in S.colpolys])
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    roww = np.diff(M.indptr)
    out = {'instance': iid, 'q': q, 'd': d, 'rows': S.nrows, 'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ,
           'nnz': int(colw.sum()), 'colw_min': int(colw.min()), 'colw_med': float(np.median(colw)), 'colw_max': int(colw.max()),
           'roww_min': int(roww.min()), 'roww_med': float(np.median(roww)), 'roww_max': int(roww.max()),
           'rows_w1': int((roww == 1).sum()), 'rows_w0': int((roww == 0).sum()), 'F_terms': len(F), 't_build': round(tb, 1)}
    print(json.dumps(out), flush=True)
    t0 = time.time()
    r = fixed_d.solve_compressed(S)
    out2 = {'instance': iid, 'd': d, 'solver': 'compressed', 'consistent': r['consistent'], 'rank': r['rank'],
            'nullity': r['nullity'], 'k': r['k'], 't_compress': r['t_compress'], 't_rref': r['t_rref'], 'total': round(time.time() - t0, 1)}
    print(json.dumps(out2), flush=True)
