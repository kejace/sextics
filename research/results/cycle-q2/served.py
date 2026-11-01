"""served_dim(D): dimension of the space of degree-12 targets (in A's character sector under Stab(D))
served by denominator D, i.e. rank([M | T]) - rank(M) where T = columns D^2 * m.  One modular LU."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu, numpy as np
from hodge import core

def served_dim(F, A, q, D, p=modlu.P22, k_extra=40, seed=1, verbose=False):
    t0 = time.time()
    S = eq.EqSystem(F, A, q, D, p=p, with_target=True, verbose=verbose)
    Mc = eq._compressed_int(S, k_extra=k_extra, seed=seed, with_rhs=False)
    rank, piv = modlu.mod_echelon(Mc, p)
    C_NQ = S.ncols - S.nT
    served = sum(1 for j in piv if j >= C_NQ)
    # is A itself served?  A is one of the target columns; check consistency of rhs separately
    res = eq.solve_fast(S if False else eq.EqSystem(F, A, q, D, p=p, verbose=False), seed=seed)
    return {'served': served, 'sector': S.nT, 'rank_NQ': rank - served, 'cols_NQ': C_NQ, 'rows': S.nrows,
            'order': S.order, 'consistent': res['consistent'], 'seconds': round(time.time() - t0, 1)}

def log_served(D, fam, note, r):
    with open('/tmp/hodge/results/cycle-q2/log.jsonl', 'a') as f:
        f.write(json.dumps({'D': core.to_str(D), 'd': core.degree(D), 'family': fam, 'e1': None, 'consistent': r['consistent'],
                            'rank': r['rank_NQ'], 'nullity': r['cols_NQ'] - r['rank_NQ'], 'seconds': r['seconds'],
                            'served_dim': r['served'], 'sector_dim': r['sector'], 'order': r['order'], 'note': note + ' [served-dim test]'}) + '\n')

if __name__ == '__main__':
    F, A, q, rec = core.load_instance('cycle-q2')
    D = core.parse(sys.argv[1])
    r = served_dim(F, A, q, D, verbose=True)
    print(r)
