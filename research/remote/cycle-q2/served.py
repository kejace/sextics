"""unserved(D): number of independent degree-12 targets in A's character sector (under Stab(D)) that
are NOT served by denominator D: rank([M | T]) - rank(M), T = columns D^2 * m.  baseline(D) = same for
D = 1 in the same sector (= dim of the nonzero-class part of the sector).  new = baseline - unserved =
number of new (nonzero) classes reached by D.  One modular LU each."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu, numpy as np
from hodge import core

def _unserved(F, A, q, D, p, Hw, k_extra=40, seed=1, verbose=False):
    S = eq.EqSystem(F, A, q, D, p=p, with_target=True, verbose=verbose, Hw=Hw)
    Mc = eq._compressed_int(S, k_extra=k_extra, seed=seed, with_rhs=False)
    rank, piv = modlu.mod_echelon(Mc, p)
    C_NQ = S.ncols - S.nT
    uns = sum(1 for j in piv if j >= C_NQ)
    return uns, S.nT, rank - uns, C_NQ, S

def served_dim(F, A, q, D, p=modlu.P22, k_extra=40, seed=1, verbose=False, with_consistency=True):
    t0 = time.time()
    Gw = eq.symmetry_group(F)
    Hw = eq.stabilizer(eq.stabilizer(Gw, D), A)
    uns, nT, rank_NQ, C_NQ, S = _unserved(F, A, q, D, p, Hw, k_extra, seed, verbose)
    base, nT1, _, _, _ = _unserved(F, A, q, core.parse('1'), p, Hw, k_extra, seed, False)
    assert nT1 == nT
    out = {'unserved': uns, 'baseline': base, 'new': base - uns, 'sector': nT, 'rank_NQ': rank_NQ, 'cols_NQ': C_NQ,
           'rows': S.nrows, 'order': S.order}
    if with_consistency:
        res = eq.solve_fast(eq.EqSystem(F, A, q, D, p=p, verbose=False, Hw=Hw), seed=seed)
        out['consistent'] = res['consistent']; out['rank'] = res['rank']; out['nullity'] = res['nullity']
    out['seconds'] = round(time.time() - t0, 1)
    return out

def log_served(D, fam, note, r, Dstr=None):
    with open('/tmp/hodge/results/cycle-q2/log.jsonl', 'a') as f:
        f.write(json.dumps({'D': Dstr or core.to_str(D), 'd': core.degree(D), 'family': fam, 'e1': None, 'consistent': r.get('consistent'),
                            'rank': r.get('rank', r['rank_NQ']), 'nullity': r.get('nullity', r['cols_NQ'] - r['rank_NQ']), 'seconds': r['seconds'],
                            'unserved': r['unserved'], 'baseline_unserved': r['baseline'], 'new_classes_served': r['new'],
                            'sector_dim': r['sector'], 'order': r['order'], 'engine': 2, 'note': note + ' [served test]'}) + '\n')

if __name__ == '__main__':
    F, A, q, rec = core.load_instance('cycle-q2')
    D = core.parse(sys.argv[1])
    print(served_dim(F, A, q, D, verbose=True))
