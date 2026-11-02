"""Full certificate test mod p for solutions of a relaxation (msolve chart output): rebuild the family exactly as
w2_relax.py did (same prime, same gauge), assemble D mod p for each F_p-rational point, and run the compressed
fixed-D consistency test mod p (fixed_d.FixedDSystem + solve_compressed).  Consistent points are candidates for a
rational lift.  Usage: python3 w2_fulltest_modp.py <d> <parspec> <relax_base> [--max-points N]"""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_relax import load_kernel, gauge_indices
from msparse import solutions

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('parspec'); ap.add_argument('base')
    ap.add_argument('--max-points', type=int, default=20); ap.add_argument('--no-kernel', action='store_true')
    a = ap.parse_args()
    meta = json.load(open(a.base + '_meta.json'))
    p = meta['p']; d = a.d
    pars = [0, 1] if a.parspec == 'm' else [int(a.parspec)]
    F, A, q, rec = load()
    reps, nums = [], []
    if not a.no_kernel:
        for par in pars:
            K, src = load_kernel(d, par, p)
            if K:
                N = particular_numerators(F, A, K, p, verbose=False)
                reps += K; nums += N
    fam = Family(F, A, d, kernel_reps=reps, kernel_nums=nums, parities=pars)
    gidx = set(gauge_indices(fam))
    keep = [k for k in range(fam.K) if k not in gidx]
    assert len(keep) == meta['n'], (len(keep), meta['n'])
    info, pts = solutions(a.base + '_chart.ms', a.base + '_chart.msout')
    print('msolve chart:', {k: v for k, v in info.items() if k not in ('elim', 'den', 'params')}, 'points:', len(pts), flush=True)
    out = []
    for pt in pts[:a.max_points]:
        if pt['x'] is None:
            out.append({'t': pt['t'], 'note': pt.get('note')}); continue
        D = {}
        for i, k in enumerate(keep):
            if pt['x'][i]:
                D = add(D, scale(reduce_mod(fam.members[k]['G'], p), pt['x'][i], p), p)
        if not D:
            out.append({'t': pt['t'], 'note': 'D = 0'}); continue
        parD = core.parity_of(D)
        rec_ = {'t': pt['t'], 'parity': parD, 'terms': len(D)}
        if parD is None:
            rec_['note'] = 'mixed parity D: full test not implemented here (use MixedSystem)'; out.append(rec_); continue
        t0 = time.time()
        Dint = {e: int(v) for e, v in D.items()}
        S = fixed_d.FixedDSystem(F, A, q, Dint, p, verbose=False)
        r = fixed_d.solve_compressed(S, verbose=False)
        rec_.update({'consistent': bool(r['consistent']), 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows,
                     'cols': S.ncols, 'seconds': round(time.time() - t0, 1)})
        print('point', rec_, flush=True)
        out.append(rec_)
    logrec = {'kind': 'relaxation-points-fulltest', 'd': d, 'parspec': a.parspec, 'p': p, 'base': a.base,
              'n_points': len(pts), 'results': out}
    w2log(logrec); print(json.dumps(logrec)[:2000], flush=True)

if __name__ == '__main__':
    main()
