"""Rank/nullity statistics per (d, sector): full sector-restricted fixed-D system at a random E1-admissible D mod p
(rows, cols, rank, nullity = gauge dimension) and the Koszul-form system (rows, cols, generic Koszul gauge).
usage: python3 fpt_ranks.py <dmin> <dmax> [p]  -> wave2-search/ranks.jsonl"""
import os, sys, json, time, random
os.environ['OPENBLAS_NUM_THREADS'] = '2'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
import fpt_msolve as fm
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
p = int(sys.argv[3]) if len(sys.argv) > 3 else 2147483647
OUT = f'{LANE}/wave2-search/ranks.jsonl'
inst = Instance()
Fp = reduce_mod(inst.F, p); Ap = reduce_mod(inst.A, p)
rng = random.Random(5)
for d in range(dmin, dmax + 1):
    for lab in sorted(inst.labels_in_degree(d)):
        t0 = time.time()
        Dm, Bt = fm.modp_family(inst, d, lab, p)
        K = len(Bt)
        rec = {'kind': 'ranks', 'd': d, 'sector': list(lab), 'p': p, 'K0': len(Dm), 'K': K, 'engine': ENGINE}
        if K == 0:
            log_line(OUT, rec); continue
        # random family point mod p
        a = [rng.randrange(1, p) for _ in range(K)]
        D = {}
        for k in range(K):
            for j in range(len(Dm)):
                if Bt[k][j]:
                    D[Dm[j]] = (D.get(Dm[j], 0) + a[k] * Bt[k][j]) % p
        D = {m: v for m, v in D.items() if v}
        try:
            S, subg = symmetry.restricted_system(Fp, Ap, inst.q, D, p, inst.G)
            r = S.solve()
            rec.update({'full_rows': S.nrows, 'full_cols': S.ncols, 'full_rank': r['rank'], 'full_nullity': r['nullity'], 'full_consistent': bool(r['consistent']),
                        'full_cols_unrestricted': S.full_cols, 'subgroup': len(subg)})
        except Exception as ex:
            rec['full_error'] = f'{ex.__class__.__name__}: {ex}'
        try:
            text, info = fm.build_system(inst, d, lab, p, verbose=False)
            rec.update({'koszul_rows': info.get('rows'), 'koszul_cols': (info.get('nC') or 0) + (info.get('nQ') or 0), 'koszul_nC': info.get('nC'), 'koszul_nQ': info.get('nQ'),
                        'koszul_gauge_generic': info.get('koszul_gauge_dim_generic'), 'e1_lift_rows': info.get('e1_rows'), 'nN': info.get('nN')})
        except Exception as ex:
            rec['koszul_error'] = f'{ex.__class__.__name__}: {ex}'
        rec['seconds'] = round(time.time() - t0, 1)
        log_line(OUT, rec)
        print(json.dumps({k: rec[k] for k in rec if k != 'kind'}), flush=True)
print('[ranks] done', flush=True)
