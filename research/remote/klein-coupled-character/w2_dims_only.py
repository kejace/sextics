"""Log sector E1 dimensions (1 prime, no representatives) for the given degrees: python3 w2_dims_only.py d1,d2,... [parD:k list]"""
import os, sys, time
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
from kcc2 import Ctx
F, A, q = L.load()
degs = [int(x) for x in sys.argv[1].split(',')]
secs = None
if len(sys.argv) > 2:
    secs = set()
    for item in sys.argv[2].split(','):
        pd, kr = item.split(':'); lo, hi = (kr.split('-') + [kr])[:2] if '-' in kr else (kr, kr)
        for kk in range(int(lo), int(hi) + 1): secs.add((int(pd), kk))
ctx = Ctx(F, A, L.P0)
for d in degs:
    for parD in (0, 1):
        for k in range(7):
            if secs is not None and (parD, k) not in secs: continue
            t = time.time()
            S = L.sector_kernel_full(ctx, d, k, parD)
            nk = len(S['kervecs']); jd = S['jd']
            rec = {'family': 'sector-dims', 'd': d, 'parD': parD, 'k': k, 'nD': len(S['Dm']), 'ring_dim_target': S['nl'], 'kernel_total': nk,
                   'J_d': jd, 'nontrivial': nk - jd, 'reps': 0, 'rep_heights': None, 'lift_ok': None, 'lift_fail': 0, 'seconds': round(time.time() - t, 1),
                   'p': L.P0, 'note': 'E1-admissible sector dimensions mod p (1 prime, no representatives)'}
            L.wlog(rec)
            print(f"d={d} parD={parD} k={k}: nD={rec['nD']} (R/J)_(d+6)={rec['ring_dim_target']} kernel={nk} J_d={jd} nontrivial={nk-jd} [{rec['seconds']}s]", flush=True)
