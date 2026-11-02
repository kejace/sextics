"""Exact (mod p) sector dimension table for d = dmin..dmax: dim S_d sector, J_d, (J:A)_d kernel, nontrivial part,
plus the float e1_space cross-check (K, syz) and sizes of the Koszul system.  One JSON line per sector."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2fast as wf, w2lib
from hodge import core
F, A, q = wf.load()
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
out = open('w2_dims.jsonl', 'a')
for d in range(dmin, dmax + 1):
    for parD in (0, 1):
        for chi in range(7):
            t0 = time.time()
            info = w2lib.sector_space(F, A, d, parD, chi, want_rational=False)
            sec = wf.Sector(F, A, d, parD, chi)
            rec = {'d': d, 'parD': parD, 'chi': chi, 'nD': info['nD'], 'J_d': info['J_d'], 'kernel_total': info['kernel_total'],
                   'nontrivial': info['nontrivial'], 'ring_dim_target': info['ring_dim_target'], 'p': info['p'],
                   'nN': sec.nN, 'nP': sec.nP, 'nQ': sec.nQ, 'rows': sec.nrows}
            if d <= 12:
                Bc, Bn, e1 = sec.e1_space()
                rec['K_float'] = e1['K']; rec['syz_float'] = e1.get('syz'); rec['e1_check'] = e1.get('e1_check')
                Fc = sec.f_multiples(); rec['F_mult'] = int(Fc.shape[1])
            rec['seconds'] = round(time.time() - t0, 1)
            out.write(json.dumps(rec) + '\n'); out.flush()
            print(rec, flush=True)
