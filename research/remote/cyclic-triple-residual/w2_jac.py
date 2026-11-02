"""Job: Fermat calibration + Jacobian-type degrees d=5,6,7 (parity-restricted only = all G-sectors at once):
mod-J relaxation and full-E2 Koszul varpro over the w-space (D = w.grad F).
usage: python3 w2_jac.py [procs] [scale]"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import time, json
from w2core import *
from w2run import *
procs = int(sys.argv[1]) if len(sys.argv) > 1 else 8
scale_ = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
JOB = 'w2jac'
F, A, q = load()
G = full_group(F, A)
t00 = time.time()
# ---- Fermat calibration (known witness D0 = X0^6+X1^6+X2^6 in sector (0,0), d=6)
F0 = fermat()
S0 = Sect(F0, A, G)
r0 = e1_space(S0, F0, A, q, 6, (0, 0))
fam0 = Family(S0, F0, A, q, 6, (0, 0), r0['basis'], verbose=True)
for relax in (True, False):
    e, cands, prom = search_batch(fam0, 'fermat-cal-d6-c00', F0, A, q, JOB, n_restarts=int(400 * scale_), procs=procs,
                                  relax=relax, promote_thresh=1e-8, best_path=f'{WD}/w2-best-fermatcal-{"relax" if relax else "full"}.json',
                                  extra={'calibration': True, 'known_witness': 'X0^6+X1^6+X2^6 (+ lambda F0 gauge)'})
    if relax:
        refine_relax_to_full(fam0, cands, 'fermat-cal-d6-c00', F0, A, q, JOB, top=10)
print(f'calibration done {time.time()-t00:.0f}s', flush=True)
# ---- lane F, parity group, d = 5, 6, 7
Sp = Sect(F, A, parity_group())
plan = {5: (int(2000 * scale_), int(2000 * scale_), None), 6: (int(3000 * scale_), int(1500 * scale_), 3600),
        7: (int(3000 * scale_), int(160 * scale_), 5400)}
for d in (5, 6, 7):
    for c in sorted(Sp.sectors(d)):
        r = e1_space(Sp, F, A, q, d, c, verbose=True)
        if r['kerA_in_S'] == 0:
            w2log({'batch': 'e1space', 'job': JOB, 'd': d, 'group': 'parity', 'sector': list(c), 'note': 'no E1-admissible D', **{k: r[k] for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J')}})
            continue
        fam = Family(Sp, F, A, q, d, c, r['basis'], verbose=True)
        n_relax, n_full, tl = plan[d]
        label = f'jac-d{d}-par{csec(c)}'
        e, cands, prom = search_batch(fam, label, F, A, q, JOB, n_restarts=n_relax, procs=procs, relax=True,
                                      best_path=f'{WD}/w2-best-{label}-relax.json', extra={'group': 'parity', 'e1space': {k: r[k] for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J')}})
        refine_relax_to_full(fam, cands, label, F, A, q, JOB, top=8)
        e, cands, prom = search_batch(fam, label, F, A, q, JOB, n_restarts=n_full, procs=procs, relax=False, time_limit=tl,
                                      best_path=f'{WD}/w2-best-{label}-full.json', extra={'group': 'parity'})
        print(f'--- {label} done, elapsed {time.time()-t00:.0f}s', flush=True)
print('done', flush=True)
