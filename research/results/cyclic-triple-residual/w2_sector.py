"""Job: sector-restricted joint search at degree d over the E1-admissible space (kernel + J_d) of every G-sector.
usage: python3 w2_sector.py d [n_full] [n_relax] [seed0] [time_limit_per_sector_s] [procs] [only_sector e.g. 04]"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import time, json
from w2core import *
from w2run import *
d = int(sys.argv[1])
n_full = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
n_relax = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
seed0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
tl = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != 'none' else None
procs = int(sys.argv[6]) if len(sys.argv) > 6 else 8
only = sys.argv[7] if len(sys.argv) > 7 else None
JOB = f'w2sector-d{d}-s{seed0}'
F, A, q = load()
G = full_group(F, A)
S = Sect(F, A, G)
t00 = time.time()
secs = sorted(S.sectors(d))
if only:
    secs = [c for c in secs if csec(c) == only]
for c in secs:
    path = f'{WD}/w2-e1space-d{d}-c{csec(c)}.json'
    if os.path.exists(path):
        r = load_space(path)
    else:
        r = e1_space(S, F, A, q, d, c, verbose=True)
        save_space(r, path)
    if r['kerA_in_S'] == 0 or not r['basis']:
        print(f'sector {c}: no E1-admissible D', flush=True)
        continue
    fam = Family(S, F, A, q, d, c, r['basis'], verbose=True)
    label = f'sec-d{d}-c{csec(c)}'
    extra = {'group': 'G18', 'e1space': {k: r[k] for k in ('dim_S', 'dim_J', 'kerA_in_S', 'ker_mod_J', 'dim_R_d', 'dim_R_t', 'forced') if k in r}}
    if n_relax:
        e, cands, prom = search_batch(fam, label, F, A, q, JOB, n_restarts=n_relax, procs=procs, relax=True, seed0=seed0,
                                      time_limit=tl, best_path=f'{WD}/w2-best-{label}-relax-s{seed0}.json', extra=extra)
        refine_relax_to_full(fam, cands, label, F, A, q, JOB, top=6)
    if n_full:
        e, cands, prom = search_batch(fam, label, F, A, q, JOB, n_restarts=n_full, procs=procs, relax=False, seed0=seed0,
                                      time_limit=tl, best_path=f'{WD}/w2-best-{label}-full-s{seed0}.json', extra=extra)
    print(f'--- {label} done, elapsed {time.time()-t00:.0f}s', flush=True)
print('done', flush=True)
