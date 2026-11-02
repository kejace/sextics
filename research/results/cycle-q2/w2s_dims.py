"""Step 1: E1-admissible denominator spaces per degree and character sector (Y-even), exact bases cached.
Usage: python3 w2s_dims.py --dmin 4 --dmax 12 [--cpus N]
"""
import sys, time, json, argparse, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import w2s_lib as L
from hodge import core

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=4); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--cpus', type=int, default=1); ap.add_argument('--job', default='dims')
args = ap.parse_args()
I = L.Inst()
print('group order', I.n, 'elements', len(I.elements), 'gens', I.gens, 'separate', I.gens_separate, flush=True)

def do(task):
    d, key = task
    try:
        kd = L.load_kernel(I, d, key, parity=0, verbose=False)
        wb = L.wspace_basis(I, d, key)
        return {k: kd[k] for k in ('d', 'key', 'nmons', 'kernel_dim', 'jac_dim', 'nonjac_dim', 'ring_dim', 'primes', 'shape', 'seconds', 'exact', 'agree', 'warning') if k in kd} | {'wdim': len(wb)}
    except Exception as ex:
        import traceback
        return {'d': d, 'key': list(key), 'error': repr(ex)[:300], 'trace': traceback.format_exc()[-800:]}

tasks = [(d, key) for d in range(args.dmin, args.dmax + 1) for key in I.keys_of_degree(d, 0)]
print(len(tasks), 'tasks', flush=True)
if args.cpus > 1:
    import multiprocessing as mp
    with mp.get_context('fork').Pool(args.cpus) as pool:
        for rec in pool.imap_unordered(do, tasks):
            rec['step'] = 'dims'; rec['job'] = args.job
            L.log_line(rec)
            print(json.dumps(rec), flush=True)
else:
    for t in tasks:
        rec = do(t); rec['step'] = 'dims'; rec['job'] = args.job
        L.log_line(rec); print(json.dumps(rec), flush=True)
print('DONE', flush=True)
