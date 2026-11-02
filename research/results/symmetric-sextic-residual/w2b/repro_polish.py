"""Reproduce one remote restart (same family, same seed) locally and polish it with Gauss-Newton."""
import os, sys, time, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
chi, d, pD, job = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
MIR = '/tmp/hodge/remote/symmetric-sextic-residual/w2b'
recs = []
for fn in glob.glob(f'{MIR}/{job}_{chi}-{d}-{pD}_w*.jsonl'):
    for l in open(fn):
        try: recs.append(json.loads(l))
        except Exception: pass
recs.sort(key=lambda r: r['rel']); best = recs[0]
print('remote best:', {k: best.get(k) for k in ('seed', 'rel', 'iters', 'nfev', 'sparse_init')}, flush=True)
F, A, q = load()
t0 = time.time(); fam, info = build_sector_family(F, A, q, chi, d, pD, label=f'{chi}/d{d}/p{pD}')
print('family built', {k: info[k] for k in ('dimK', 'dimJ', 'ker_beyond_J', 'dimF', 'dimP', 'rows', 'nN', 'nQ')}, f'{time.time()-t0:.0f}s', flush=True)
if best.get('sparse_init'):
    print('NOTE: remote restart used a sparse init (worker rng) - cannot reproduce exactly; using plain seed init')
t0 = time.time(); c = fam.search(n_restarts=1, seed=best['seed'], maxiter=400)[0]
print('local re-run: rel', f"{c['rel']:.3e}", 'iters', c['iters'], f'{time.time()-t0:.0f}s', flush=True)
t0 = time.time(); a, x, rel, hist = polish(fam, c['a'], iters=10)
print('polish history:', [f'{h:.3e}' for h in hist], f'{time.time()-t0:.0f}s', flush=True)
# continue L-BFGS from the polished point with more iterations
c2 = fam.search(n_restarts=1, seed=0, maxiter=1500, inits=(a,))[0]
print('L-BFGS 1500 iters from polished point: rel', f"{c2['rel']:.3e}", 'iters', c2['iters'], flush=True)
a2, x2, rel2, hist2 = polish(fam, c2['a'], iters=10)
print('second polish:', [f'{h:.3e}' for h in hist2], flush=True)
j = int(np.argmax(np.abs(a2))); an = a2 / a2[j]
print('dominant coordinates:', [(fam.reps[k], round(float(an[k]), 5)) for k in np.argsort(-np.abs(an))[:10]])
wlog({'job': 'local-polish', 'family': 'polish-check', 'chi': chi, 'd': d, 'pD': pD, 'seed': best['seed'], 'remote_rel': best['rel'],
      'local_rerun_rel': c['rel'], 'polished_rel': rel, 'lbfgs1500_rel': c2['rel'], 'polished2_rel': rel2, 'dimP': fam.dimP})
