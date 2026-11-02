"""Job: exact (mod p) sweep of structured candidate denominators with the reduced Koszul E2 system.
usage: python3 w2_sweep.py candidates.json out_tag [--dmax D] [--e1only-above D]"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, time, argparse, traceback
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import w2lib as L
from hodge import core
ap = argparse.ArgumentParser()
ap.add_argument('cands'); ap.add_argument('tag'); ap.add_argument('--dmax', type=int, default=10)
ap.add_argument('--e1only-above', type=int, default=10); ap.add_argument('--only', default=None)
args = ap.parse_args()
out = f'w2-sweep-{args.tag}.jsonl'
C = json.load(open(args.cands))
done = set()
try:
    for line in open(out):
        r = json.loads(line); done.add(r.get('label'))
except FileNotFoundError:
    pass
for fam, label, expr, frame in C:
    if label in done: continue
    if args.only and args.only not in label: continue
    t0 = time.time()
    rec = {'experiment': 'structured-sweep', 'family': fam, 'label': label, 'frame': frame}
    try:
        D = core.parse(expr)
        if frame == 'uv':
            D = L.uv_to_xy(D)
        d = core.degree(D)
        rec['d'] = d
        if d > args.dmax:
            rec['note'] = f'skipped (d > {args.dmax})'
        else:
            r = L.exact_fixed_D(D, p=L.P0, verbose=True, lift=True, e1_only=(d > args.e1only_above))
            rec.update(r)
            rec['D_xy'] = core.to_str(D)[:3000]
    except Exception as ex:
        rec['error'] = f'{ex.__class__.__name__}: {ex}'; rec['traceback'] = traceback.format_exc()[-1500:]
    rec['seconds'] = round(time.time() - t0, 1)
    L.log_line(out, rec)
    print(json.dumps({k: rec.get(k) for k in ('family', 'label', 'd', 'e1', 'consistent', 'rank', 'nullity', 'note', 'error', 'seconds')}), flush=True)
    if rec.get('witness'):
        json.dump(rec['witness'], open('witness-candidate.json', 'w'), indent=1)
        print('WITNESS FOUND', label, flush=True)
print('sweep done', flush=True)
