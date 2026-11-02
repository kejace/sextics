"""Job: restricted E1 kernels of [A] on structured families (vert / coneXM / coneLY), degrees given on argv.
Usage: python3 w2c_vert.py <out.jsonl> d1 d2 ...   (both parities each)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2c_lib import *
import flint
try: flint.ctx.threads = int(os.environ.get('FLINT_THREADS', '8'))
except Exception: pass

out_path = sys.argv[1]
degs = [int(x) for x in sys.argv[2:]]
F, A, q, rec = core.load_instance(IID)
p = core.PRIMES[0]
for d in degs:
    for par in (0, 1):
        for fname, fn in (('vert', fam_vert), ('coneXM', fam_coneXM), ('coneLY', fam_coneLY)):
            fam = fn(d, par)
            if not fam: continue
            t0 = time.time()
            vecs, info = restricted_kernel(F, A, fam, d, p)
            an = analyze_kernel(F, vecs, fam, d, p)
            rec_out = {'exp': 'restricted-E1-kernel', 'family': fname, 'd': d, 'parD': par, 'prime': p,
                       'fam_dim': info['fam_dim'], 'rows': info['rows'], 'cols': info['cols'],
                       'kernel_dim': info['kernel_dim'], 'kernel': an, 'seconds': round(time.time() - t0, 1),
                       'job': 'construct-spp-vert', 'result': ('kernel=%d' % info['kernel_dim'])}
            with open(out_path, 'a') as f: f.write(json.dumps(rec_out) + '\n')
            print('DONE', fname, d, par, 'kernel', info['kernel_dim'], [ (a['F_divides'], a['in_J']) for a in an ], flush=True)
print('ALL DONE', flush=True)
