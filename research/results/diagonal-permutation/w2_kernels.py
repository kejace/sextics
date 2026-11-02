"""Wave 2 job (remote): E1 kernel (mult_kernel) at the degrees given on the command line, both Y-parities; bases saved mod p.
usage: python3 w2_kernels.py 10 11 12   -> appends to w2-kernels.jsonl, writes w2-kernel-d<d>-par<par>.json"""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
LANE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(LANE, 'w2-kernels.jsonl')
def log(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(OUT, 'a') as f: f.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:600], flush=True)
F, A, q, rec = core.load_instance('diagonal-permutation')
for d in [int(x) for x in sys.argv[1:]]:
    for par in (0, 1):
        t0 = time.time()
        K = jacobian.mult_kernel(F, A, d, parity=par, verbose=True)
        log({'exp': 'mult_kernel', 'd': d, 'parity': par, 'dim_mod_J': K['dim'], 'kernel_dim': K['kernel_dim'], 'jacobian_dim': K['jacobian_dim'],
             'ring_dim': K['ring_dim'], 'rows': K['rows'], 'cols': K['cols'], 'seconds': round(time.time() - t0, 1)})
        if K['dim']:
            json.dump({'d': d, 'parity': par, 'p': K['p'], 'basis_modp': [{','.join(map(str, e)): v for e, v in b.items()} for b in K['basis']]},
                      open(os.path.join(LANE, f'w2-kernel-d{d}-par{par}.json'), 'w'))
print('KERNEL JOB DONE', flush=True)
