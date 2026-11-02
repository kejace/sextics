"""Job: E1-admissible denominator spaces per degree d=4..12 and parity 0/1 (mult_kernel + J(F)_d dims)."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from hodge import core, jacobian
import w2lib as L
F, A, q = L.load()
out_path = 'w2-kernel-table.jsonl'
dmin = int(sys.argv[1]) if len(sys.argv) > 1 else 4
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
for d in range(dmin, dmax + 1):
    for par in (0, 1):
        t0 = time.time()
        res = jacobian.mult_kernel(F, A, d, p=L.P0, parity=par, verbose=True)
        rec = {'experiment': 'e1-kernel-table', 'd': d, 'parity': par, 'kernel_dim_mod_J': res['dim'],
               'kernel_dim_in_R': res['kernel_dim'], 'J_d_dim': res['jacobian_dim'], 'RJ_d_dim': res['ring_dim'],
               'R_d_dim': res['jacobian_dim'] + res['ring_dim'], 'admissible_dim': res['kernel_dim'],
               'rows': res['rows'], 'cols': res['cols'], 'p': res['p'], 'seconds': round(time.time() - t0, 1)}
        L.log_line(out_path, rec)
        print(json.dumps(rec), flush=True)
        if res['dim'] and d >= 10:
            # save kernel representatives (standard-monomial normal forms, mod p) for later structured probes
            with open(f'w2-kernel-basis-d{d}-p{par}.json', 'w') as f:
                json.dump({'d': d, 'parity': par, 'p': res['p'], 'basis': [{core.to_str({e: 1}): v for e, v in b.items()} for b in res['basis']]}, f)
print('done', flush=True)
