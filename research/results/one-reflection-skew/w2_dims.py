"""Wave-2 task (1): E1-admissible denominator spaces per degree and character sector.

5-variable reduction (exact, lane-relevant): for d = 4..12 and each (rho, sigma) sector of D, the dimension of
{D : A D in J_f}, of J_f,d, of the f-multiples, and the number of non-Jacobian directions; the float kernel is
cross-checked against the mod-p dimensions inside w2lib.e1_kernel.  Also the 6-variable mult_kernel of the engine
for Y-even D at d = 4..dmax6 (cross-check with the dossier table).  Writes w2_dims.json and log lines.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import w2lib as W
from w2lib import core

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=4)
ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--dmax6', type=int, default=9)
ap.add_argument('--inst', default='one-reflection-skew')
a = ap.parse_args()
prob = W.Problem.from_instance(a.inst)
out = {'inst': a.inst, 'five_var': [], 'six_var': [], 'ts': W.now()}
path = f'{W.WD}/w2_dims.json'
for d in range(a.dmin, a.dmax + 1):
    for pD in (0, 1):
        for sD in ((0, 1) if prob.use_sigma else (None,)):
            t0 = time.time()
            try:
                Kb, info = W.e1_kernel(prob, d, pD, sD)
                # Jacobian-type? (kernel == J_d)   and the target-side numbers
                info['jacobian_type'] = (info['ker_mod_J'] == 0)
                info['float_kernel_dim'] = int(Kb.shape[1])
            except Exception as ex:
                info = {'d': d, 'pD': pD, 'sD': sD, 'error': f'{ex.__class__.__name__}: {ex}'}
            info['secs'] = round(time.time() - t0, 1)
            out['five_var'].append(info)
            print(json.dumps(info), flush=True)
            W.log_event({'event': 'e1-dims-5var', **info})
            with open(path, 'w') as fh:
                json.dump(out, fh, indent=1)
# six-variable cross-check with the engine's mult_kernel (Y-even D), modest degrees only
from hodge import jacobian
F, A, q, rec = core.load_instance(a.inst)
for d in range(a.dmin, a.dmax6 + 1):
    t0 = time.time()
    try:
        r = jacobian.mult_kernel(F, A, d, parity=0)
        info = {'d': d, 'parity': 0, 'kernel_dim': r['kernel_dim'], 'ker_mod_J': r['dim'], 'jacobian_dim': r['jacobian_dim'],
                'ring_dim': r['ring_dim'], 'secs': round(time.time() - t0, 1)}
    except Exception as ex:
        info = {'d': d, 'parity': 0, 'error': f'{ex.__class__.__name__}: {ex}', 'secs': round(time.time() - t0, 1)}
    out['six_var'].append(info)
    print(json.dumps(info), flush=True)
    W.log_event({'event': 'e1-dims-6var', **info})
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1)
print('DONE')
