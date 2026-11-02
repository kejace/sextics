"""Precompute (and cache) the exact E1-kernel bases for the search families: all degrees/sectors/targets."""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import w2lib as W
from w2lib import core
inst = sys.argv[1] if len(sys.argv) > 1 else 'one-reflection-skew'
degs = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [7, 8, 9, 10, 11, 12]
prob = W.Problem.from_instance(inst)
targets = [('A', prob.A)]
if inst == 'one-reflection-skew':
    targets.append(('Aplus', core.parse('X0**3*(Y0**3+Y1**3)/2')))
for d in degs:
    for (pD, sD) in ((0, 0), (0, 1), (1, 0), (1, 1)):
        for tn, T in targets:
            if tn == 'Aplus' and (pD, sD) not in ((0, 0), (0, 1)):
                continue
            t0 = time.time()
            try:
                B, info = W.e1_kernel_exact(prob, d, pD, sD, T=T)
                rec = {'event': 'exact-kernel', 'inst': inst, 'd': d, 'rho': pD, 'sigma': sD, 'target': tn, 'kerD': info['kerD'],
                       'primes_used': info.get('primes_used'), 'max_height_bits': info.get('max_height_bits'), 'cached': info.get('cached', False),
                       'secs': round(time.time() - t0, 1)}
            except Exception as ex:
                rec = {'event': 'exact-kernel', 'inst': inst, 'd': d, 'rho': pD, 'sigma': sD, 'target': tn, 'error': f'{ex.__class__.__name__}: {ex}',
                       'secs': round(time.time() - t0, 1)}
            W.log_event(rec)
            print(json.dumps(rec), flush=True)
print('DONE', flush=True)
