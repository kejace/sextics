"""Full fixed-D certificate test of the degree-8 E1 class K8 (never completed in wave 1).
(1) A*K8 in J mod new primes; (2) engine fast compressed solver (PRIMES[0]); (3) dense character-restricted solver at a new prime."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
from common import *
from hodge import fixed_d, symmetry
F, A, q, rec = load()
K8 = core.parse(json.load(open(f'{RES}/kernel_d8_X03Y03.json'))[0])
expr = core.to_str(K8)
for p in NEWPRIMES[:2]:
    t0 = time.time()
    ok, info = jacobian.in_jacobian(F, mul(A, K8), p)
    log({'exp': 'k8-e1-check', 'D': 'K8', 'p': p, 'AK8_in_J': ok, 'info': info, 'seconds': round(time.time() - t0, 1)})
which = sys.argv[1] if len(sys.argv) > 1 else 'both'
if which in ('fast', 'both'):
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, K8, verbose=True)
    out2 = {k: v for k, v in out.items() if k not in ('N', 'Q')}
    log({'exp': 'k8-fixedD-fast', 'D': 'K8', 'd': 8, 'p': core.PRIMES[0], 'consistent': out['consistent'],
         'rank': out['rank'], 'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
         'exact': out.get('exact'), 'seconds': round(time.time() - t0, 1), 'rigor': 'modular (compressed rows, k_extra=40)'})
    if out['consistent'] and out.get('exact'):
        okv, infov = core.verify_witness(F, A, q, K8, out['N'], out['Q'])
        log({'exp': 'k8-fixedD-verify', 'verified': okv})
        if okv:
            json.dump({'D': expr, 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])},
                      open(f'{PK}/WITNESS_K8.json', 'w'), indent=1)
if which in ('sym', 'both'):
    t0 = time.time()
    st, data = symmetry.solve_fixed_D_symmetric(F, A, q, K8, primes=[NEWPRIMES[0], NEWPRIMES[1]], verbose=True)
    data2 = {k: v for k, v in data.items() if k not in ('N', 'Q')}
    log({'exp': 'k8-fixedD-symmetric-dense', 'D': 'K8', 'd': 8, 'status': st, 'data': data2,
         'seconds': round(time.time() - t0, 1), 'rigor': 'modular (dense rref, character-restricted system)'})
