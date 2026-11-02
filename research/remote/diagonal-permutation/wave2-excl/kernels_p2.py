"""Independent re-verification of the E1 kernel dims (attack1 §2.1) with the ENGINE implementation (jacobian.mult_kernel)
at a DIFFERENT prime (default core.PRIMES[1] = 2147483629; wave 1 used 2147483647 and its own e1kernel.py).
usage: python3 kernels_p2.py [p] [dmin] [dmax]"""
import sys, time
from common import *
from hodge import jacobian
F, A, q, rec = core.load_instance(IID)
p = int(sys.argv[1]) if len(sys.argv) > 1 else core.PRIMES[1]
dmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 9
for d in range(dmin, dmax + 1):
    for par in (0, 1):
        t0 = time.time()
        K = jacobian.mult_kernel(F, A, d, p=p, parity=par, verbose=False)
        e = {'exp': 'e1-kernel-mult_kernel', 'd': d, 'parity': par, 'p': p, 'primitive_kernel_dim': K['dim'],
             'kernel_dim': K['kernel_dim'], 'jacobian_dim': K['jacobian_dim'], 'ring_dim': K['ring_dim'],
             'seconds': round(time.time() - t0, 1)}
        log(e)
        if K['dim']:
            receipt(f'e1kernel_d{d}_par{par}_p{p}.json', {'d': d, 'parity': par, 'p': p,
                    'basis_modp': [{','.join(map(str, m)): int(v) for m, v in b.items()} for b in K['basis']]})
