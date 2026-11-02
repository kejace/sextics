"""(a) re-verification 3b: XY-frame cross-check of the E1 table with the engine's jacobian.mult_kernel at two
new primes (independent code path: dense Jacobian generators in the original coordinates)."""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '.')
from hodge import core, jacobian
from logutil import log
F, A, q, rec = core.load_instance('polarized-k3-remainder')
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for p in (1000000007, 65521):
    for d in range(0, dmax + 1):
        for par in (0, 1):
            t0 = time.time()
            r = jacobian.mult_kernel(F, A, d, p, parity=par)
            rec_ = dict(exp='a3b_e1_xy', prime=p, d=d, parity=par, ring_dim=r['ring_dim'], jacobian_dim=r['jacobian_dim'],
                        kernel_mod_J=r['dim'], kernel_dim_raw=r['kernel_dim'], seconds=round(time.time() - t0, 1),
                        rigor='mod p; full rank mod p is a lower bound for the rank over Q')
            print(rec_, flush=True); log(**rec_)
