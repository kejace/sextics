"""(a) Re-verify the 4-variable core E1 kernels ker(a : R(H)_k -> R(H)_{k+6}) at NEW primes with J-rank certification.
Writes corekernel-verify.json; logs one line per (k,par,p)."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
primes = [PRIMES[1], PRIMES[2]]   # 2147483629, 2147483587 (wave 1 used PRIMES[0] and, for k<=8, PRIMES[1])
res = {}
for k in range(0, kmax + 1):
    for par in (0, 1):
        for p in primes:
            r = core_kernel(k, par, p, want_basis=False)
            res[f'{k},{par},{p}'] = r
            log(dict(exp='a1-corekernel-verify', **r))
json.dump(res, open('corekernel-verify.json', 'w'), indent=1)
