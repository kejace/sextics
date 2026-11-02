"""Rational reconstruction of the E1-kernel classes (canonical normal forms w.r.t. rref of J_d^chi) from several primes.
Usage: python3 e1recon.py <d> <chi ab> <nprimes>   -> e1reps_d<d>_chi<ab>_Qk.json (list, None where reconstruction failed)"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
d = int(sys.argv[1]); c = tuple(int(x) for x in sys.argv[2]); k = int(sys.argv[3])
primes = L.PRIMES[2:2 + k]
per_prime = []
dims = []
for p in primes:
    t0 = time.time()
    r = L.e1_kernel(d, c, p, verbose=False)
    per_prime.append((p, r['reps'])); dims.append(r['kernel_mod_J'])
    print(f'p={p}: kernel mod J = {r["kernel_mod_J"]} ({time.time()-t0:.0f}s)', flush=True)
assert len(set(dims)) == 1, dims
n = dims[0]
out = []; nok = 0
for idx in range(n):
    ok, Pq = L.ratrecon_polys([(p, reps[idx]) for p, reps in per_prime])
    out.append(core.to_str({e: v for e, v in Pq.items() if v is not None}) if ok else None)
    nok += ok
fn = f'{L.OUT}/e1reps_d{d}_chi{c[0]}{c[1]}_Qk.json'
json.dump(out, open(fn, 'w'))
L.receipt('e1recon', exp='e1-kernel-ratrecon-k', d=d, chi=c, primes=primes, n_classes=n, reconstructed=nok, file=fn)
print('E1RECON DONE', nok, 'of', n, flush=True)
