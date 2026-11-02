"""Exact rational lift of the primitive E1 kernel at degree d (parity par): run jacobian.mult_kernel at several primes,
CRT the canonical rref representatives (standard-monomial basis is canonical for a fixed monomial order) and rationally
reconstruct; verify A*D in J(F) exactly mod two further primes.  Output: kernel-d<d>-par<par>-exact.json {'D': [str,...]}.
usage: python3 kernel_lift.py <d> <par> [nprimes=4]"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from fractions import Fraction
LANE = '/tmp/hodge/results/diagonal-permutation'
d = int(sys.argv[1]); par = int(sys.argv[2]); npr = int(sys.argv[3]) if len(sys.argv) > 3 else 4
F, A, q, rec = core.load_instance('diagonal-permutation')
def log(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); e['lens'] = 'search'
    with open(f'{LANE}/wave2-search.log.jsonl', 'a') as f: f.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:400], flush=True)
bases = []
for p in core.PRIMES[:npr]:
    t0 = time.time()
    K = jacobian.mult_kernel(F, A, d, p=p, parity=par, verbose=True)
    log({'exp': 'mult_kernel', 'd': d, 'parity': par, 'p': p, 'dim_mod_J': K['dim'], 'kernel_dim': K['kernel_dim'],
         'jacobian_dim': K['jacobian_dim'], 'ring_dim': K['ring_dim'], 'seconds': round(time.time() - t0, 1)})
    bases.append((p, K['basis']))
    json.dump({'d': d, 'parity': par, 'p': p, 'basis_modp': [{','.join(map(str, e)): v for e, v in b.items()} for b in K['basis']]},
              open(f'{LANE}/w2search/kernel-d{d}-par{par}-p{p}.json', 'w'))
dims = {len(b) for _, b in bases}
assert len(dims) == 1, dims
n = dims.pop()
# CRT per basis element / monomial
out = []
ok_all = True
for i in range(n):
    supp = set()
    for _, b in bases: supp |= set(b[i].keys())
    D = {}
    for e in supp:
        a, m = 0, 1
        for p, b in bases:
            v = b[i].get(e, 0)
            a, m = core.crt_pair(a, m, v, p), m * p
        fr = core.ratrecon(a, m)
        if fr is None:
            ok_all = False; fr = None
            break
        if fr != 0: D[e] = fr
    if fr is None:
        log({'exp': 'kernel-lift', 'd': d, 'parity': par, 'element': i, 'status': 'ratrecon failed with %d primes' % npr})
        continue
    # verify at two fresh primes
    good = True
    for p in core.PRIMES[npr:npr + 2]:
        okj, info = jacobian.in_jacobian(F, core.mul(A, D), p)
        good &= bool(okj)
    out.append(core.to_str(D))
    log({'exp': 'kernel-lift', 'd': d, 'parity': par, 'element': i, 'terms': len(D), 'maxden': max(abs(v.denominator) for v in D.values()),
         'maxnum': max(abs(v.numerator) for v in D.values()), 'verified_2_more_primes': good})
    ok_all &= good
json.dump({'d': d, 'parity': par, 'n': len(out), 'all_ok': ok_all, 'D': out}, open(f'{LANE}/w2search/kernel-d{d}-par{par}-exact.json', 'w'))
log({'exp': 'kernel-lift-done', 'd': d, 'parity': par, 'n': len(out), 'all_ok': ok_all})
