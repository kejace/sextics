"""E1-kernel re-verification: ker([A]: (S/J)_d -> (S/J)_{d+6}) for both Y-parities and the Klein characters, one prime.
Usage: python3 e1ker.py <prime> <dmin> <dmax> [parities e.g. 01]
Writes kernel bases (mod p) to wave2-excl/e1ker_p<prime>_d<d>_par<par>.json and logs one JSON line per (d, parity)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
from common import *
p = int(sys.argv[1]); dmin, dmax = int(sys.argv[2]), int(sys.argv[3])
pars = [int(c) for c in (sys.argv[4] if len(sys.argv) > 4 else '01')]
F, A, q, rec = load()
Fp = reduce_mod(F, p)
for d in range(dmin, dmax + 1):
    for par in pars:
        t0 = time.time()
        out = jacobian.mult_kernel(F, A, d, p, parity=par, verbose=False)
        parT = (par + core.parity_of(A)) % 2
        # smoothness / Hilbert-function check at the target degree d+6 and at d (ranks of J must equal the smooth values)
        Jt, _ = jacobian.jacobian_matrix(F, d + 6, p, parT)
        rankT = Jt.rank() if Jt.nrows() else 0
        nT = len(monomials(d + 6, parT))
        expT = expected_ring_dim(d + 6, parT); expD = expected_ring_dim(d, par)
        okT = (nT - rankT == expT); okD = (out['ring_dim'] == expD)
        # Klein character split of the kernel basis (each basis element is reduced; split by monomial character)
        chars = {}
        for b in out['basis']:
            cs = set(klein(e) for e in b)
            key = str(sorted(cs))
            chars[key] = chars.get(key, 0) + 1
        rec_ = {'exp': 'e1-kernel', 'd': d, 'parity': par, 'p': p, 'ker_dim_modJ': out['dim'],
                'ker_dim_raw': out['kernel_dim'], 'ring_dim_d': out['ring_dim'], 'jac_dim_d': out['jacobian_dim'],
                'ring_dim_target': nT - rankT, 'hilbert_ok_target': okT, 'hilbert_ok_d': okD,
                'basis_chars': chars, 'seconds': round(time.time() - t0, 1),
                'rigor': 'modular kernel dim is an upper bound for the characteristic-0 kernel dim when hilbert_ok_target (rank_p J = rank_Q J)'}
        log(rec_)
        json.dump({'p': p, 'd': d, 'parity': par, 'basis': [core.to_str(b) for b in out['basis']]},
                  open(f'{PK}/e1ker_p{p}_d{d}_par{par}.json', 'w'), indent=1)
