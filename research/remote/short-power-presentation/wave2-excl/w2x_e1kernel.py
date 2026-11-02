"""Independent re-verification of the E1 kernel K_d^par / J_d^par = {D in S_d^par : A D in J} / J_d^par,
using hodge.jacobian.mult_kernel (engine code path, dense nullspace) with prime core.PRIMES[1],
plus the parity-graded Hilbert-function check that makes the mod-p answer an upper bound for the
characteristic-zero dimension (and hence a proof when it is 0).
Usage: python3 w2x_e1kernel.py 1e 1o 2e ... (d + parity letter); writes e1kernel_d{d}_p{par}.json."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave2-excl')
from w2x_lib import *

def run(F, A, d, par, p=P):
    t0 = time.time()
    out = jacobian.mult_kernel(F, A, d, p=p, parity=par, verbose=True)
    t1 = time.time()
    # Hilbert checks: J_d^par and J_{d+6}^{par+1}
    rJd = out['jacobian_dim']; eJd = expected_rank_J(d, par)
    rJ6 = out['cols'] - len(core.monomials(d, par))  # number of Jacobian generator columns, not rank
    rk6, n6 = jac_rank(F, d + 6, (par + 1) % 2, p)
    e6 = expected_rank_J(d + 6, (par + 1) % 2)
    rec = {'exp': 'e1kernel', 'd': d, 'parD': par, 'p': p,
           'dim S_d': len(core.monomials(d, par)), 'dim J_d (mod p)': rJd, 'dim J_d (expected smooth)': eJd,
           'rank J_{d+6} (mod p)': rk6, 'rank J_{d+6} (expected)': e6, 'dim S_{d+6}': n6,
           'kernel incl J': out['kernel_dim'], 'kernel mod J': out['dim'], 'ring_dim': out['ring_dim'],
           'hilbert_ok': (rJd == eJd and rk6 == e6),
           'rows': out['rows'], 'cols': out['cols'], 'seconds': round(time.time() - t0, 1)}
    basis = [core.to_str(b) for b in out['basis']]
    with open(f'{PK}/e1kernel_d{d}_p{par}.json', 'w') as f:
        json.dump({'rec': rec, 'basis_modp': basis}, f, indent=1)
    print(json.dumps(rec), flush=True)
    log(rec)
    return rec, out

if __name__ == '__main__':
    F, A, q, rec = load()
    for j in sys.argv[1:]:
        d = int(j[:-1]); par = 0 if j[-1] == 'e' else 1
        run(F, A, d, par)
