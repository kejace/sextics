"""(a) adversarial re-verification of the E1-kernel ladder with a DIFFERENT prime and a DIFFERENT code path
(hodge.jacobian.mult_kernel instead of lane e1kernel.py), plus the J-block ranks needed for the over-Q argument."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2lib import *
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2147483587
ds = [int(a) for a in sys.argv[2:]] or list(range(1, 10))
for d in ds:
    for par in (0, 1):
        t0 = time.time()
        out = jacobian.mult_kernel(F, A, d, p=p, parity=par)
        # rank of the J_{d+6} block in the target parity (needed: must equal CI value for the rigorous over-Q argument)
        tpar = (1 + par) % 2
        Jm, mons = jacobian.jacobian_matrix(F, d + 6, p, tpar)
        rJ_target = Jm.rank()
        ci_target = len(mons) - hilbert_CI_parity(d + 6, tpar)
        ci_d = len(monomials(d, par)) - hilbert_CI_parity(d, par)
        rec = {'exp': 'recheck-e1-kernel', 'method': 'jacobian.mult_kernel', 'p': p, 'd': d, 'parity': par,
               'ker_mod_J': out['dim'], 'kernel_dim': out['kernel_dim'], 'J_d_rank_p': out['jacobian_dim'],
               'J_d_CI': ci_d, 'J_target_rank_p': rJ_target, 'J_target_CI': ci_target,
               'ring_dim_p': out['ring_dim'], 'R_d_CI': hilbert_CI_parity(d, par),
               'rows': out['rows'], 'cols': out['cols'], 'seconds': round(time.time() - t0, 1),
               'rigorous_over_Q': bool(out['dim'] == 0 and out['jacobian_dim'] == ci_d and rJ_target == ci_target)}
        print(json.dumps(rec), flush=True)
        log(rec)
        if out['dim'] and d >= 10:
            import pickle
            with open(f'{PK}/kernel_basis_d{d}_par{par}_p{p}.pkl', 'wb') as f:
                pickle.dump({'p': p, 'basis': out['basis']}, f)
