"""Referee replay of the symsweep method at d = 7 and 8: regenerate the agent's deterministic kernel reps (small_kernel_reps),
re-verify them exactly, then test them in the character-sector-restricted full system with the referee's own code at p = 10^9+9."""
import sys, time; sys.path.insert(0, '/tmp/hodge/work/referee42'); sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import reflib as R
from fractions import Fraction
LOG = '/tmp/hodge/work/referee42/out/ref_t5.jsonl'
PNEW = 1000000009
F, A = R.load_instance()
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
import w2lib as W
for k, par, idxs in ((7, 0, (0,)), (7, 1, (0,)), (8, 0, (0,)), (8, 1, (0,))):
    t0 = time.time()
    r = W.small_kernel_reps(W.H, W.CORE, W.A, k, par)
    reps = r['reps_modJ']
    for i in idxs:
        D = {e: Fraction(c) for e, c in reps[i].items()}
        Jt = R.jac_cols(H, R.CORE, k + 6, (par + 1) % 2); rJ = R.rank_exact(Jt)
        inker = R.rank_exact(Jt + [R.pmul(A, D)]) == rJ
        Jk = R.jac_cols(H, R.CORE, k, par); notinJ = R.rank_exact(Jk + [D]) == R.rank_exact(Jk) + 1
        allowed = R.sector_filter(F, A, D)
        res = R.fixed_D_test(F, A, D, PNEW, allowed=allowed)
        res.update({'exp': 'symsweep-replay', 'D': f'ker{k} rep #{i} par{par}', 'n_reps': len(reps), 'A*D_in_J(H)_exact': inker, 'D_not_in_J(H)_k_exact': notinJ,
                    'gen_seconds': round(time.time() - t0, 1)})
        R.log(LOG, res)
