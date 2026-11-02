"""Job A (excl): re-verify the E1-kernel structure with new primes and EXACT rational computations.
 (1) regular-sequence (smoothness) checks of H (core) and B (binary) mod P1,P2,P3 and over Q;
 (2) core multiplication kernels ker(a. : R(H)_k -> R(H)_{k+6}) for k=0..16, both parities, mod P1,P2,P3 and
     EXACT over Q for k<=12 (basis reps saved, LLL-reduced);
 (3) direct 6-variable E1 kernels {D : A D in J(F)} for d<=8, both parities, mod P2, compared with the
     Thom-Sebastiani prediction  ker_d(F) = sum_k ker_{d-k}(H) (x) R(B)_k.
Outputs: remote/<job>.jsonl, remote/<job>-corereps.json"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log, save_json
from hodge import core
from hodge.core import monomials, mul, add, shift, pderiv, reduce_mod
from fractions import Fraction

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-kernels'

# (1) regular sequences
for name, P, idx in (('B', W.B, W.BIN), ('H', W.H, W.CORE)):
    for p in (W.P1, W.P2, W.P3, None):
        t0 = time.time()
        ok, det = W.regular_sequence_check(P, idx, p=p)
        log(JOB, {'exp': 'regular-sequence', 'poly': name, 'prime': p, 'result': 'regular' if ok else 'NOT regular',
                  'hilbert': det['hilbert'], 'seconds': round(time.time() - t0, 1),
                  'note': 'Artinian Jacobian ring with CI Hilbert function => partials form a regular sequence (=> V smooth); '
                          'exact over Q when prime is null'})

# (2) core kernels
dims = {}
reps_out = {}
for k in range(0, 17):
    for par in (0, 1):
        row = {}
        for p in (W.P1, W.P2, W.P3):
            t0 = time.time()
            r = W.mult_kernel(W.H, W.CORE, W.A, k, par, p=p)
            row[str(p)] = (r['dim_full'], r['dim_J'], r['dim_modJ'])
            log(JOB, {'exp': 'core-kernel-dim', 'd': k, 'parity': par, 'prime': p, 'n_mons': r['n_mons'],
                      'dim_full': r['dim_full'], 'dim_J': r['dim_J'], 'dim_modJ': r['dim_modJ'],
                      'seconds': round(time.time() - t0, 2)})
        if k <= 12:
            t0 = time.time()
            r = W.mult_kernel(W.H, W.CORE, W.A, k, par, p=None)
            row['Q'] = (r['dim_full'], r['dim_J'], r['dim_modJ'])
            reps = W.small_height_reps(r['reps']) if r['reps'] else []
            # sanity: LLL-reduced reps still independent mod J and still in the kernel (exact)
            Hi = [pderiv(W.H, i) for i in W.CORE]
            Jk = []
            if k >= 5:
                for j, i in enumerate(W.CORE):
                    pj = 0 if i < 3 else 1
                    for m in W.mons(k - 5, W.CORE, (par - pj) % 2):
                        Jk.append(shift(Hi[j], m, 1))
            rJ = W.rank_exact(Jk) if Jk else 0
            indep = (W.rank_exact(Jk + reps) if (Jk or reps) else 0) == rJ + len(reps)
            # kernel membership: a*rep in J(H)_{k+6}
            tpar = (par + 1) % 2
            J6 = []
            for j, i in enumerate(W.CORE):
                pj = 0 if i < 3 else 1
                for m in W.mons(k + 1, W.CORE, (tpar - pj) % 2):
                    J6.append(shift(Hi[j], m, 1))
            rJ6 = W.rank_exact(J6) if J6 else 0
            inker = all(W.rank_exact(J6 + [mul(W.A, rp)]) == rJ6 for rp in reps)
            maxh = max((abs(int(c)) for rp in reps for c in rp.values()), default=0)
            reps_out[f'{k}-{par}'] = [W.poly_to_json(rp) for rp in reps]
            log(JOB, {'exp': 'core-kernel-exact', 'd': k, 'parity': par, 'prime': None, 'n_mons': r['n_mons'],
                      'dim_full': r['dim_full'], 'dim_J': r['dim_J'], 'dim_modJ': r['dim_modJ'], 'n_reps': len(reps),
                      'reps_independent_modJ': indep, 'reps_in_kernel': inker, 'max_abs_coeff': maxh,
                      'seconds': round(time.time() - t0, 2),
                      'note': 'EXACT over Q (fmpz nullspace, LLL-reduced representatives modulo J(H)_k)'})
        dims[f'{k}-{par}'] = row
    save_json(JOB, 'corereps', {'dims': dims, 'reps': reps_out,
                                'format': 'reps[k-par] = list of polys {"e0 e1 e2 e3 e4 e5": "coeff"}; exact rational kernel reps modulo J(H)_k'})

# (3) direct 6-variable kernels vs Thom-Sebastiani
RB = [1, 2, 3, 4, 5, 4, 3, 2, 1]   # dim R(B)_k
def ts_pred(d, par):
    tot = 0
    for kb in range(0, min(d, 8) + 1):
        kc = d - kb
        if kc > 16: continue
        # core parity = par - (Y2-degree parity of the binary monomial); R(B)_kb splits by Y2-parity:
        # standard monomials X2^a Y2^b with a+b=kb, a,b<=4 (Fermat-type basis mod J(B) has the same parity split)
        for a_ in range(0, kb + 1):
            b_ = kb - a_
            if a_ > 4 or b_ > 4: continue
            cp = (par - b_) % 2
            tot += dims[f'{kc}-{cp}'][str(W.P2)][2]
    return tot
for d in range(1, 9):
    for par in (0, 1):
        t0 = time.time()
        r = W.mult_kernel(W.F, tuple(range(6)), W.A, d, par, p=W.P2)
        log(JOB, {'exp': '6var-kernel-dim', 'd': d, 'parity': par, 'prime': W.P2, 'n_mons': r['n_mons'],
                  'dim_full': r['dim_full'], 'dim_J': r['dim_J'], 'dim_modJ': r['dim_modJ'],
                  'ts_prediction_modJ': ts_pred(d, par), 'agree': ts_pred(d, par) == r['dim_modJ'],
                  'seconds': round(time.time() - t0, 2),
                  'note': 'direct 6-variable {D: A D in J(F)} modulo J(F)_d vs Thom-Sebastiani sum_k ker_{d-k}(H) x R(B)_k'})
log(JOB, {'exp': 'done'})
