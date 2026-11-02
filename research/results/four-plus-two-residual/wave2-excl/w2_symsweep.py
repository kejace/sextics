"""Job E (excl): core-only E1-kernel elements of degrees 7..12 (exact rational reps modulo J(H)_k, both
(Y0,Y1)-parities) in the FULL 6-variable fixed-D system restricted to the character sector forced by the
diagonal group G (order 24; lossless for G-eigenvectors, which all core-only pure-parity D are).
For each D: mod-p ranks (P2) + EXACT rational ranks of M and [M|b] (fmpz_mat) when the restricted system is
small enough => [EXACT] exclusion over C; else exact rank of M + modular augmented rank => [RANK].
Also random small-integer combos of the reps per (k, parity).  Multiprocessing."""
import sys, time, json, random, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core, symmetry
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from multiprocessing import Pool

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-symsweep'
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 32
KMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 7
KMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 12
EXACT_COLS_MAX = int(sys.argv[5]) if len(sys.argv) > 5 else 2500

def task(args):
    tag, family, Dj, p = args
    D = {tuple(int(x) for x in k.split()): Fraction(v) for k, v in Dj.items()}
    t0 = time.time()
    try:
        d = core.degree(D)
        allowed, sub_ = W.sector_allowed(D)
        S = W.build_fixed_D(D, p, allowed=allowed)
        r = W.modp_consistency(S)
        out = {'exp': 'sym-fixedD', 'D': tag, 'family': family, 'd': d, 'prime': p, 'stabilizer': len(sub_),
               'rows': r['rows'], 'cols': r['cols'], 'rank_M_modp': r['rank_M'], 'rank_Mb_modp': r['rank_Mb'],
               'consistent_modp': r['consistent'], 'nullity_modp': r['nullity']}
        if r['consistent']:
            out['result'] = 'CONSISTENT mod p -- LIFT NEEDED'
        elif d <= 8:
            # [RANK]: exact nullity of M from the Koszul parametrisation (restricted to the sector), compared with the modular nullity
            kz = W.koszul_nullity_exact(D, allowed=allowed)
            out['nullity_exact_Q'] = kz['nullity_total']
            out['koszul'] = {k_: v_ for k_, v_ in kz.items() if k_ != 'seconds'}
            out['result'] = 'inconsistent over C [RANK]' if kz['nullity_total'] == r['nullity'] else 'inconsistent mod p (nullity mismatch: inconclusive)'
        elif r['cols'] <= EXACT_COLS_MAX:
            Sx = W.build_fixed_D(D, None, allowed=allowed)
            rM = W.rank_exact(Sx['colpolys'])
            out['rank_M_Q'] = rM
            out['result'] = 'inconsistent over C [RANK]' if rM == r['rank_M'] else 'inconsistent mod p (rank_Q(M) > rank_p(M): inconclusive)'
        else:
            out['result'] = 'inconsistent mod p [MODP]'
        out['seconds'] = round(time.time() - t0, 1)
        return out
    except Exception as ex:
        return {'exp': 'sym-fixedD', 'D': tag, 'family': family, 'prime': p, 'result': f'ERROR {ex!r}', 'seconds': round(time.time() - t0, 1)}

if __name__ == '__main__':
    t0 = time.time()
    rng = random.Random(77)
    tasks = []
    for k in range(KMIN, KMAX + 1):
        for par in (0, 1):
            r = W.small_kernel_reps(W.H, W.CORE, W.A, k, par)
            reps = r['reps_modJ']
            log(JOB, {'exp': 'sym-exact-reps', 'd': k, 'parity': par, 'n_reps': len(reps), 'dim_modJ': len(reps), 'dim_J': r['dim_J'], 'verified_exact': r['verified'],
                      'max_abs_coeff': max((abs(int(c)) for rp in reps for c in rp.values()), default=0), 'seconds': round(time.time() - t0, 1)})
            for i, c in enumerate(reps):
                tasks.append((f'ker{k} rep #{i} par{par}', f'ker{k}-rep-par{par}', W.poly_to_json(c), W.P2))
            for t in range(3):
                D = {}
                co = [rng.randint(-3, 3) for _ in reps]
                for c, x in zip(reps, co):
                    if x: D = add(D, scale(c, Fraction(x)))
                if D:
                    tasks.append((f'random combo of ker{k} reps par{par} #{t}', f'ker{k}-random-par{par}', W.poly_to_json(D), W.P2))
    log(JOB, {'exp': 'sym-plan', 'n_tasks': len(tasks), 'workers': NW, 'kmin': KMIN, 'kmax': KMAX})
    # schedule large degrees first (longest tasks) for better packing
    tasks.sort(key=lambda t: -int(t[1].split('-')[0][3:]))
    with Pool(NW) as pool:
        for res in pool.imap_unordered(task, tasks):
            log(JOB, res)
    log(JOB, {'exp': 'done', 'seconds': round(time.time() - t0, 1)})
