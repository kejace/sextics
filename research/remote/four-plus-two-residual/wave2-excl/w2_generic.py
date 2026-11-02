"""Job G (excl): generic members of the full E1-admissible family at degrees 7..10, FULL unrestricted 6-variable
system.  Thom-Sebastiani description of {D : A D in J(F)}_d :
   D = sum_i w_i F_i  (w_i of degree d-5)  +  sum_{kb} sum_j c_{kb,j}(x) b_{kb,j}(X2,Y2),
   c in ker_{d-kb}(H) exact reps, b random binary forms of degree kb.
Random small-integer coefficients => rational D => [RANK] via exact Koszul nullity for d <= 8; [MODP] (dense for d<=8,
compressed for d = 9, 10) otherwise.  Pure Y-parities and mixed parity (parity-free system)."""
import sys, time, json, random, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from multiprocessing import Pool

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-generic'
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
DMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 7
DMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 10
NTRIALS = int(sys.argv[5]) if len(sys.argv) > 5 else 2
F, A = W.F, W.A
Fi = [pderiv(F, i) for i in range(NV)]

def task(args):
    tag, family, Dj, p, pf, compressed = args
    D = {tuple(int(x) for x in k.split()): Fraction(v) for k, v in Dj.items()}
    t0 = time.time()
    try:
        d = core.degree(D)
        r = W.fixed_D_test(D, p, parity_free=pf, compressed=compressed, exact_nullity=(d <= 8))
        return {'exp': 'generic-fixedD', 'D': tag, 'family': family, 'prime': p, **r}
    except Exception as ex:
        return {'exp': 'generic-fixedD', 'D': tag, 'family': family, 'prime': p, 'result': f'ERROR {ex!r}', 'seconds': round(time.time() - t0, 1)}

if __name__ == '__main__':
    t0 = time.time()
    rng = random.Random(4242)
    reps = {}
    for k in range(6, DMAX + 1):
        for par in (0, 1):
            r = W.small_kernel_reps(W.H, W.CORE, W.A, k, par)
            reps[(k, par)] = r['reps_modJ']
    log(JOB, {'exp': 'generic-reps', 'dims': {f'{k}-{par}': len(v) for (k, par), v in reps.items()}, 'seconds': round(time.time() - t0, 1)})
    def rand_bin(kb, y2par):
        """random binary form of degree kb in (X2,Y2) with Y2-degree parity y2par (None: any)"""
        out = {}
        for a_ in range(kb + 1):
            b_ = kb - a_
            if y2par is not None and b_ % 2 != y2par: continue
            x = rng.randint(-2, 2)
            if x: out[(0, 0, a_, 0, 0, b_)] = Fraction(x)
        return out
    def member(d, par):
        """random E1-admissible D of degree d and Y-parity par (None = mixed)"""
        D = {}
        for i in range(NV):
            wpar = None if par is None else (par - (0 if i < 3 else 1)) % 2
            w = {m: Fraction(rng.randint(-1, 1)) for m in monomials(d - 5, wpar) if rng.random() < 0.5}
            D = add(D, mul(w, Fi[i]))
        for kb in range(0, min(d - 6, 8) + 1):
            kc = d - kb
            for cpar in (0, 1):
                if (kc, cpar) not in reps: continue
                y2par = None if par is None else (par - cpar) % 2
                for c in reps[(kc, cpar)]:
                    if rng.random() < 0.6:
                        b = rand_bin(kb, y2par)
                        if b: D = add(D, mul(c, b))
        return D
    tasks = []
    for d in range(DMIN, DMAX + 1):
        comp = d >= 9
        for par in (0, 1, None):
            for t in range(NTRIALS):
                D = member(d, par)
                if not D: continue
                # sanity: parity
                assert (core.parity_of(D) == par) if par is not None else True
                tag = f'random E1-admissible member d={d} parity={par} #{t}'
                fam = f'E1-family-random-d{d}-' + ('mixed' if par is None else f'par{par}')
                tasks.append((tag, fam, W.poly_to_json(D), W.P2, par is None, comp or par is None and d >= 8))
    log(JOB, {'exp': 'generic-plan', 'n_tasks': len(tasks), 'workers': NW})
    tasks.sort(key=lambda t: -int(t[1].split('-d')[1].split('-')[0]))
    with Pool(NW) as pool:
        for res in pool.imap_unordered(task, tasks):
            log(JOB, res)
    log(JOB, {'exp': 'done', 'seconds': round(time.time() - t0, 1)})
