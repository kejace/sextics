"""Job C (excl): degree-6 denominators, FULL 6-variable fixed-D system (no sector restriction), two new primes,
exact-nullity [RANK] argument for every rational D.
E1-admissible D of degree 6 (Thom-Sebastiani):  D = sum_i w_i F_i (w_i linear)  +  c(x),  c in ker_6(H) (20 even + 13 odd reps).
Tested: the 33 exact reps; H, B and other structural D; random small-integer combos of the reps (per parity);
random members of the full family J(F)_6 + ker (per parity); mixed-parity members (parity-free system);
pencils c + lambda*B.  Multiprocessing."""
import sys, time, json, random, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from multiprocessing import Pool

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-d6'
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 16
F, A = W.F, W.A
Fi = [pderiv(F, i) for i in range(NV)]
P = core.parse

def task(args):
    tag, family, Dj, p, pf = args
    D = {tuple(int(x) for x in k.split()): Fraction(v) for k, v in Dj.items()}
    try:
        if not pf:
            ok, info = W.e1_test(D, p)
            if not ok:
                return {'exp': 'd6-fixedD', 'D': tag, 'family': family, 'd': 6, 'prime': p, 'e1': False,
                        'result': 'E1 fails (A D not in J(F)) [MODP]', 'rank': info['rank'], 'seconds': 0.0}
        r = W.fixed_D_test(D, p, parity_free=pf)
        return {'exp': 'd6-fixedD', 'D': tag, 'family': family, 'prime': p, 'e1': True, **r}
    except Exception as ex:
        return {'exp': 'd6-fixedD', 'D': tag, 'family': family, 'prime': p, 'result': f'ERROR {ex!r}'}

if __name__ == '__main__':
    t0 = time.time()
    reps = {}
    for par in (0, 1):
        r = W.mult_kernel(W.H, W.CORE, W.A, 6, par, p=None)
        reps[par] = W.small_height_reps(r['reps'])
        log(JOB, {'exp': 'd6-exact-reps', 'parity': par, 'n_reps': len(reps[par]), 'dim_modJ': r['dim_modJ'],
                  'max_abs_coeff': max(abs(int(c)) for rp in reps[par] for c in rp.values()), 'seconds': round(time.time() - t0, 1)})
    rng = random.Random(2026)
    tasks = []
    def J(D): return W.poly_to_json(D)
    # structural
    f = {e: c for e, c in W.H.items() if e[3] == 0 and e[4] == 0}
    g = {e: c for e, c in W.H.items() if e[0] == 0 and e[1] == 0}
    K = sub(sub(W.H, f), g)
    struct = [('H', W.H), ('B', W.B), ('f', f), ('g', g), ('K', K), ('f+g', add(f, g)), ('H-B', sub(W.H, W.B)),
              ('X0^3*Y0^3 (=A)', W.A), ('X0^6+Y0^6', P('X0**6+Y0**6')), ('X0*X1*Y0*Y1*(X0^2+2*Y0^2)', P('X0*X1*Y0*Y1*(X0**2+2*Y0**2)')),
              ('X0^2 Y0^2 (X0^2+Y0^2)', P('X0**2*Y0**2*(X0**2+Y0**2)')), ('(X0,Y0)-radial dF', add(mul(P('X0'), Fi[0]), mul(P('Y0'), Fi[3]))),
              ('(X2,Y2)-radial dF = 6B', add(mul(P('X2'), Fi[2]), mul(P('Y2'), Fi[5])))]
    for name, D in struct:
        for p in (W.P2, W.P3):
            tasks.append((name, 'structural', J(D), p, core.parity_of(D) is None))
    # the 33 reps
    for par in (0, 1):
        for i, c in enumerate(reps[par]):
            for p in (W.P2, W.P3):
                tasks.append((f'ker6 rep #{i} par{par}', f'ker6-rep-par{par}', J(c), p, False))
    # random small-integer combos of reps (rational D)
    for par in (0, 1):
        for t in range(4):
            D = {}
            co = [rng.randint(-3, 3) for _ in reps[par]]
            for c, x in zip(reps[par], co):
                if x: D = add(D, scale(c, Fraction(x)))
            for p in (W.P2, W.P3):
                tasks.append((f'random combo of ker6 reps par{par} #{t} coeffs={co}', f'ker6-random-par{par}', J(D), p, False))
    # random members of the full family J(F)_6 + ker (per parity): sum_i w_i F_i + sum c_j rep_j
    def rand_lin(par):
        return {m: Fraction(rng.randint(-2, 2)) for m in monomials(1, par)}
    for par in (0, 1):
        for t in range(4):
            D = {}
            for i in range(NV):
                wpar = (par - (0 if i < 3 else 1)) % 2
                D = add(D, mul(rand_lin(wpar), Fi[i]))
            for c in reps[par]:
                x = rng.randint(-2, 2)
                if x: D = add(D, scale(c, Fraction(x)))
            for p in (W.P2, W.P3):
                tasks.append((f'random member of J(F)_6+ker6 par{par} #{t}', f'E1-family-random-par{par}', J(D), p, False))
    # mixed parity members: D_even + D_odd (parity-free system)
    for t in range(3):
        D = {}
        for par in (0, 1):
            for c in reps[par]:
                x = rng.randint(-2, 2)
                if x: D = add(D, scale(c, Fraction(x)))
        for i in range(NV):
            D = add(D, mul({m: Fraction(rng.randint(-1, 1)) for m in monomials(1)}, Fi[i]))
        for p in (W.P2, W.P3):
            tasks.append((f'random MIXED-parity member of J(F)_6+ker6 #{t}', 'E1-family-random-mixed', J(D), p, True))
    # pencils c + lambda B for the first 3 even reps and lambda in {1,-1,3}
    for i, c in enumerate(reps[0][:3]):
        for lam in (1, -1, 3):
            D = add(c, scale(W.B, Fraction(lam)))
            tasks.append((f'ker6 rep #{i} par0 + ({lam})*B', 'pencil-rep+B', J(D), W.P2, False))
    # c + beta with beta in J(B)_6 (non-eigenvector): Y2*B_X2 or X2*B_Y2
    BX, BY = pderiv(W.B, 2), pderiv(W.B, 5)
    for i, c in enumerate(reps[0][:2]):
        tasks.append((f'ker6 rep #{i} par0 + Y2*dB/dX2', 'rep+J(B)', J(add(c, mul(P('Y2'), BX))), W.P2, True))
        tasks.append((f'ker6 rep #{i} par0 + X2*dB/dY2', 'rep+J(B)', J(add(c, mul(P('X2'), BY))), W.P2, True))
    log(JOB, {'exp': 'd6-plan', 'n_tasks': len(tasks), 'workers': NW})
    with Pool(NW) as pool:
        for res in pool.imap_unordered(task, tasks):
            log(JOB, res)
    log(JOB, {'exp': 'done', 'seconds': round(time.time() - t0, 1)})
