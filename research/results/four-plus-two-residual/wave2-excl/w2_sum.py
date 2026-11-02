"""Job D (excl): served-target sectors and the sum principle (q=1): served(D1 D2) = served(D1) + served(D2) mod J(F).
Computes (mod P2, character-sector restricted, lossless for eigenvector D) the space served(D) of targets A' in the
G-character sector of A (40 Y-odd core sextic monomials) for: the binary coordinate lines X2, Y2, the six linear
factors of B over a prime where B splits, random lines X2 - rho Y2, the core kernel reps (k=6..10, both parities),
structural core D's; then tests whether A lies in served(c) + served(l) + J for every pair (c core, l binary),
which by the sum principle excludes every product denominator c(x) * (product of tested binary lines).
Consistency check of the sum principle: served(B) directly vs sum over its six linear factors."""
import sys, time, json, random, os, itertools
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log, save_json
from hodge import core, symmetry
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV, reduce_mod
from fractions import Fraction
from flint import nmod_poly, nmod_mat
from multiprocessing import Pool

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-sum'
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 24
P = core.parse

def splitting_prime(start=1000000007):
    """A prime p > start such that b(r,1) = r^3 + 15/7 r^2 + 10/7 r + 1 has three roots in F_p that are all squares
    => B = prod (X2 - s_i Y2)(X2 + s_i Y2) splits into six F_p-rational lines."""
    from sympy import nextprime
    p = start
    while True:
        p = int(nextprime(p))
        inv7 = pow(7, -1, p)
        f = nmod_poly([1, 10 * inv7 % p, 15 * inv7 % p, 1], p)   # coefficients low->high: 1 + 10/7 r + 15/7 r^2 + r^3
        fac = f.factor()[1]
        roots = [(-int(g.coeffs()[0]) * pow(int(g.coeffs()[1]), -1, p)) % p for g, e in fac if g.degree() == 1]
        if len(roots) == 3:
            sq = []
            for r in roots:
                s = pow(r, (p + 1) // 4, p) if p % 4 == 3 else None
                if s is None:
                    # Tonelli-Shanks via sympy
                    from sympy.ntheory import sqrt_mod
                    s = sqrt_mod(r, p)
                if s is not None and (s * s) % p == r:
                    sq.append(s)
            if len(sq) == 3:
                return p, roots, sq

def served(args):
    tag, Dj, p, pf = args
    D = {tuple(int(x) for x in k.split()): (int(v) if pf == 'modp' else Fraction(v)) for k, v in Dj.items()}
    t0 = time.time()
    try:
        basis, S = W.served_sector(D, p)
        return {'exp': 'served-sector', 'D': tag, 'd': core.degree(D), 'prime': p, 'dim_served_sector': len(basis),
                'cols': len(S['cols']), 'nT': S['nT'], 'seconds': round(time.time() - t0, 1), 'basis': [W.poly_to_json(b) for b in basis]}
    except Exception as ex:
        return {'exp': 'served-sector', 'D': tag, 'prime': p, 'result': f'ERROR {ex!r}'}

if __name__ == '__main__':
    t0 = time.time()
    p, roots, sq = splitting_prime()
    log(JOB, {'exp': 'splitting-prime', 'prime': p, 'roots_of_b': roots, 'sqrt_roots': sq, 'seconds': round(time.time() - t0, 1),
              'note': 'B(X2,Y2) = prod_i (X2^2 - r_i Y2^2) = prod (X2 -/+ s_i Y2) over F_p'})
    # verify factorisation mod p
    Bp = reduce_mod(W.B, p)
    prod = {(0,) * NV: 1}
    lines = []
    for s in sq:
        for sg in (1, -1):
            l = {(0, 0, 1, 0, 0, 0): 1, (0, 0, 0, 0, 0, 1): (-sg * s) % p}
            lines.append(l); prod = mul(prod, l, p)
    log(JOB, {'exp': 'B-splits-check', 'prime': p, 'result': 'B == product of 6 lines mod p' if prod == Bp else 'MISMATCH'})
    assert prod == Bp
    Dlist = []
    def addD(tag, D, pf='exact'):
        Dlist.append((tag, {' '.join(map(str, e)): str(c) for e, c in D.items()}, p, pf))
    addD('X2', P('X2')); addD('Y2', P('Y2'))
    for i, l in enumerate(lines): addD(f'line_{i}: X2 - ({l[(0,0,0,0,0,1)]})*Y2 [factor of B]', l, 'modp')
    rng = random.Random(3)
    for t in range(3):
        rho = rng.randrange(2, p)
        addD(f'random line X2 - {rho} Y2', {(0, 0, 1, 0, 0, 0): 1, (0, 0, 0, 0, 0, 1): (-rho) % p}, 'modp')
    addD('B', W.B); addD('X2^2 - r0 Y2^2', {(0, 0, 2, 0, 0, 0): 1, (0, 0, 0, 0, 0, 2): (-roots[0]) % p}, 'modp')
    f = {e: c for e, c in W.H.items() if e[3] == 0 and e[4] == 0}
    g = {e: c for e, c in W.H.items() if e[0] == 0 and e[1] == 0}
    K = sub(sub(W.H, f), g)
    for name, D in (('H', W.H), ('f', f), ('g', g), ('K', K), ('f+g', add(f, g)), ('X0', P('X0')), ('Y0', P('Y0')), ('X0*X1*Y0*Y1', P('X0*X1*Y0*Y1')),
                    ('X0^2+Y0^2', P('X0**2+Y0**2')), ('X0*X1+Y0*Y1', P('X0*X1+Y0*Y1')), ('X0^3+Y0^3', P('X0**3+Y0**3')), ('X0^6+Y0^6', P('X0**6+Y0**6')),
                    ('dH/dX0', pderiv(W.H, 0)), ('dH/dY0', pderiv(W.H, 3))):
        addD(name, D)
    reps_all = {}
    for k in range(6, 11):
        for par in (0, 1):
            r = W.small_kernel_reps(W.H, W.CORE, W.A, k, par)
            reps = r['reps_modJ']
            reps_all[(k, par)] = reps
            for i, c in enumerate(reps):
                addD(f'ker{k} rep #{i} par{par}', c)
    log(JOB, {'exp': 'sum-plan', 'n_D': len(Dlist), 'workers': NW})
    results = {}
    with Pool(NW) as pool:
        for res in pool.imap_unordered(served, Dlist):
            results[res['D']] = res
            log(JOB, {k: v for k, v in res.items() if k != 'basis'})
    save_json(JOB, 'served', results)
    # ---- membership tests (mod p): A in served(c) + served(l) + J_sector ?
    Jcols = [c for c in W.jacobian_target_cols(p) if all(e in set(W.TARGET_SECTOR) for e in c)]
    Ap = reduce_mod(W.A, p)
    def basis_of(tag):
        r = results.get(tag, {})
        return [{tuple(int(x) for x in k.split()): int(v) for k, v in b.items()} for b in r.get('basis', [])]
    def span_contains_A(tags):
        polys = list(Jcols)
        for t in tags: polys += basis_of(t)
        r0 = W.rank_modp(polys, p) if polys else 0
        r1 = W.rank_modp(polys + [Ap], p)
        return r1 == r0, r0
    bin_tags = ['X2', 'Y2', 'B'] + [t for t in results if t.startswith('line_') or t.startswith('random line') or t.startswith('X2^2')]
    core_tags = [t for t in results if t not in bin_tags]
    # sum principle check: served(B) vs sum of its six lines
    inB, rB = span_contains_A(['B']); lines_tags = [t for t in results if t.startswith('line_')]
    polysB = list(Jcols) + basis_of('B'); polysL = list(Jcols)
    for t in lines_tags: polysL += basis_of(t)
    rB_ = W.rank_modp(polysB, p); rL = W.rank_modp(polysL, p); rBL = W.rank_modp(polysB + polysL, p)
    log(JOB, {'exp': 'sum-principle-check', 'prime': p, 'dim_served_B_plus_J': rB_, 'dim_sum_lines_plus_J': rL, 'dim_join': rBL,
              'equal': (rB_ == rL == rBL), 'note': 'served(B) should equal sum_i served(line_i) modulo J (q=1 sum principle)'})
    # all binary subsets (union of all binary lines) + each core D
    allbin = [t for t in results if t.startswith('line_') or t in ('X2', 'Y2') or t.startswith('random line')]
    hits = []
    for ct in core_tags:
        ok_single, r_single = span_contains_A([ct])
        ok_all, r_all = span_contains_A([ct] + allbin)
        log(JOB, {'exp': 'sum-membership', 'core_D': ct, 'A_in_served_core_plus_J': ok_single, 'dim_core_plus_J': r_single,
                  'A_in_served_core_plus_ALLbinary_plus_J': ok_all, 'dim_all': r_all, 'prime': p})
        if ok_single or ok_all: hits.append(ct)
    ok_bin, r_bin = span_contains_A(allbin)
    log(JOB, {'exp': 'sum-membership', 'core_D': '(none)', 'A_in_served_ALLbinary_plus_J': ok_bin, 'dim': r_bin, 'prime': p})
    # all core reps together + all binary
    ok_tot, r_tot = span_contains_A(core_tags + allbin)
    log(JOB, {'exp': 'sum-membership-total', 'A_in_sum_of_everything_plus_J': ok_tot, 'dim': r_tot, 'n_core': len(core_tags), 'n_bin': len(allbin), 'prime': p,
              'hits': hits, 'note': 'if False, NO product of the tested denominators (any multiplicities) serves A [MODP, sum principle]'})
    log(JOB, {'exp': 'done', 'seconds': round(time.time() - t0, 1)})
