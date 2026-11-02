"""Sample F_p-rational points of the trivial-sector relaxation surface V0 (D = sum a_i z_i F_i, deg 6) and run the full
E1+E2 sector solve at each (screening of the generic member of V0; a CONSISTENT point would be a witness candidate).
usage: python3 wave2-excl/v0_sample.py [--n 6] [--prime P] [--seed S] [--farkas]
"""
import sys, json, time, argparse, os, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hodge import core
from hodge.core import mul, add, scale, pderiv
from flint import nmod_poly, nmod

E01 = (3,1,0,1,1,0); E12 = (0,3,1,0,1,1); E20 = (1,0,3,1,0,1)

def quad_coeffs(e, a, unknown, p):
    """coefficients (c2, c1, c0) of Q_e(a) as a polynomial in a[unknown] (other entries fixed ints), mod p."""
    assert all(e[i] == 0 for i in range(6) if a[i] is None and i != unknown), 'unknown ordering'
    s0 = sum(e[i] * a[i] for i in range(6) if i != unknown and e[i]); t0 = sum(e[i] * a[i] * a[i] for i in range(6) if i != unknown and e[i])
    eu = e[unknown]
    # (s0 + eu x)^2 - 6 (t0 + eu x^2) = (eu^2 - 6 eu) x^2 + 2 s0 eu x + (s0^2 - 6 t0)
    return ((eu * eu - 6 * eu) % p, (2 * s0 * eu) % p, (s0 * s0 - 6 * t0) % p)

def roots(coeffs, p):
    c2, c1, c0 = coeffs
    f = nmod_poly([c0, c1, c2], p)
    if f.degree() < 0: return None   # identically zero
    return [int(r) for r, m in f.roots()]

def eliminant_a5(a0, a1, a2, p):
    """lex Groebner basis (Singular) of the three quadrics in (a3,a4,a5) with a0,a1,a2 fixed; returns the univariate
    polynomial in a5 (coefficient list, low to high) or None."""
    import subprocess, re
    script = (f'ring r = {p}, (a3,a4,a5), lp;\n'
              f'ideal I = (3*{a0}+{a1}+a3+a4)^2 - 6*(3*{a0}^2+{a1}^2+a3^2+a4^2), (3*{a1}+{a2}+a4+a5)^2 - 6*(3*{a1}^2+{a2}^2+a4^2+a5^2), (3*{a2}+{a0}+a5+a3)^2 - 6*(3*{a2}^2+{a0}^2+a5^2+a3^2);\n'
              'ideal G = std(I); int i; for (i = 1; i <= size(G); i++) { if (size(variables(G[i])) == 1) { if (variables(G[i])[1] == a5) { "ELIM"; G[i]; } } }\nquit;\n')
    r = subprocess.run(['Singular', '-q'], input=script, capture_output=True, text=True, timeout=120)
    m = re.search(r'ELIM\n(.*)\n', r.stdout)
    if not m: return None
    expr = m.group(1).replace('^', '**')
    # evaluate coefficients via sympy-free parsing: use python eval on nmod_poly? simpler: parse terms
    coeffs = {}
    for term in re.findall(r'[+-]?[^+-]+', expr.replace(' ', '')):
        if not term: continue
        sign = -1 if term.startswith('-') else 1
        term = term.lstrip('+-')
        if 'a5' in term:
            c, _, rest = term.partition('a5')
            c = int(c.rstrip('*')) if c.rstrip('*') else 1
            d = int(rest.lstrip('**')) if rest else 1
        else:
            c = int(term); d = 0
        coeffs[d] = (coeffs.get(d, 0) + sign * c) % p
    deg = max(coeffs)
    return [coeffs.get(i, 0) for i in range(deg + 1)]

def points_on_V0(p, rng, want):
    pts = []; tries = 0
    while len(pts) < want and tries < 60:
        tries += 1
        a0, a1, a2 = rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p)
        cf = eliminant_a5(a0, a1, a2, p)
        if not cf: continue
        f = nmod_poly(cf, p)
        for r5, mult in f.roots():
            a = [a0, a1, a2, None, None, int(r5)]
            # q12 -> a4 (quadratic, a3 absent)
            for x4 in roots(quad_coeffs(E12, a, 4, p), p) or []:
                a[4] = x4
                r3a = set(roots(quad_coeffs(E01, a, 3, p), p) or []); r3b = set(roots(quad_coeffs(E20, a, 3, p), p) or [])
                for x3 in r3a & r3b:
                    a[3] = x3
                    pts.append(list(a))
            a[4] = None
    return pts, tries

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=6); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--farkas', action='store_true')
    a_ = ap.parse_args(); p = a_.prime
    rng = random.Random(a_.seed)
    F, A, q, G = w2lib.load()
    Z = [core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
    ZF = [mul(Z[i], pderiv(F, i)) for i in range(6)]
    pts, tries = points_on_V0(p, rng, a_.n)
    print(f'{len(pts)} points on V0(F_p) after {tries} fibres', flush=True)
    for a in pts:
        t0 = time.time()
        D = {}
        for i in range(6):
            D = add(D, scale(ZF[i], a[i]))
        Dp = {e: int(v) % p for e, v in D.items()}
        Dp = {e: v for e, v in Dp.items() if v}
        Dsym = {e: (v if v <= p // 2 else v - p) for e, v in Dp.items()}
        if a_.farkas:
            fr = w2lib.farkas(F, A, q, Dsym, p, G)
            res = fr['status']; rk = f"{fr.get('rank_M')}/{fr.get('rank_aug')} of {fr['cols']}"
        else:
            out, S, r = w2lib.solve_modp(F, A, q, Dsym, p, G)
            res = 'consistent-mod-p' if out['consistent'] else 'inconsistent-mod-p'; rk = f"{out['rank']}/{out['cols']}"
        rec = {'exp': 'V0-sample', 'd': 6, 'prime': p, 'a': a, 'result': res, 'rank': rk, 'seconds': round(time.time() - t0, 1), 'job': os.environ.get('W2JOB', 'v0_sample'),
               'note': 'random F_p point of the trivial-sector relaxation surface V0; full sector E1+E2 solve'}
        w2lib.log(**rec)
        print(rec, flush=True)

if __name__ == '__main__':
    main()
