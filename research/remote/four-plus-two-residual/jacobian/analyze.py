"""Reconstruct L-polynomials of T_B, E', Prym P from twisted sums; run splitting tests."""
import json, glob, math
from sympy import symbols, Poly, factor_list, sqrt, QQ, discriminant, factorint, cyclotomic_poly, resultant, gcd, expand, Rational, nroots, I, re as sre, im as sim
T, S, x, zz = symbols('T S x z')
# Z[zeta] arithmetic: (a,b) = a + b zeta, zeta^2 = -1 - zeta
def zmul(u, v):
    a, b = u; c, d = v
    return (a*c - b*d, a*d + b*c - b*d)
def zadd(u, v): return (u[0]+v[0], u[1]+v[1])
def zneg(u): return (-u[0], -u[1])
def zconj(u): return (u[0]-u[1], -u[1])
def znorm(u): return u[0]**2 - u[0]*u[1] + u[1]**2
def zdiv(u, v):
    """exact division in Z[zeta]; returns None if not exact"""
    n = znorm(v); w = zmul(u, zconj(v))
    if w[0] % n or w[1] % n: return None
    return (w[0]//n, w[1]//n)
def zscal(c, u): return (c*u[0], c*u[1])
def ztoexpr(u):  # as sympy expr with zeta = (-1+sqrt(-3))/2
    return u[0] + u[1]*(Rational(-1,2) + sqrt(-3)/2)
def zint(u):
    assert u[1] == 0; return u[0]
def norm_p_elements(p):
    out = []
    for a in range(-2*int(math.isqrt(p))-2, 2*int(math.isqrt(p))+3):
        for b in range(-2*int(math.isqrt(p))-2, 2*int(math.isqrt(p))+3):
            if a*a - a*b + b*b == p: out.append((a, b))
    return out
def vp(u, pi, p):
    """valuation of u at the prime (pi) of Z[zeta] (pi of norm p)."""
    v = 0
    while u != (0, 0):
        w = zdiv(u, pi)
        if w is None: break
        u = w; v += 1
    return v

recs = {}
for fn in sorted(glob.glob('sums_*.jsonl')):
    for line in open(fn):
        r = json.loads(line); recs[r['p']] = r

def lpoly_from_cubic(b, p):
    """L_+(T) = 1 + b1 T + b2 T^2 + b3 T^3 over Z[zeta]; return full L_P(T) in Z[T] (degree 6) as sympy Poly."""
    b1, b2, b3 = b
    Lp = 1 + ztoexpr(b1)*T + ztoexpr(b2)*T**2 + ztoexpr(b3)*T**3
    Lm = 1 + ztoexpr(zconj(b1))*T + ztoexpr(zconj(b2))*T**2 + ztoexpr(zconj(b3))*T**3
    L = Poly(expand(Lp*Lm), T)
    assert all(c.is_integer for c in L.all_coeffs()), L
    return L

def rootunity_ratios(L, p, mmax=200):
    """Check whether any ratio of two distinct roots of the Weil polynomial (reverse of L) is a root of unity."""
    f = Poly(L.all_coeffs()[::-1], x)   # f(x) = x^6 L(1/x)*... careful: L(T)=prod(1-alpha T) -> reverse gives prod(x - alpha)
    f = Poly(expand(x**L.degree() * L.as_expr().subs(T, 1/x)), x)
    g = Poly(expand(f.as_expr().subs(x, zz*x)), x)   # f(z x)
    R = Poly(resultant(f.as_expr(), g.as_expr(), x), zz)   # = prod_{i,j} (z alpha_i - alpha_j) up to const
    # remove the trivial factor (z-1)^deg from i=j
    hits = []
    for m in range(2, mmax+1):
        c = Poly(cyclotomic_poly(m, zz), zz)
        if R.rem(c).is_zero: hits.append(m)
    return hits


def zprime_gens(ell):
    """generators of primes of Z[zeta] above rational prime ell"""
    if ell == 3: return [(1, -1)]            # 1 - zeta
    if ell % 3 == 2: return [(ell, 0)]
    g = [x for x in norm_p_elements(ell)]
    g0 = g[0]; return [g0, zconj(g0)]
def disc_class(b, p):
    b1, b2, b3 = b
    # disc of x^3 + b1 x^2 + b2 x + b3
    D = zadd(zadd(zmul(zmul(b1,b1), zmul(b2,b2)), zscal(-4, zmul(b2, zmul(b2,b2)))),
             zadd(zadd(zscal(-4, zmul(zmul(b1,b1), zmul(b1,b3))), zscal(-27, zmul(b3,b3))), zscal(18, zmul(b1, zmul(b2,b3)))))
    N = znorm(D)
    if N == 0: return {'D': D, 'square_class': 'zero'}
    odd = []; rest = D
    for ell, e in factorint(N).items():
        for g in zprime_gens(ell):
            v = vp(D, g, p)
            if v % 2: odd.append((ell, g))
            for _ in range(v): rest = zdiv(rest, g)
    # rest is a unit; its class mod squares: zeta^k are squares, -zeta^k are not
    unit_sign = 1 if rest in [(1,0),(0,1),(-1,-1)] else -1
    assert rest in [(1,0),(0,1),(-1,-1),(-1,0),(0,-1),(1,1)], rest
    return {'D': D, 'N(D)': factorint(N), 'odd_primes': odd, 'unit_sign': unit_sign}
def rootunity_numeric(L, mmax=2000):
    import mpmath
    mpmath.mp.dps = 60
    cs = [mpmath.mpf(int(c)) for c in L.all_coeffs()[::-1]]   # reversed: roots alpha
    try:
        rts = mpmath.polyroots(cs, maxsteps=500, extraprec=200); tol = mpmath.mpf(10)**-40
    except mpmath.libmp.libhyper.NoConvergence:
        import numpy as np
        rts = [mpmath.mpc(complex(r)) for r in np.roots([float(c) for c in L.all_coeffs()])]; tol = mpmath.mpf(10)**-7
    hits = set()
    for i in range(len(rts)):
        for j in range(len(rts)):
            if i == j: continue
            ang = mpmath.arg(rts[i]/rts[j]) / (2*mpmath.pi)
            for m in range(1, mmax+1):
                if abs(ang*m - mpmath.nint(ang*m)) < tol:
                    hits.add(m); break
    return sorted(hits)

results = []
for p in sorted(recs):
    r = recs[p]
    res = {'p': p, 'class': 'split' if p % 3 == 1 else 'inert'}
    if p % 3 == 1:
        S1B, S1E = tuple(r['TB_k1']['S']), tuple(r['E_k1']['S'])
        S2B, S2E = tuple(r['TB_k2']['S']), tuple(r['E_k2']['S'])
        om = r['TB_k1']['omega']; assert om == r['E_k1']['omega'] == r['TB_k2']['omega']
        piE = zneg(S1E)                         # Frobenius eigenvalue of E' on H_+
        assert znorm(piE) == p, (p, piE)
        assert zneg(S2E) == zmul(piE, piE), (p, S2E, piE)  # consistency of E' over F_{p^2}
        aE = zint(zadd(piE, zconj(piE)))
        assert r['E_k1']['npts'] == p + 1 - aE
        res['piE'] = piE; res['aE'] = aE; res['aE2-4p'] = aE*aE - 4*p
        # Prym traces
        t1 = zneg(zadd(S1B, zneg(S1E))); t2 = zneg(zadd(S2B, zneg(S2E)))
        b1 = zneg(t1); b2 = zdiv(zadd(zmul(t1, t1), zneg(t2)), (2, 0))
        assert b2 is not None
        res['t1'] = t1; res['t2'] = t2; res['b1'] = b1; res['b2'] = b2
        # b3 = p*b2/conj(b1)
        if b1 != (0, 0):
            b3 = zdiv(zscal(p, b2), zconj(b1))
            res['b3'] = b3
            res['b3_ok'] = (b3 is not None) and znorm(b3) == p**3
        else:
            b3 = None; res['b3'] = None; res['b3_ok'] = None
        if f'TB_k3' in r:
            S3B, S3E = tuple(r['TB_k3']['S']), tuple(r['E_k3']['S'])
            assert zneg(S3E) == zmul(piE, zmul(piE, piE))
            t3 = zneg(zadd(S3B, zneg(S3E)))
            # Newton: p3 = e1 p2 - e2 p1 + 3 e3, e1=-b1, e2=b2, e3=-b3
            e1, e2 = zneg(b1), b2
            e3 = zdiv(zadd(t3, zadd(zneg(zmul(e1, t2)), zmul(e2, t1))), (3, 0))
            b3_direct = zneg(e3)
            res['b3_from_Fp3'] = b3_direct
            res['b3_crosscheck'] = (b3 == b3_direct) if b3 is not None else None
            if b3 is None: b3 = b3_direct
        if b3 is None:
            results.append(res); continue
        # p-adic calibration: prime pp = (p, zeta - omega): generator pi with pi = a + b zeta, a + b*omega = 0 mod p
        gens = [g for g in norm_p_elements(p) if (g[0] + g[1]*om) % p == 0]
        pi = gens[0]
        res['v_pp(piE)'] = vp(piE, pi, p)
        res['newton_pp(L+P)'] = [vp(b1, pi, p), vp(b2, pi, p), vp(b3, pi, p)]
        res['newton_ppbar(L+P)'] = [vp(b1, zconj(pi), p), vp(b2, zconj(pi), p), vp(b3, zconj(pi), p)]
        # elliptic CM factor test: roots of L_+ in Z[zeta] of norm p  (1 - pi'' T divides L_+)
        Lplus = lambda u: zadd(zadd((1, 0), zmul(b1, u)), zadd(zmul(b2, zmul(u, u)), zmul(b3, zmul(u, zmul(u, u)))))
        cm_roots = []
        for g in norm_p_elements(p):
            # test L_+(1/g) = 0  <=>  g^3 + b1 g^2 + b2 g + b3 = 0
            val = zadd(zadd(zmul(g, zmul(g, g)), zmul(b1, zmul(g, g))), zadd(zmul(b2, g), b3))
            if val == (0, 0): cm_roots.append(g)
        res['cm_elliptic_roots_in_L+'] = cm_roots
        L = lpoly_from_cubic((b1, b2, b3), p)
        res['L_P'] = L.all_coeffs()[::-1]
        # Weil check numerically
        import numpy as np
        rts = np.roots([float(c) for c in L.all_coeffs()])
        res['weil_ok'] = bool(all(abs(abs(rt) - 1/math.sqrt(p)) < 1e-6 for rt in rts))
        fl = factor_list(L.as_expr(), T)
        res['factor_Q'] = [(str(f), e) for f, e in fl[1]]
        fl3 = factor_list(L.as_expr(), T, extension=sqrt(-3))
        res['factor_Qzeta'] = [(str(f), e) for f, e in fl3[1]]
        res['rootunity_hits'] = rootunity_ratios(L, p)
        res['rootunity_numeric'] = rootunity_numeric(L)
        res['disc_class'] = disc_class((b1, b2, b3), p)
        # Frobenius field invariant: squarefree part of disc of the reversed sextic (if irreducible)
        if len(fl[1]) == 1 and fl[1][0][1] == 1 and Poly(fl[1][0][0], T).degree() == 6:
            d = discriminant(Poly(L.all_coeffs()[::-1], x).as_expr(), x)
            fi = factorint(d)
            sf = 1
            for q, e in fi.items():
                if e % 2: sf *= q
            res['disc_sqfree'] = sf
        # full T_B count checks
        nB = r['TB_k1']['npts']; SB = S1B
        assert nB == p + 1 + 2*SB[0] - SB[1]
        res['#TB(Fp)'] = nB; res['#E(Fp)'] = r['E_k1']['npts']
    else:
        S2B, S2E = tuple(r['TB_k2']['S']), tuple(r['E_k2']['S'])
        res['#E(Fp2)'] = r['E_k2']['npts']; res['#TB(Fp2)'] = r['TB_k2']['npts']
        res['E_supersingular_check'] = (r['E_k2']['npts'] == p*p + 1 + 2*p) and (zneg(S2E) == (-p, 0))
        t2 = zneg(zadd(S2B, zneg(S2E)))
        res['t2'] = t2
        res['t2_real'] = (t2[1] == 0)
        b = -t2[0]; c = b - p
        res['M(S)'] = f'(1 + {p} S)(1 + ({c}) S + {p*p} S^2)'
        res['c'] = c; res['c_over_p'] = Rational(c, p)
        res['second_factor_j0_type'] = c in {0, p, -p, 2*p, -2*p}
        res['second_factor_ordinary'] = (c % p != 0)
        if 'TB_k4' in r:
            n4 = r['TB_k4']['npts']
            pred = p**4 + 1 - 2*(p - b)**2
            res['Fp4_check'] = (n4, pred, n4 == pred)
    results.append(res)

json.dump(results, open('analysis.json', 'w'), indent=1, default=str)
for res in results:
    print(json.dumps(res, default=str))
