#!/usr/bin/env python3
"""Structured candidate denominators suggested by the lane's geometry: E1 screen, then exact
sector-restricted fixed-D solve (mod p) through the wave-1 sector solver (common.test_D), logged to
log.jsonl AND summarized into wave2-search.log.jsonl.

  python3 w2_struct.py <degree> [filter]

Families (all built from the block structure S0 = X0^2+Y0^2, P0 = X0Y0, delta = X0^2-Y0^2, the contour
quartics 2G_S +- G_P, the r1-analogue contour (F_X1^2 - F_Y1^2)/(X1^2-Y1^2), block radial derivatives,
Hessian minors, polars of the fixed points of the reflections, and products with invariant quadrics/quartics).
Mixed-parity candidates are replaced by D*tau(D).  Non-eigen (character-mixed) pure-parity candidates go
through the generic engine (fixed_d.solve_fixed_D_fast) when d <= 10.
"""
import sys, os, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from common import *          # F, A, q, chi, chi_of, test_D, log (log.jsonl), X0..Y2, Z, lin, dir_deriv, exact_div
import w2lib as W
import sympy as sp

d_want = int(sys.argv[1])
only = sys.argv[2] if len(sys.argv) > 2 else None

def sq(P): return mul(P, P)
def tau(P):
    return {e: (c if (e[3] + e[4] + e[5]) % 2 == 0 else -c) for e, c in P.items()}

S0 = add(mul(X0, X0), mul(Y0, Y0)); P0 = mul(X0, Y0); delta = sub(mul(X0, X0), mul(Y0, Y0))
S1 = add(mul(X1, X1), mul(Y1, Y1)); P1 = mul(X1, Y1); d1 = sub(mul(X1, X1), mul(Y1, Y1))
S2 = add(mul(X2, X2), mul(Y2, Y2)); P2 = mul(X2, Y2); d2 = sub(mul(X2, X2), mul(Y2, Y2))
Q1 = add(add(sq(P1), mul(P1, P2)), sq(P2))
Q2 = mul(P1, add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
GS = add(add(scale(sq(S0), 3), scale(sq(P0), 2)), scale(Q1, 5))
GP = add(scale(mul(P0, S0), 4), scale(Q2, Fraction(5, 7)))
D8 = sub(scale(sq(GS), 4), sq(GP))
Fd = [pderiv(F, i) for i in range(NV)]
Hess = [[pderiv(Fd[i], j) for j in range(NV)] for i in range(NV)]
# r1-analogue contour: (F_X1^2 - F_Y1^2)/(X1^2 - Y1^2)
D8r = exact_div(sub(sq(Fd[1]), sq(Fd[4])), d1)
# X2/Y2 "pseudo-contour" (no exchange symmetry, but same algebra): (F_X2^2 - F_Y2^2) may not be divisible
try:
    D8s = exact_div(sub(sq(Fd[2]), sq(Fd[5])), d2)
except AssertionError:
    D8s = None
rad = [add(mul(Z[i], Fd[i]), mul(Z[i + 3], Fd[i + 3])) for i in range(3)]   # block radial derivatives (deg 6)
# polar quartics at the s0-fixed points p+- = [1:0:0:+-1:0:0] and the r1-fixed points [0:1:0:0:+-1:0]
pol_s0 = [exact_div(add(Fd[0], Fd[3]), add(X0, Y0)), exact_div(sub(Fd[0], Fd[3]), sub(X0, Y0))]
pol_r1 = [exact_div(add(Fd[1], Fd[4]), add(X1, Y1)), exact_div(sub(Fd[1], Fd[4]), sub(X1, Y1))]
# Hessian minors
def minor(idx):
    M = sp.Matrix(len(idx), len(idx), lambda a, b: 0)
    xs = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
    def to_sp(P):
        return sum(sp.Rational(c.numerator, c.denominator) * sp.prod([x ** k for x, k in zip(xs, e)]) for e, c in P.items())
    M = sp.Matrix([[to_sp(Hess[i][j]) for j in idx] for i in idx])
    det = sp.Poly(sp.expand(M.det()), *xs)
    out = {}
    for mon, c in det.terms():
        out[tuple(int(k) for k in mon)] = Fraction(int(c.p), int(c.q))
    return out

cands = {}
quads = {'S0': S0, 'P0^2(deg4)': None, 'X0^2': mul(X0, X0), 'Y0^2': mul(Y0, Y0), 'S1': S1, 'S2': S2, 'X1^2': mul(X1, X1),
         'Y1^2': mul(Y1, Y1), 'X2^2': mul(X2, X2), 'Y2^2': mul(Y2, Y2), 'P1': P1, 'P2': P2, 'delta': delta, 'd1': d1, 'd2': d2}
quads = {k: v for k, v in quads.items() if v is not None}
quartics = {'S0^2': sq(S0), 'P0^2': sq(P0), 'delta^2': sq(delta), 'GS': GS, 'GP': GP, 'Q1': Q1, 'Q2': Q2, 'S0S1': mul(S0, S1),
            'S0S2': mul(S0, S2), 'S1S2': mul(S1, S2), 'P1^2': sq(P1), 'P2^2': sq(P2), 'P0P1': mul(P0, P1), 'P0P2': mul(P0, P2),
            'X0^4+Y0^4': add(power(X0, 4), power(Y0, 4)), 'X1^4+Y1^4': add(power(X1, 4), power(Y1, 4)), 'X2^4+Y2^4': add(power(X2, 4), power(Y2, 4)),
            'pol_s0+': pol_s0[0], 'pol_s0-': pol_s0[1], 'pol_r1+': pol_r1[0], 'pol_r1-': pol_r1[1]}
if d_want == 8:
    cands.update({'D8r (r1-contour)': D8r, 'D8s (X2Y2 pseudo-contour)': D8s, 'pol_r1+ * pol_r1-': mul(pol_r1[0], pol_r1[1]),
                  'pol_s0+ * pol_r1+': mul(pol_s0[0], pol_r1[0]), 'pol_s0+ * pol_r1-': mul(pol_s0[0], pol_r1[1]),
                  'pol_s0+^2': sq(pol_s0[0]), 'pol_r1+^2': sq(pol_r1[0]),
                  'minor(X0,Y0)': minor([0, 3]), 'minor(X1,Y1)': minor([1, 4]), 'minor(X2,Y2)': minor([2, 5]),
                  'minor(X0,X1)': minor([0, 1]), 'minor(X0,Y1)': minor([0, 4]), 'minor(Y0,X1)': minor([3, 1]),
                  'minor(X0,X2)': minor([0, 2]), 'minor(X1,X2)': minor([1, 2]),
                  'rad0*S0': mul(rad[0], S0), 'rad0*delta': mul(rad[0], delta), 'rad0*X0^2': mul(rad[0], mul(X0, X0)),
                  'rad1*S1': mul(rad[1], S1), 'rad1*d1': mul(rad[1], d1), 'rad2*S2': mul(rad[2], S2), 'rad1*S0': mul(rad[1], S0),
                  'rad0*S1': mul(rad[0], S1), 'rad0*S2': mul(rad[0], S2), 'rad0*P1': mul(rad[0], P1), 'rad1*P0^?': None,
                  'GS*S0^2': mul(GS, sq(S0)), 'GS*delta^2': mul(GS, sq(delta)), 'GS*P0^2': mul(GS, sq(P0)), 'GP*P0*S0': mul(GP, mul(P0, S0)),
                  'GS*Q1': mul(GS, Q1), 'GS*GP?': None, 'D8*': D8})
    for qa, qb in itertools.combinations(sorted(quartics), 2):
        cands[f'{qa}*{qb}'] = mul(quartics[qa], quartics[qb])
elif d_want == 10:
    for qn, qp in quads.items():
        cands[f'D8* * {qn}'] = mul(D8, qp)
        cands[f'D8r * {qn}'] = mul(D8r, qp)
    cands.update({'rad0*GS': mul(rad[0], GS), 'rad0*GP': mul(rad[0], GP), 'rad1*GS': mul(rad[1], GS), 'rad2*GS': mul(rad[2], GS),
                  'rad0*S0^2': mul(rad[0], sq(S0)), 'rad0*delta^2': mul(rad[0], sq(delta)), 'rad0*Q1': mul(rad[0], Q1),
                  'rad0*pol_s0+': mul(rad[0], pol_s0[0]), 'rad1*pol_r1+': mul(rad[1], pol_r1[0]),
                  'minor(X0,Y0)*S0': mul(minor([0, 3]), S0), 'minor(X0,Y0)*delta': mul(minor([0, 3]), delta),
                  'minor(X1,Y1)*S1': mul(minor([1, 4]), S1), 'minor(X0,Y0)*S1': mul(minor([0, 3]), S1),
                  'minor3(X0,Y0,X1)*?': None, 'FX0^2+FY0^2': add(sq(Fd[0]), sq(Fd[3])), 'FX1^2+FY1^2': add(sq(Fd[1]), sq(Fd[4])),
                  'FX2^2+FY2^2': add(sq(Fd[2]), sq(Fd[5])), 'FX0^2-FY0^2 (=delta*D8*)': sub(sq(Fd[0]), sq(Fd[3])),
                  'sum Fi^2': add(add(add(sq(Fd[0]), sq(Fd[3])), add(sq(Fd[1]), sq(Fd[4]))), add(sq(Fd[2]), sq(Fd[5]))),
                  'FX0FY0?': mul(Fd[0], Fd[3]), 'FX1FY1': mul(Fd[1], Fd[4]), 'FX0FX1+FY0FY1': add(mul(Fd[0], Fd[1]), mul(Fd[3], Fd[4]))})
elif d_want == 12:
    cands.update({'D8* * S0^2': mul(D8, sq(S0)), 'D8* * P0^2': mul(D8, sq(P0)), 'D8* * delta^2': mul(D8, sq(delta)),
                  'D8* * GS': mul(D8, GS), 'D8* * GP': mul(D8, GP), 'D8* * Q1': mul(D8, Q1), 'D8* * S1^2': mul(D8, sq(S1)),
                  'D8* * D8r (deg16!)': None, 'D8r * S1^2': mul(D8r, sq(S1)), 'D8r * GS': mul(D8r, GS), 'D8r * P1^2': mul(D8r, sq(P1)),
                  'D8r * d1^2': mul(D8r, sq(d1)), 'rad0^2': sq(rad[0]), 'rad0*rad1': mul(rad[0], rad[1]), 'rad1^2': sq(rad[1]),
                  'rad0*rad2': mul(rad[0], rad[2]), 'minor3(X0,Y0,X1)': None, 'minor(X0,Y0)*minor(X1,Y1)': mul(minor([0, 3]), minor([1, 4])),
                  'minor(X0,Y0)^... S0^2': mul(minor([0, 3]), sq(S0)), 'minor(X0,Y0)*GS': mul(minor([0, 3]), GS),
                  'det Hess(X0,Y0,X1) (deg12)': minor([0, 3, 1]), 'det Hess(X0,Y0,Y1)': minor([0, 3, 4]), 'det Hess(X0,Y0,X2)': minor([0, 3, 2]),
                  'det Hess(X1,Y1,X2)': minor([1, 4, 2]), 'det Hess(X0,X1,X2)': minor([0, 1, 2]), 'det Hess(X0,Y1,Y2)': minor([0, 4, 5]),
                  '(FX0^2+FY0^2)*S0': mul(add(sq(Fd[0]), sq(Fd[3])), S0), '(FX0^2-FY0^2)*delta': mul(sub(sq(Fd[0]), sq(Fd[3])), delta),
                  '(FX1^2-FY1^2)*d1': mul(sub(sq(Fd[1]), sq(Fd[4])), d1), 'pol_s0+^3': power(pol_s0[0], 3),
                  'pol_s0+^2*pol_s0-': mul(sq(pol_s0[0]), pol_s0[1]), 'pol_r1+^2*pol_r1-': mul(sq(pol_r1[0]), pol_r1[1])})
cands = {k: v for k, v in cands.items() if v is not None}

done = 0
for name, D in cands.items():
    if only and only not in name:
        continue
    if not D:
        continue
    dd = core.degree(D)
    if dd != d_want or not core.is_homogeneous(D):
        print('skip (deg/hom)', name, dd, file=sys.stderr); continue
    if core.parity_of(D) is None:
        D = mul(D, tau(D)); name = name + ' *tau'
        dd = core.degree(D)
        if dd > 12:
            print('skip (tau-doubled degree too big)', name, file=sys.stderr); continue
    t0 = time.time()
    r = test_D(D, f'w2struct-d{dd}', note=name, verbose=True)
    done += 1
    W.log_line({'task': 'struct', 'd': dd, 'name': name, 'chi': None if chi_of(D) is None else W.chi_str(chi_of(D)),
                'e1': r.get('e1'), 'consistent': r.get('consistent'), 'rank': r.get('rank'), 'nullity': r.get('nullity'),
                'rows': r.get('rows'), 'cols': r.get('cols'), 'seconds': round(time.time() - t0, 1),
                'WITNESS': bool(r.get('consistent'))})
    print(name, '->', {k: v for k, v in r.items() if k not in ('system', 'res', 'out', 'D')}, f'({time.time()-t0:.0f}s)', flush=True)
    if r.get('consistent'):
        print('!!! CONSISTENT candidate', name, '-> exact re-solve', flush=True)
        status, data = W.exact_resolve(D, verbose=True)
        W.log_line({'task': 'struct-exact', 'd': dd, 'name': name, 'status': status,
                    'verified': data.get('verified'), 'coprime': data.get('coprime')})
        if status == 'solved' and data.get('verified'):
            W.save_witness_candidate(D, data['N'], data['Q'], extra={'coprime': data.get('coprime'), 'tag': name})
            print('WITNESS FOUND', name, flush=True)
print('done', done, 'candidates at degree', d_want)
