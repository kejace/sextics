"""Structured Gamma-invariant (or G_aff-eigen) candidates; E1 filter then full sector solve. Usage: families.py <degree>"""
import sys, json, time
from common import *
import sympy as sp

d_want = int(sys.argv[1]) if len(sys.argv) > 1 else 10
only = sys.argv[2] if len(sys.argv) > 2 else None

# building blocks (exact rationals)
S = add(mul(X0, X0), mul(Y0, Y0)); Pm = mul(X0, Y0); delta = sub(mul(X0, X0), mul(Y0, Y0))
Q1 = add(add(mul(mul(X1, X1), mul(Y1, Y1)), mul(mul(X1, X2), mul(Y1, Y2))), mul(mul(X2, X2), mul(Y2, Y2)))
Q2 = mul(mul(X1, Y1), add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
R = add(add(power(X1, 6), power(X2, 6)), add(power(Y1, 6), power(Y2, 6)))
f0 = add(power(S, 3), scale(mul(mul(Pm, Pm), S), 2))
GS = add(add(scale(mul(S, S), 3), scale(mul(Pm, Pm), 2)), scale(Q1, 5))
GP = add(scale(mul(Pm, S), 4), scale(Q2, Fraction(5, 7)))
D8 = sub(scale(mul(GS, GS), 4), mul(GP, GP))
D8a = sub(scale(GS, 2), GP); D8b = add(scale(GS, 2), GP)
Fd = [pderiv(F, i) for i in range(NV)]
Hess = [[pderiv(Fd[i], j) for j in range(NV)] for i in range(NV)]
def minor2(i, j):
    return sub(mul(Hess[i][i], Hess[j][j]), mul(Hess[i][j], Hess[i][j]))
def sq(P): return mul(P, P)
def sumsq(*Ps):
    out = {}
    for P in Ps: out = add(out, mul(P, P))
    return out
S2 = add(mul(X2, X2), mul(Y2, Y2)); S1 = add(mul(X1, X1), mul(Y1, Y1))
T0 = add(power(X0, 4), power(Y0, 4)); T2 = add(power(X2, 4), power(Y2, 4)); P1 = mul(X1, Y1); P2 = mul(X2, Y2)
rad0 = add(mul(X0, Fd[0]), mul(Y0, Fd[3])); rad1 = add(mul(X1, Fd[1]), mul(Y1, Fd[4])); rad2 = add(mul(X2, Fd[2]), mul(Y2, Fd[5]))
quads_inv = {'X0^2': mul(X0, X0), 'Y0^2': mul(Y0, Y0), 'X2^2': mul(X2, X2), 'Y2^2': mul(Y2, Y2), 'S0': S, 'S2': S2, 'X0^2+X2^2': add(mul(X0, X0), mul(X2, X2)), 'S0+S2': add(S, S2), 'S0-S2': sub(S, S2)}
quartics_inv = {'S0^2': sq(S), 'P0^2': sq(Pm), 'T0': T0, 'S2^2': sq(S2), 'T2': T2, 'P2^2': sq(P2), 'P1^2': sq(P1), 'P1P2': mul(P1, P2), 'P0P1': mul(Pm, P1), 'P0P2': mul(Pm, P2), 'Q1': Q1, 'S0S2': mul(S, S2), 'delta^2': sq(delta), 'GS': GS, 'GP': GP}

cands = {}
if d_want == 8:
    cands.update({'D8*': D8, 'D8a': D8a, 'D8b': D8b, 'G_S': GS, 'G_P': GP, 'minor2(X0,Y0)': minor2(0, 3), 'minor2(X1,Y1)': minor2(1, 4),
                  'minor2(X2,Y2)': minor2(2, 5), 'FX0^2+FY0^2 (deg10!)': None, 'f0*X0^2': None,
                  'GS^2': sq(GS), 'GP^2': sq(GP), 'f0*S0': mul(f0, S), 'S0^4': power(S, 4), 'P0^4': power(Pm, 4), 'rad0*S0': mul(rad0, S), 'rad0*X0^2': mul(rad0, mul(X0, X0)),
                  'rad1*S0': mul(rad1, S), 'rad2*S0': mul(rad2, S), 'rad0*S2': mul(rad0, S2), 'f0*S2': mul(f0, S2), 'R*S0': mul(R, S), '(X2^6+Y2^6)*S0': mul(add(power(X2, 6), power(Y2, 6)), S),
                  'S0^2*S2^2': mul(sq(S), sq(S2)), 'f0^2/S0^2=(S0^2+2P0^2)^2': sq(add(sq(S), scale(sq(Pm), 2)))})
    cands = {k: v for k, v in cands.items() if v is not None}
elif d_want == 10:
    for qn, qp in quads_inv.items():
        cands[f'D8* * {qn}'] = mul(D8, qp)
    cands.update({'FX0^2+FY0^2': sumsq(Fd[0], Fd[3]), 'FX1^2+FY1^2': sumsq(Fd[1], Fd[4]), 'FX2^2+FY2^2': sumsq(Fd[2], Fd[5]),
                  'sum all Fi^2': sumsq(*Fd), 'FX0^2+FY0^2 - (FX2^2+FY2^2)': sub(sumsq(Fd[0], Fd[3]), sumsq(Fd[2], Fd[5])),
                  'X0^5FX0+Y0^5FY0': add(mul(power(X0, 5), Fd[0]), mul(power(Y0, 5), Fd[3])),
                  'X2^5FX2+Y2^5FY2': add(mul(power(X2, 5), Fd[2]), mul(power(Y2, 5), Fd[5])),
                  'rad0*S0^2': mul(rad0, sq(S)), 'rad0*T0': mul(rad0, T0), 'rad0*P0^2': mul(rad0, sq(Pm)), 'rad0*Q1': mul(rad0, Q1), 'rad0*S2^2': mul(rad0, sq(S2)),
                  'rad1*S0^2': mul(rad1, sq(S)), 'rad2*S0^2': mul(rad2, sq(S)), 'rad0*GS': mul(rad0, GS), 'rad0*GP': mul(rad0, GP), 'rad1*GS': mul(rad1, GS), 'rad2*GS': mul(rad2, GS),
                  'D8a*D8b*S0 (=D8*S0)': None, 'minor2(X0,Y0)*S0': mul(minor2(0, 3), S), 'minor2(X0,Y0)*S2': mul(minor2(0, 3), S2),
                  'GS*rad0': None, 'GS*(X0FX0-Y0FY0)??': None})
    cands = {k: v for k, v in cands.items() if v is not None}
elif d_want == 12:
    cands.update({'D8*^2': None, 'D8* * S0^2': mul(D8, sq(S)), 'D8* * P0^2': mul(D8, sq(Pm)), 'D8* * T0': mul(D8, T0), 'D8* * S2^2': mul(D8, sq(S2)), 'D8* * Q1': mul(D8, Q1),
                  'D8* * GS': mul(D8, GS), 'D8* * GP': mul(D8, GP), 'D8* * P1^2': mul(D8, sq(P1)), 'D8* * delta^2': mul(D8, sq(delta)),
                  '(FX0^2+FY0^2)*S0': mul(sumsq(Fd[0], Fd[3]), S), 'rad0^2': sq(rad0), 'rad0*rad1': mul(rad0, rad1), 'rad0*rad2': mul(rad0, rad2),
                  'FX0^2FY0^2?? (deg20)': None, 'f0*D8a? mixed': None})
    cands = {k: v for k, v in cands.items() if v is not None}

for name, D in cands.items():
    if only and only not in name:
        continue
    d = core.degree(D)
    if d != d_want:
        print('skip (deg)', name, d, file=sys.stderr); continue
    if not core.is_homogeneous(D):
        print('skip (inhom)', name, file=sys.stderr); continue
    t0 = time.time()
    r = test_D(D, f'family-d{d}', note=name, verbose=True)
    print(name, '->', {k: v for k, v in r.items() if k not in ('system', 'res', 'out')}, f'({time.time()-t0:.0f}s)', flush=True)
