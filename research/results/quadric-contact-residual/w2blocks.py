"""Exact building blocks for quadric-contact-residual (copied from families.py without the test loop)."""
from common import *
S0 = add(mul(X0, X0), mul(Y0, Y0)); P0 = mul(X0, Y0); delta = sub(mul(X0, X0), mul(Y0, Y0))
S1 = add(mul(X1, X1), mul(Y1, Y1)); P1 = mul(X1, Y1); S2 = add(mul(X2, X2), mul(Y2, Y2)); P2 = mul(X2, Y2)
Q1 = add(add(mul(P1, P1), mul(P1, P2)), mul(P2, P2))          # B = P1^2 + P1 P2 + P2^2
Q2 = mul(P1, add(mul(X2, X2), scale(mul(Y2, Y2), 2)))          # P1 (X2^2 + 2 Y2^2)
R16 = add(add(power(X1, 6), power(X2, 6)), add(power(Y1, 6), power(Y2, 6)))
f0 = add(power(S0, 3), scale(mul(mul(P0, P0), S0), 2))         # X0^6+5X0^4Y0^2+5X0^2Y0^4+Y0^6
GS = add(add(scale(mul(S0, S0), 3), scale(mul(P0, P0), 2)), scale(Q1, 5))
GP = add(scale(mul(P0, S0), 4), scale(Q2, Fraction(5, 7)))
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
T0 = add(power(X0, 4), power(Y0, 4)); T2 = add(power(X2, 4), power(Y2, 4)); T1 = add(power(X1, 4), power(Y1, 4))
rad0 = add(mul(X0, Fd[0]), mul(Y0, Fd[3])); rad1 = add(mul(X1, Fd[1]), mul(Y1, Fd[4])); rad2 = add(mul(X2, Fd[2]), mul(Y2, Fd[5]))
