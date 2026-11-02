"""Usage: python3 batch.py <batchname>  -- runs a named list of denominators."""
import sys, time
import lib
from lib import *
name = sys.argv[1]
B = {}
# structural degree-6 pieces of F (expressions are strings or (poly, expr) pairs)
B['A'] = [
    (f, 'f(X)', 'structural', 'X-block; strict split denominator at t=0'),
    (g, 'g(Y)', 'structural', 'Y-block'),
    (RX, '6f+4tP (eps_X F)', 'structural', 'X-radial derivative of F'),
    (RY, '6g+2tP (eps_Y F)', 'structural', 'Y-radial derivative of F'),
    (core.sub(F, f), 'F-f = g+tP', 'structural', 'F minus X block'),
    (core.sub(F, g), 'F-g = f+tP', 'structural', 'F minus Y block'),
    (core.add(f, g), 'f+g', 'structural', 'split part of F'),
    (Pc, 'P = L^2 M', 'structural', 'coupling block'),
    (core.mul(core.mul(L, M), M), 'L M^2', 'structural', 'mirror coupling (2,4)'),
    (core.sub(f, g), 'f-g', 'structural', ''),
    (core.sub(RX, RY), 'eps_X F - eps_Y F = 6f-6g+2tP', 'structural', ''),
]
B['B'] = [
    (P('X0'), 'X0', 'cone-linear', ''), (P('Y0'), 'Y0', 'cone-linear', ''),
    (P('X1'), 'X1', 'cone-linear', ''), (P('X2'), 'X2', 'cone-linear', ''),
    (P('Y1'), 'Y1', 'cone-linear', ''), (P('Y2'), 'Y2', 'cone-linear', ''),
    (L, 'L', 'coupling-form', ''), (M, 'M', 'coupling-form', ''),
    (P('X0**2'), 'X0^2', 'cone', ''), (P('Y0**2'), 'Y0^2', 'cone', ''),
    (P('X0*Y0'), 'X0*Y0', 'cone', 'Y-odd D'),
    (P('X0**2+Y0**2'), 'X0^2+Y0^2', 'cone-quadric', ''),
    (P('X0**2-Y0**2'), 'X0^2-Y0^2', 'cone-quadric', ''),
    (core.mul(L, M), 'L*M', 'coupling-form', ''),
    (core.mul(L, L), 'L^2', 'coupling-form', ''), (core.mul(M, M), 'M^2', 'coupling-form', ''),
    (core.add(L, M), 'L+M', 'coupling-form', ''), (core.sub(L, M), 'L-M', 'coupling-form', ''),
    (P('X0**3'), 'X0^3', 'cone', ''), (P('Y0**3'), 'Y0^3', 'cone', ''),
    (P('X0**2*Y0'), 'X0^2*Y0', 'cone', 'Y-odd D'),
    (P('X0**3+Y0**3'), 'X0^3+Y0^3', 'cone-fermat', 'Y-odd D'),
    (P('X0**4+Y0**4'), 'X0^4+Y0^4', 'cone-fermat', ''),
    (P('X0**3*Y0**3'), 'A itself', 'cone', 'Y-odd D = A'),
    (P('X0**2*Y0**2'), 'X0^2 Y0^2', 'cone', ''),
    (core.mul(L, core.mul(L, L)), 'L^3', 'coupling-form', ''), (core.mul(M, core.mul(M, M)), 'M^3', 'coupling-form', ''),
    (core.mul(core.mul(L, L), core.mul(M, M)), 'L^2 M^2', 'coupling-form', 'd=8', ),
]
for D, expr, fam, note in B[name]:
    r = run_D(D, expr, fam, note)
    show(r)
