#!/usr/bin/env python3
"""Degree-5 closure test, any parity: if the E1 kernel mod J vanishes in degree 5 (both parities) and [A] is injective on R_4
(E1 kernel 0 at d=4), every E1-admissible degree-5 D is D = w(F) with w constant and E2 mod J forces w^T Hess(F) w = 0 as a
quartic.  If the 21 second partials F_ij are linearly independent, this forces all w_i w_j = 0, i.e. w = 0: no degree-5
denominator of any parity.  Prints the rank of the 21 quartics (exact over Q via fmpq_mat) per lane."""
import sys, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from flint import fmpq_mat, fmpq
from fractions import Fraction
lanes = sys.argv[1:] or ['cycle-q2','cyclic-triple-residual','dense-rational','diagonal-permutation','eight-plane-configuration','four-plus-two-residual','klein-coupled-character','off-paired','one-reflection-skew','polarized-k3-remainder','quadratic-cubic-pullback','quadric-contact-residual','rank-one-biquartic','reflection-coupled-residual','short-power-presentation','symmetric-sextic-residual']
for ln in lanes:
    F, A, q, rec = core.load_instance(ln)
    Fij = []
    for i in range(6):
        Fi = core.pderiv(F, i)
        for j in range(i, 6):
            Fij.append(core.pderiv(Fi, j))
    mons = sorted({e for P in Fij for e in P})
    idx = {e: t for t, e in enumerate(mons)}
    M = fmpq_mat(len(Fij), len(mons))
    for r, P in enumerate(Fij):
        for e, c in P.items():
            c = Fraction(c); M[r, idx[e]] = fmpq(c.numerator, c.denominator)
    print(json.dumps({'lane': ln, 'rank_second_partials': M.rank(), 'n': 21}))
