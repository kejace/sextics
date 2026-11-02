#!/usr/bin/env python3
"""Small exact re-verifications (local, seconds):
 (1) A = X0^3 Y0^3 not in J(F)_6  (new prime 2147483587): rank[J_6^- | A] > rank[J_6^-].
 (2) Degree-5 Hessian kernel: the linear map Sym^2(C^6) -> S_4,  S -> sum_{i,j} S_ij F_ij, is injective (exact over Q).
     Consequence (with (J:A)_4 = 0 in both parities from the kernel sweep): a degree-5 denominator D = v.grad F
     (the only E1-admissible ones, any direction v in C^6, mixed parity allowed) needs A v^T Hess(F) v in J_10,
     i.e. v^T Hess(F) v = 0 as a polynomial, i.e. v v^T in the kernel => v = 0. So deg D = 5 is impossible for every v,
     not only the six coordinate directions."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
from fractions import Fraction
import flint
F, A, q = load()
p = 2147483587
t0 = time.time()
ok, info = fixed_d.e1_solvable(F, A, q, core.parse('1'), p=p)
rec = log2(exp='A-not-in-J6', prime=p, result='A in J' if ok else 'A NOT in J_6 (rigorous: modular rank of [J_6^-|A] exceeds rank of J_6^-, and rank J_6^- = 18 attains the Koszul count)', **{k: v for k, v in info.items() if isinstance(v, (int, bool))}, seconds=round(time.time()-t0, 2))
print(rec)
t0 = time.time()
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
mons4 = monomials(4); idx = {m: j for j, m in enumerate(mons4)}
pairs = [(i, j) for i in range(NV) for j in range(i, NV)]
M = flint.fmpq_mat(len(mons4), len(pairs))
for c, (i, j) in enumerate(pairs):
    P_ = Fij[i][j] if i == j else scale(Fij[i][j], 2)
    for m, v in P_.items():
        fr = Fraction(v); M[idx[m], c] = flint.fmpq(fr.numerator, fr.denominator)
rk = M.rank()
rec = log2(exp='degree5-hessian-kernel', method='exact fmpq_mat rank', rows=len(mons4), cols=len(pairs), rank=rk,
           result='injective (kernel 0): no degree-5 denominator for any direction v' if rk == len(pairs) else f'KERNEL of dim {len(pairs)-rk}',
           seconds=round(time.time()-t0, 2))
print(rec)
