#!/usr/bin/env python3
"""Exact replays: (a) the 56 cubics in Q_0..Q_5 meet J_6^even exactly in C.F; (b) projective diagonal symmetry group of F
via the Smith normal form of the exponent-difference lattice (own code, sympy)."""
import sys, os, itertools, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from hodge import core
import referee_replay as R
t0 = time.time()
# (a)
cub = []
for e in itertools.combinations_with_replacement(range(6), 3):
    P_ = {(0, 0, 0, 0, 0, 0): Fraction(1)}
    for i in e: P_ = core.mul(P_, R.Qs[i])
    cub.append(P_)
cols = R.monos(6, 0)
gens = [core.shift(R.Fi[i], m, 1) for i in range(6) for m in R.monos(1, (0 - R.PARZ[i]) % 2)]
rC = R.qmat(cub, cols).rank(); rJ = R.qmat(gens, cols).rank(); rCJ = R.qmat(cub + gens, cols).rank()
Fin = R.qmat(cub + [R.F], cols).rank()
R.log(dict(check='cubics_in_Q_meet_J6_exact', exact=True, n_cubics=56, rank_cubics=rC, rank_J6even=rJ, rank_sum=rCJ,
           intersection_dim=rC + rJ - rCJ, F_in_span_of_cubics=(Fin == rC), conclusion='intersection = C.F' if rC + rJ - rCJ == 1 else 'unexpected',
           seconds=round(time.time() - t0, 1)))
# (b) diagonal symmetries: t in (C*)^6 with t^e constant over the support of F  <=>  characters trivial on the lattice L of differences
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
exps = list(R.F.keys()); e0 = exps[0]
L = Matrix([[a - b for a, b in zip(e, e0)] for e in exps[1:]])
S = smith_normal_form(L, domain=ZZ)
diag = [S[i, i] for i in range(min(S.shape))]
nz = [abs(int(d)) for d in diag if d != 0]
torsion = [d for d in nz if d > 1]
rank = len(nz)
R.log(dict(check='diagonal_symmetry_group_SNF', exact=True, lattice_rank=rank, smith_invariants=nz,
           note='group of diagonal symmetries modulo scalars = Hom(Z^6/L, C*) / scalars; rank 5 lattice => finite group, order = product of invariants > 1',
           finite_group_order=(1 if not torsion else eval('*'.join(map(str, torsion)))), conclusion=('Z/2' if torsion == [2] and rank == 5 else 'see invariants')))
