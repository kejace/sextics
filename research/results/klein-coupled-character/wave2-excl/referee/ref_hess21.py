"""Exact rank of the 21 second partials of 7F (degree-5 closure input) and exact dim J_4 = 0 (nothing of degree 4 in J)."""
import sys
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from fractions import Fraction
F, Fi, H = L.Fpolys(None)
mons = L.monomials(4); idx = {e: r for r, e in enumerate(mons)}
ent = {}; c = 0; pairs = []
for i in range(6):
    for j in range(i, 6):
        for e, v in H[i][j].items(): ent[(idx[e], c)] = Fraction(v)
        pairs.append((i, j)); c += 1
M = L.mat(len(mons), 21, ent, None)
print('exact rank of the 21 quartics d_i d_j (7F):', M.rank(), 'of 21')
print('HESS21_DONE')
