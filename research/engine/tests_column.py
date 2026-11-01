import sys
sys.path.insert(0, '.')
from hodge import core, fixed_d
from fractions import Fraction
p = core.PRIMES[0]
F, A, q, rec = core.load_instance('contact46006')
d, D, N, Q = core.load_witness('contact46006')
S = fixed_d.FixedDSystem(F, A, q, D, p)
# pick a few N columns and compare with core.defect of the single-monomial field
bad = 0
for j in [0, 5, 100, 777, 1500]:
    kind, i, m = S.cols[j]
    if kind != 'N': continue
    Nj = [{} for _ in range(6)]; Nj[i][m] = 1
    E = core.defect(core.reduce_mod(F,p), core.reduce_mod(A,p), core.reduce_mod(D,p), q, Nj, p)
    # subtract A D^2 contribution
    AD2 = core.mul(core.reduce_mod(A,p), core.power(core.reduce_mod(D,p),2,p), p)
    E = core.sub(E, AD2, p)
    col = S.colpolys[j]
    diff = core.sub(E, col, p)
    print(j, kind, i, m, 'terms', len(E), len(col), 'diff terms', len(diff))
    bad += bool(diff)
print('BAD' if bad else 'columns OK')
