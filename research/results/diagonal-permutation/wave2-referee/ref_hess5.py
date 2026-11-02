"""Degree-5 step: exact rational rank of the 21 second partials F_ij (i<=j) of F.
If rank = 21, then sum_i,j c_i c_j F_ij = 0 (as a quartic) forces c_i c_j = 0 for all i,j, i.e. c = 0."""
import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
from flint import fmpq_mat, fmpq
F, A, q = R.load_instance()
t0 = time.time()
parts = [R.pderiv(F, i) for i in range(6)]
second = {(i, j): R.pderiv(parts[i], j) for i in range(6) for j in range(i, 6)}
mons = R.monomials(4); idx = {m: k for k, m in enumerate(mons)}
rows = []
for (i, j), P in second.items():
    row = [fmpq(0)] * len(mons)
    for e, c in P.items():
        row[idx[e]] = fmpq(c.numerator, c.denominator)
    rows.append(row)
M = fmpq_mat(rows)
rk = M.rank()
# also: exact rank of the 6 first partials in degree 5 (sanity) 
mons5 = R.monomials(5); idx5 = {m: k for k, m in enumerate(mons5)}
rows5 = []
for P in parts:
    row = [fmpq(0)] * len(mons5)
    for e, c in P.items(): row[idx5[e]] = fmpq(c.numerator, c.denominator)
    rows5.append(row)
rk5 = fmpq_mat(rows5).rank()
rec = {'exp': 'second_partials_exact_rank', 'n_second_partials': len(second), 'rank_Q': int(rk), 'first_partials_rank_Q': int(rk5),
       'conclusion': ('rank 21 = all 21 quadratic monomials c_i c_j independent => w^T Hess(F) w = 0 forces w = 0' if rk == 21 else 'RANK DEFICIENT'),
       'seconds': round(time.time() - t0, 2), 'arithmetic': 'exact rational (python-flint fmpq_mat)'}
print(json.dumps(rec)); open('ref_hess5.json', 'w').write(json.dumps(rec, indent=1))
