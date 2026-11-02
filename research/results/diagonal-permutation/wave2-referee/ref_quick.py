"""Quick fresh-prime checks: A not in J_6 (degree-0 denominators impossible); F has 24 terms; A's two terms."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
from flint import nmod_mat
p = 999999937
F, A, q = R.load_instance()
plus, minus = R.koszul_series(12)
parts = [R.pderiv(F, i) for i in range(6)]
partsp = [{e: R.modp(c, p) for e, c in P.items()} for P in parts]
Ap = {e: R.modp(c, p) for e, c in A.items()}
rows = R.monomials(6, 1); idx = {m: k for k, m in enumerate(rows)}
cols = []
for i in range(6):
    mpar = 1 if i < 3 else 0
    for m in R.monomials(1, mpar): cols.append(R.mul({m: 1}, partsp[i]))
def rank(cs):
    M = nmod_mat(len(rows), len(cs), p)
    for j, P in enumerate(cs):
        for e, c in P.items(): M[idx[e], j] = c % p
    return M.rank()
rJ = rank(cols); rA = rank(cols + [Ap])
rec = {'p': p, 'dim_S6_minus': len(rows), 'rank_J6_minus': rJ, 'ring_dim_6_minus': len(rows) - rJ, 'koszul': minus[6], 'A_in_J': rA == rJ, 'F_terms': len(F), 'A_terms': len(A)}
print(json.dumps(rec)); json.dump(rec, open('ref_quick.json', 'w'), indent=1)
