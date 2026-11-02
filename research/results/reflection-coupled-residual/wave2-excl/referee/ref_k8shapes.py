"""Referee: does the affine family K8 + J_8 contain a member of a monomial-defined shape U?  Exact over Q (fmpq_mat):
K8 in U + J_8  <=>  rank[J_8^(00) gens; U monomials; K8] == rank[J_8^(00) gens; U monomials]."""
import sys, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
from flint import fmpq_mat, fmpq
from fractions import Fraction
F, A, q = R.load_instance()
K8 = R.parse(json.load(open('/tmp/hodge/results/reflection-coupled-residual/wave2-excl/K8_E1_exact.json'))['D'])
mons = R.monomials(8, (0, 0)); idx = {e: i for i, e in enumerate(mons)}; C = len(mons)
Fi = [F.deriv(i) for i in range(6)]
Jrows = []
for i in range(6):
    for m in R.monomials(3, R.cadd((0, 0), R.CHI_Z[i])):
        Jrows.append({tuple(a + b for a, b in zip(e, m)): c for e, c in Fi[i].t.items()})
def rank_of(rows):
    M = fmpq_mat(len(rows), C)
    for r, P in enumerate(rows):
        for e, c in P.items(): M[r, idx[e]] = fmpq(c.numerator, c.denominator)
    return M.rank()
rJ = rank_of(Jrows)
shapes = {
    'all Y-exponents even (Y0,Y1,Y2 even; X arbitrary)': lambda e: e[3] % 2 == 0 and e[4] % 2 == 0 and e[5] % 2 == 0,
    'polynomial in squares (all exponents even)': lambda e: all(x % 2 == 0 for x in e),
    'independent of X0,Y2 (cone with line vertex)': lambda e: e[0] == 0 and e[5] == 0,
    'multiple of X0': lambda e: e[0] >= 1,
    'Y-degree <= 4': lambda e: e[3] + e[4] + e[5] <= 4,
    'X-degree <= 4': lambda e: e[0] + e[1] + e[2] <= 4,
}
out = {'rank_J8_00': rJ}
for name, pred in shapes.items():
    U = [{e: Fraction(1)} for e in mons if pred(e)]
    r1 = rank_of(Jrows + U); r2 = rank_of(Jrows + U + [K8.t])
    out[name] = dict(dimU=len(U), rank_J_plus_U=r1, K8_in_U_plus_J=(r2 == r1))
    print(name, out[name], flush=True)
R.receipt('ref_k8shapes.jsonl', dict(exp='ref-k8-shapes-exact-Q', **{k: v for k, v in out.items()}))
