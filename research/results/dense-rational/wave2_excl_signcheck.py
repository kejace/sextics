"""Adversarial check of the relaxation identity used for T4/T5: for random quadratic w (sector even), D = w(F), N = A w:
   B_D(N) := D div N - N(D)  satisfies  B_D(N) + A w^T Hess(F) w  in J_14  (and with the opposite sign it does NOT).
Also checks E1: N(F) = A D, and that the cleared defect equals F * B_D(N)."""
import sys, random, time
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
p = 2147483587
F, A, q = load(); Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
random.seed(3)
w = []
for i in range(NV):
    w.append({m: random.randrange(1, p) for m in monomials(2, par_i(i))})
D = {}
for i in range(NV): D = add(D, mul(w[i], Fi[i], p), p)
N = [mul(Ap, w[i], p) for i in range(NV)]
# E1
lhs = {}
for i in range(NV): lhs = add(lhs, mul(N[i], Fi[i], p), p)
assert not sub(lhs, mul(Ap, D, p), p), 'E1 fails'
divN = {}
for i in range(NV): divN = add(divN, pderiv(N[i], i, p), p)
ND = {}
for i in range(NV): ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
B = sub(mul(D, divN, p), ND, p)
E = core.defect(Fp, Ap, D, 1, N, p)
assert not sub(E, mul(Fp, B, p), p), 'defect != F*B_D(N)'
wHw = {}
for i in range(NV):
    for j in range(NV): wHw = add(wHw, mul(mul(w[i], w[j], p), H[i][j], p), p)
AwHw = mul(Ap, wHw, p)
def in_J(P_):
    deg = core.degree(P_); par = core.parity_of(P_)
    cols = []
    for i in range(NV):
        for m in monomials(deg - 5, (par - par_i(i)) % 2): cols.append(shift(Fi[i], m, 1, p))
    rows = {}
    for c in cols + [P_]:
        for e in c: rows.setdefault(e, len(rows))
    M = nmod_mat(len(rows), len(cols) + 1, p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = v
    for e, v in P_.items(): M[rows[e], len(cols)] = v
    M0 = nmod_mat(len(rows), len(cols), p)
    for j, c in enumerate(cols):
        for e, v in c.items(): M0[rows[e], j] = v
    return M.rank() == M0.rank()
t0 = time.time()
plus = in_J(add(B, AwHw, p)); minus = in_J(sub(B, AwHw, p)); raw = in_J(AwHw)
rec = log2(exp='relaxation-sign-check', prime=p, d=7, sector='even', B_plus_AwHw_in_J=plus, B_minus_AwHw_in_J=minus, AwHw_in_J=raw,
           result='identity B_D(Aw) = -A w^T Hess w (mod J) confirmed' if (plus and not minus) else 'SIGN PROBLEM', seconds=round(time.time()-t0, 1))
print(rec)
