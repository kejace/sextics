"""First-order deformation of the Fermat subset-Euler certificates toward dense-rational.
F_eps = F0 + eps F1, F1 = 9973 (F - F0).  D0 = D_I = sum_{i in I} Z_i^6, N0 = A E_I / 6 (canonical), Q0 = 0.
First-order equation (linear in D1, N1, Q1):
  A D0 D1 + F0 [D1 div N0 - N0(D1)] - D0 N1(F0) + F0 [D0 div N1 - N1(D0)] - F0^2 Q1 = D0 N0(F1)   (RHS)
Report consistency (mod p) with D1 in S_6^+ (or D1 = 0), N1_i in S_7 with parities, Q1 in S_6^-.
"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
F, A, q = load()
p = P
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
F1 = scale(sub(F, F0), 9973)
assert all(Fraction(c).denominator == 1 for c in F1.values())
F0p = core.reduce_mod(F0, p); F1p = core.reduce_mod(F1, p); Ap = core.reduce_mod(A, p)
Z = ['X0','X1','X2','Y0','Y1','Y2']
def euler_I(I, Pl):
    out = {}
    for i in I:
        e = [0]*6; e[i] = 1
        out = add(out, shift(pderiv(Pl, i, p), tuple(e), 1, p), p)
    return out
def run(I, allow_D1=True):
    D0 = core.reduce_mod({tuple(6 if k == i else 0 for k in range(6)): 1 for i in I}, p)
    # N0 = A E_I / 6 : components N0_i = A Z_i /6 for i in I
    inv6 = pow(6, -1, p)
    N0 = [{} for _ in range(6)]
    for i in I:
        e = [0]*6; e[i] = 1
        N0[i] = shift(Ap, tuple(e), inv6, p)
    # check order 0 identity
    E0 = core.defect(F0p, Ap, D0, 1, N0, p)
    assert not E0, 'order-0 identity fails'
    divN0 = {}
    for i in range(6): divN0 = add(divN0, pderiv(N0[i], i, p), p)
    F0i = [pderiv(F0p, i, p) for i in range(6)]
    D0i = [pderiv(D0, i, p) for i in range(6)]
    F02 = power(F0p, 2, p)
    cols = []; names = []
    # D1 columns (even sextic monomials)
    if allow_D1:
        for m in monomials(6, 0):
            mm = {m: 1}
            col = mul(mul(Ap, D0, p), mm, p)
            # F0 [D1 div N0 - N0(D1)]
            t = mul(mm, divN0, p)
            for i in I:
                t = sub(t, mul(N0[i], pderiv(mm, i, p), p), p)
            col = add(col, mul(F0p, t, p), p)
            cols.append(col); names.append(('D1', m))
    # N1 columns: -D0 N1_i F0_i + F0 (D0 d_i N1_i - N1_i D0_i)
    pars = fixed_d.numerator_parities(D0, Ap)
    for i in range(6):
        for m in monomials(7, pars[i]):
            mm = {m: 1}
            col = scale(mul(D0, mul(mm, F0i[i], p), p), p-1, p)
            col = add(col, mul(F0p, sub(mul(D0, pderiv(mm, i, p), p), mul(mm, D0i[i], p), p), p), p)
            cols.append(col); names.append(('N1', i, m))
    for m in monomials(6, 1):
        cols.append(scale(shift(F02, m, 1, p), p-1, p)); names.append(('Q1', m))
    rhs = {}
    for i in I:
        rhs = add(rhs, mul(N0[i], pderiv(F1p, i, p), p), p)
    rhs = mul(D0, rhs, p)
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*(C+1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e]*(C+1)+j] = v
    for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = v
    M = nmod_mat(R, C+1, flat, p)
    r_aug = M.rank()
    M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    r0 = M0.rank()
    return {'I': [Z[i] for i in I], 'allow_D1': allow_D1, 'rows': R, 'cols': C, 'rank': r0, 'rank_aug': r_aug, 'consistent': r0 == r_aug}
if __name__ == '__main__':
    out = []
    subsets = [I for I in itertools.combinations(range(6), 3) if (0 in I) != (3 in I)]
    for I in subsets:
        t0 = time.time()
        r = run(I, allow_D1=True)
        r['seconds'] = round(time.time()-t0, 1)
        print(r, flush=True)
        out.append(r)
    with open('/tmp/hodge/results/dense-rational/deform_first_order.json', 'w') as f:
        json.dump(out, f, indent=1)
