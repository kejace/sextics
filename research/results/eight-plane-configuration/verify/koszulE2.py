"""Independent Koszul-reduced E2 system for Jacobian-type denominators D = w(F) (q = 1), exact over Q.
Derivation (own): N(F) = A D + mu F  =>  N = A w + mu E/6 + K_c, K_c = sum_{i<j} c_ij (F_j d_i - F_i d_j), c_ij in S_{d-4};
E = 0  <=>  D div N - N(D) - D mu = F Q, and mu cancels, so a certificate exists iff exist c_ij (right parity) and Q in S_{2d-6}
with   [D div(Aw) - (Aw)(D)] + sum_{i<j} [ D ((d_i c_ij) F_j - (d_j c_ij) F_i) - c_ij (F_j D_i - F_i D_j) ]  =  F Q.
Here w is a polynomial vector field with D = w(F) = sum_i w_i F_i (deg w_i = d-5).
Unknowns: c_ij coefficients + Q coefficients; equation in degree 2d (parity of A, i.e. odd).  Exact rank over Q via flint fmpq_mat,
or mod p via nmod_mat."""
import sys, json, time; sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/verify')
from v_lib import *
from flint import fmpq_mat, fmpq, nmod_mat

def reduced_system(F, A, w, p=None):
    """w: list of 6 polys (the vector field). Returns (cols, rhs, info)."""
    if p:
        F = reduce_mod(F, p); A = reduce_mod(A, p); w = [reduce_mod(x, p) if x else {} for x in w]
    Fi = [pderiv(F, i, p) for i in range(6)]
    D = {}
    for i in range(6):
        if w[i]: D = add(D, mul(w[i], Fi[i], p), p)
    d = deg(D); pD = parity_of(D); pA = parity_of(A)
    assert pD is not None, 'mixed-parity D: use the general builder'
    Di = [pderiv(D, i, p) for i in range(6)]
    # particular part: D div(Aw) - (Aw)(D)
    Aw = [mul(A, x, p) if x else {} for x in w]
    divAw = {}
    for i in range(6):
        if Aw[i]: divAw = add(divAw, pderiv(Aw[i], i, p), p)
    AwD = {}
    for i in range(6):
        if Aw[i]: AwD = add(AwD, mul(Aw[i], Di[i], p), p)
    part = add(mul(D, divAw, p), scale(AwD, -1, p), p)
    # Koszul unknowns: c_ij of degree d-4 with parity so that component N_i = c F_j has the right parity
    parsN = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(6)]
    cols = []; info = []
    degc = d - 4
    for i in range(6):
        for j in range(i+1, 6):
            # N_i = c F_j : parity(c) + parity(F_j) = parsN[i]
            pc = (parsN[i] - (0 if j < 3 else 1)) % 2
            # consistency with N_j = -c F_i
            assert pc == (parsN[j] - (0 if i < 3 else 1)) % 2
            for m in monomials(degc, pc):
                c = {m: 1}
                dic = pderiv(c, i, p); djc = pderiv(c, j, p)
                t = add(mul(dic, Fi[j], p) if dic else {}, scale(mul(djc, Fi[i], p), -1, p) if djc else {}, p)
                t = mul(D, t, p) if t else {}
                t2 = add(mul(Fi[j], Di[i], p), scale(mul(Fi[i], Di[j], p), -1, p), p)
                t = add(t, scale(shift(t2, m), -1, p), p)
                cols.append(t); info.append(('c', i, j, m))
    for m in monomials(2*d - 6, pA):
        cols.append(shift(F, m, -1, p)); info.append(('Q', m))
    rhs = scale(part, -1, p)   # cols . x = -part
    return cols, rhs, info, d

def rank_test(cols, rhs, p=None):
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    if p:
        flat = [0]*(R*(C+1))
        for j, c in enumerate(cols):
            for e, v in c.items(): flat[rows[e]*(C+1)+j] = v % p
        for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = v % p
        Ma = nmod_mat(R, C+1, flat, p); ra = Ma.rank()
        M = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p); r = M.rank()
    else:
        Ma = fmpq_mat(R, C+1)
        for j, c in enumerate(cols):
            for e, v in c.items(): Ma[rows[e], j] = fmpq(v.numerator, v.denominator)
        for e, v in rhs.items(): Ma[rows[e], C] = fmpq(v.numerator, v.denominator)
        ra = Ma.rank()
        M = fmpq_mat(R, C)
        for j, c in enumerate(cols):
            for e, v in c.items(): M[rows[e], j] = fmpq(v.numerator, v.denominator)
        r = M.rank()
    return {'rows': R, 'cols': C, 'rank': r, 'rank_aug': ra, 'consistent': r == ra}

if __name__ == '__main__':
    F, A, q, _ = load_instance()
    X = [parse_poly(v) for v in VARS]
    one = {(0,)*6: Fraction(1)}
    def unit(i, poly=one): w = [{}]*6; w = list(w); w[i] = poly; return w
    fields = {
        'dF/dX0': unit(0), 'dF/dX1': unit(1), 'dF/dX2': unit(2), 'dF/dY0': unit(3), 'dF/dY1': unit(4), 'dF/dY2': unit(5),
        'E0(F)': [X[0], {}, {}, X[3], {}, {}], 'E1(F)': [{}, X[1], {}, {}, X[4], {}], 'E2(F)': [{}, {}, X[2], {}, {}, X[5]],
        'epsX(F)': [X[0], X[1], X[2], {}, {}, {}], 'X0*dF/dX0': [X[0], {}, {}, {}, {}, {}],
        'R0(F)': [X[3], {}, {}, scale(X[0], -1), {}, {}], 'Y0*dF/dY0': [{}, {}, {}, X[3], {}, {}],
        'X0*dF/dY0': [{}, {}, {}, X[0], {}, {}], 'Y0*dF/dX0': [X[3], {}, {}, {}, {}, {}],
    }
    p = None
    args = sys.argv[1:]
    if args and args[0].startswith('p='): p = int(args[0][2:]); args = args[1:]
    for name in args:
        t0 = time.time()
        cols, rhs, info, d = reduced_system(F, A, fields[name], p)
        out = rank_test(cols, rhs, p); out.update({'D': name, 'd': d, 'p': p or 'exact-Q', 'seconds': round(time.time()-t0, 1), 'method': 'own Koszul-reduced E2 (koszulE2.py)'})
        print(json.dumps(out), flush=True)
        log({'check': f'koszulE2 {name}', **out})
        with open('/tmp/hodge/results/eight-plane-configuration/verify/koszulE2.jsonl', 'a') as f: f.write(json.dumps(out)+'\n')
