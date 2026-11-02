"""Independent fixed-D certificate system (q=1), built from scratch.
Unknowns: N_i (deg d+1, parity forced by D and A), Q (deg 2d-6, parity of A).  Equation: L(N,Q) = -A D^2, with
L = -q D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i D_i) - F^2 Q.
Consistency mod p decided by dense rank of [L | b] vs L (no compression).  Also returns the exact
sparse system over Q (rows scaled to integers) for certificate work."""
import sys, json, time; sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/verify')
from v_lib import *
from flint import nmod_mat
P_DEFAULT = 4294967291

def parities(D, A):
    pD, pA = parity_of(D), parity_of(A)
    assert pD is not None
    return [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(6)], pA

def build_system(F, A, q, D, p=None):
    """Return (cols: list of polys, rhs poly, colinfo) over Q (p=None) or mod p."""
    if p:
        F = reduce_mod(F, p); A = reduce_mod(A, p); D = reduce_mod(D, p)
    d = deg(D); parsN, pA = parities(D, A)
    Fi = [pderiv(F, i, p) for i in range(6)]; Di = [pderiv(D, i, p) for i in range(6)]
    FD = mul(F, D, p); F2 = mul(F, F, p)
    G = [add(scale(mul(D, Fi[i], p), -q, p), scale(mul(F, Di[i], p), -1, p), p) for i in range(6)]  # -(q D F_i + F D_i)
    cols = []; info = []
    for i in range(6):
        for m in monomials(d+1, parsN[i]):
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(col); info.append(('N', i, m))
    degQ = 2*d - 6
    if degQ >= 0:
        for m in monomials(degQ, pA):
            cols.append(shift(F2, m, -1, p)); info.append(('Q', -1, m))
    rhs = scale(mul(A, mul(D, D, p), p), -1, p)
    return cols, rhs, info

def consistency_mod_p(F, A, q, D, p=P_DEFAULT):
    t0 = time.time()
    cols, rhs, info = build_system(F, A, q, D, p)
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*(C+1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e]*(C+1)+j] = v
    for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = v
    Maug = nmod_mat(R, C+1, flat, p)
    t1 = time.time()
    r_aug = Maug.rank()
    M = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    r = M.rank()
    t2 = time.time()
    return {'d': deg(D), 'p': p, 'rows': R, 'cols': C, 'nN': sum(1 for i in info if i[0]=='N'), 'nQ': sum(1 for i in info if i[0]=='Q'),
            'rank': r, 'rank_aug': r_aug, 'consistent_mod_p': r_aug == r, 'nullity': C - r, 't_build': round(t1-t0,1), 't_rank': round(t2-t1,1),
            'method': 'own code fixedD.py dense rank, no compression'}

if __name__ == '__main__':
    F, A, q, _ = load_instance()
    p = P_DEFAULT
    args = sys.argv[1:]
    if args and args[0].startswith('p='): p = int(args[0][2:]); args = args[1:]
    Fi = [pderiv(F, i) for i in range(6)]
    named = {'dF/dX0': Fi[0], 'dF/dX1': Fi[1], 'dF/dX2': Fi[2], 'dF/dY0': Fi[3], 'dF/dY1': Fi[4], 'dF/dY2': Fi[5]}
    X = [parse_poly(v) for v in VARS]
    named['E0(F)'] = add(mul(X[0], Fi[0]), mul(X[3], Fi[3])); named['E1(F)'] = add(mul(X[1], Fi[1]), mul(X[4], Fi[4])); named['E2(F)'] = add(mul(X[2], Fi[2]), mul(X[5], Fi[5]))
    named['epsX(F)'] = add(add(mul(X[0], Fi[0]), mul(X[1], Fi[1])), mul(X[2], Fi[2]))
    named['X0*dF/dX0'] = mul(X[0], Fi[0])
    named['R0(F)'] = add(mul(X[3], Fi[0]), scale(mul(X[0], Fi[3]), -1))
    for name in args:
        D = named[name] if name in named else parse_poly(name)
        out = consistency_mod_p(F, A, q, D, p); out['D'] = name
        print(json.dumps(out), flush=True)
        log({'check': f'fixedD {name}', **out})
        with open('/tmp/hodge/results/eight-plane-configuration/verify/fixedD.jsonl', 'a') as f: f.write(json.dumps(out)+'\n')
