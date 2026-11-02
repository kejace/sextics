"""Independent Hessian-relaxation check for Jacobian denominators D = w(F), w a LINEAR vector field (d = 6).
Necessary condition (own derivation): Hess_F(w,w) in J_6.  Sectors: 'even' (w_+ : X_j d_{X_i}, Y_j d_{Y_i}; 18 params, minus Euler),
'odd' (w_- : Y_j d_{X_i}, X_j d_{Y_i}; 18 params), 'mixed' (all 36 minus Euler).
We build the quadratic forms q_k(w) = coordinates of NF_{J_6}(Hess_F(w,w)) (normal form via rref of J_6 basis mod p),
then the Macaulay matrix of degree k (q_k times monomials of degree k-2) and test whether it has full column rank
(=> (S/I)_k = 0 => V(I) = {0} over F_p-bar, hence over C provided rank_p J_6 = dim_Q J_6)."""
import sys, json, time; sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/verify')
from v_lib import *
from flint import nmod_mat
from itertools import combinations_with_replacement
P_DEFAULT = 4294967291

def setup(p):
    F, A, q, _ = load_instance()
    Fp = reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    Hess = [[pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
    X = [parse_poly(v) for v in VARS]
    # linear fields: basis list of (i, j) meaning Z_j d_i
    even = [(i, j) for i in range(6) for j in range(6) if (i < 3) == (j < 3)]   # 18
    odd = [(i, j) for i in range(6) for j in range(6) if (i < 3) != (j < 3)]    # 18
    return Fp, Fi, Hess, even, odd

def nf_setup(Fi, p, par):
    """rref basis of J_6^par (generators Z_j F_i with right parity); returns function reducing a sextic to NF coordinates."""
    gens = []
    for i in range(6):
        pm = (par - (0 if i < 3 else 1)) % 2
        for j in range(6):
            if ypar(tuple(1 if k == j else 0 for k in range(6))) == pm:
                gens.append(shift(Fi[i], tuple(1 if k == j else 0 for k in range(6)), 1, p))
    mons = monomials(6, par); idx = {m: k for k, m in enumerate(mons)}
    G = nmod_mat(len(gens), len(mons), [0]*(len(gens)*len(mons)), p)
    for r, g in enumerate(gens):
        for e, v in g.items(): G[r, idx[e]] = v
    Gr, rk = G.rref()
    piv = []
    for r in range(rk):
        c = 0
        while c < len(mons) and int(Gr[r, c]) == 0: c += 1
        piv.append(c)
    nonpiv = [c for c in range(len(mons)) if c not in set(piv)]
    def nf(poly):
        v = [0]*len(mons)
        for e, c in poly.items(): v[idx[e]] = c % p
        for r, c in enumerate(piv):
            if v[c]:
                f = v[c]
                for cc in range(c, len(mons)):
                    g = int(Gr[r, cc])
                    if g: v[cc] = (v[cc] - f*g) % p
        return [v[c] for c in nonpiv]
    return nf, rk, len(nonpiv)

def hess_bilinear(Hess, f1, f2, p):
    """Hess_F(w1, w2) for basis fields w1 = Z_j1 d_i1, w2 = Z_j2 d_i2 : Z_j1 Z_j2 F_{i1 i2}."""
    (i1, j1), (i2, j2) = f1, f2
    m = [0]*6; m[j1] += 1; m[j2] += 1
    return shift(Hess[i1][i2], tuple(m), 1, p)

def quadrics(sector, p):
    Fp, Fi, Hess, even, odd = setup(p)
    if sector == 'even': basis = even
    elif sector == 'odd': basis = odd
    else: basis = even + odd
    # remove Euler direction from the even part by dropping the coordinate X0 d_X0 (any w can be shifted by h*Euler)
    if sector in ('even', 'mixed'):
        basis = [b for b in basis if b != (0, 0)]
    n = len(basis)
    nf0, rk0, n0 = nf_setup(Fi, p, 0); nf1, rk1, n1 = nf_setup(Fi, p, 1)
    exactJ6 = dimJ_exact(6, 0)
    # quadratic forms: for each pair (a<=b) of basis fields, the sextic Hess(w_a, w_b) (times 2 if a != b -- symmetric form; we
    # store coefficient of monomial w_a w_b, which equals (a==b ? Hess(a,a) : 2 Hess(a,b)))
    pairs = list(combinations_with_replacement(range(n), 2))
    # Q[k][pair] = coordinate k of NF(...) ; coordinates: even NF coords (n0) then odd NF coords (n1)
    cols = {}
    for (a, b) in pairs:
        h = hess_bilinear(Hess, basis[a], basis[b], p)
        if a != b: h = scale(h, 2, p)
        hp = parity_of(h) if h else None
        if hp is None and h: raise RuntimeError('mixed parity sextic?')
        if not h: continue
        vec = nf0(h) if hp == 0 else nf1(h)
        off = 0 if hp == 0 else n0
        for k, v in enumerate(vec):
            if v: cols.setdefault(off + k, {})[(a, b)] = v
    # the quadrics are the rows k (each is a quadratic form in n vars)
    qlist = [cols[k] for k in sorted(cols)]
    return basis, n, qlist, {'rank_p_J6_even': rk0, 'rank_p_J6_odd': rk1, 'exact_J6': exactJ6, 'nquadrics_raw': len(qlist)}

def macaulay_rank(n, qlist, k, p):
    """rank of degree-k Macaulay matrix: rows = q * (monomial of degree k-2), columns = monomials of degree k in n vars."""
    def monos(deg):
        return list(combinations_with_replacement(range(n), deg))
    colmons = monos(k); cidx = {m: i for i, m in enumerate(colmons)}
    mult = monos(k-2)
    rows = []
    for qd in qlist:
        for mm in mult:
            row = {}
            for (a, b), v in qd.items():
                m = tuple(sorted((a, b) + mm))
                row[cidx[m]] = (row.get(cidx[m], 0) + v) % p
            rows.append(row)
    R, C = len(rows), len(colmons)
    flat = [0]*(R*C)
    for r, row in enumerate(rows):
        for c, v in row.items(): flat[r*C + c] = v
    M = nmod_mat(R, C, flat, p)
    return M.rank(), R, C

if __name__ == '__main__':
    p = P_DEFAULT
    args = sys.argv[1:]
    if args and args[0].startswith('p='): p = int(args[0][2:]); args = args[1:]
    sector = args[0]; ks = [int(a) for a in args[1:]] or [2, 3]
    t0 = time.time()
    basis, n, qlist, meta = quadrics(sector, p)
    # linearization rank (degree 2): rank of quadrics as vectors in Sym^2
    out = {'sector': sector, 'p': p, 'nvars': n, **meta, 't_setup': round(time.time()-t0, 1)}
    for k in ks:
        t1 = time.time()
        rk, R, C = macaulay_rank(n, qlist, k, p)
        out[f'mac{k}'] = {'rank': rk, 'rows': R, 'cols': C, 'full_column_rank': rk == C, 'dim_quotient_deg_k': C - rk, 't': round(time.time()-t1, 1)}
        print(json.dumps(out), flush=True)
    log({'check': f'hessian-relaxation d=6 sector={sector}', **out, 'method': 'own code hessmac.py (NF mod J_6 + Macaulay rank)'})
    with open('/tmp/hodge/results/eight-plane-configuration/verify/hessmac.jsonl', 'a') as f: f.write(json.dumps(out)+'\n')
