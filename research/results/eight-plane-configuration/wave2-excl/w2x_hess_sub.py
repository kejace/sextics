"""Hessian relaxation restricted to a structured subfamily of vector fields (d = 7, 8, 9).

Given a list of fields (6-tuples of polys mod p, coefficients in S_{d-5}, correct parity) spanning a subspace U of W_d^par,
remove U cap T_d (T_d = {h*Euler}) and write the quadric system  w' -> [A * Hess_F(w',w')] in R_{2d}^{odd}  (d >= 8)
or  [Hess_F(w',w')] in R_{2d-6}^{even}  (d <= 7) as msolve input.  V' = {0} (msolve quotient dimension finite)
proves rigorously that no coprime D = w(F) with w in U admits a certificate.

Subfamilies (name -> fields):
  arr      : arrangement-preserving fields  sum_i h_i R_i + sum_i g_i E_i,  R_i = Y_i d_{X_i} - X_i d_{Y_i},  E_i = X_i d_{X_i} + Y_i d_{Y_i}
  ham:<G>  : Hamiltonian fields  w = C grad(G)  for all antisymmetric constant C of the right parity, G from a named list
  grad     : gradient fields  w = grad(G), G in S_{d-5+1} of the right parity (all G)  -> D = grad F . grad G
  blockX   : fields with components only along d_{X_i} (w_{Y_i} = 0)
  blockY   : fields with components only along d_{Y_i}
  coefX    : fields whose coefficients lie in C[X0,X1,X2]
  coefB0   : fields whose coefficients lie in C[X0,Y0]
  rand:<k> : k random fields (sanity: should give V' = {0} quickly)
Usage: python3 w2x_hess_sub.py d par family [p] [tag]
"""
import sys, time, json, random, itertools
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from w2x_hess import field_basis, euler_mult_basis

def named_forms(p):
    K = [reduce_mod(core.parse(s), p) for s in json.load(open(f'{WD}/K_decomp.json'))['K']]
    H = [reduce_mod(core.parse(s), p) for s in json.load(open(f'{WD}/K_decomp.json'))['H']]
    Q = [reduce_mod(core.parse(f'X{i}**2+Y{i}**2'), p) for i in range(3)]
    diag = [reduce_mod(core.parse(f'X{i}**4-X{i}**2*Y{i}**2+Y{i}**4'), p) for i in range(3)]
    forms = {}
    for i in range(3):
        forms[f'Q{i}'] = Q[i]; forms[f'K{i}'] = K[i]; forms[f'H{i}'] = H[i]; forms[f'diag{i}'] = diag[i]; forms[f'Q{i}^2'] = mul(Q[i], Q[i], p)
    for i in range(3):
        for j in range(i + 1, 3): forms[f'Q{i}Q{j}'] = mul(Q[i], Q[j], p)
    forms['Q0Q1Q2'] = mul(mul(Q[0], Q[1], p), Q[2], p)
    forms['sumQ'] = add(add(Q[0], Q[1], p), Q[2], p)
    forms['F'] = reduce_mod(F, p)
    Fp, Fi = partials(p)
    for i in range(6): forms[f'F_{core.VARS[i]}'] = Fi[i]
    return forms

def to_vec(field, B, bidx, p):
    """field: list of 6 polys -> coefficient vector over field_basis B (dict k -> coeff). Raises if outside W_d^par."""
    v = {}
    for i in range(NV):
        for m, c in field[i].items():
            k = bidx.get((i, m))
            if k is None: raise ValueError(f'field component {i} monomial {m} not in W_d^par')
            v[k] = (v.get(k, 0) + c) % p
    return v

def subfamily(name, d, par, p):
    deg = d - 5
    fields = []
    X = [tuple(1 if k == i else 0 for k in range(NV)) for i in range(NV)]
    if name == 'arr':
        # h_i R_i: components (X_i: Y_i h, Y_i: -X_i h); parity of h: w_{X_i} = Y_i h must have parity par -> h parity par+1
        for i in range(3):
            for h in monomials(deg - 1, (par + 1) % 2):
                w = [{} for _ in range(NV)]
                w[i] = {tuple(a + b for a, b in zip(h, X[i + 3])): 1}
                w[i + 3] = {tuple(a + b for a, b in zip(h, X[i])): p - 1}
                fields.append(w)
            for g in monomials(deg - 1, par):
                w = [{} for _ in range(NV)]
                w[i] = {tuple(a + b for a, b in zip(g, X[i])): 1}
                w[i + 3] = {tuple(a + b for a, b in zip(g, X[i + 3])): 1}
                fields.append(w)
    elif name.startswith('ham:'):
        G = named_forms(p)[name[4:]]
        assert degree(G) == deg + 1, (degree(G), deg + 1)
        Gi = [pderiv(G, i, p) for i in range(NV)]
        pG = parity_of(G)
        for i in range(NV):
            for j in range(i + 1, NV):
                # w = G_j d_i - G_i d_j ; parity(w_{X_i}) must be par: parity(G_j) = pG + (1 if j>=3 else 0)
                ok_i = ((pG + (1 if j >= 3 else 0)) % 2) == ((par + (1 if i >= 3 else 0)) % 2)
                ok_j = ((pG + (1 if i >= 3 else 0)) % 2) == ((par + (1 if j >= 3 else 0)) % 2)
                if ok_i and ok_j:
                    w = [{} for _ in range(NV)]; w[i] = Gi[j]; w[j] = scale(Gi[i], p - 1, p)
                    if any(w): fields.append(w)
    elif name == 'grad':
        for G in monomials(deg + 1, par):
            w = [pderiv({G: 1}, i, p) for i in range(NV)]
            if any(w): fields.append(w)
    elif name in ('blockX', 'blockY'):
        rng = range(3) if name == 'blockX' else range(3, 6)
        for i in rng:
            for m in monomials(deg, (par + (0 if i < 3 else 1)) % 2):
                w = [{} for _ in range(NV)]; w[i] = {m: 1}; fields.append(w)
    elif name in ('coefX', 'coefB0'):
        allowed = (lambda m: m[3] == m[4] == m[5] == 0) if name == 'coefX' else (lambda m: m[1] == m[2] == m[4] == m[5] == 0)
        for i in range(NV):
            for m in monomials(deg, (par + (0 if i < 3 else 1)) % 2):
                if allowed(m):
                    w = [{} for _ in range(NV)]; w[i] = {m: 1}; fields.append(w)
    elif name.startswith('comp:'):
        i = int(name[5:])
        for m in monomials(deg, (par + (0 if i < 3 else 1)) % 2):
            w = [{} for _ in range(NV)]; w[i] = {m: 1}; fields.append(w)
    elif name.startswith('hamrand:'):
        seed = int(name[8:]); random.seed(seed)
        pG = seed % 2
        G = {m: random.randint(1, p - 1) for m in monomials(deg + 1, pG)}
        Gi = [pderiv(G, i, p) for i in range(NV)]
        for i in range(NV):
            for j in range(i + 1, NV):
                ok_i = ((pG + (1 if j >= 3 else 0)) % 2) == ((par + (1 if i >= 3 else 0)) % 2)
                ok_j = ((pG + (1 if i >= 3 else 0)) % 2) == ((par + (1 if j >= 3 else 0)) % 2)
                if ok_i and ok_j:
                    w = [{} for _ in range(NV)]; w[i] = Gi[j]; w[j] = scale(Gi[i], p - 1, p)
                    if any(w): fields.append(w)
    elif name.startswith('rand:'):
        k = int(name[5:]); random.seed(k)
        B = field_basis(d, par)
        for r in range(k):
            w = [{} for _ in range(NV)]
            for (i, m) in B:
                c = random.randint(0, p - 1)
                w[i][m] = c
            fields.append(w)
    else:
        raise ValueError(name)
    return fields

def build_sub(d, par, p, fields, verbose=True):
    t0 = time.time()
    Fp, Fi = partials(p)
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    B = field_basis(d, par); bidx = {b: k for k, b in enumerate(B)}; n = len(B)
    T = euler_mult_basis(d, par)
    U = [to_vec(w, B, bidx, p) for w in fields]
    # complement of (span U cap span T) inside span U: rref [T; U] and keep U rows that increase the rank
    def mat(vs):
        flat = [0] * (len(vs) * n)
        for r, v in enumerate(vs):
            for k, c in v.items(): flat[r * n + k] = c
        return nmod_mat(len(vs), n, flat, p)
    rT = mat(T).rank() if T else 0
    rU = mat(U).rank()
    keep = []; cur = rT; basisrows = list(T)
    for v in U:
        r = mat(basisrows + [v]).rank()
        if r > cur: keep.append(v); basisrows.append(v); cur = r
    m = len(keep)
    info = {'d': d, 'parD': par, 'p': p, 'dimW': n, 'dimT': rT, 'dimU': rU, 'dimU_cap_T': rU - m, 'dimUprime': m}
    # fields of the complement basis as 6-tuples
    comp_fields = []
    for v in keep:
        w = [{} for _ in range(NV)]
        for k, c in v.items():
            i, mo = B[k]; w[i][mo] = c
        comp_fields.append(w)
    # quadric coefficients: Hess(e_a, e_b) = sum_ij (e_a)_i (e_b)_j F_ij, reduced
    deg = 2 * d - 6
    red = Reducer(Fi, deg, 0, p)
    pairs = [(a, b) for a in range(m) for b in range(a, m)]
    C = len(red.mons)
    flat = [0] * (len(pairs) * C)
    for r, (a, b) in enumerate(pairs):
        h = {}
        for i in range(NV):
            if not comp_fields[a][i]: continue
            for j in range(NV):
                if not comp_fields[b][j]: continue
                h = add(h, mul(mul(comp_fields[a][i], comp_fields[b][j], p), Fij[i][j], p), p)
        if a != b: h = scale(h, 2, p)
        for e, v in h.items(): flat[r * C + red.idx[e]] = v
    V = nmod_mat(len(pairs), C, flat, p); del flat
    def sub(M, rows, cols):
        return nmod_mat(len(rows), len(cols), [int(M[i, j]) for i in rows for j in cols], p)
    Vstd = sub(V, range(len(pairs)), red.standard); Vpiv = sub(V, range(len(pairs)), red.pivots); del V
    RS = nmod_mat(red.rank, len(red.standard), [c for row in red.rows_std for c in row], p)
    NF = Vstd - Vpiv * RS
    info.update({'deg_hess': deg, 'rank_p_J_deg': red.rank, 'exact_J_deg': J_dim(deg, 0), 'dim_R_deg': len(red.standard)})
    if d >= 8:
        red2 = Reducer(Fi, 2 * d, 1, p)
        nstd = len(red.standard); C2 = len(red2.mons)
        Ap = reduce_mod(A, p)
        flat = [0] * (nstd * C2)
        for r, j in enumerate(red.standard):
            for e, v in shift(Ap, red.mons[j], 1, p).items(): flat[r * C2 + red2.idx[e]] = v
        W = nmod_mat(nstd, C2, flat, p); del flat
        Wstd = sub(W, range(nstd), red2.standard); Wpiv = sub(W, range(nstd), red2.pivots); del W
        RS2 = nmod_mat(red2.rank, len(red2.standard), [c for row in red2.rows_std for c in row], p)
        MA = Wstd - Wpiv * RS2
        NF = NF * MA
        info.update({'rank_p_J_2d': red2.rank, 'exact_J_2d': J_dim(2 * d, 1), 'dim_R_2d_odd': len(red2.standard), 'rank_multA': MA.rank(), 'dim_AnnA_deg': nstd - MA.rank()})
    rk = NF.rank()
    info.update({'dimSym2': len(pairs), 'rank_linearization': rk, 'linearization_proves_empty': rk == len(pairs), 't_build': round(time.time() - t0, 1)})
    return info, NF, pairs, m

def write_ms(NF, pairs, m, p, path):
    NT = NF.transpose(); NTr, rk = NT.rref()
    with open(path, 'w') as f:
        f.write(','.join(f'w{a}' for a in range(m)) + '\n'); f.write(f'{p}\n')
        polys = []
        for i in range(rk):
            terms = []
            for r, (a, b) in enumerate(pairs):
                c = int(NTr[i, r])
                if c: terms.append(f'{c}*w{a}*w{b}' if a != b else f'{c}*w{a}^2')
            polys.append('+'.join(terms))
        f.write(',\n'.join(polys) + '\n')
    return rk

if __name__ == '__main__':
    d = int(sys.argv[1]); par = int(sys.argv[2]); fam = sys.argv[3]
    p = int(sys.argv[4]) if len(sys.argv) > 4 else P_NEW
    tag = sys.argv[5] if len(sys.argv) > 5 else 'local'
    fields = subfamily(fam, d, par, p)
    info, NF, pairs, m = build_sub(d, par, p, fields)
    info['family'] = fam
    safe = fam.replace(':', '-').replace('^', 'p')
    if m and not info['linearization_proves_empty']:
        path = f'{PK}/hsub_d{d}_par{par}_{safe}_p{p}.ms'
        info['msolve_input'] = path; info['n_quadrics_written'] = write_ms(NF, pairs, m, p, path)
    print(json.dumps(info), flush=True)
    with open(f'{PK}/hess_sub' + ('' if tag == 'local' else f'.{tag}') + '.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
    wlog({'exp': 'b-hess-subfamily', 'family': f'{fam} subfamily of J_{d}, parity {par}', 'd': d, 'prime': p,
          'result': ('empty: dimU\'=0 (subfamily inside T, i.e. F-multiples only)' if m == 0 else 'V\'={0} by linearization (rigorous)' if info['linearization_proves_empty'] else f"linearization rank {info['rank_linearization']} < {info['dimSym2']}: msolve input written"),
          'rank': info['rank_linearization'], 'seconds': info['t_build'], 'job': tag,
          'note': f"dimU={info['dimU']}, dimU cap T={info['dimU_cap_T']}, dimU'={m}"},
         fn=None if tag == 'local' else f'{PK}/wave2-excl-remote-{tag}.jsonl')
