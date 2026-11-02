"""Hessian relaxation for Jacobian denominators D = w(F), deg D = d, w a vector field with coefficients in S_{d-5}.

Facts (exact over Q, proved in the report):
  * for d <= 9 every E1-admissible D of degree d lies in J(F), so D = w(F) (w of the parity matching D);
  * gcd(D,F) = 1  <=>  w not in T_d := {h*Euler : h in S_{d-6}}  (for d <= 9);
  * E2 implies the necessary condition  A * Hess_F(w,w) in J_{2d};  for d <= 7 (mult. by A injective on R_{2d-6})
    this is  Hess_F(w,w) in J_{2d-6};  the condition only depends on w modulo T_d.
This script builds, for a parity sector and a complement W' of T_d, the quadrics  w' -> NF_J(Hess_F(w',w'))
(or their pairing with A*R_{...} for d = 8,9), and
  d = 5 : checks Hess_F(w,w) = 0 has only w = 0 (exact, sympy Groebner over Q),
  d = 6 : linearization: rank of Sym^2(W') -> R_6^+ ; rank = dim Sym^2(W') proves V' = {0} (modular rank is a lower bound),
  d >= 7: writes the quadric system (mod p) for msolve and reports the linearization rank + dim of the quadric span.
"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *

def field_basis(d, par):
    """Basis of W_d^par: vector fields w = sum_i w_i d_i, w_i in S_{d-5}, with w_{X_i} parity par, w_{Y_i} parity par+1.
    Returns list of (i, m) meaning monomial m at component i."""
    out = []
    for i in range(NV):
        pm = (par + (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, pm):
            out.append((i, m))
    return out

def euler_mult_basis(d, par):
    """T_d^par = { h * Euler : h in S_{d-6}, parity(h) = par } as coefficient vectors over field_basis."""
    if d < 6: return []
    B = field_basis(d, par); bidx = {b: k for k, b in enumerate(B)}
    T = []
    for h in monomials(d - 6, par):
        v = {}
        for i in range(NV):
            m = list(h); m[i] += 1; m = tuple(m)
            v[bidx[(i, m)]] = 1
        T.append(v)
    return T

def complement_basis(B, T, p):
    """choose coordinates: a complement W' of span(T) inside W: drop one coordinate per T-vector (pivot after rref)."""
    if not T: return list(range(len(B)))
    n = len(B)
    flat = [0] * (len(T) * n)
    for r, v in enumerate(T):
        for k, c in v.items(): flat[r * n + k] = c % p
    M = nmod_mat(len(T), n, flat, p)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < n and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    assert rk == len(T)
    return [k for k in range(n) if k not in set(piv)]

def hess_bilinear(Fij, B, a, b, p):
    """Hess_F(e_a, e_b) = sum_ij (e_a)_i (e_b)_j F_ij for basis monomial fields e_a=(i,m), e_b=(j,m')."""
    i, m = B[a]; j, m2 = B[b]
    mm = tuple(x + y for x, y in zip(m, m2))
    return shift(Fij[i][j], mm, 1, p)

def build(d, par, p, verbose=True):
    t0 = time.time()
    Fp, Fi = partials(p)
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    B = field_basis(d, par); T = euler_mult_basis(d, par)
    comp = complement_basis(B, T, p)
    n = len(comp)
    deg = 2 * d - 6
    # target space: R_deg^{even}; for d >= 8 we must use Ann(A): condition A*Hess in J_{2d}  <=>  Hess in Ann_R(A)_{deg}
    red = Reducer(Fi, deg, 0, p)
    std = red.standard; nstd = len(std)
    # quadrics: for each pair (a<=b) of complement basis vectors, NF of Hess(e_a,e_b) (symmetrized: coefficient of w_a w_b is 2*Hess(e_a,e_b) for a<b)
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    nP = len(pairs)
    C = len(red.mons)
    flat = [0] * (nP * C)
    for r, (a, b) in enumerate(pairs):
        h = hess_bilinear(Fij, B, comp[a], comp[b], p)
        if a != b: h = scale(h, 2, p)
        for e, v in h.items(): flat[r * C + red.idx[e]] = v
    V = nmod_mat(nP, C, flat, p); del flat
    def sub(M, rows, cols):
        return nmod_mat(len(rows), len(cols), [int(M[i, j]) for i in rows for j in cols], p)
    Vstd = sub(V, range(nP), std); Vpiv = sub(V, range(nP), red.pivots); del V
    RS = nmod_mat(red.rank, nstd, [c for row in red.rows_std for c in row], p)
    NF = Vstd - Vpiv * RS           # nP x nstd : row (a,b) = coefficients of NF(Hess(e_a,e_b)) in standard monomials
    t1 = time.time()
    info = {'d': d, 'parD': par, 'p': p, 'dimW': len(B), 'dimT': len(T), 'dimWprime': n, 'dimSym2': nP,
            'deg_hess': deg, 'rank_p_J_deg': red.rank, 'exact_J_deg': J_dim(deg, 0), 'dim_R_deg': nstd}
    if d >= 8:
        # multiply by A and reduce in degree 2d (odd): only the image matters.  Build the map R_deg -> R_{2d}^{odd}: std monomial s -> NF(A*s)
        red2 = Reducer(Fi, 2 * d, 1, p)
        C2 = len(red2.mons)
        flat = [0] * (nstd * C2)
        Ap = reduce_mod(A, p)
        for r, j in enumerate(std):
            s = red.mons[j]
            for e, v in shift(Ap, s, 1, p).items(): flat[r * C2 + red2.idx[e]] = v
        W = nmod_mat(nstd, C2, flat, p); del flat
        Wstd = sub(W, range(nstd), red2.standard); Wpiv = sub(W, range(nstd), red2.pivots); del W
        RS2 = nmod_mat(red2.rank, len(red2.standard), [c for row in red2.rows_std for c in row], p)
        MA = Wstd - Wpiv * RS2     # nstd x dim R_{2d}^{odd}
        NF = NF * MA               # nP x dim R_{2d}^odd : quadric coefficients after multiplying by A
        info.update({'rank_p_J_2d': red2.rank, 'exact_J_2d': J_dim(2 * d, 1), 'dim_R_2d_odd': len(red2.standard), 'rank_multA': MA.rank()})
    # linearization rank and quadric span
    rank_lin = NF.rank()
    info['rank_linearization'] = rank_lin
    info['linearization_proves_empty'] = (rank_lin == nP)
    info['n_quadrics_independent'] = rank_lin   # = dim of span of the quadrics in Sym^2(W'^*) (same matrix, transposed view)
    info['t_build'] = round(t1 - t0, 1)
    return info, NF, pairs, comp, B, red

def write_msolve(d, par, p, NF, pairs, n, path):
    """write the independent quadrics (columns of NF = linear functionals on Sym^2) as polynomials in w1..wn for msolve."""
    # columns of NF: each column k gives quadric sum_{(a,b)} NF[(a,b),k] * w_a w_b ; take a basis of the column space via rref of NF^T
    NT = NF.transpose()
    NTr, rk = NT.rref()
    names = [f'w{a}' for a in range(n)]
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n')
        f.write(f'{p}\n')
        polys = []
        for i in range(rk):
            terms = []
            for r, (a, b) in enumerate(pairs):
                c = int(NTr[i, r])
                if c:
                    terms.append(f'{c}*w{a}*w{b}' if a != b else f'{c}*w{a}^2')
            polys.append('+'.join(terms))
        f.write(',\n'.join(polys) + '\n')
    return rk

if __name__ == '__main__':
    d = int(sys.argv[1]); par = int(sys.argv[2]); p = int(sys.argv[3]) if len(sys.argv) > 3 else P_NEW
    tag = sys.argv[4] if len(sys.argv) > 4 else 'local'
    info, NF, pairs, comp, B, red = build(d, par, p)
    n = len(comp)
    if d >= 6:
        path = f'{PK}/hess_d{d}_par{par}_p{p}.ms'
        rk = write_msolve(d, par, p, NF, pairs, n, path)
        info['msolve_input'] = path; info['n_quadrics_written'] = rk
    print(json.dumps(info), flush=True)
    with open(f'{PK}/hess.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
    wlog({'exp': 'b-hessian-relaxation', 'family': f'J_{d} (D = w(F)), parity {par}', 'd': d, 'prime': p,
          'result': ('V\'={0} by linearization (rigorous: modular rank = dim Sym^2)' if info['linearization_proves_empty'] else f"linearization rank {info['rank_linearization']} < {info['dimSym2']}; quadric system written for msolve"),
          'rank': info['rank_linearization'], 'seconds': info['t_build'], 'job': tag,
          'note': f"dim W={info['dimW']}, dim T={info['dimT']}, dim W'={info['dimWprime']}, Hess target R_{info['deg_hess']}^+ dim {info['dim_R_deg']}; rank_p J_{info['deg_hess']}={info['rank_p_J_deg']} exact {info['exact_J_deg']}"},
         fn=None if tag == 'local' else f'{PK}/wave2-excl-remote-{tag}.jsonl')
