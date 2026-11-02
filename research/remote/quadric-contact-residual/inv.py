"""Gamma-invariant E1 kernel and quadratic E2 closure for quadric-contact-residual.

Gamma = G_aff (order 48, abelian, characters chi) x| <s0 = (0<->3), r1 = (1<->4)>.
All computations mod p = P0 (31-bit prime), polynomials as dicts {exp6: int}.
"""
import sys, time, json, itertools
from common import *          # F, A, q, chi, chi_add, chi_neg, CHI_A, CHI_Z, monomials, mul, add, sub, scale, pderiv, shift, P0, log, SectorSystem, sector_e1
from flint import nmod_mat

import os
p = int(os.environ.get("HODGE_P", P0))
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
CHI0 = (0, 0, 0, 0)
INV4 = pow(4, p - 2, p)

# ---------------------------------------------------------------- group <s0, r1>
def perm_poly(P, pi):
    out = {}
    for e, c in P.items():
        e2 = tuple(e[pi[j]] for j in range(NV))
        out[e2] = c
    return out

ID = (0, 1, 2, 3, 4, 5); S0 = (3, 1, 2, 0, 4, 5); R1 = (0, 4, 2, 3, 1, 5); S0R1 = (3, 4, 2, 0, 1, 5)
GRP = [ID, S0, R1, S0R1]

def sym_poly(P):
    out = {}
    for g in GRP:
        out = add(out, perm_poly(P, g), p)
    return scale(out, INV4, p)

def perm_field(N, g):
    # (gN)_i = g(N_{g(i)})
    return [perm_poly(N[g[i]], g) for i in range(NV)]

def sym_field(N):
    out = [{} for _ in range(NV)]
    for g in GRP:
        gN = perm_field(N, g)
        out = [add(out[i], gN[i], p) for i in range(NV)]
    return [scale(c, INV4, p) for c in out]

def orbit(e):
    return sorted({tuple(e[g[j]] for j in range(NV)) for g in GRP})

def inv_basis(d):
    """Basis of Gamma-invariant polynomials of degree d: orbit sums of chi=0 monomials (coeff 1 each)."""
    seen, basis = set(), []
    for e in monomials(d):
        if chi(e) != CHI0 or e in seen:
            continue
        orb = orbit(e)
        seen.update(orb)
        basis.append({m: 1 for m in orb})
    return basis

# ---------------------------------------------------------------- sector linear algebra
class Sector:
    """Monomials of degree m with character c; normal forms modulo J_m[c]."""
    def __init__(self, m, c, verbose=False):
        t0 = time.time()
        self.m, self.c = m, c
        self.mons = [e for e in monomials(m) if chi(e) == c]
        self.idx = {e: i for i, e in enumerate(self.mons)}
        self.n = len(self.mons)
        rows = []
        self.gens = []    # (i, m') generators
        if m >= 5:
            for mm in monomials(m - 5):
                for i in range(NV):
                    if chi_add(chi(mm), chi_neg(CHI_Z[i])) == c:
                        rows.append(shift(Fi[i], mm, 1, p))
                        self.gens.append((i, mm))
        self.nJ = len(rows)
        flat = [0] * (self.nJ * self.n)
        for r, P in enumerate(rows):
            for e, v in P.items():
                flat[r * self.n + self.idx[e]] = v
        M = nmod_mat(self.nJ, self.n, flat, p) if self.nJ else None
        if M is not None:
            Mr, rk = M.rref()
            self.rank = rk
            piv = []
            for i in range(rk):
                j = 0
                while int(Mr[i, j]) == 0:
                    j += 1
                piv.append(j)
            self.piv = piv
            # reduced basis rows (rk x n) as nmod_mat
            self.Rb = nmod_mat(rk, self.n, [int(Mr[i, j]) for i in range(rk) for j in range(self.n)], p) if rk else None
            self.Mgen = M
        else:
            self.rank, self.piv, self.Rb, self.Mgen = 0, [], None, None
        self.nonpiv = [j for j in range(self.n) if j not in set(self.piv)]
        self.dimR = self.n - self.rank
        if verbose:
            print(f'[Sector m={m} c={c}] n={self.n} Jrows={self.nJ} rank={self.rank} dimR={self.dimR} ({time.time()-t0:.1f}s)', file=sys.stderr)

    def vec(self, P):
        v = [0] * self.n
        for e, c in P.items():
            v[self.idx[e]] = c % p
        return v

    def nf_matrix(self, vecs):
        """Normal forms (as nmod_mat rows restricted to non-pivot columns) of a list of vectors."""
        k = len(vecs)
        V = nmod_mat(k, self.n, [x for v in vecs for x in v], p)
        if self.rank:
            Pm = nmod_mat(k, self.rank, [v[j] for v in vecs for j in self.piv], p)
            V = V - Pm * self.Rb
        # restrict to nonpivot columns
        return nmod_mat(k, len(self.nonpiv), [int(V[i, j]) for i in range(k) for j in self.nonpiv], p)

    def in_J(self, P):
        N = self.nf_matrix([self.vec(P)])
        return all(int(N[0, j]) == 0 for j in range(N.ncols()))

    def lift(self, P):
        """Express P (in J) as sum over gens: returns list of coefficient per generator, or None."""
        v = self.vec(P)
        # solve Mgen^T x = v ; augmented rref
        n, nJ = self.n, self.nJ
        flat = []
        for j in range(n):
            flat.extend(int(self.Mgen[r, j]) for r in range(nJ))
            flat.append(v[j])
        Aug = nmod_mat(n, nJ + 1, flat, p)
        R, rk = Aug.rref()
        x = [0] * nJ
        for i in range(rk):
            j = 0
            while int(R[i, j]) == 0:
                j += 1
            if j == nJ:
                return None
            x[j] = int(R[i, nJ])
        return x

# ---------------------------------------------------------------- E1 kernel in the invariant sector
def e1_kernel(d, verbose=True, c=CHI0, sym=True):
    t0 = time.time()
    if sym:
        assert c == CHI0
        basis = inv_basis(d)                 # Gamma-invariant S_d
    else:
        basis = [{e: 1} for e in monomials(d) if chi(e) == c]
    nb = len(basis)
    S6 = Sector(d + 6, chi_add(c, CHI_A), verbose=verbose)
    vecs = [S6.vec(mul(Ap, D, p)) for D in basis]
    NF = S6.nf_matrix(vecs)                  # nb x dimR
    # kernel: c with c . NF = 0  -> nullspace of NF^T
    X, nul = NF.transpose().nullspace()
    ker = [[int(X[i, j]) for i in range(nb)] for j in range(nul)]
    # Jacobian part: symmetrized m F_i, with C vectors
    jac = []       # list of (Dpoly, Cfield); the first n_gauge entries are the gauge directions s*F
    n_gauge = 0
    INV6 = pow(6, p - 2, p)
    if d >= 6:
        gauge_s = inv_basis(d - 6) if sym else [{e: 1} for e in monomials(d - 6) if chi(e) == c]
        for s in gauge_s:
            C = [scale(mul(s, {tuple(1 if j == i else 0 for j in range(NV)): 1}, p), INV6, p) for i in range(NV)]
            Dp = mul(s, Fp, p)
            jac.append((Dp, C)); n_gauge += 1
    if d >= 5:
        for mm in monomials(d - 5):
            for i in range(NV):
                if chi_add(chi(mm), chi_neg(CHI_Z[i])) == c:
                    C = [{} for _ in range(NV)]; C[i] = {mm: 1}
                    Cs = sym_field(C) if sym else C
                    Dp = {}
                    for j in range(NV):
                        Dp = add(Dp, mul(Cs[j], Fi[j], p), p)
                    if Dp:
                        jac.append((Dp, Cs))
    # express jac elements in the invariant basis and find independent ones
    bidx = {}
    for k, B in enumerate(basis):
        for e in B:
            bidx[e] = k
    def coords(P):
        v = [0] * nb
        for e, c in P.items():
            v[bidx[e]] = c % p
        return v
    jvecs = [coords(Dp) for Dp, _ in jac]
    if jvecs:
        JM = nmod_mat(len(jvecs), nb, [x for v in jvecs for x in v], p)
        JR, jrank = JM.rref()
        # choose independent subset greedily
        indep = []
        cur = None
        for k, v in enumerate(jvecs):
            test = nmod_mat(len(indep) + 1, nb, [x for u in [jvecs[t] for t in indep] + [v] for x in u], p)
            if test.rank() == len(indep) + 1:
                indep.append(k)
            if len(indep) == jrank:
                break
    else:
        jrank, indep = 0, []
    # extras: kernel vectors modulo the Jacobian span
    extras = []
    if nul > jrank:
        span = [jvecs[t] for t in indep]
        for v in ker:
            test = nmod_mat(len(span) + 1, nb, [x for u in span + [v] for x in u], p)
            if test.rank() == len(span) + 1:
                span.append(v); extras.append(v)
    info = {'d': d, 'sector': str(c), 'sym': sym, 'n_gauge': n_gauge, 'dim_inv_Sd': nb, 'dim_kernel': nul, 'dim_J_inv': jrank, 'n_extra': len(extras),
            'sector_rows': S6.nJ, 'sector_cols': S6.n, 'dimR_d6': S6.dimR, 'seconds': round(time.time() - t0, 1)}
    if verbose:
        print('[e1_kernel]', info, file=sys.stderr)
    def topoly(v):
        P = {}
        for k, c in enumerate(v):
            if c % p:
                P = add(P, scale(basis[k], c, p), p)
        return P
    extra_polys = [topoly(v) for v in extras]
    jac_sel = [jac[t] for t in indep]
    return info, basis, jac_sel, extra_polys, S6

def lift_extra(Dp, d, S6, sym=True):
    """Solve N(F) = A*Dp for N (deg d+1), optionally symmetrize over <s0,r1>; returns N list."""
    x = S6.lift(mul(Ap, Dp, p))
    assert x is not None, 'extra not in J?!'
    N = [{} for _ in range(NV)]
    for (i, mm), c in zip(S6.gens, x):
        if c:
            N[i] = add(N[i], {mm: c}, p)
    if sym:
        N = sym_field(N)
    # check
    chk = {}
    for i in range(NV):
        chk = add(chk, mul(N[i], Fi[i], p), p)
    assert sub(chk, mul(Ap, Dp, p), p) == {}, 'lift check failed'
    return N

# ---------------------------------------------------------------- quadratic E2 closure
def B_matrix(N, d, Sd0, S2d):
    """Matrix of D -> B_D(N) = D div N - N(D) from S_d[0] (monomials Sd0.mons) to S_{2d}[chi_A]."""
    divN = {}
    for i in range(NV):
        divN = add(divN, pderiv(N[i], i, p), p)
    rows = []
    for e in Sd0.mons:
        P = shift(divN, e, 1, p)
        for i in range(NV):
            if e[i]:
                e2 = list(e); e2[i] -= 1
                P = sub(P, shift(N[i], tuple(e2), e[i], p), p)
        rows.append(S2d.vec(P))
    return nmod_mat(len(rows), S2d.n, [x for r in rows for x in r], p)

def quadratic_closure(d, verbose=True, out_prefix=None, c=CHI0, sym=True):
    t0 = time.time()
    info, basis, jac_sel, extras, S6 = e1_kernel(d, verbose=verbose, c=c, sym=sym)
    Ds = [Dp for Dp, _ in jac_sel] + extras
    Ns = [[mul(Ap, cc, p) for cc in Cs] for _, Cs in jac_sel]
    for Dp in extras:
        Ns.append(lift_extra(Dp, d, S6, sym=sym))
    n = len(Ds)
    if n == 0:
        return {'d': d, 'sector': str(c), 'n_params': 0}, [], [], []
    Sd0 = Sector(d, c, verbose=False)          # only need monomial list
    S2d = Sector(2 * d, chi_add(CHI_A, chi_add(c, c)), verbose=verbose)
    Dvecs = nmod_mat(n, Sd0.n, [x for Dp in Ds for x in Sd0.vec(Dp)], p)
    # B_{kl} = D_k . Bmat(N_l) ; collect all n*n then normal form
    allrows = []
    for l in range(n):
        Bm = B_matrix(Ns[l], d, Sd0, S2d)
        Prod = Dvecs * Bm                      # n x S2d.n ; row k = B_{D_k}(N_l)
        for k in range(n):
            allrows.append([int(Prod[k, j]) for j in range(S2d.n)])
        if verbose and (l % 5 == 0 or l == n - 1):
            print(f'  B matrices {l+1}/{n} ({time.time()-t0:.0f}s)', file=sys.stderr)
    NF = S2d.nf_matrix(allrows)                # (n*n) x dimR
    m = NF.ncols()
    # symmetric quadratic forms: Q_j(c) = sum_{k<=l} q_{kl,j} c_k c_l with q = NF[k,l] + NF[l,k] (k<l), NF[k,k]
    quads = []
    for j in range(m):
        coeffs = {}
        for k in range(n):
            for l in range(k, n):
                v = int(NF[k * n + l, j]) + (int(NF[l * n + k, j]) if l != k else 0)
                v %= p
                if v:
                    coeffs[(k, l)] = v
        if coeffs:
            quads.append(coeffs)
    # rank of the quadric system as linear forms in the monomials c_k c_l
    pairs = [(k, l) for k in range(n) for l in range(k, n)]
    pidx = {pr: i for i, pr in enumerate(pairs)}
    QM = nmod_mat(len(quads), len(pairs), [0] * (len(quads) * len(pairs)), p) if quads else None
    if quads:
        flat = [0] * (len(quads) * len(pairs))
        for r, cf in enumerate(quads):
            for pr, v in cf.items():
                flat[r * len(pairs) + pidx[pr]] = v
        QM = nmod_mat(len(quads), len(pairs), flat, p)
        qrank = QM.rank()
    else:
        qrank = 0
    # square-column test: c_k forced to 0 if the pure square column (k,k) is independent of all others
    forced = []
    if quads:
        QR, _ = QM.rref()
        # column (k,k) independent of all other columns <=> rank drops when removing it ... use: rank(QM) vs rank(QM without column)
        full = qrank
        for k in range(n):
            j = pidx[(k, k)]
            cols = [i for i in range(len(pairs)) if i != j]
            sub_ = nmod_mat(len(quads), len(cols), [int(QM[r, i]) for r in range(len(quads)) for i in cols], p)
            if sub_.rank() < full:
                forced.append(k)
    res = {'d': d, 'sector': str(c), 'sym': sym, 'n_params': n, 'n_gauge': info['n_gauge'], 'n_jac': len(jac_sel), 'n_extra': len(extras), 'dimR_2d_sector': m,
           'n_quadrics': len(quads), 'quadric_rank': qrank, 'n_pair_monomials': len(pairs),
           'forced_zero_by_square_columns': forced, 'seconds': round(time.time() - t0, 1)}
    res.update({'e1_' + k: v for k, v in info.items() if k != 'd'})
    if verbose:
        print('[quadratic_closure]', res, file=sys.stderr)
    if out_prefix:
        # Singular input for the quadric ideal
        tag = f'_d{d}' + ('' if sym else '_c' + ''.join(str(x) for x in c))
        ng = info['n_gauge']
        with open(out_prefix + f'_quads{tag}.sing', 'w') as f:
            f.write(f'// gauge fixed: c(1..{ng}) = 0 (directions s*F)\nring r = {p}, (c({ng+1}..{n})), dp;\nideal I;\n')
            for cf in quads:
                terms = '+'.join(f'{v}*c({k+1})*c({l+1})' for (k, l), v in cf.items() if k >= ng and l >= ng)
                if terms:
                    f.write(f'I = I, {terms};\n')
            f.write('I = simplify(I,2);\nideal G = std(I);\nprint("dim = " + string(dim(G)) + "  (0 means only c=0 over the algebraic closure, projectively empty)");\n')
            f.write('print("vdim = " + string(vdim(G)));\nif (dim(G) > 0) { print("--- components (minAssGTZ) ---"); LIB "primdec.lib"; list L = minAssGTZ(G); int i; for (i=1;i<=size(L);i++){ print(L[i]); print("----"); } }\nquit;\n')
        with open(out_prefix + f'_quads{tag}.ms', 'w') as f:
            f.write(','.join(f'c{k+1}' for k in range(ng, n)) + '\n' + str(p) + '\n')
            lines = []
            for cf in quads:
                terms = '+'.join(f'{v}*c{k+1}*c{l+1}' for (k, l), v in cf.items() if k >= ng and l >= ng)
                if terms:
                    lines.append(terms)
            f.write(',\n'.join(lines) + '\n')
        with open(out_prefix + f'_family{tag}.json', 'w') as f:
            json.dump({'res': res, 'D': [core.to_str({e: c for e, c in Dp.items()}) for Dp in Ds]}, f)
    return res, Ds, Ns, quads
