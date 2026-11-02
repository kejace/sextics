"""Wave-2 EXCLUSION library (lane rank-one-biquartic).

Rigorous low-degree exclusion by the Jacobian-type reduction:
  * Euler gauge: any certificate (D,N,Q) with D coprime to F can be normalized to N(F) = A D exactly
    (replace N by N - (h/6) Euler, which leaves the cleared identity invariant).
  * Then N = N0 + K(c) with K(c)_i = sum_j c_ij F_j, c antisymmetric (syzygies of the regular sequence
    F_0..F_5 are Koszul), deg c_ij = d-4, and the identity reduces to  E2:  D div N - N(D) = F Q.
  * D ranges over the E1-kernel V_d = {D in S_d : A D in J(F)} (all parities together: a mixed-parity D is
    allowed by the contract), parametrized by an exact Q-basis (D_k, N_k) with N_k(F) = A D_k.
  * E2 is quadratic in (a, c) (a = coordinates of D) and linear in Q: weight-2 homogeneous with weights
    a:1, c:1, Q:2.  The spurious locus Sp = F * S_{d-6} (not coprime) is a union of basis directions.
  * Exclusion certificate: for each coordinate u of V_d / Sp, an exact rational combination of the E2
    equations equal to u^2 (weight 2) proves V(I) subset {u = 0}, i.e. no certificate of degree d.
Everything exact: Fractions for polynomials, python-flint nmod_mat for modular ranks, fmpq_mat for exact solves.
"""
import sys, time, json, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, yparity
from flint import nmod_mat, fmpq_mat, fmpz_mat, fmpq, fmpz

LANE = 'rank-one-biquartic'
WORK = f'/tmp/hodge/results/{LANE}'
PK = f'{WORK}/wave2-excl'
F, A, q, rec = core.load_instance(LANE)
assert q == 1
Fi = [pderiv(F, i) for i in range(NV)]
P = core.parse
t = Fraction(5, 7)
f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
L = P('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2')
M = P('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
L2M = mul(mul(L, L), M)
assert sub(F, add(add(f, g), scale(L2M, t))) == {}
EULER = [{tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)} for i in range(NV)]  # unit monomials z_i

def mono(e, c=1):
    return {tuple(e): Fraction(c)}

def divergence(N, p=None):
    out = {}
    for i in range(NV):
        if N[i]:
            out = add(out, pderiv(N[i], i, p), p)
    return out

def apply_field(N, D, p=None):
    """N(D) = sum_i N_i dD/dz_i"""
    out = {}
    for i in range(NV):
        if N[i]:
            out = add(out, mul(N[i], pderiv(D, i, p), p), p)
    return out

def B_op(D, N, p=None):
    """E2 operator  B_D(N) = D div N - N(D)."""
    return sub(mul(D, divergence(N, p), p), apply_field(N, D, p), p)

def dF(G, p=None):
    """polar dF(G) = sum_i G_i F_i"""
    out = {}
    for i in range(NV):
        if G[i]:
            out = add(out, mul(G[i], Fi[i] if p is None else core.reduce_mod(Fi[i], p), p), p)
    return out

def koszul_basis(d):
    """Basis of Koszul fields of degree d+1 (= deg N): for i<j and monomial m of degree d-4:
    K_i = m F_j, K_j = -m F_i.  Returns list of (label, N-vector)."""
    out = []
    if d - 4 < 0:
        return out
    for i in range(NV):
        for j in range(i + 1, NV):
            for m in monomials(d - 4):
                N = [{} for _ in range(NV)]
                N[i] = shift(Fi[j], m, 1)
                N[j] = shift(Fi[i], m, -1)
                out.append((('c', i, j, m), N))
    return out

def polar_basis(d, parity=None):
    """Basis of J_d: (label, D=dF(G), N=A*G) for G = m e_i, deg m = d-5 (all parities unless given)."""
    out = []
    if d - 5 < 0:
        return out
    for i in range(NV):
        for m in monomials(d - 5):
            G = [{} for _ in range(NV)]
            G[i] = mono(m)
            D = dF(G)
            if parity is not None and core.parity_of(D) != parity:
                continue
            N = [mul(A, G[k]) if G[k] else {} for k in range(NV)]
            out.append((('J', i, m), D, N))
    return out

class FamilySystem:
    """Weight-2 Macaulay matrix of the E2 system for the family D = sum a_k D_k, N = sum a_k N_k + K(c), Q.

    basis: list of (label, D_k, N_k) exact over Q with N_k(F) = A D_k.
    Columns: ('aa',k,l) k<=l ; ('ac',k,j) ; ('Q',m).  Rows: monomials of degree 2d (all parities).
    """
    def __init__(self, d, basis, parity_D=None, verbose=True, check=True):
        t0 = time.time()
        self.d = d
        self.basis = basis
        K = len(basis)
        self.K = K
        if check:
            for lab, D, N in basis:
                assert core.degree(D) == d, (lab, core.degree(D))
                lhs = apply_field(N, F)
                assert sub(lhs, mul(A, D)) == {}, f'E1 fails for basis element {lab}'
        self.kos = koszul_basis(d)
        self.Qmons = monomials(2 * d - 6) if 2 * d - 6 >= 0 else []
        # columns as sparse polys of degree 2d
        cols = []
        colpolys = []
        # a_k a_l terms: P_kl = D_k div N_l - N_l(D_k); coefficient of a_k a_l (k<l) is P_kl + P_lk
        Pm = {}
        for k in range(K):
            for l in range(K):
                Pm[(k, l)] = B_op(basis[k][1], basis[l][2])
        for k in range(K):
            for l in range(k, K):
                poly = Pm[(k, k)] if k == l else add(Pm[(k, l)], Pm[(l, k)])
                cols.append(('aa', k, l)); colpolys.append(poly)
        t1 = time.time()
        for k in range(K):
            Dk = basis[k][1]
            for j, (lab, Kf) in enumerate(self.kos):
                cols.append(('ac', k, j)); colpolys.append(B_op(Dk, Kf))
        t2 = time.time()
        for m in self.Qmons:
            cols.append(('Q', m)); colpolys.append(shift(F, m, -1))
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows:
                    rows[e] = len(rows)
        self.cols, self.colpolys, self.rows = cols, colpolys, rows
        self.colidx = {c: j for j, c in enumerate(cols)}
        self.nrows, self.ncols = len(rows), len(cols)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[family d={d}] K={K} koszul={len(self.kos)} Q={len(self.Qmons)} rows={self.nrows} cols={self.ncols} '
                  f'nnz={sum(len(c) for c in colpolys)} (aa {t1-t0:.0f}s, ac {t2-t1:.0f}s, total {self.build_time:.0f}s)', flush=True)

    def nmod(self, p):
        R, C = self.nrows, self.ncols
        flat = [0] * (R * C)
        rows = self.rows
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                v = Fraction(v)
                flat[rows[e] * C + j] = (v.numerator * pow(v.denominator, -1, p)) % p
        return nmod_mat(R, C, flat, p)

    def rowspace_analysis(self, p, verbose=True):
        """rref mod p; returns dict with rank, pivots, and for every 'aa' column whether its unit vector lies in the
        row space (i.e. the monomial a_k a_l is in I_2 mod p), plus the rref restricted to aa-columns of the rows whose
        pivots are aa-columns (the quadratic forms in a that lie in I_2 mod p, modulo the others)."""
        t0 = time.time()
        Mx = self.nmod(p)
        t1 = time.time()
        Mr, rank = Mx.rref()
        t2 = time.time()
        C = self.ncols
        pivots = []
        for i in range(rank):
            j = 0
            while j < C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        pivset = set(pivots)
        aa_cols = [j for j, c in enumerate(self.cols) if c[0] == 'aa']
        unit_in = {}
        # rows whose pivot is an aa column and which are supported only on aa columns -> quadratic forms in a inside I_2
        quad_forms = []
        for i, j in enumerate(pivots):
            if self.cols[j][0] != 'aa':
                continue
            row = {k: int(Mr[i, k]) for k in range(j, C) if int(Mr[i, k])}
            if all(self.cols[k][0] == 'aa' for k in row):
                quad_forms.append({self.cols[k]: v for k, v in row.items()})
                if len(row) == 1:
                    unit_in[self.cols[j]] = True
        info = dict(p=p, rank=rank, nrows=self.nrows, ncols=self.ncols, n_aa=len(aa_cols),
                    n_unit_aa_in_I2=len(unit_in), n_quad_forms_in_I2=len(quad_forms),
                    t_build=round(t1 - t0, 1), t_rref=round(t2 - t1, 1))
        if verbose:
            print('[rowspace]', json.dumps(info), flush=True)
        return info, pivots, unit_in, quad_forms, Mr

def exact_certificate(FS, pivots, Mr_unused, targets, p, verbose=True, max_rows=None):
    """Exact rational certificates.  targets: list of (name, tvec) with tvec a dict {column index: Fraction} describing
    a weight-2 polynomial in the unknowns (e.g. u^2 for a coordinate u).  Finds lambda in Q^rows with lambda^T E = tvec
    (so the polynomial lies in the ideal generated by the E2 equations), using the pivot rows/cols found mod p
    (square submatrix invertible mod p, hence over Q), an exact fmpq solve, and an exact verification on all columns.
    Returns dict name -> (lambda as sparse dict row->Fraction) or None."""
    t0 = time.time()
    R, C = FS.nrows, FS.ncols
    Mx = FS.nmod(p)
    MT = Mx.transpose()
    MTr, r2 = MT.rref()
    prow = []
    for i in range(r2):
        j = 0
        while j < R and int(MTr[i, j]) == 0:
            j += 1
        prow.append(j)
    r = len(pivots)
    assert r2 == r, (r2, r)
    del Mx, MT, MTr
    S = fmpq_mat(r, r)
    rowpos = {ri: a for a, ri in enumerate(prow)}
    colpos = {cj: b for b, cj in enumerate(pivots)}
    rows = FS.rows
    for b, cj in enumerate(pivots):
        for e, v in FS.colpolys[cj].items():
            ri = rows[e]
            if ri in rowpos:
                v = Fraction(v)
                S[rowpos[ri], b] = fmpq(v.numerator, v.denominator)
    T = fmpq_mat(r, len(targets))
    for tcol, (name, tvec) in enumerate(targets):
        for j0, v in tvec.items():
            if j0 in colpos:
                v = Fraction(v)
                T[colpos[j0], tcol] = fmpq(v.numerator, v.denominator)
            # a target component on a non-pivot column is fine: it will be checked in the verification
    t1 = time.time()
    try:
        MU = S.transpose().solve(T)
    except Exception as ex:
        print('[cert] solve failed', ex, flush=True)
        return None
    t2 = time.time()
    out = {}
    # verification: lambda^T E as an fmpq_mat product (1 x r) * (r x C) restricted to prow rows -> use sparse python for small r,
    # flint for large.
    Esub = None
    for tcol, (name, tvec) in enumerate(targets):
        lam = {}
        for a_, ri in enumerate(prow):
            v = MU[a_, tcol]
            if v != 0:
                lam[ri] = Fraction(int(v.p), int(v.q))
        ok = True
        for cj, col in enumerate(FS.colpolys):
            s = Fraction(0)
            for e, v in col.items():
                ri = rows[e]
                if ri in lam:
                    s += lam[ri] * Fraction(v)
            want = Fraction(tvec.get(cj, 0))
            if s != want:
                ok = False
                break
        out[name] = lam if ok else None
        if verbose:
            h = max((max(abs(x.numerator), x.denominator) for x in lam.values()), default=1)
            print(f'[cert] target {name} verified={ok} support={len(lam)} maxheight_digits={len(str(h))}', flush=True)
    if verbose:
        print(f'[cert] submatrix {r}x{r} build {t1-t0:.0f}s solve {t2-t1:.0f}s verify {time.time()-t2:.0f}s', flush=True)
    return out

def reduce_by_rref(Mr, pivots, tvec, p, C):
    """Is the vector tvec (dict col->Fraction) in the row space of the rref Mr (mod p)?"""
    v = [0] * C
    for j, x in tvec.items():
        x = Fraction(x); v[j] = (x.numerator * pow(x.denominator, -1, p)) % p
    for i, j in enumerate(pivots):
        c = v[j]
        if c:
            for k in range(j, C):
                m = int(Mr[i, k])
                if m:
                    v[k] = (v[k] - c * m) % p
    return not any(v)
