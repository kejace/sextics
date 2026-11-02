"""Wave-2 shared machinery for lane rank-one-biquartic (q = 1, A = X0^3 Y0^3).

Normalized E2 ("complete correction criterion", book 5.2).  For an E1-admissible D (A D in J(F))
with one lift N0 (N0(F) = A D, q = 1) a certificate exists iff there are Koszul coefficients
c_ij (deg d-4, parity restricted) and C in S_{2d-6}^{par A} with

    B_D(N0) + T_D(c) + F C = 0,      B_D(N) = D div N - N(D),   T_D(c) = B_D(K(c)),
    K(c) = sum_{i<j} c_ij (F_j d_i - F_i d_j),
    T_D(c) = sum_{i<j} D (d_i c_ij F_j - d_j c_ij F_i) - c_ij (F_j D_i - F_i D_j).

Float version (variable-projection search over a linear family D = sum a_m D_m with lifts
N0 = sum a_m N_m) and exact mod-p helpers (single denominators, relaxation ranks).
Variables X0,X1,X2,Y0,Y1,Y2; Y-parity of a monomial = (e3+e4+e5) % 2.
"""
import sys, os, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, power

IID = 'rank-one-biquartic'
LANE = f'/tmp/hodge/results/{IID}'
OUT = f'{LANE}/wave2/out'
os.makedirs(OUT, exist_ok=True)


def par_idx(i):
    return 0 if i < 3 else 1


def parity_of(P):
    return core.parity_of(P)


def log_line(path, rec):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(path, 'a') as fh:
        fh.write(json.dumps(rec, default=str) + '\n')


# ----------------------------------------------------------------------------
# monomial spaces and dense vectors
# ----------------------------------------------------------------------------
class Space:
    """Monomials of one degree and Y-parity (par None = all), with dense-vector helpers."""

    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par) if deg >= 0 else []
        self.idx = {m: k for k, m in enumerate(self.mons)}
        self.n = len(self.mons)
        self._keys = None

    def vec(self, P):
        v = np.zeros(self.n)
        for e, c in P.items():
            v[self.idx[e]] = float(c)
        return v

    def vec_mod(self, P, p):
        v = [0] * self.n
        for e, c in P.items():
            v[self.idx[e]] = int(c) % p
        return v

    def poly(self, v, tol=0.0):
        return {m: float(v[k]) for k, m in enumerate(self.mons) if abs(v[k]) > tol}

    def keys(self, base):
        if self._keys is None or self._keys[0] != base:
            E = np.array(self.mons, dtype=np.int64).reshape(-1, 6)
            self._keys = (base, (E * (base ** np.arange(6, dtype=np.int64))).sum(1))
        return self._keys[1]


def mult_op(P, src, dst):
    """Sparse dst.n x src.n matrix of multiplication by the fixed polynomial P (float coefficients)."""
    rows, cols, vals = [], [], []
    for k, m in enumerate(src.mons):
        for e, c in P.items():
            s = (m[0]+e[0], m[1]+e[1], m[2]+e[2], m[3]+e[3], m[4]+e[4], m[5]+e[5])
            rows.append(dst.idx[s]); cols.append(k); vals.append(float(c))
    return sp.csr_matrix((vals, (rows, cols)), shape=(dst.n, src.n))


def diff_op(i, src, dst):
    rows, cols, vals = [], [], []
    for k, m in enumerate(src.mons):
        if m[i]:
            e = list(m); e[i] -= 1
            rows.append(dst.idx[tuple(e)]); cols.append(k); vals.append(float(m[i]))
    return sp.csr_matrix((vals, (rows, cols)), shape=(dst.n, src.n))


def lookup_index(space, keys_query, base):
    """Map integer monomial keys (same base) to indices in `space` (-1 if absent)."""
    ks = space.keys(base)
    order = np.argsort(ks)
    sk = ks[order]
    pos = np.searchsorted(sk, keys_query)
    pos = np.clip(pos, 0, len(sk) - 1)
    ok = sk[pos] == keys_query
    out = np.where(ok, order[pos], -1)
    return out


class ProdTable:
    """Bilinear product of dense vectors  S_a x S_b -> S_c  via a precomputed index table."""

    def __init__(self, sa, sb, sc):
        self.sa, self.sb, self.sc = sa, sb, sc
        base = sc.deg + 1
        ka, kb = sa.keys(base), sb.keys(base)
        S = (ka[:, None] + kb[None, :]).ravel()
        idx = lookup_index(sc, S, base)
        assert (idx >= 0).all(), 'product table: monomial outside target space (parity/degree mismatch)'
        self.idx = idx.astype(np.int64).reshape(sa.n, sb.n)
        self.flat = self.idx.ravel()

    def prod(self, u, v):
        return np.bincount(self.flat, weights=np.multiply.outer(u, v).ravel(), minlength=self.sc.n)

    def corr_a(self, w, v):
        """vector over S_a:  R[k] = sum_l w[idx[k,l]] v[l]  (adjoint of u -> u*v)."""
        return w[self.idx] @ v

    def corr_b(self, w, u):
        return w[self.idx].T @ u


# ----------------------------------------------------------------------------
# instance structure
# ----------------------------------------------------------------------------
def load(iid=None):
    F, A, q, rec = core.load_instance(iid or os.environ.get('HODGE_IID', IID))
    assert q == 1
    return F, A, q


def structure(F):
    P = core.parse
    f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
    g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
    L = P('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2')
    M = P('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
    t = Fraction(5, 7)
    Pc = mul(mul(L, L), M)
    assert sub(F, add(add(f, g), scale(Pc, t))) == {}
    return dict(f=f, g=g, L=L, M=M, t=t, P=Pc)


def euler_part(P, idxs):
    out = {}
    for i in idxs:
        out = add(out, mul({tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}, pderiv(P, i)))
    return out


# ----------------------------------------------------------------------------
# exact (dict) helpers
# ----------------------------------------------------------------------------
def B_of(D, N, p=None):
    """B_D(N) = D div N - N(D) (exact / mod p)."""
    divN = {}
    ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)


def koszul_field(F, c, p=None):
    """K(c) for c = {(i,j): poly}, i<j:  component i += c_ij F_j, component j -= c_ij F_i."""
    Fi = [pderiv(F, i, p) for i in range(NV)]
    N = [{} for _ in range(NV)]
    for (i, j), cij in c.items():
        if not cij:
            continue
        N[i] = add(N[i], mul(cij, Fi[j], p), p)
        N[j] = sub(N[j], mul(cij, Fi[i], p), p)
    return N


def normalize_lift(F, A, D, N, q=1):
    """Given a cleared certificate (D, N) with N(F) = A D + F h, return N' = N - (h/6) Z with N'(F) = A D exactly."""
    NF = {}
    for i in range(NV):
        NF = add(NF, mul(N[i], pderiv(F, i)))
    res = sub(scale(NF, q), mul(A, D))
    if not res:
        return N
    ok, h = jacobian.divisible_by(F, res)  # mod p -> need exact: do exact division instead
    # exact division by F via linear algebra over Q (small): use fixed-degree solve
    dG = core.degree(res) - 6
    par = (core.parity_of(res) - core.parity_of(F)) % 2
    gm = monomials(dG, par)
    cols = [shift(F, m, 1) for m in gm]
    rows = {}
    for col in cols + [res]:
        for e in col:
            rows.setdefault(e, len(rows))
    from flint import fmpq_mat, fmpq
    R, C = len(rows), len(cols)
    Mq = fmpq_mat(R, C)
    bq = fmpq_mat(R, 1)
    for j, col in enumerate(cols):
        for e, v in col.items():
            v = Fraction(v); Mq[rows[e], j] = fmpq(v.numerator, v.denominator)
    for e, v in res.items():
        v = Fraction(v); bq[rows[e], 0] = fmpq(v.numerator, v.denominator)
    # least squares via normal equations exact
    Mt = Mq.transpose()
    x = (Mt * Mq).solve(Mt * bq)
    hpoly = {}
    for j, m in enumerate(gm):
        v = x[j, 0]
        v = Fraction(int(v.p), int(v.q))
        if v:
            hpoly[m] = v
    assert sub(mul(F, hpoly), res) == {}, 'residual not divisible by F'
    hz = scale(hpoly, Fraction(1, 6 * q))
    N2 = []
    for i in range(NV):
        zi = {tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}
        N2.append(sub(N[i], mul(hz, zi)))
    return N2


# ----------------------------------------------------------------------------
# normalized E2 problem (float), q = 1
# ----------------------------------------------------------------------------
class E2Float:
    """Assembles the normalized E2 system for denominators in S_d^{pD}; q = 1.

    Per evaluation (given D as a dense vector over S_d^{pD} and N0 over the N-spaces):
        B = D div N0 - N0(D)               (dense vector over Srow = S_{2d}^{pA})
        M(D) = Koszul column block          (Srow.n x n_c)
        F-reduction rho: Srow -> standard monomials (not divisible by X0^6)
    """

    def __init__(self, F, A, d, pD, verbose=True, want_reduction=True):
        t0 = time.time()
        self.F, self.A, self.d, self.pD = F, A, d, pD
        pA = core.parity_of(A); self.pA = pA
        self.pN = [(pA - pD - par_idx(i)) % 2 for i in range(NV)]
        self.pairs = [(i, j) for i in range(NV) for j in range(i + 1, NV)]
        self.pc = {(i, j): (pA - pD + par_idx(i) + par_idx(j)) % 2 for (i, j) in self.pairs}
        self.SD = Space(d, pD)
        self.SN = [Space(d + 1, self.pN[i]) for i in range(NV)]
        self.nN = sum(s.n for s in self.SN)
        self.Noff = np.cumsum([0] + [s.n for s in self.SN])
        self.Srow = Space(2 * d, pA)
        self.SC = Space(2 * d - 6, pA)
        self.Sc = {ij: Space(d - 4, self.pc[ij]) for ij in self.pairs}
        self.coff = {}
        off = 0
        for ij in self.pairs:
            self.coff[ij] = off; off += self.Sc[ij].n
        self.n_c = off
        self.Fi = [pderiv(F, i) for i in range(NV)]
        # derivative spaces of D and N
        self.SDi = [Space(d - 1, (pD + par_idx(i)) % 2) for i in range(NV)]
        self.SdivN = Space(d, (pA - pD) % 2)  # parity(d_i N_i) = pN_i + par_i = pA - pD
        # operators on D
        self.dD = [diff_op(i, self.SD, self.SDi[i]) for i in range(NV)]
        self.SDF = [Space(d + 5, (pD + par_idx(j)) % 2) for j in range(NV)]
        self.mDF = [mult_op(self.Fi[j], self.SD, self.SDF[j]) for j in range(NV)]
        self.SP = {}
        self.mP = {}
        for (i, j) in self.pairs:
            spc = Space(d + 4, (pD + par_idx(i) + par_idx(j)) % 2)
            self.SP[(i, j)] = spc
            op = mult_op(self.Fi[j], self.SDi[i], spc) @ self.dD[i] - mult_op(self.Fi[i], self.SDi[j], spc) @ self.dD[j]
            self.mP[(i, j)] = op.tocsr()
        # operators on N: divergence pieces and Koszul map c -> N
        self.dN = [diff_op(i, self.SN[i], self.SdivN) for i in range(NV)]
        kos_rows, kos_cols, kos_vals = [], [], []
        for (i, j) in self.pairs:
            sc = self.Sc[(i, j)]
            # component i += c_ij F_j ; component j -= c_ij F_i
            opi = mult_op(self.Fi[j], sc, self.SN[i]).tocoo()
            opj = mult_op(self.Fi[i], sc, self.SN[j]).tocoo()
            kos_rows += list(opi.row + self.Noff[i]); kos_cols += list(opi.col + self.coff[(i, j)]); kos_vals += list(opi.data)
            kos_rows += list(opj.row + self.Noff[j]); kos_cols += list(opj.col + self.coff[(i, j)]); kos_vals += list(-opj.data)
        self.Kmap = sp.csr_matrix((kos_vals, (kos_rows, kos_cols)), shape=(self.nN, self.n_c))
        # product tables for B
        self.T_DdivN = ProdTable(self.SD, self.SdivN, self.Srow)
        self.T_NiDi = [ProdTable(self.SN[i], self.SDi[i], self.Srow) for i in range(NV)]
        # Koszul column assembly data
        self._build_columns()
        # F reduction
        self.rho_std, self.rho_red_rows, self.R = None, None, None
        if want_reduction:
            self._build_reduction()
        self.build_time = time.time() - t0
        if verbose:
            print(f'[E2Float d={d} pD={pD}] rows={self.Srow.n} n_c={self.n_c} n_C={self.SC.n} nN={self.nN} '
                  f'n_red={self.n_red} col_nnz={len(self.col_rows)} built {self.build_time:.1f}s', file=sys.stderr, flush=True)

    def _build_columns(self):
        d = self.d
        base = 2 * d + 1
        rows_all, cols_all, src_all, coef_all = [], [], [], []
        # source vector layout: DF_0..DF_5 then P_ij
        self.src_off = {}
        off = 0
        for j in range(NV):
            self.src_off[('DF', j)] = off; off += self.SDF[j].n
        for ij in self.pairs:
            self.src_off[('P', ij)] = off; off += self.SP[ij].n
        self.n_src = off
        rowkeys = self.Srow.keys(base)
        for (i, j) in self.pairs:
            sc = self.Sc[(i, j)]
            kDFj = self.SDF[j].keys(base); kDFi = self.SDF[i].keys(base); kP = self.SP[(i, j)].keys(base)
            for kc, mc in enumerate(sc.mons):
                col = self.coff[(i, j)] + kc
                mkey = int(np.dot(np.array(mc, dtype=np.int64), base ** np.arange(6, dtype=np.int64)))
                ei = base ** i; ej = base ** j
                # term1: + mc[i] * shift_{mc - e_i}(DF_j)
                if mc[i]:
                    keys = kDFj + (mkey - ei)
                    r = lookup_index(self.Srow, keys, base); assert (r >= 0).all()
                    rows_all.append(r); cols_all.append(np.full(len(r), col, dtype=np.int64))
                    src_all.append(np.arange(self.SDF[j].n) + self.src_off[('DF', j)]); coef_all.append(np.full(len(r), float(mc[i])))
                if mc[j]:
                    keys = kDFi + (mkey - ej)
                    r = lookup_index(self.Srow, keys, base); assert (r >= 0).all()
                    rows_all.append(r); cols_all.append(np.full(len(r), col, dtype=np.int64))
                    src_all.append(np.arange(self.SDF[i].n) + self.src_off[('DF', i)]); coef_all.append(np.full(len(r), -float(mc[j])))
                keys = kP + mkey
                r = lookup_index(self.Srow, keys, base); assert (r >= 0).all()
                rows_all.append(r); cols_all.append(np.full(len(r), col, dtype=np.int64))
                src_all.append(np.arange(self.SP[(i, j)].n) + self.src_off[('P', (i, j))]); coef_all.append(np.full(len(r), -1.0))
        if rows_all:
            self.col_rows = np.concatenate(rows_all); cols = np.concatenate(cols_all)
            self.col_src = np.concatenate(src_all); self.col_coef = np.concatenate(coef_all)
            self.col_flat = self.col_rows * self.n_c + cols
        else:
            self.col_rows = np.zeros(0, dtype=np.int64); self.col_src = np.zeros(0, dtype=np.int64)
            self.col_coef = np.zeros(0); self.col_flat = np.zeros(0, dtype=np.int64)

    def _build_reduction(self):
        """Normal form modulo F w.r.t. lex (X0 > X1 > ... > Y2): rows -> monomials not divisible by LM(F)."""
        Srow = self.Srow
        LM = max(self.F)  # lex-largest exponent tuple
        lc = float(self.F[LM])
        Fother = {e: -float(c) / lc for e, c in self.F.items() if e != LM}
        self.LM = LM

        def divisible(m):
            return all(m[t] >= LM[t] for t in range(6))
        std = [k for k, m in enumerate(Srow.mons) if not divisible(m)]
        red = [k for k, m in enumerate(Srow.mons) if divisible(m)]
        self.n_red = len(std)
        std_pos = {k: t for t, k in enumerate(std)}
        self.rho_std = np.array(std, dtype=np.int64)
        self.rho_red_rows = np.array(red, dtype=np.int64)
        if not red:
            self.R = np.zeros((self.n_red, 0)); self.rho_nnz = 0; return
        nf_cache = {}
        order = sorted(red, key=lambda k: Srow.mons[k])  # lex ascending: t*m' < LM*m' is processed first
        for k in order:
            m = Srow.mons[k]
            mp = tuple(m[t] - LM[t] for t in range(6))
            acc = {}
            for e, c in Fother.items():
                s = (mp[0]+e[0], mp[1]+e[1], mp[2]+e[2], mp[3]+e[3], mp[4]+e[4], mp[5]+e[5])
                ks = Srow.idx[s]
                if not divisible(s):
                    acc[std_pos[ks]] = acc.get(std_pos[ks], 0.0) + c
                else:
                    for kk, v in nf_cache[ks].items():
                        acc[kk] = acc.get(kk, 0.0) + c * v
            nf_cache[k] = acc
        R = np.zeros((self.n_red, len(red)))
        for t, k in enumerate(red):
            for kk, v in nf_cache[k].items():
                R[kk, t] = v
        self.R = R
        self.rho_nnz = int((R != 0).sum())

    # -- evaluation pieces -------------------------------------------------------------------
    def sources(self, Dv):
        src = np.zeros(self.n_src)
        for j in range(NV):
            o = self.src_off[('DF', j)]
            src[o:o + self.SDF[j].n] = self.mDF[j] @ Dv
        for ij in self.pairs:
            o = self.src_off[('P', ij)]
            src[o:o + self.SP[ij].n] = self.mP[ij] @ Dv
        return src

    def koszul_matrix(self, Dv):
        """Dense Srow.n x n_c matrix of c -> T_D(c)."""
        src = self.sources(Dv)
        vals = self.col_coef * src[self.col_src]
        M = np.bincount(self.col_flat, weights=vals, minlength=self.Srow.n * self.n_c)
        return M.reshape(self.Srow.n, self.n_c)

    def Bvec(self, Dv, Nv):
        """B_D(N) as a dense vector over Srow."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ Nv[self.Noff[i]:self.Noff[i + 1]]
        out = self.T_DdivN.prod(Dv, divN)
        for i in range(NV):
            Di = self.dD[i] @ Dv
            out -= self.T_NiDi[i].prod(Nv[self.Noff[i]:self.Noff[i + 1]], Di)
        return out

    def Bpieces_norm(self, Dv, Nv):
        """Sum of norms of the reduced pieces rho(D div N), rho(N_i D_i): a cancellation-free scale for residuals."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ Nv[self.Noff[i]:self.Noff[i + 1]]
        tot = np.linalg.norm(self.reduce(self.T_DdivN.prod(Dv, divN)))
        for i in range(NV):
            Di = self.dD[i] @ Dv
            tot += np.linalg.norm(self.reduce(self.T_NiDi[i].prod(Nv[self.Noff[i]:self.Noff[i + 1]], Di)))
        return float(tot)

    def reduce(self, X):
        """Apply the F-reduction rho to a vector or matrix (rows over Srow)."""
        if X.ndim == 1:
            return X[self.rho_std] + self.R @ X[self.rho_red_rows]
        return X[self.rho_std, :] + self.R @ X[self.rho_red_rows, :]

    def reduce_T(self, Y):
        """Transpose of rho applied to a reduced vector: returns vector over Srow."""
        out = np.zeros(self.Srow.n)
        out[self.rho_std] = Y
        out[self.rho_red_rows] = self.R.T @ Y
        return out

    def adjoint_B(self, w, Dv, Nv):
        """Gradients of <w, B(D, N)> with respect to D (vector over SD) and N (vector over N-spaces)."""
        divN = np.zeros(self.SdivN.n)
        for i in range(NV):
            divN += self.dN[i] @ Nv[self.Noff[i]:self.Noff[i + 1]]
        gD = self.T_DdivN.corr_a(w, divN)
        gN = np.zeros(self.nN)
        cD = self.T_DdivN.corr_b(w, Dv)  # over SdivN
        for i in range(NV):
            gN[self.Noff[i]:self.Noff[i + 1]] += self.dN[i].T @ cD
            Di = self.dD[i] @ Dv
            Ni = Nv[self.Noff[i]:self.Noff[i + 1]]
            gD -= self.dD[i].T @ self.T_NiDi[i].corr_b(w, Ni)
            gN[self.Noff[i]:self.Noff[i + 1]] -= self.T_NiDi[i].corr_a(w, Di)
        return gD, gN

    def Nvec(self, N):
        v = np.zeros(self.nN)
        for i in range(NV):
            v[self.Noff[i]:self.Noff[i + 1]] = self.SN[i].vec(N[i])
        return v

    def Cmatrix(self):
        """Dense Srow.n x SC.n block of C -> F C (for un-reduced solves)."""
        return mult_op(self.F, self.SC, self.Srow).toarray()

    def residual_full(self, Dv, Nv, cond=1e-10):
        """Exact-norm least squares with explicit C block: returns (relres, c, C)."""
        B = self.Bvec(Dv, Nv)
        M = self.koszul_matrix(Dv)
        MC = np.hstack([M, self.Cmatrix()])
        x, *_ = sla.lstsq(MC, -B, cond=cond, lapack_driver='gelsd', check_finite=False)
        r = B + MC @ x
        return float(np.linalg.norm(r) / max(np.linalg.norm(B), 1e-300)), x[:self.n_c], x[self.n_c:]


# ----------------------------------------------------------------------------
# E1-admissible family in floats (J_d part + kernel part) with lifts
# ----------------------------------------------------------------------------
class E1Family:
    """Orthonormal basis of {D in S_d^{pD} : A D in J(F)} with lifts N (N(F) = A D), computed in floats.

    Uses the SVD of the Jacobian-generator matrix in degree d+6 (parity pA+pD): the admissible space is the
    kernel of  Q_perp^T (A . m)  where Q_perp spans the orthogonal complement of the column space.
    """

    def __init__(self, prob: E2Float, rcond=1e-9, verbose=True):
        t0 = time.time()
        F, A, d, pD = prob.F, prob.A, prob.d, prob.pD
        self.prob = prob
        pT = (prob.pA + pD) % 2
        ST = Space(d + 6, pT)
        self.ST = ST
        # Jacobian generator matrix: columns N_i monomials -> F_i m
        blocks = [mult_op(prob.Fi[i], prob.SN[i], ST) for i in range(NV)]
        Jm = sp.hstack(blocks).toarray()
        self.Jm_shape = Jm.shape
        U, s, Vt = sla.svd(Jm, full_matrices=True, check_finite=False, lapack_driver='gesdd')
        tol = rcond * s[0]
        rank = int((s > tol).sum())
        self.J_rank = rank
        self.J_gap = (float(s[rank - 1]), float(s[rank]) if rank < len(s) else 0.0)
        Qperp = U[:, rank:]                       # ST.n x (ST.n - rank)
        self.ring_dim = ST.n - rank
        # A . m for m in SD
        AM = mult_op(A, prob.SD, ST).toarray()    # ST.n x SD.n
        C = Qperp.T @ AM                          # ring_dim x SD.n
        ns = sla.null_space(C, rcond=rcond)       # SD.n x K
        self.K = ns.shape[1]
        self.Dbasis = ns.T                        # K x SD.n (orthonormal rows)
        # singular values of C near the cut (diagnostics)
        sC = sla.svdvals(C)
        self.C_gap = (float(sC[-self.K - 1]) if self.K < len(sC) else None, float(sC[-self.K]) if self.K > 0 else None)
        # lifts: solve Jm N = A D  (min norm) via the SVD pieces
        rhs = AM @ self.Dbasis.T                  # ST.n x K
        Ur = U[:, :rank]; sr = s[:rank]; Vr = Vt[:rank, :]
        self.Nbasis = (Vr.T @ ((Ur.T @ rhs) / sr[:, None])).T   # K x nN
        e1res = Jm @ self.Nbasis.T - rhs
        self.e1_rel = float(np.linalg.norm(e1res) / max(np.linalg.norm(rhs), 1e-300))
        # J_d part identification: dimension of the family inside J_d (E1 trivial)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[E1Family d={d} pD={pD}] Jm {Jm.shape} rank={rank} ring_dim={self.ring_dim} gap={self.J_gap} '
                  f'K={self.K} C_gap={self.C_gap} e1_rel={self.e1_rel:.2e} ({self.build_time:.1f}s)', file=sys.stderr, flush=True)

    def spurious_subspace(self):
        """Coordinates (in the K-basis) of F * S_{d-6}^{pD}; returns orthonormal K x n_sp matrix (columns)."""
        prob = self.prob
        d = prob.d
        if d < 6:
            return np.zeros((self.K, 0))
        SG = Space(d - 6, (prob.pD - core.parity_of(prob.F)) % 2)
        if SG.n == 0:
            return np.zeros((self.K, 0))
        FG = mult_op(prob.F, SG, prob.SD).toarray()   # SD.n x SG.n
        coords = self.Dbasis @ FG                     # K x SG.n
        # check these lie in the family
        recon = self.Dbasis.T @ coords
        self.spurious_fit = float(np.linalg.norm(recon - FG) / np.linalg.norm(FG))
        Q, _ = np.linalg.qr(coords)
        return Q

    def J_part_dim(self):
        """Dimension of J_d^{pD} (numerically) and of its intersection with the family."""
        prob = self.prob
        d = prob.d
        if d < 5:
            return 0
        gens = []
        for i in range(NV):
            par = (prob.pD - par_idx(i)) % 2
            Sm = Space(d - 5, par)
            if Sm.n:
                gens.append(mult_op(prob.Fi[i], Sm, prob.SD).toarray())
        G = np.hstack(gens)
        r = np.linalg.matrix_rank(G, tol=1e-9 * np.linalg.norm(G, 2))
        return int(r)


# ----------------------------------------------------------------------------
# restricted sub-families of the admissible space (float)
# ----------------------------------------------------------------------------
def ybase_polys(F, d, par):
    """Y-base-type polynomials P(X) M(Y)^a g(Y)^b of degree d and parity par (the packet's compact-Y class)."""
    st = structure(F)
    M, g = st['M'], st['g']
    out = []
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0:
                continue
            W = power(M, a) if a else {(0,) * 6: Fraction(1)}
            if b:
                W = mul(W, g)
            for m in monomials(k):
                if m[3] or m[4] or m[5]:
                    continue
                P_ = mul({m: Fraction(1)}, W)
                if core.parity_of(P_) == par:
                    out.append(P_)
    return out


def restricted_family(prob, fam, kind):
    """Sub-family (Dbasis, Nbasis, info) of the admissible family `fam`:
       'kernel'       : complement of J_d inside the admissible space (pure E1-kernel directions)
       'kernel-mixed' : kernel directions orthogonal to the Y-base-type span P(X) M^a g^b as well
    """
    d, pD = prob.d, prob.pD
    gens = []
    for i in range(NV):
        par = (pD - par_idx(i)) % 2
        Sm = Space(d - 5, par)
        if Sm.n:
            gens.append(mult_op(prob.Fi[i], Sm, prob.SD).toarray())
    Q_list = []
    rJ = 0
    if gens:
        G = np.hstack(gens)
        QJ, sJ, _ = sla.svd(G, full_matrices=False)
        rJ = int((sJ > 1e-9 * sJ[0]).sum())
        Q_list.append(QJ[:, :rJ])
    rU = 0
    if kind == 'kernel-mixed':
        U = np.array([prob.SD.vec(P_) for P_ in ybase_polys(prob.F, d, pD)]).T
        QU, sU, _ = sla.svd(U, full_matrices=False)
        rU = int((sU > 1e-9 * sU[0]).sum())
        Q_list.append(QU[:, :rU])
    Pperp = fam.Dbasis.T.copy()
    for Qx in Q_list:
        Pperp = Pperp - Qx @ (Qx.T @ Pperp)
    Uf, s, Vt = sla.svd(Pperp, full_matrices=False)
    rk = int((s > 1e-7 * s[0]).sum()) if s.size else 0
    Dbasis = Uf[:, :rk].T
    coords = fam.Dbasis @ Dbasis.T
    Nbasis = (fam.Nbasis.T @ coords).T
    info = {'K': int(rk), 'K_full': fam.K, 'dim_J_d': int(rJ), 'dim_ybase_span': int(rU),
            'sv_gap': [float(s[rk - 1]) if rk else None, float(s[rk]) if rk < len(s) else None]}
    return Dbasis, Nbasis, info
