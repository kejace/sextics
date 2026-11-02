"""wave-2 core for lane cyclic-triple-residual: sector-restricted E1-admissible spaces (exact over Q) and a
Koszul-reduced variable-projection joint search (numerical) with rational reconstruction + exact confirmation.

Certificate split form (q = 1 here, kept general):  E1: q sum_i N_i F_i = A D ;  E2: D div N - N(D) = F^q Q.
Every E1 lift is N = N0 + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j), deg c_ij = d + 6q - 10.
For a linear family D(a) = sum_k a_k D_k with exact lifts N0_k (q N0_k(F) = A D_k) the E2 residual is
    r(a, c, Q) = Bil(D(a), N0(a) + K(c)) - F^q Q,   Bil(D, N) = D div N - N(D)   (bilinear),
linear in (c, Q) for fixed a: variable projection f(a) = min_{c,Q} |r|^2 on the unit sphere in a.
Mod-J ("Hessian") relaxation: f_J(a) = |P_perp^T Bil(D(a), N0(a))|^2 with P_perp an orthonormal basis of the
orthogonal complement of J(F)_{2d} in the E2 row space (a necessary condition, since im T and F^q S lie in J).

Character bookkeeping (sect.py conventions, all mod 6):  chi(N_i) = chiA + c + chiE_i,  chi(Q) = chiA + 2c,
chi(c_ij) = chiA + c + chiE_i + chiE_j,  E2 rows chi = chiA + 2c,  J_d^c generators m F_i with chi(m) = c + chiE_i.
"""
import os, sys, time, json, math, itertools
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from flint import fmpz_mat, nmod_mat
from hodge import core, fixed_d
from hodge.core import NV, mul, add, sub, scale, pderiv, shift
import sect

IID = 'cyclic-triple-residual'
WD = '/tmp/hodge/results/cyclic-triple-residual'
W2LOG = f'{WD}/wave2-search.log.jsonl'
PRIMES = core.PRIMES
E6 = sect.E6


def w2log(entry, path=W2LOG):
    entry = dict(entry)
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('engine', 2)
    entry.setdefault('lane', IID)
    with open(path, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')


def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q


def fermat():
    return core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')


# ----------------------------------------------------------------------------------------------------------
# sectors
# ----------------------------------------------------------------------------------------------------------
class Sect:
    """Character machinery for a subgroup `sub` (list of weight tuples mod 6, last entry 0) fixing F (char 0)."""

    def __init__(self, F, A, sub):
        self.sub = list(sub)
        self.gens = sect.generators(self.sub)
        self.chi = sect.make_chi(self.gens)
        self.chiA = sect.chi_of(sect.integral(A), self.chi)
        self.chiF = sect.chi_of(sect.integral(F), self.chi)
        self.chiE = [self.chi(E6[i]) for i in range(6)]
        self.order = len(self.sub)
        self._cache = {}

    def mons(self, deg, c):
        key = (deg, tuple(c))
        if key not in self._cache:
            self._cache[key] = [m for m in sect.mons(deg) if self.chi(m) == tuple(c)] if deg >= 0 else []
        return self._cache[key]

    def sectors(self, deg):
        out = {}
        for m in sect.mons(deg):
            out.setdefault(self.chi(m), 0)
            out[self.chi(m)] += 1
        return out

    def chiN(self, c, i):
        return sect.cadd(self.chiA, tuple(c), self.chiE[i])

    def chiQ(self, c):
        return sect.cadd(self.chiA, tuple(c), tuple(c))

    def chiC(self, c, i, j):
        return sect.cadd(self.chiA, tuple(c), self.chiE[i], self.chiE[j])

    def chiJgen(self, c, i):
        return sect.cadd(tuple(c), self.chiE[i])


def full_group(F, A):
    Fz, Az = sect.integral(F), sect.integral(A)
    return sect.semi_invariant_subgroup(sect.elems(Fz), [Az])


def parity_group():
    return [(0, 0, 0, 0, 0, 0), (3, 3, 3, 0, 0, 0)]


def csec(c):
    return ''.join(str(x) for x in c)


# ----------------------------------------------------------------------------------------------------------
# polynomial helpers (exact dicts)
# ----------------------------------------------------------------------------------------------------------
def poly_norm(P):
    return math.sqrt(sum(float(v) ** 2 for v in P.values()))


def to_float(P):
    return {e: float(v) for e, v in P.items()}


def padd(*Ps):
    out = {}
    for P in Ps:
        for e, v in P.items():
            out[e] = out.get(e, 0) + v
    return {e: v for e, v in out.items() if v}


def e1_check(F, A, q, D, N):
    """Exact E1 residual q sum N_i F_i - A D."""
    acc = {}
    for i in range(6):
        if N[i]:
            acc = add(acc, mul(N[i], pderiv(F, i)))
    return sub(scale(acc, q), mul(A, D))


def bil(D, N):
    """Bil(D, N) = D div N - N(D), exact."""
    divN, ND = {}, {}
    for i in range(6):
        if N[i]:
            divN = add(divN, pderiv(N[i], i))
            ND = add(ND, mul(N[i], pderiv(D, i)))
    return sub(mul(D, divN), ND)


# ----------------------------------------------------------------------------------------------------------
# E1-admissible space in degree d, sector c (exact over Q via fmpz nullspace, dims mod p)
# ----------------------------------------------------------------------------------------------------------
def _flat_matrix(cols, rows):
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    return flat, R, C


def rank_mod_p(polys, p=PRIMES[0]):
    """Rank (mod p) of a list of integer/rational polys (as vectors)."""
    polys = [sect.integral(P) if any(isinstance(v, Fraction) for v in P.values()) else P for P in polys]
    rows = {}
    for P in polys:
        for e in P:
            rows.setdefault(e, len(rows))
    if not rows or not polys:
        return 0
    flat, R, C = _flat_matrix(polys, rows)
    return nmod_mat(R, C, [x % p for x in flat], p).rank()


def pivot_subset_mod_p(polys, p=PRIMES[0]):
    """Indices of a maximal independent subset of the given polys (greedy by rref pivots mod p)."""
    polys_i = [sect.integral(P) if any(isinstance(v, Fraction) for v in P.values()) else P for P in polys]
    rows = {}
    for P in polys_i:
        for e in P:
            rows.setdefault(e, len(rows))
    if not polys_i:
        return []
    flat, R, C = _flat_matrix(polys_i, rows)   # columns = polys
    M = nmod_mat(R, C, [x % p for x in flat], p)
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j < C and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    return piv


def e1_space(S, F, A, q, d, c, p=PRIMES[0], verbose=False, exact=True):
    """E1-admissible denominators {D in S_d^c : A D in J(F)} with exact lifts.

    Returns dict: Dm, basis (list of (D, N) exact), dims {dim_S, dim_J, kerA_in_S, ker_mod_J, dim_R_d, dim_R_t, forced}.
    The basis is: a maximal independent subset of J_d^c generators m F_i (lift A m e_i / q), plus nullspace
    vectors of [A*m | m' F_i] whose D-parts are independent modulo J (lift from the nullspace vector)."""
    t0 = time.time()
    c = tuple(c)
    Fz = sect.integral(F)                     # 7 F for the lane's F
    LF = Fraction(Fz[next(iter(Fz))]) / Fraction(F[next(iter(Fz))])   # Fz = LF * F
    Az = sect.integral(A)
    assert Az == {e: int(v) for e, v in A.items()}, 'A must be integral monic-ish for this routine'
    dFz = [sect.ideriv(Fz, i) for i in range(6)]
    Dm = S.mons(d, c)
    nD = len(Dm)
    out = {'d': d, 'sector': list(c), 'dim_S': nD}
    if nD == 0:
        out.update({'dim_J': 0, 'kerA_in_S': 0, 'ker_mod_J': 0, 'basis': [], 'Dm': Dm, 'seconds': 0.0})
        return out
    # J_d^c generators
    Jgens = []
    if d >= 5:
        for i in range(6):
            for m in S.mons(d - 5, S.chiJgen(c, i)):
                Jgens.append((i, m, sect.imul({m: 1}, dFz[i])))
    # E1 target: degree d+6, chi = chiA + c ; columns A*m and m' Fz_i (deg m' = d+1, chi(m') = chiN(c,i))
    tchi = sect.cadd(S.chiA, c)
    rows = {}
    for e in S.mons(d + 6, tchi):
        rows[e] = len(rows)
    colsA = [sect.imul(Az, {m: 1}) for m in Dm]
    colsJ, metaJ = [], []
    for i in range(6):
        for mp in S.mons(d + 1, S.chiN(c, i)):
            colsJ.append(sect.imul({mp: 1}, dFz[i])); metaJ.append((i, mp))
    allcols = colsA + colsJ
    for col in allcols:
        for e in col:
            assert e in rows, 'row set mismatch'
    flat, R, C = _flat_matrix(allcols, rows)
    # ranks mod p (fast): J-block rank -> dim R_{d+6}^{tchi}; full rank -> kerA_in_S = C - rank
    Mp = nmod_mat(R, C, [x % p for x in flat], p)
    rank_full = Mp.rank()
    flatJ, RJ, CJ = _flat_matrix(colsJ, rows)
    rank_J6 = nmod_mat(RJ, CJ, [x % p for x in flatJ], p).rank() if CJ else 0
    nullity_p = C - rank_full
    kerA_in_S_p = nullity_p - (CJ - rank_J6)          # dim of projection of the kernel to the D coordinates
    dim_J = rank_mod_p([g[2] for g in Jgens], p) if Jgens else 0
    out.update({'dim_J': dim_J, 'dim_R_d': nD - dim_J, 'dim_R_t': R - rank_J6,
                'kerA_in_S': kerA_in_S_p, 'ker_mod_J': kerA_in_S_p - dim_J,
                'forced': max(0, (nD - dim_J) - (R - rank_J6)), 'rows': R, 'cols': C})
    if verbose:
        print(f'[e1_space] d={d} c={c} S={nD} J={dim_J} kerA={kerA_in_S_p} ker_mod_J={out["ker_mod_J"]} '
              f'R_d={out["dim_R_d"]} R_t={out["dim_R_t"]} forced={out["forced"]} ({R}x{C})', file=sys.stderr, flush=True)
    # basis: independent J generators
    basis = []
    if Jgens:
        piv = pivot_subset_mod_p([g[2] for g in Jgens], p)
        for j in piv:
            i, m, poly = Jgens[j]
            D = {e: Fraction(v) for e, v in sect.imul({m: 1}, [pderiv(F, k) for k in range(6)][i]).items()}
            N = [{} for _ in range(6)]
            N[i] = {e: Fraction(v, q) for e, v in mul(A, {m: Fraction(1)}).items()}
            basis.append((D, N))
    if out['ker_mod_J'] > 0 and exact:
        # exact nullspace of the full E1 matrix, pick vectors independent modulo J
        t1 = time.time()
        M = fmpz_mat(R, C, flat)
        X, nullity = M.nullspace()
        assert nullity == nullity_p, (nullity, nullity_p)
        # D-parts mod p -> choose columns whose D-parts are independent modulo J_d^c
        Jpolys = [basis_k[0] for basis_k in basis]
        cur = list(Jpolys)
        cur_rank = rank_mod_p(cur, p) if cur else 0
        chosen = []
        for k in range(nullity):
            a = [int(X[j, k]) for j in range(nD)]
            if not any(a):
                continue
            Dk = {Dm[j]: Fraction(a[j]) for j in range(nD) if a[j]}
            r = rank_mod_p(cur + [Dk], p)
            if r > cur_rank:
                cur.append(Dk); cur_rank = r
                b = [int(X[nD + j, k]) for j in range(CJ)]
                N = [{} for _ in range(6)]
                for (i, mp), v in zip(metaJ, b):
                    if v:
                        N[i][mp] = N[i].get(mp, 0) - Fraction(v) * LF / q
                chosen.append((Dk, N))
                if len(chosen) == out['ker_mod_J']:
                    break
        assert len(chosen) == out['ker_mod_J'], (len(chosen), out['ker_mod_J'])
        for Dk, N in chosen:
            assert not e1_check(F, A, q, Dk, N), 'exact E1 lift failed'
        basis.extend(chosen)
        out['t_nullspace'] = round(time.time() - t1, 1)
    out['basis'] = basis
    out['Dm'] = Dm
    out['seconds'] = round(time.time() - t0, 1)
    return out


def save_space(res, path):
    ser = {k: v for k, v in res.items() if k not in ('basis', 'Dm')}
    ser['basis'] = [{'D': core.to_str(D), 'N': [core.to_str(n) for n in N]} for D, N in res['basis']]
    with open(path, 'w') as f:
        json.dump(ser, f)


def load_space(path):
    with open(path) as f:
        ser = json.load(f)
    ser['basis'] = [(core.parse(b['D']), [core.parse(s) for s in b['N']]) for b in ser['basis']]
    return ser


# ----------------------------------------------------------------------------------------------------------
# numerical family (Koszul-reduced variable projection, sector restricted)
# ----------------------------------------------------------------------------------------------------------
class _Prod:
    """Bilinear product of coefficient vectors over monomial lists U, V into the row index `rows`."""

    def __init__(self, U, V, rows, nrows):
        I, J, Rk = [], [], []
        for i, u in enumerate(U):
            for j, v in enumerate(V):
                e = (u[0]+v[0], u[1]+v[1], u[2]+v[2], u[3]+v[3], u[4]+v[4], u[5]+v[5])
                r = rows.get(e)
                if r is not None:
                    I.append(i); J.append(j); Rk.append(r)
        self.I, self.J, self.R = np.array(I, dtype=np.int64), np.array(J, dtype=np.int64), np.array(Rk, dtype=np.int64)
        self.nU, self.nV, self.nR = len(U), len(V), nrows

    def __call__(self, u, v):
        if not len(self.I):
            return np.zeros(self.nR)
        return np.bincount(self.R, weights=u[self.I] * v[self.J], minlength=self.nR)

    def adj_u(self, r, v):
        if not len(self.I):
            return np.zeros(self.nU)
        return np.bincount(self.I, weights=r[self.R] * v[self.J], minlength=self.nU)

    def adj_v(self, r, u):
        if not len(self.I):
            return np.zeros(self.nV)
        return np.bincount(self.J, weights=r[self.R] * u[self.I], minlength=self.nV)


def _deriv_matrix(src, dst, i):
    """Sparse matrix of d/dz_i from monomial list src to monomial list dst."""
    idx = {m: k for k, m in enumerate(dst)}
    ri, ci, vals = [], [], []
    for j, m in enumerate(src):
        if m[i]:
            e = list(m); e[i] -= 1; e = tuple(e)
            k = idx.get(e)
            if k is not None:
                ri.append(k); ci.append(j); vals.append(float(m[i]))
    return sp.csr_matrix((vals, (ri, ci)), shape=(len(dst), len(src)))


class Family:
    """Numerical Koszul-reduced variable-projection family for D(a) = sum a_k D_k in sector c, degree d."""

    def __init__(self, S, F, A, q, d, c, basis, exclude_F=True, verbose=False, cond=1e-10, orthonormal=True):
        t0 = time.time()
        self.S, self.F, self.A, self.q, self.d, self.c = S, F, A, q, d, tuple(c)
        self.cond = cond
        c = self.c
        K = len(basis)
        self.K = K
        self.basis = basis
        Ff = to_float(F)
        # monomial sets
        self.Dm = S.mons(d, c); self.nD = len(self.Dm); Didx = {m: k for k, m in enumerate(self.Dm)}
        self.Nm = [S.mons(d + 1, S.chiN(c, i)) for i in range(6)]
        self.Noff = np.cumsum([0] + [len(x) for x in self.Nm]); self.nN = int(self.Noff[-1])
        Nidx = [{m: k for k, m in enumerate(self.Nm[i])} for i in range(6)]
        self.divm = S.mons(d, sect.cadd(S.chiA, c)); ndiv = len(self.divm)
        self.gradm = [S.mons(d - 1, sect.cadd(c, tuple((-x) % 6 for x in S.chiE[i]))) for i in range(6)]
        self.rowsm = S.mons(2 * d, S.chiQ(c)); self.nrows = len(self.rowsm)
        rows = {m: k for k, m in enumerate(self.rowsm)}; self.rows = rows
        degc = d + 6 * q - 10
        self.cm = []   # list of (i, j, monomial)
        for i in range(6):
            for j in range(i + 1, 6):
                for m in (S.mons(degc, S.chiC(c, i, j)) if degc >= 0 else []):
                    self.cm.append((i, j, m))
        self.nC = len(self.cm)
        self.Qm = S.mons(2 * d - 6, S.chiQ(c)) if 2 * d - 6 >= 0 else []
        self.nQ = len(self.Qm)
        # D and N0 matrices (float), normalized per basis element
        Dmat = np.zeros((self.nD, K)); N0mat = np.zeros((self.nN, K))
        self.scales = []
        for k, (D, N) in enumerate(basis):
            nrm = poly_norm(D)
            self.scales.append(nrm)
            for e, v in D.items():
                Dmat[Didx[e], k] = float(v) / nrm
            for i in range(6):
                for e, v in N[i].items():
                    N0mat[self.Noff[i] + Nidx[i][e], k] = float(v) / nrm
        self.Dmat_raw = Dmat
        if orthonormal and K > 1:
            Qo, Ro = np.linalg.qr(Dmat)
            Dmat = Qo
            N0mat = np.linalg.solve(Ro.T, N0mat.T).T   # N0 R^{-1}
            self.Rfac = Ro
        else:
            self.Rfac = np.eye(K)
        self.Dmat, self.N0mat = Dmat, N0mat
        # derivative operators
        Fi = [pderiv(F, i) for i in range(6)]
        self.Div = sp.hstack([_deriv_matrix(self.Nm[i], self.divm, i) for i in range(6)]).tocsr()
        self.Grad = [_deriv_matrix(self.Dm, self.gradm[i], i) for i in range(6)]
        self.P1 = _Prod(self.Dm, self.divm, rows, self.nrows)
        self.P2 = [_Prod(self.Nm[i], self.gradm[i], rows, self.nrows) for i in range(6)]
        # Koszul map K(c): N-vector per c-column, and U: vec(Bil(e_m, K e_c)) over (row, c) by D-monomial m
        Fif = [to_float(P) for P in Fi]
        Kcols = []   # list of per-column N polys (list of 6 float dicts)
        for (i, j, m) in self.cm:
            N = [{} for _ in range(6)]
            N[i] = {tuple(a + b for a, b in zip(m, e)): v for e, v in Fif[j].items()}
            N[j] = {tuple(a + b for a, b in zip(m, e)): -v for e, v in Fif[i].items()}
            Kcols.append(N)
        self._build_U(Kcols)
        # Q block: -F^q * m
        Fq = core.power(F, q)
        Fqf = to_float(Fq)
        ri, ci, vals = [], [], []
        for jq, m in enumerate(self.Qm):
            for e, v in Fqf.items():
                r = rows[tuple(a + b for a, b in zip(m, e))]
                ri.append(r); ci.append(jq); vals.append(-v)
        self.MQ = sp.csr_matrix((vals, (ri, ci)), shape=(self.nrows, self.nQ)).toarray() if self.nQ else np.zeros((self.nrows, 0))
        self.ncols = self.nC + self.nQ
        # F-multiples exclusion
        self.P = np.eye(K); self.excluded_dim = 0
        if exclude_F and d >= 6:
            self._exclude_F_multiples()
        self.dimP = self.P.shape[1]
        # relaxation projector (lazy)
        self._Pperp = None
        self.build_time = round(time.time() - t0, 1)
        self._nfev = 0
        if verbose:
            print(f'[family] d={d} c={c} K={K} (dimP {self.dimP}) rows={self.nrows} cols={self.ncols} '
                  f'(c {self.nC}, Q {self.nQ}) nN={self.nN} nnz(U)={self.U.nnz} built {self.build_time}s',
                  file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    def _build_U(self, Kcols):
        """U: ((nrows * nC) x nD) sparse with U[(r, c), m] = coefficient of row r in Bil(e_m, K e_c)."""
        rows, nC = self.rows, self.nC
        ri, ci, vals = [], [], []
        # per column: div(K e_c) poly and the 6 component polys
        for cidx, N in enumerate(Kcols):
            divK = {}
            for i in range(6):
                for e, v in N[i].items():
                    if e[i]:
                        e2 = list(e); e2[i] -= 1; e2 = tuple(e2)
                        divK[e2] = divK.get(e2, 0.0) + v * e[i]
            divK = {e: v for e, v in divK.items() if v}
            for mi, m in enumerate(self.Dm):
                # term 1: e_m * divK
                for e, v in divK.items():
                    r = rows.get((m[0]+e[0], m[1]+e[1], m[2]+e[2], m[3]+e[3], m[4]+e[4], m[5]+e[5]))
                    if r is not None:
                        ri.append(r * nC + cidx); ci.append(mi); vals.append(v)
                # term 2: - sum_l K_l * d_l e_m = - sum_l m_l e_{m - e_l} * K_l
                for l in range(6):
                    if m[l] and N[l]:
                        ml = m[l]
                        for e, v in N[l].items():
                            ee = (m[0]+e[0], m[1]+e[1], m[2]+e[2], m[3]+e[3], m[4]+e[4], m[5]+e[5])
                            ee = ee[:l] + (ee[l] - 1,) + ee[l+1:]
                            r = rows.get(ee)
                            if r is not None:
                                ri.append(r * nC + cidx); ci.append(mi); vals.append(-v * ml)
        self.U = sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                               shape=(self.nrows * nC, self.nD))
        self.U.sum_duplicates()

    def _exclude_F_multiples(self):
        d = self.d
        Gm = self.S.mons(d - 6, self.c)
        if not Gm:
            return
        Ff = to_float(self.F)
        Didx = {m: k for k, m in enumerate(self.Dm)}
        Fm = np.zeros((self.nD, len(Gm)))
        for j, m in enumerate(Gm):
            for e, v in Ff.items():
                Fm[Didx[tuple(a + b for a, b in zip(m, e))], j] = -v
        Mb = np.hstack([self.Dmat, Fm])
        ns = sla.null_space(Mb)
        Sa = ns[:self.K, :]
        if Sa.size == 0:
            return
        Uo, s, _ = np.linalg.svd(Sa, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        self.P = Uo[:, r:]
        self.excluded_dim = r

    # ------------------------------------------------------------------
    def bilinear(self, Dv, Nv):
        """Bil(D, N) for coefficient vectors."""
        out = self.P1(Dv, self.Div @ Nv)
        for i in range(6):
            gi = self.Grad[i] @ Dv
            if gi.size and self.Nm[i]:
                out -= self.P2[i](Nv[self.Noff[i]:self.Noff[i+1]], gi)
        return out

    def bil_adj_D(self, r, Nv):
        """w with w . D = r . Bil(D, N) for all D."""
        divN = self.Div @ Nv
        w = self.P1.adj_u(r, divN)
        for i in range(6):
            if self.gradm[i] and self.Nm[i]:
                w -= self.Grad[i].T @ self.P2[i].adj_v(r, Nv[self.Noff[i]:self.Noff[i+1]])
        return w

    def bil_adj_N(self, r, Dv):
        """y with y . N = r . Bil(D, N) for all N."""
        y = self.Div.T @ self.P1.adj_v(r, Dv)
        for i in range(6):
            gi = self.Grad[i] @ Dv
            if gi.size and self.Nm[i]:
                y[self.Noff[i]:self.Noff[i+1]] -= self.P2[i].adj_u(r, gi)
        return y

    def assemble(self, a):
        Dv = self.Dmat @ a
        N0v = self.N0mat @ a
        T = (self.U @ Dv).reshape(self.nrows, self.nC) if self.nC else np.zeros((self.nrows, 0))
        b = -self.bilinear(Dv, N0v)
        M = np.hstack([T, self.MQ]) if self.nQ else T
        return M, b, Dv, N0v

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        M, b, Dv, N0v = self.assemble(a)
        if M.shape[1]:
            x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            r = M @ x - b
        else:
            x = np.zeros(0); r = -b
        self._nfev += 1
        return float(r @ r), x, r, b

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b, Dv, N0v = self.assemble(a)
        if M.shape[1]:
            x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            r = M @ x - b
        else:
            x = np.zeros(0); r = -b
        f = float(r @ r)
        cvec = x[:self.nC]
        # d/da_k of r.(M x): r . T_k c = (U^T vec(r c^T)) . Dmat[:, k]
        if self.nC:
            z = self.U.T @ np.outer(r, cvec).ravel()
            gT = self.Dmat.T @ z
        else:
            gT = np.zeros(self.K)
        # d/da_k of r.(-b) = r . Bil(D_k, N0(a)) + r . Bil(D(a), N0_k)
        wD = self.bil_adj_D(r, N0v)
        yN = self.bil_adj_N(r, Dv)
        gb = self.Dmat.T @ wD + self.N0mat.T @ yN
        g = 2.0 * (gT + gb)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        return float(np.linalg.norm(self.bilinear(self.Dmat @ a, self.N0mat @ a)))

    # ---------------- mod-J relaxation -------------------------------
    def Pperp(self):
        """Orthonormal basis of the orthogonal complement of J(F)_{2d}^{chiQ} inside the E2 row space."""
        if self._Pperp is None:
            S, d, c = self.S, self.d, self.c
            Ff = [to_float(pderiv(self.F, i)) for i in range(6)]
            cols = []
            for i in range(6):
                for m in S.mons(2 * d - 5, sect.cadd(S.chiQ(c), S.chiE[i])):
                    cols.append((i, m))
            Jm = np.zeros((self.nrows, len(cols)))
            for j, (i, m) in enumerate(cols):
                for e, v in Ff[i].items():
                    Jm[self.rows[tuple(a + b for a, b in zip(m, e))], j] = v
            if Jm.shape[1] == 0:
                self._Pperp = np.eye(self.nrows); self.rank_J2d = 0
            else:
                Uo, s, _ = np.linalg.svd(Jm, full_matrices=True)
                rk = int(np.sum(s > 1e-9 * s[0]))
                self.rank_J2d = rk
                self._Pperp = Uo[:, rk:]
        return self._Pperp

    def relax_value_and_grad(self, u):
        Pp = self.Pperp()
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        Dv = self.Dmat @ a; N0v = self.N0mat @ a
        b = self.bilinear(Dv, N0v)
        pb = Pp.T @ b
        f = float(pb @ pb)
        rho = Pp @ pb
        wD = self.bil_adj_D(rho, N0v)
        yN = self.bil_adj_N(rho, Dv)
        g = 2.0 * (self.Dmat.T @ wD + self.N0mat.T @ yN)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def relax_value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        b = self.bilinear(self.Dmat @ a, self.N0mat @ a)
        pb = self.Pperp().T @ b
        return float(pb @ pb), float(b @ b)

    # ---------------- searches ----------------------------------------
    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, relax=False, verbose=False, stop_rel=1e-13,
               init=None, time_limit=None):
        rng = np.random.default_rng(seed)
        fun = self.relax_value_and_grad if relax else self.value_and_grad
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(fun, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            if relax:
                f, bn2 = self.relax_value(a)
                rel = f / max(bn2, 1e-300)
            else:
                f, x, r, b = self.value(a)
                bn2 = float(b @ b)
                rel = f / max(bn2, 1e-300)
            spur, srel = self.is_spurious(a)
            out.append({'value': f, 'rel': rel, 'bnorm2': bn2, 'a': a, 'iters': int(res.nit), 'restart': it,
                        'spurious': spur, 'div_rel': srel, 'seed': seed, 'relax': relax})
            if verbose:
                print(f'[search{"-relax" if relax else ""}] seed {seed} restart {it}: f={f:.3e} rel={rel:.3e} '
                      f'iters={res.nit} spur={spur} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
            if rel < stop_rel and not spur:
                break
            if time_limit and time.time() - t0 > time_limit:
                break
        out.sort(key=lambda z: (z['spurious'], z['rel']))
        return out

    def is_spurious(self, a, tol=1e-8):
        """Is D(a) numerically divisible by F?"""
        d = self.d
        if d < 6:
            return False, 1.0
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Dv = self.Dmat @ a
        Gm = self.S.mons(d - 6, self.c)
        if not Gm:
            return False, 1.0
        Ff = to_float(self.F)
        Didx = {m: k for k, m in enumerate(self.Dm)}
        Fm = np.zeros((self.nD, len(Gm)))
        for j, m in enumerate(Gm):
            for e, v in Ff.items():
                Fm[Didx[tuple(x + y for x, y in zip(m, e))], j] = v
        x, *_ = sla.lstsq(Fm, Dv, check_finite=False)
        rel = float(np.linalg.norm(Fm @ x - Dv) / max(np.linalg.norm(Dv), 1e-300))
        return rel < tol, rel

    # ---------------- reconstruction ----------------------------------
    def D_float(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Dv = self.Dmat @ a
        return {m: Dv[k] for k, m in enumerate(self.Dm) if abs(Dv[k]) > 1e-14 * np.max(np.abs(Dv))}

    def reconstruct(self, a, dens=(1, 6, 12, 60, 420, 2520, 10**4, 10**5, 10**6, 10**7), tol=None):
        """Rational candidates for D(a): coefficient-wise limit_denominator after normalizing the largest to 1,
        and (second route) coordinates w.r.t. the raw rational basis. Returns list of (D exact, descr, err)."""
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Dv = self.Dmat @ a
        j = int(np.argmax(np.abs(Dv)))
        Dn = Dv / Dv[j]
        cands = []
        seen = set()
        for den in dens:
            fr = {self.Dm[k]: Fraction(float(Dn[k])).limit_denominator(den) for k in range(self.nD)}
            fr = {m: v for m, v in fr.items() if v}
            if not fr:
                continue
            err = math.sqrt(sum((float(fr.get(self.Dm[k], 0)) - Dn[k]) ** 2 for k in range(self.nD))) / np.linalg.norm(Dn)
            key = core.to_str(fr)
            if key in seen:
                continue
            seen.add(key)
            cands.append((fr, f'coef-limit_den({den})', err))
        # coordinates w.r.t. raw basis: a_raw = R^{-1} a (since Dmat = Draw R^{-1} ... Draw = Q R => Q a = Draw R^{-1} a)
        try:
            araw = np.linalg.solve(self.Rfac, a)       # coefficients on the normalized raw columns
            araw = araw / np.array(self.scales)        # on the exact basis polys
            j2 = int(np.argmax(np.abs(araw)))
            an = araw / araw[j2]
            for den in dens:
                fr = [Fraction(float(v)).limit_denominator(den) for v in an]
                D = {}
                for k, v in enumerate(fr):
                    if v:
                        D = add(D, scale(self.basis[k][0], v))
                if not D:
                    continue
                key = core.to_str(D)
                if key in seen:
                    continue
                seen.add(key)
                err = max(abs(float(f) - v) for f, v in zip(fr, an))
                cands.append((D, f'basis-limit_den({den})', err))
        except Exception:
            pass
        return cands


# ----------------------------------------------------------------------------------------------------------
# exact confirmation
# ----------------------------------------------------------------------------------------------------------
def confirm_exact(F, A, q, D, exact=None, verbose=False):
    """E1 (mod p) then the wave-1 sector solver (exact ranks over Q when the system is small, else 2 primes)."""
    t0 = time.time()
    d = core.degree(D)
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    out = {'d': d, 'e1': bool(ok), 'e1_info': {k: info[k] for k in ('rank', 'rank_aug', 'rows', 'cols')}}
    if not ok:
        out['consistent'] = False; out['seconds'] = round(time.time() - t0, 1)
        return out
    if exact is None:
        exact = d <= 9
    res = sect.test(F, A, q, D, exact=exact, verbose=verbose, lift=True)
    out.update({k: v for k, v in res.items() if k not in ('N', 'Q')})
    if res.get('verified'):
        out['witness'] = {'d': d, 'D': core.to_str(D), 'N': res['N'], 'Q': res['Q']}
    out['seconds'] = round(time.time() - t0, 1)
    return out


def promote(fam, cand, F, A, q, tag, job, extra=None, max_cands=6):
    """Rational reconstruction + exact confirmation of a numerical candidate; logs everything."""
    results = []
    for D, descr, err in fam.reconstruct(cand['a'])[:max_cands]:
        if core.degree(D) != fam.d:
            continue
        try:
            divF, _ = __import__('hodge.jacobian', fromlist=['jacobian']).divisible_by(F, D)
        except Exception:
            divF = None
        rec = {'batch': 'promote', 'job': job, 'tag': tag, 'd': fam.d, 'sector': list(fam.c), 'descr': descr,
               'recon_err': err, 'num_rel': cand['rel'], 'F_divides_D': divF, 'D': core.to_str(D)[:600], 'terms': len(D)}
        if divF:
            rec['note'] = 'reconstructed D divisible by F (spurious)'
            w2log(rec); results.append(rec); continue
        conf = confirm_exact(F, A, q, D)
        rec.update({k: v for k, v in conf.items() if k != 'witness'})
        if conf.get('witness'):
            rec['WITNESS'] = True
            path = f'{WD}/witness-candidate.json'
            with open(path, 'w') as f:
                json.dump({'instanceId': IID, 'lane': IID, 'job': job, 'tag': tag, **conf['witness']}, f, indent=1)
            rec['witness_path'] = path
            print('\n*** WITNESS FOUND ***', core.to_str(D)[:200], '\n', flush=True)
        if extra:
            rec.update(extra)
        w2log(rec)
        results.append(rec)
        if conf.get('witness'):
            break
    return results


# ----------------------------------------------------------------------------------------------------------
# self-tests
# ----------------------------------------------------------------------------------------------------------
def check_direct(fam, a=None, seed=0):
    """Compare the vectorized residual with a direct float polynomial evaluation of E2 for random (a, c, Q)."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(fam.K) if a is None else np.asarray(a, float)
    a /= np.linalg.norm(a)
    x = rng.standard_normal(fam.ncols)
    M, b, Dv, N0v = fam.assemble(a)
    r_vec = M @ x - b
    # direct: D(a), N = N0(a) + K(c), Q
    Dp = {m: Dv[k] for k, m in enumerate(fam.Dm)}
    N = [{} for _ in range(6)]
    for i in range(6):
        for k, m in enumerate(fam.Nm[i]):
            N[i][m] = N0v[fam.Noff[i] + k]
    Fi = [to_float(pderiv(fam.F, i)) for i in range(6)]
    for cidx, (i, j, m) in enumerate(fam.cm):
        cval = x[cidx]
        for e, v in Fi[j].items():
            ee = tuple(p + q_ for p, q_ in zip(m, e)); N[i][ee] = N[i].get(ee, 0.0) + cval * v
        for e, v in Fi[i].items():
            ee = tuple(p + q_ for p, q_ in zip(m, e)); N[j][ee] = N[j].get(ee, 0.0) - cval * v
    Qp = {m: x[fam.nC + k] for k, m in enumerate(fam.Qm)}
    divN, ND = {}, {}
    for i in range(6):
        divN = padd(divN, {tuple(e[:i] + (e[i]-1,) + e[i+1:]): v * e[i] for e, v in N[i].items() if e[i]})
        ND = padd(ND, mul(N[i], {tuple(e[:i] + (e[i]-1,) + e[i+1:]): v * e[i] for e, v in Dp.items() if e[i]}))
    E2 = padd(mul(Dp, divN), {e: -v for e, v in ND.items()}, {e: -v for e, v in mul(to_float(core.power(fam.F, fam.q)), Qp).items()})
    r_dir = np.zeros(fam.nrows)
    extra = 0.0
    for e, v in E2.items():
        if e in fam.rows:
            r_dir[fam.rows[e]] += v
        else:
            extra += v * v
    err = np.linalg.norm(r_dir - r_vec) / max(np.linalg.norm(r_vec), 1e-300)
    return err, math.sqrt(extra)


def check_grad(fam, seed=0, relax=False, h=1e-6):
    rng = np.random.default_rng(seed)
    u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
    fun = fam.relax_value_and_grad if relax else fam.value_and_grad
    f0, g = fun(u)
    errs = []
    for t in range(min(5, fam.dimP)):
        e = np.zeros(fam.dimP); e[rng.integers(fam.dimP)] = 1.0
        fp, _ = fun(u + h * e); fm, _ = fun(u - h * e)
        fd = (fp - fm) / (2 * h)
        errs.append(abs(fd - g @ e) / max(abs(fd), abs(g @ e), 1e-12))
    return max(errs), f0
