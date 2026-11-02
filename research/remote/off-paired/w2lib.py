"""Wave-2 search library for lane off-paired (q = 1, A = X0^3 Y0^3).

Sector-restricted, *reduced* variable projection: the denominator ranges over a rational
basis of an E1-admissible space (J:A)_{d,chi} (Jacobian part J_{d,chi} / non-Jacobian kernel /
both), the Euler gauge F*R_{d-6} removed.  For D = sum_k a_k D_k we write
    N = sum_k a_k N_k + sum_j s_j S_j,       q*sum_i N_k,i F_i = A D_k,   S_j = Koszul syzygies,
so that E1 holds identically and the certificate reduces to
    B_D(N) - F Q = 0,      B_D(N) = D div N - N . grad D      (degree 2d, character chi_A + 2 chi_D).
For fixed a this is linear in (s, Q): f(a) = min_{s,Q} || M(a)(s;Q) - b(a) ||^2 with M(a) linear and
b(a) quadratic in a -- exactly the structure of hodge.varpro.VarProFamily, whose optimizer we reuse.
"Hessian relaxation" = the same family with no Koszul unknowns (N = sum a_k N_k).

Everything numeric is float64; only the exact re-solve (symmetry.solve_fixed_D_symmetric on the
rationally reconstructed D) counts as evidence.  Engine v2.
"""
import sys, time, json, math, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import flint
from hodge import core, fixed_d, jacobian, symmetry, varpro
from hodge.core import NV, monomials, mul, add, sub, scale, shift, pderiv, power

LANE = 'off-paired'
LANE_DIR = f'/tmp/hodge/results/{LANE}'
LOG = f'{LANE_DIR}/wave2-search.log.jsonl'
P0 = 1000003
P1 = 2147483647

# diagonal group Z/6 x Z/6 of off-paired (dossier 1.3); verified in check_group()
GW = [(0, 2, 0, 1, 3, 3), (0, 0, 1, 0, 0, 5)]
NG = 6
WZ = [(GW[0][i], GW[1][i]) for i in range(NV)]   # character of the coordinate z_i


def chi(m):
    return (sum(a * b for a, b in zip(m, GW[0])) % NG, sum(a * b for a, b in zip(m, GW[1])) % NG)


def cadd(*cs):
    a = sum(c[0] for c in cs) % NG
    b = sum(c[1] for c in cs) % NG
    return (a, b)


def cneg(c):
    return ((-c[0]) % NG, (-c[1]) % NG)


def poly_chi(Pl):
    cs = {chi(m) for m in Pl}
    return cs.pop() if len(cs) == 1 else None


def mons_chi(d, c):
    if d < 0:
        return []
    return [m for m in monomials(d) if chi(m) == c]


def all_chars_parity(par):
    """Characters (a,b) whose g-weight parity equals the Y-parity par."""
    return [(a, b) for a in range(NG) for b in range(NG) if a % 2 == par]


def load():
    F, A, q, _rec = core.load_instance(LANE)
    return F, A, q


def check_group(F):
    G = symmetry.diagonal_group(F)
    assert G.order == 6 and len(G.elements) == 36, (G.order, len(G.elements))
    els = set(G.elements)
    for w in GW:
        # the engine normalizes w_5 = 0; compare projectively (shift by a scalar weight)
        shifts = [tuple((w[i] - s) % 6 for i in range(6)) for s in range(6)]
        assert any(sh in els for sh in shifts), w
    return G


# ----------------------------------------------------------------------------------------------
def fmpz_rows_from_polys(polys, idx, n, den):
    """Integer matrix (len(polys) x n) of den*coefficients (den clears the denominators)."""
    rows = []
    for Pl in polys:
        r = [0] * n
        for e, c in Pl.items():
            v = Fraction(c) * den
            assert v.denominator == 1, (c, den)
            r[idx[e]] = int(v)
        rows.append(r)
    return rows


def lcm_den(polys):
    L = 1
    for Pl in polys:
        for c in Pl.values():
            dn = Fraction(c).denominator
            L = L * dn // math.gcd(L, dn)
    return L


# ----------------------------------------------------------------------------------------------
# E1-admissible spaces in a character sector (exact over Q)
class SectorSpaces:
    """Rational bases of, in degree d and D-character c:
         gauge  : F * R_{d-6, c}
         jac    : J(F)_{d,c} = span{ m F_i }           (complement of gauge, basis = selected m F_i)
         ker    : complement of (gauge + jac) inside (J:A)_{d,c}   (non-Jacobian E1 kernel)
       and dimensions.  Basis polys have Fraction coefficients.
    """

    def __init__(self, F, A, q, d, c, p=P0, verbose=False):
        t0 = time.time()
        self.F, self.A, self.q, self.d, self.c, self.p = F, A, q, d, c, p
        self.cA = poly_chi(A)
        Dm = mons_chi(d, c)
        self.Dmons = Dm
        self.dimR = len(Dm)
        K = len(Dm)
        Fi = [pderiv(F, i) for i in range(NV)]
        self.gauge = [shift(F, m, 1) for m in mons_chi(d - 6, c)] if d >= 6 else []
        jac_gens = []
        for i in range(NV):
            for m in mons_chi(d - 5, cadd(c, WZ[i])):
                jac_gens.append(((i, m), shift(Fi[i], m, 1)))
        self.jac_gens = jac_gens
        self.E1 = []
        if K:
            cAD = cadd(self.cA, c)
            Ncols = []
            for i in range(NV):
                for m in mons_chi(d + 6 - 5, cadd(cAD, WZ[i])):
                    Ncols.append(shift(Fi[i], m, 1))
            Acols = [mul(A, {m: Fraction(1)}) for m in Dm]
            ridx = {}
            for Pl in Acols + Ncols:
                for e in Pl:
                    ridx.setdefault(e, len(ridx))
            R = len(ridx)
            den = lcm_den(Acols + Ncols)
            colsI = fmpz_rows_from_polys(Acols + [scale(Pl, -1) for Pl in Ncols], ridx, R, den)
            nC = len(colsI)
            M = flint.fmpz_mat(R, nC, [colsI[j][r] for r in range(R) for j in range(nC)])
            X, nul = M.nullspace()
            E1 = []
            for j in range(nul):
                colv = [int(X[i, j]) for i in range(K)]
                if any(colv):
                    g = 0
                    for v in colv:
                        g = math.gcd(g, abs(v))
                    E1.append([v // g for v in colv])
            self.e1_raw = E1
            self.E1 = self._independent(E1, K)
        self.dimE1 = len(self.E1)
        gauge_vecs = [self._vec(Pl, Dm) for Pl in self.gauge]
        self.dim_gauge = self._rank_rows([self._modrow(v) for v in gauge_vecs], K)
        jac_vecs = [self._vec(Pl, Dm) for _, Pl in jac_gens]
        sel_jac = self._extend(gauge_vecs, jac_vecs, K)
        self.jac_basis = [jac_gens[j][1] for j in sel_jac]
        self.jac_tags = [jac_gens[j][0] for j in sel_jac]
        self.dim_jac_total = self._rank_rows([self._modrow(v) for v in gauge_vecs + jac_vecs], K)
        base = gauge_vecs + [jac_vecs[j] for j in sel_jac]
        e1_vecs = [[Fraction(v) for v in row] for row in self.E1]
        sel_ker = self._extend(base, e1_vecs, K)
        self.ker_basis = [self._poly(self.E1[j], Dm) for j in sel_ker]
        self.dim_ker = len(self.ker_basis)
        self.dim_jac = len(self.jac_basis)
        self.time = time.time() - t0
        if verbose:
            print(f'[spaces] d={d} chi={c}: R={self.dimR} E1={self.dimE1} gauge={self.dim_gauge} '
                  f'jac(mod gauge)={self.dim_jac} ker={self.dim_ker} ({self.time:.1f}s)', file=sys.stderr)

    def _vec(self, Pl, Dm):
        idx = {m: i for i, m in enumerate(Dm)}
        v = [Fraction(0)] * len(Dm)
        for e, cf in Pl.items():
            v[idx[e]] = Fraction(cf)
        return v

    def _poly(self, row, Dm):
        return {m: Fraction(int(v)) for m, v in zip(Dm, row) if v}

    def _modrow(self, v):
        p = self.p
        out = []
        for x in v:
            fx = Fraction(x)
            out.append(int(fx.numerator % p * pow(fx.denominator % p, -1, p) % p))
        return out

    def _rank_rows(self, rows, K):
        if not rows or K == 0:
            return 0
        M = flint.nmod_mat(len(rows), K, [x for r in rows for x in r], self.p)
        return M.rank()

    def _extend(self, base, cands, K):
        """Greedy: indices of cands that extend span(base) (mod p)."""
        rows = [self._modrow(v) for v in base]
        r = self._rank_rows(rows, K)
        sel = []
        for j, v in enumerate(cands):
            rows2 = rows + [self._modrow(v)]
            r2 = self._rank_rows(rows2, K)
            if r2 > r:
                rows, r = rows2, r2
                sel.append(j)
        return sel

    def _independent(self, introws, K):
        if not introws:
            return []
        sel = self._extend([], [[Fraction(x) for x in r] for r in introws], K)
        return [introws[j] for j in sel]

    def basis(self, kind):
        if kind == 'jac':
            return list(self.jac_basis)
        if kind == 'ker':
            return list(self.ker_basis)
        if kind == 'full':
            return list(self.jac_basis) + list(self.ker_basis)
        raise ValueError(kind)

    def summary(self):
        return {'d': self.d, 'chi': list(self.c), 'dimR': self.dimR, 'dimE1': self.dimE1,
                'gauge': self.dim_gauge, 'jac_mod_gauge': self.dim_jac, 'ker': self.dim_ker,
                'jac_total_with_gauge': self.dim_jac_total, 'seconds': round(self.time, 2)}


# ----------------------------------------------------------------------------------------------
def koszul_generators(F, d, c, cA=(3, 0)):
    """Koszul syzygies S = e (F_j e_i - F_i e_j), deg e = d-4, in the N-character sector of (d, c).
    Returns list of (tag, [6 polys])."""
    Fi = [pderiv(F, i) for i in range(NV)]
    out = []
    if d - 4 < 0:
        return out
    for i in range(NV):
        for j in range(i + 1, NV):
            ce = cadd(cA, c, WZ[i], WZ[j])
            for e in mons_chi(d - 4, ce):
                N = [dict() for _ in range(NV)]
                N[i] = shift(Fi[j], e, 1)
                N[j] = shift(Fi[i], e, -1)
                out.append(((i, j, e), N))
    return out


def n_lift_exact_jac(A, q, tag):
    """For D = m F_i (tag=(i,m)) the exact E1 lift N = (1/q) A m e_i."""
    i, m = tag
    N = [dict() for _ in range(NV)]
    N[i] = scale(shift(A, m, 1), Fraction(1, q))
    return N


class E1Lifter:
    """Float least-squares E1 lifts: q sum_i N_i F_i = A D in the sector of (d, c)."""

    def __init__(self, F, A, q, d, c):
        self.F, self.A, self.q = F, A, q
        Fi = [pderiv(F, i) for i in range(NV)]
        cA = poly_chi(A)
        cAD = cadd(cA, c)
        cols = []
        for i in range(NV):
            for m in mons_chi(d + 1, cadd(cAD, WZ[i])):
                cols.append(((i, m), scale(shift(Fi[i], m, 1), q)))
        self.cols = cols
        ridx = {}
        for _, Pl in cols:
            for e in Pl:
                ridx.setdefault(e, len(ridx))
        self.ridx = ridx
        Mtx = np.zeros((len(ridx), len(cols)))
        for j, (_, Pl) in enumerate(cols):
            for e, v in Pl.items():
                Mtx[ridx[e], j] = float(v)
        self.M = Mtx

    def lift(self, D):
        AD = mul(self.A, D)
        b = np.zeros(len(self.ridx))
        for e, v in AD.items():
            if e not in self.ridx:
                raise ValueError('A D outside the column span support (not E1-admissible?)')
            b[self.ridx[e]] = float(v)
        x, *_ = sla.lstsq(self.M, b, cond=1e-12, lapack_driver='gelsd', check_finite=False)
        res = np.linalg.norm(self.M @ x - b) / max(np.linalg.norm(b), 1e-300)
        N = [dict() for _ in range(NV)]
        for (i, m), v in zip([t for t, _ in self.cols], x):
            if v != 0.0:
                N[i][m] = v
        return N, res


def B_bilinear(D, N):
    """B_D(N) = D div N - N . grad D  (float or Fraction coefficients)."""
    divN = {}
    for i in range(NV):
        divN = add(divN, pderiv(N[i], i))
    out = mul(D, divN)
    for i in range(NV):
        out = sub(out, mul(N[i], pderiv(D, i)))
    return out


# ----------------------------------------------------------------------------------------------
class ReducedVarPro(varpro.VarProFamily):
    """Reduced sector-restricted variable projection (see module docstring).

    basis   : list of D polys (Fraction coefficients), all of degree d and character c
    nlifts  : list of N (6 polys, float/Fraction) with q sum N_i F_i = A D_k
    koszul  : list of N (6 polys) with sum S_i F_i = 0 (may be empty: Hessian relaxation)
    """

    def __init__(self, F, A, q, basis, nlifts, koszul, c, verbose=False, cond=1e-9, label=''):
        t0 = time.time()
        self.cond = cond
        self.F, self.A, self.q = F, A, q
        self.basis = [dict(B) for B in basis]
        self.label = label
        K = len(basis); self.K = K
        degs = {core.degree(B) for B in basis}
        assert len(degs) == 1, degs
        d = degs.pop(); self.d = d
        self.pD = core.parity_of(basis[0])
        self.c = c
        cA = poly_chi(A)
        cQ = cadd(cA, c, c)
        self.degN, self.degQ = d + 6 * q - 5, 2 * d - 6
        self.Qmons = mons_chi(self.degQ, cQ) if self.degQ >= 0 else []
        self.nS, self.nQ = len(koszul), len(self.Qmons)
        self.nN = self.nS   # the parent's x[:nN] slot holds the Koszul coefficients
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        Dk_i = [[pderiv(B, i) for i in range(NV)] for B in self.basis]
        divS = []
        for S in koszul:
            dv = {}
            for i in range(NV):
                dv = add(dv, pderiv(S[i], i))
            divS.append(dv)
        blocks = []
        for k, B in enumerate(self.basis):
            ri, ci, vals = [], [], []
            for j, S in enumerate(koszul):
                col = mul(B, divS[j])
                for i in range(NV):
                    if S[i] and Dk_i[k][i]:
                        col = sub(col, mul(S[i], Dk_i[k][i]))
                for e, v in col.items():
                    ri.append(ridx(e)); ci.append(j); vals.append(float(v))
            blocks.append((ri, ci, vals))
        qri, qci, qvals = [], [], []
        negF = scale(F, -1)
        for j, m in enumerate(self.Qmons):
            col = shift(negF, m, 1)
            for e, v in col.items():
                qri.append(ridx(e)); qci.append(j); qvals.append(float(v))
        divN = []
        for N in nlifts:
            dv = {}
            for i in range(NV):
                dv = add(dv, pderiv(N[i], i))
            divN.append(dv)
        pair = {}
        for k, B in enumerate(self.basis):
            for l, N in enumerate(nlifts):
                v = mul(B, divN[l])
                for i in range(NV):
                    if N[i] and Dk_i[k][i]:
                        v = sub(v, mul(N[i], Dk_i[k][i]))
                v = scale(v, -1)
                pair[(k, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows); self.rows = R; self.row_index = rows
        self.Mk = [sp.csr_matrix((np.array(vals, dtype=float), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nS)) for ri, ci, vals in blocks]
        self.MQ = sp.csr_matrix((np.array(qvals, dtype=float), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nQ)).toarray() if self.nQ else np.zeros((R, 0))
        self.Vk = []
        for k in range(K):
            ri, ci, vals = [], [], []
            for l in range(K):
                v = add(pair[(k, l)], pair[(l, k)])
                for e, cf in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(0.5 * float(cf))
            self.Vk.append(sp.csr_matrix((np.array(vals, dtype=float), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, K)))
        self.ncols = self.nS + self.nQ
        self.P = np.eye(K); self.excluded_dim = 0; self.dimP = K
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[rvp {label}] d={d} chi={c} K={K} rows={R} cols={self.ncols} (S {self.nS}, Q {self.nQ}) '
                  f'built in {self.build_time:.1f}s', file=sys.stderr)

    def confirm(self, frac, verbose=False, group=None):
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, self.q, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'e1': info, 'D': D}
        status, data = symmetry.solve_fixed_D_symmetric(self.F, self.A, self.q, D, verbose=verbose, group=group)
        data['D'] = D
        return status, data

    def info(self):
        return {'K': self.K, 'rows': self.rows, 'nS': self.nS, 'nQ': self.nQ, 'd': self.d, 'chi': list(self.c),
                'build_seconds': round(self.build_time, 1)}


def build_family(F, A, q, d, c, kind, hessian=False, verbose=False, spaces=None):
    """kind in {'jac','ker','full'}; hessian=True drops the Koszul unknowns."""
    sp_ = spaces or SectorSpaces(F, A, q, d, c, verbose=verbose)
    basis, lifts = [], []
    if kind in ('jac', 'full'):
        for tag, B in zip(sp_.jac_tags, sp_.jac_basis):
            basis.append(B); lifts.append(n_lift_exact_jac(A, q, tag))
    if kind in ('ker', 'full') and sp_.ker_basis:
        L = E1Lifter(F, A, q, d, c)
        for B in sp_.ker_basis:
            N, res = L.lift(B)
            if res > 1e-8:
                raise RuntimeError(f'E1 lift failed for kernel vector (res {res:.2e})')
            basis.append(B); lifts.append(N)
    if not basis:
        return None, sp_
    kos = [] if hessian else [N for _, N in koszul_generators(F, d, c, poly_chi(A))]
    fam = ReducedVarPro(F, A, q, basis, lifts, kos, c, verbose=verbose,
                        label=f'{kind}{"-hess" if hessian else ""}')
    return fam, sp_


# ----------------------------------------------------------------------------------------------
def log_batch(rec, path=LOG):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('lane', LANE)
    rec.setdefault('engine', 2)
    with open(path, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')


def D_str(Pl):
    return core.to_str(Pl)
