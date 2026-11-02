"""Shared helpers for the wave-2 exclusion lens, lane four-plus-two-residual.
Everything mod p uses python-flint nmod_mat (dense). Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5.
F = H(X0,X1,Y0,Y1) + B(X2,Y2); CORE = [0,1,3,4]; a = A = X0^3 Y0^3 is core-only.
"""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, mul, add, scale, pderiv, reduce_mod, ratrecon, crt_pair, degree
from flint import nmod_mat
from fractions import Fraction

LANE = '/tmp/hodge/results/four-plus-two-residual'
LOG = LANE + '/wave2-excl.log.jsonl'
F, A, q, rec = core.load_instance('four-plus-two-residual')
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
B = {e: c for e, c in F.items() if e[2] or e[5]}
assert add(H, B) == F
CORE = [0, 1, 3, 4]
PRIMES = core.PRIMES

def log(ent):
    ent = dict(ent); ent.setdefault('t', time.strftime('%Y-%m-%d %H:%M:%S'))
    with open(LOG, 'a') as f: f.write(json.dumps(ent) + '\n')
    print(json.dumps(ent), flush=True)

def cmon(deg, par=None):
    """core monomials (no X2, Y2) of given degree and Y-parity"""
    return [e for e in monomials(deg, par) if e[2] == 0 and e[5] == 0]

def build(cols, p, rows=None):
    """columns = list of polys (dict mono->coeff); returns (nmod_mat, row index dict)"""
    rows = {} if rows is None else rows
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p), rows

def build_rows(polys, mons, p):
    """polys as rows in the basis `mons` (list of exponent tuples)"""
    idx = {m: i for i, m in enumerate(mons)}
    R, C = len(polys), len(mons)
    flat = [0] * (R * C)
    for i, P in enumerate(polys):
        for e, v in P.items(): flat[i * C + idx[e]] = v % p
    return nmod_mat(R, C, flat, p)

# ---------- Hilbert functions (exact, characteristic 0) ----------
def poly_mul(a, b):
    out = {}
    for (i, s), x in a.items():
        for (j, t), y in b.items():
            k = (i + j, (s + t) % 2); out[k] = out.get(k, 0) + x * y
    return out
def hilbert_bigraded(nvars_even, nvars_odd, gens_even, gens_odd, maxdeg=40):
    """bigraded (degree, Y-parity) Hilbert function of R/(regular sequence of quintics):
    R has nvars_even Y-even variables and nvars_odd Y-odd ones; the sequence has gens_even Y-even quintics and gens_odd Y-odd ones.
    Series: prod (1 - t^5 s^{p_i}) / prod (1 - t s^{p_j}). Returns dict (deg, par) -> dim."""
    # numerator and denominator expansions truncated at maxdeg
    def series_inv(par):  # 1/(1 - t s^par) = sum t^k s^{k par}
        return {(k, (k * par) % 2): 1 for k in range(maxdeg + 1)}
    S = {(0, 0): 1}
    for _ in range(nvars_even): S = poly_mul(S, series_inv(0))
    for _ in range(nvars_odd): S = poly_mul(S, series_inv(1))
    S = {k: v for k, v in S.items() if k[0] <= maxdeg}
    for par, n in ((0, gens_even), (1, gens_odd)):
        for _ in range(n):
            S = poly_mul(S, {(0, 0): 1, (5, par): -1})
            S = {k: v for k, v in S.items() if k[0] <= maxdeg}
    return S
HF_RH = hilbert_bigraded(2, 2, 2, 2)      # R(H) = C[X0,X1,Y0,Y1]/J(H)
HF_RF = hilbert_bigraded(3, 3, 3, 3)      # R(F)
HF_RB = hilbert_bigraded(1, 1, 1, 1)      # R(B)
def dim_RH(k, par): return HF_RH.get((k, par), 0)
def dim_RF(k, par): return HF_RF.get((k, par), 0)
def dim_RB(k, par): return HF_RB.get((k, par), 0)
def dim_Jcore(k, par):  # dim J(H)_k^{par} in C[X0,X1,Y0,Y1]
    return len(cmon(k, par)) - dim_RH(k, par)
def dim_JF(k, par):
    return len(monomials(k, par)) - dim_RF(k, par)

# ---------- 4-variable core kernel with J-rank certification ----------
def core_kernel(k, par, p, want_basis=True):
    """{D core, deg k, parity par : a D in J(H)} mod p.  Returns dict with kernel dims and certification flags:
    rank_J_{k+6} and rank_J_k mod p compared with the exact (char 0) values from the Hilbert function.
    If both match, dim(kernel_Q / J(H)_k) <= dim(kernel_p) - dim J_k  (modular = upper bound)."""
    t0 = time.time()
    Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
    Hi = [pderiv(Hp, i, p) for i in CORE]
    mons = cmon(k, par)
    tpar = (1 + par) % 2
    Jcols = []
    for idx, i in enumerate(CORE):
        for m2 in cmon(k + 1, (tpar - (0 if i < 3 else 1)) % 2):
            Jcols.append(shift(Hi[idx], m2, 1, p))
    Acols = [shift(ap, m, 1, p) for m in mons]
    MJ, _ = build(Jcols, p); rJ = MJ.rank()
    MAJ, _ = build(Acols + Jcols, p); rAJ = MAJ.rank()
    dimKer = len(mons) - (rAJ - rJ)
    Jk = []
    if k >= 5:
        for idx, i in enumerate(CORE):
            for m3 in cmon(k - 5, (par - (0 if i < 3 else 1)) % 2):
                Jk.append(shift(Hi[idx], m3, 1, p))
    rJk = build(Jk, p)[0].rank() if Jk else 0
    out = dict(k=k, par=par, p=p, n_mons=len(mons), kernel_dim_p=dimKer, rank_J_k6_p=rJ, dim_J_k6_exact=dim_Jcore(k + 6, tpar),
               rank_J_k_p=rJk, dim_J_k_exact=dim_Jcore(k, par), kernel_mod_J_p=dimKer - rJk,
               ring_dim_exact=dim_RH(k, par), seconds=round(time.time() - t0, 2))
    out['certified_upper_bound'] = (rJ == out['dim_J_k6_exact']) and (rJk == out['dim_J_k_exact'])
    if want_basis and dimKer - rJk > 0:
        X, nullity = MAJ.nullspace()
        nA = len(mons)
        vecs = [[int(X[j, kk]) for j in range(nA)] for kk in range(nullity)]
        vecs = [v for v in vecs if any(v)]
        Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
        Mr, rk = Mv.rref()
        basis = [{mons[j]: int(Mr[i, j]) for j in range(nA) if int(Mr[i, j])} for i in range(rk)]
        # independent modulo J(H)_k
        indep = []; cur = list(Jk); currank = rJk
        for b in basis:
            Mt, _ = build(cur + [b], p)
            if Mt.rank() > currank:
                indep.append(b); cur.append(b); currank += 1
        out['basis'] = indep
    return out

def jacobian_rref(deg, p, par=None, Fpoly=None):
    """rref of J(F)_deg (parity par or all) as rows in the monomial basis of degree deg; returns (rows matrix, mons, pivots)"""
    Fpoly = F if Fpoly is None else Fpoly
    Fp = reduce_mod(Fpoly, p); Fi = [pderiv(Fp, i, p) for i in range(6)]
    mons = monomials(deg, par)
    gens = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(deg - 5, None if par is None else (par - pi) % 2):
            gens.append(shift(Fi[i], m, 1, p))
    M = build_rows(gens, mons, p)
    R, rk = M.rref()
    piv = []
    j = 0
    for i in range(rk):
        while int(R[i, j]) == 0: j += 1
        piv.append(j); j += 1
    return R, rk, mons, piv

class NormalFormer:
    """normal form modulo J(F)_deg via rref rows: NF(v) = v - sum_i v[piv_i] * row_i"""
    def __init__(self, deg, p, par=None, Fpoly=None):
        self.p = p; self.deg = deg
        R, rk, mons, piv = jacobian_rref(deg, p, par, Fpoly)
        self.R = nmod_mat(rk, len(mons), [int(R[i, j]) for i in range(rk) for j in range(len(mons))], p) if rk else None
        self.rk = rk; self.mons = mons; self.piv = piv; self.idx = {m: i for i, m in enumerate(mons)}
        self.std = [j for j in range(len(mons)) if j not in set(piv)]
    def nf_matrix(self, M):
        """M: nmod_mat rows = polys in basis self.mons. Returns NF rows restricted to standard monomials (columns self.std)."""
        p = self.p
        if self.R is None: return M
        n = M.nrows()
        sub = nmod_mat(n, self.rk, [int(M[i, j]) for i in range(n) for j in self.piv], p)
        NF = M - sub * self.R
        return nmod_mat(n, len(self.std), [int(NF[i, j]) for i in range(n) for j in self.std], p)

def particular_solution(Mcols, rhs, p):
    """solve sum_j x_j col_j = rhs mod p; returns vector x (list) or None if inconsistent"""
    M, rows = build(Mcols + [rhs], p)
    R, rk = M.rref()
    n = len(Mcols)
    x = [0] * n
    j = 0
    for i in range(rk):
        while int(R[i, j]) == 0: j += 1
        if j == n: return None   # pivot in rhs column -> inconsistent
        x[j] = int(R[i, n]); j += 1
    return x

def e1_particular_N(D, p, Fpoly=None, Apoly=None):
    """N with N(F) = A D (mod p), deg N_i = deg D + 1; parity-restricted columns. Returns list of 6 polys."""
    Fpoly = F if Fpoly is None else Fpoly; Apoly = A if Apoly is None else Apoly
    Fp = reduce_mod(Fpoly, p); Ap = reduce_mod(Apoly, p); Dp = reduce_mod(D, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    d = degree(D); parD = core.parity_of(D); parA = core.parity_of(Apoly)
    tpar = (parA + parD) % 2
    cols = []; info = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (tpar - pi) % 2):
            cols.append(shift(Fi[i], m, 1, p)); info.append((i, m))
    rhs = mul(Ap, Dp, p)
    x = particular_solution(cols, rhs, p)
    if x is None: return None
    N = [dict() for _ in range(6)]
    for (i, m), v in zip(info, x):
        if v: N[i][m] = v
    return N

def vecfield_apply(N, P, p):
    """N(P) = sum N_i dP/dz_i"""
    out = {}
    for i in range(6):
        if N[i]: out = add(out, mul(N[i], pderiv(P, i, p), p), p)
    return out
def divergence(N, p):
    out = {}
    for i in range(6):
        if N[i]: out = add(out, pderiv(N[i], i, p), p)
    return out
def e2_expr(D, N, p):
    """D div N - N(D)"""
    return add(mul(D, divergence(N, p), p), scale(vecfield_apply(N, D, p), p - 1, p), p)
