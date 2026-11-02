"""w2lib.py -- wave-2 joint-search library for quadric-contact-residual (q = 1, A = X0^3 Y0^3).

Everything is restricted to G_aff character sectors (common.chi: (e0+e3 mod 2, e2+e5 mod 2, e1-e4 mod 6,
e3+e4+e5 mod 2); the last entry is the Y-parity).  F has character 0, A has CHI_A = (0,0,0,1).
For a denominator D with characters in a set CD, the numerators/quotient/rows are forced into
    N_i : CHI_A + c + CHI_Z[i]        Q, rows : CHI_A + c + c'        (c, c' in CD),
and restricting the unknowns to these sectors is lossless (averaging over G_aff).

Classes
  SectorVarPro   full-N variable projection over a linear family of D's (float basis), sector-restricted.
  KoszulVarPro   Jacobian-type family D = sum C_i F_i, N = A C + K (Koszul), identity divided by F.
  HessRelax      quadratic relaxation A sum C_i C_j F_ij = 0 mod J_2d for Jacobian-type families.
Functions
  colon_family(d, c)      -> exact Jacobian basis of J_d[c] (independent m F_i) + numerical extras
  jac_C_basis(d, CD)      -> C-space basis [(i, m)] for Jacobian-type families
  search_parallel(...)    -> multiprocessing random restarts (fork; family object inherited)
  reconstruct_D / confirm_exact / process_candidates -> rational reconstruction + exact re-solve (engine)
  wlog(**rec)             -> append a JSON line to wave2-search.log.jsonl
"""
import sys, os, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from common import F, A, q, chi, chi_add, chi_neg, chi_of, CHI_A, CHI_Z, P0, RES, SectorSystem, sector_e1
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
from flint import nmod_mat

W2LOG = RES + '/wave2-search.log.jsonl'
CAND = RES + '/wave2-candidates.jsonl'
p = P0
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(F, i) for i in range(NV)]
Fip = [pderiv(Fp, i, p) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
CHI0 = (0, 0, 0, 0)
ALLCHARS = [(a, b, c, y) for a in (0, 1) for b in (0, 1) for c in range(6) for y in (0, 1)]


def wlog(**rec):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    rec.setdefault('engine', 2)
    with open(W2LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec


def cstr(c):
    return ''.join(str(x) for x in c)


def chars_of_degree(d):
    return sorted({chi(e) for e in monomials(d)})


def mons_chi(d, c):
    return [e for e in monomials(d) if chi(e) == c]


def mons_chis(d, cs):
    cs = set(cs)
    return [e for e in monomials(d) if chi(e) in cs]


def chars_of_poly(P):
    return sorted({chi(e) for e in P})


def to_float(P):
    return {e: float(v) for e, v in P.items()}


def fadd(P, Q, s=1.0):
    out = dict(P)
    for e, v in Q.items():
        out[e] = out.get(e, 0.0) + s * v
    return {e: v for e, v in out.items() if v != 0.0}


def fmul(P, Q):
    out = {}
    for e1, v1 in P.items():
        for e2, v2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            out[e] = out.get(e, 0.0) + v1 * v2
    return {e: v for e, v in out.items() if v != 0.0}


def fshift(P, m, c=1.0):
    return {tuple(a + b for a, b in zip(e, m)): c * v for e, v in P.items()}


def fpderiv(P, i):
    out = {}
    for e, v in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            out[tuple(e2)] = v * e[i]
    return out


def fscale(P, s):
    return {e: s * v for e, v in P.items()}


# ---------------------------------------------------------------- exact helpers (mod p)
def modp(v):
    fr = Fraction(v)
    return (fr.numerator % p) * pow(fr.denominator % p, p - 2, p) % p


def jac_gens(deg, cs):
    """Exact generators m*F_i of J_deg with character in cs: list of (i, m, poly(Fraction))."""
    cs = set(cs)
    out = []
    if deg < 5:
        return out
    for m in monomials(deg - 5):
        cm = chi(m)
        for i in range(NV):
            if chi_add(cm, chi_neg(CHI_Z[i])) in cs:
                out.append((i, m, shift(Fi[i], m, 1)))
    return out


def indep_subset_modp(polys, mons):
    """Indices of a maximal independent subset of the given polys (mod p) over the monomial list."""
    if not polys:
        return []
    idx = {e: j for j, e in enumerate(mons)}
    n = len(mons)
    flat = [0] * (len(polys) * n)
    for r, P in enumerate(polys):
        for e, v in P.items():
            flat[r * n + idx[e]] = modp(v)
    M = nmod_mat(len(polys), n, flat, p)
    T, rk = M.transpose().rref()
    piv = []
    for i in range(rk):
        j = 0
        while int(T[i, j]) == 0:
            j += 1
        piv.append(j)
    return piv


def null_space_float(M, rtol=1e-9):
    """Orthonormal basis of the nullspace of a dense float matrix (SVD)."""
    if M.shape[0] == 0:
        return np.eye(M.shape[1])
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    tol = rtol * (s[0] if s.size else 1.0)
    r = int(np.sum(s > tol))
    return Vt[r:].T


def colon_family(d, c, verbose=False, want_extras=True):
    """E1-admissible family in S_d[c]: exact independent Jacobian generators + numerical extras.

    Returns dict(d, c, mons, jac: list of exact polys (Fraction), extras: list of float polys,
                 dimS, dimJ, dim_colon (numerical), n_extra, seconds)."""
    t0 = time.time()
    mons = mons_chi(d, c)
    n = len(mons)
    gens = jac_gens(d, [c])
    gpolys = [g[2] for g in gens]
    piv = indep_subset_modp(gpolys, mons) if gens else []
    jac = [gpolys[j] for j in piv]
    dimJ = len(jac)
    extras = []
    dim_colon = dimJ
    if want_extras and n > dimJ:
        sigma = chi_add(c, CHI_A)
        tmons = mons_chi(d + 6, sigma)
        tidx = {e: j for j, e in enumerate(tmons)}
        tg = [g[2] for g in jac_gens(d + 6, [sigma])]
        cols = [shift(A, m, 1) for m in mons] + tg
        M = np.zeros((len(tmons), len(cols)))
        for j, P in enumerate(cols):
            for e, v in P.items():
                M[tidx[e], j] = float(v)
        NS = null_space_float(M)
        Dpart = NS[:n, :]
        U, s, Vt = np.linalg.svd(Dpart, full_matrices=False)
        r = int(np.sum(s > 1e-9 * (s[0] if s.size else 1.0)))
        colon = U[:, :r]
        dim_colon = r
        if dimJ:
            JM = np.zeros((n, dimJ))
            midx = {e: j for j, e in enumerate(mons)}
            for k, P in enumerate(jac):
                for e, v in P.items():
                    JM[midx[e], k] = float(v)
            QJ, _ = np.linalg.qr(JM)
            E = colon - QJ @ (QJ.T @ colon)
        else:
            E = colon
        U2, s2, _ = np.linalg.svd(E, full_matrices=False)
        r2 = int(np.sum(s2 > 1e-7))
        for k in range(r2):
            extras.append({mons[j]: float(U2[j, k]) for j in range(n) if abs(U2[j, k]) > 1e-14})
    out = {'d': d, 'c': c, 'mons': mons, 'jac': jac, 'extras': extras, 'dimS': n, 'dimJ': dimJ,
           'dim_colon': dim_colon, 'n_extra': len(extras), 'seconds': round(time.time() - t0, 1)}
    if verbose:
        print(f'[colon d={d} c={cstr(c)}] dimS={n} dimJ={dimJ} colon={dim_colon} extras={len(extras)} ({out["seconds"]}s)', file=sys.stderr)
    return out


def jac_C_basis(d, CD):
    """C-space basis for Jacobian-type D = sum C_i F_i with chi(D) in CD: list of (i, m), deg m = d-5."""
    CD = set(CD)
    out = []
    if d < 5:
        return out
    for m in monomials(d - 5):
        cm = chi(m)
        for i in range(NV):
            if chi_add(cm, chi_neg(CHI_Z[i])) in CD:
                out.append((i, m))
    return out


# ---------------------------------------------------------------- generic variable-projection core
class _VarProBase:
    """Shared optimizer: M(a) = sum_k a_k M_k | M_const, b(a) = sum_kl a_k a_l v_kl, a on the unit
    sphere inside the column space of P (gauge complement)."""

    def _finish(self, verbose, name):
        K = self.K
        self.ncols = self.nA + self.nC
        self.dimP = self.P.shape[1]
        self._nfev = 0
        self.name = name
        if verbose:
            print(f'[{name}] d={self.d} K={K} rows={self.rows} cols={self.ncols} (a-dependent {self.nA}, const {self.nC}) '
                  f'gauge dim={self.excluded_dim} dimP={self.dimP} built {self.build_time:.1f}s', file=sys.stderr)

    def assemble(self, a):
        M = self.Mk[0] * a[0]
        for k in range(1, self.K):
            M = M + self.Mk[k] * a[k]
        M = M.toarray()
        if self.nC:
            M = np.hstack([M, self.MC])
        b = np.zeros(self.rows)
        for k in range(self.K):
            b += a[k] * (self.Vk[k] @ a)
        return M, b

    def value(self, a):
        a = np.asarray(a, dtype=float)
        a = a / np.linalg.norm(a)
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r

    def value_and_grad(self, u):
        u = np.asarray(u, dtype=float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        f = float(r @ r)
        xA = x[:self.nA]
        g = np.zeros(self.K)
        for k in range(self.K):
            dM = self.Mk[k] @ xA
            db = 2.0 * (self.Vk[k] @ a)
            g[k] = 2.0 * float(r @ (dM - db))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, b = self.assemble(a)
        return float(np.linalg.norm(b))

    def D_float(self, a):
        D = {}
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * v
        return {e: v for e, v in D.items() if v != 0.0}

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-12, verbose=False, stop_rel=1e-13, init=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r = self.value(a)
            bn = self.bnorm(a)
            rel = f / max(bn * bn, 1e-300)
            cand = {'value': f, 'rel': rel, 'a': a, 'iters': int(res.nit), 'restart': it, 'seed': seed}
            out.append(cand)
            if verbose:
                print(f'[search] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel:
                break
        out.sort(key=lambda c: c['rel'])
        return out


def _gauge_complement(basis_float, d, CD):
    """Orthonormal basis of the complement (in a-space) of {a : D(a) in F * S_{d-6}}; returns (P, excluded_dim)."""
    K = len(basis_float)
    if d < 6:
        return np.eye(K), 0
    Gm = mons_chis(d - 6, CD)
    Fcols = [to_float(shift(F, m, 1)) for m in Gm]
    mons = {}
    for B in list(basis_float) + Fcols:
        for e in B:
            mons.setdefault(e, len(mons))
    Mb = np.zeros((len(mons), K + len(Fcols)))
    for k, B in enumerate(basis_float):
        for e, v in B.items():
            Mb[mons[e], k] = v
    for j, C in enumerate(Fcols):
        for e, v in C.items():
            Mb[mons[e], K + j] = -v
    ns = null_space_float(Mb)
    S = ns[:K, :] if ns.size else np.zeros((K, 0))
    if S.shape[1] == 0:
        return np.eye(K), 0
    U, s, _ = np.linalg.svd(S, full_matrices=True)
    r = int(np.sum(s > 1e-9 * max(1.0, s[0])))
    return U[:, r:], r


class SectorVarPro(_VarProBase):
    """Full-N variable projection over D(a) = sum a_k B_k (float polys of degree d), sector-restricted."""

    def __init__(self, basis, verbose=False, cond=1e-11, name='sector', exclude_gauge=True):
        t0 = time.time()
        self.cond = cond
        self.basis = [to_float(B) for B in basis]
        K = len(self.basis); self.K = K
        degs = {core.degree(B) for B in basis}
        assert len(degs) == 1, 'family members must share one degree'
        d = degs.pop(); self.d = d
        CD = set()
        for B in self.basis:
            CD.update(chars_of_poly(B))
        self.CD = sorted(CD)
        NCH = [{chi_add(CHI_A, chi_add(c, CHI_Z[i])) for c in CD} for i in range(NV)]
        RCH = {chi_add(CHI_A, chi_add(c, c2)) for c in CD for c2 in CD}
        self.NCH, self.RCH = NCH, RCH
        degN = d + 1; degQ = 2 * d - 6
        Ncols = [(i, m) for i in range(NV) for m in monomials(degN) if chi(m) in NCH[i]]
        Qcols = mons_chis(degQ, RCH) if degQ >= 0 else []
        self.Ncols, self.Qcols = Ncols, Qcols
        self.nA, self.nC = len(Ncols), len(Qcols)
        Ff = to_float(F); Af = to_float(A)
        Fif = [to_float(P) for P in Fi]
        rows = {}
        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        blocks = []
        for Bk in self.basis:
            FBk = fmul(Ff, Bk)
            G = [fscale(fadd(fmul(Bk, Fif[i]), fmul(Ff, fpderiv(Bk, i))), -1.0) for i in range(NV)]
            ri, ci, vals = [], [], []
            for j, (i, m) in enumerate(Ncols):
                col = fshift(G[i], m)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = fadd(col, fshift(FBk, tuple(m2), float(m[i])))
                for e, v in col.items():
                    ri.append(ridx(e)); ci.append(j); vals.append(v)
            blocks.append((ri, ci, vals))
        F2 = fscale(fmul(Ff, Ff), -1.0)
        qri, qci, qvals = [], [], []
        for j, m in enumerate(Qcols):
            for e, v in fshift(F2, m).items():
                qri.append(ridx(e)); qci.append(j); qvals.append(v)
        AB = [fmul(Af, Bk) for Bk in self.basis]
        pair = {}
        for k in range(K):
            for l in range(k, K):
                v = fscale(fmul(AB[k], self.basis[l]), -1.0)
                pair[(k, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows); self.rows = R; self.row_index = rows
        self.Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nA)) for ri, ci, vals in blocks]
        self.MC = sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nC)).toarray() if self.nC else np.zeros((R, 0))
        self.Vk = []
        for k in range(K):
            ri, ci, vals = [], [], []
            for l in range(K):
                v = pair[(min(k, l), max(k, l))]
                for e, cv in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(cv)
            self.Vk.append(sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, K)))
        if exclude_gauge:
            self.P, self.excluded_dim = _gauge_complement(self.basis, d, CD)
        else:
            self.P, self.excluded_dim = np.eye(K), 0
        self.build_time = time.time() - t0
        self._finish(verbose, name)

    def unpack(self, x):
        N = [{} for _ in range(NV)]; Q = {}
        for (i, m), v in zip(self.Ncols, x[:self.nA]):
            if v: N[i][m] = v
        for m, v in zip(self.Qcols, x[self.nA:]):
            if v: Q[m] = v
        return N, Q


class KoszulVarPro(_VarProBase):
    """Jacobian-type family: D(a) = sum_k a_k D_k with D_k = m F_i for C-basis (i, m); N = A C(a) + K(z),
    K(z) = sum_j z_j Q_j (F_s e_r - F_r e_s).  Identity (divided by F): B_D(A C) + B_D(K) = F Q."""

    def __init__(self, d, CD, verbose=False, cond=1e-11, name='koszul', Cbasis=None):
        t0 = time.time()
        self.cond = cond
        self.d = d
        CD = set(CD); self.CD = sorted(CD)
        Cb = Cbasis if Cbasis is not None else jac_C_basis(d, CD)
        self.Cbasis = Cb
        K = len(Cb); self.K = K
        Ff = to_float(F); Af = to_float(A)
        Fif = [to_float(P) for P in Fi]
        self.basis = [fshift(Fif[i], m) for (i, m) in Cb]
        AC = [[{} for _ in range(NV)] for _ in Cb]
        for k, (i, m) in enumerate(Cb):
            AC[k][i] = fshift(Af, m)
        Kb = []
        if d - 4 >= 0:
            KCH = {}
            for c in CD:
                for r in range(NV):
                    for s in range(r + 1, NV):
                        KCH.setdefault((r, s), set()).add(chi_add(CHI_A, chi_add(c, chi_add(CHI_Z[r], CHI_Z[s]))))
            for m in monomials(d - 4):
                cm = chi(m)
                for r in range(NV):
                    for s in range(r + 1, NV):
                        if cm in KCH[(r, s)]:
                            Kb.append((r, s, m))
        self.Kbasis = Kb
        self.nA = len(Kb)
        RCH = {chi_add(CHI_A, chi_add(c, c2)) for c in CD for c2 in CD}
        degQ = 2 * d - 6
        Qcols = mons_chis(degQ, RCH) if degQ >= 0 else []
        self.Qcols = Qcols; self.nC = len(Qcols)
        rows = {}
        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        Kf = []
        for (r, s, m) in Kb:
            V = [{} for _ in range(NV)]
            V[r] = fshift(Fif[s], m); V[s] = fscale(fshift(Fif[r], m), -1.0)
            divV = fadd(fpderiv(V[r], r), fpderiv(V[s], s))
            Kf.append((V, divV))
        Dki_all = [[fpderiv(Dk, i) for i in range(NV)] for Dk in self.basis]
        blocks = []
        for k in range(K):
            Dk = self.basis[k]
            Dki = Dki_all[k]
            ri, ci, vals = [], [], []
            for j, (V, divV) in enumerate(Kf):
                col = fmul(Dk, divV)
                for i in range(NV):
                    if V[i] and Dki[i]:
                        col = fadd(col, fmul(V[i], Dki[i]), -1.0)
                for e, v in col.items():
                    ri.append(ridx(e)); ci.append(j); vals.append(v)
            blocks.append((ri, ci, vals))
        qri, qci, qvals = [], [], []
        Fneg = fscale(Ff, -1.0)
        for j, m in enumerate(Qcols):
            for e, v in fshift(Fneg, m).items():
                qri.append(ridx(e)); qci.append(j); qvals.append(v)
        divAC = []
        for k in range(K):
            dv = {}
            for i in range(NV):
                if AC[k][i]:
                    dv = fadd(dv, fpderiv(AC[k][i], i))
            divAC.append(dv)
        def Bkl(k, l):
            out = fmul(self.basis[k], divAC[l])
            for i in range(NV):
                if AC[l][i] and Dki_all[k][i]:
                    out = fadd(out, fmul(AC[l][i], Dki_all[k][i]), -1.0)
            return out
        pair = {}
        for k in range(K):
            for l in range(k, K):
                if k == l:
                    v = fscale(Bkl(k, k), -1.0)
                else:
                    v = fscale(fadd(Bkl(k, l), Bkl(l, k)), -0.5)
                pair[(k, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows); self.rows = R; self.row_index = rows
        self.Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nA)) for ri, ci, vals in blocks]
        self.MC = sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nC)).toarray() if self.nC else np.zeros((R, 0))
        self.Vk = []
        for k in range(K):
            ri, ci, vals = [], [], []
            for l in range(K):
                v = pair[(min(k, l), max(k, l))]
                for e, cv in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(cv)
            self.Vk.append(sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, K)))
        self.P, self.excluded_dim = _gauge_complement(self.basis, d, CD)
        self.build_time = time.time() - t0
        self._finish(verbose, name)


class HessRelax:
    """Quadratic relaxation for Jacobian-type D = sum C_i F_i:  A sum_ij C_i C_j F_ij in J_2d.
    f(c) = || proj_{J_2d^perp} ( sum_kl c_k c_l h_kl ) ||^2 on the unit sphere of the gauge complement."""

    def __init__(self, d, CD, verbose=False, Cbasis=None):
        t0 = time.time()
        self.d = d
        CD = set(CD); self.CD = sorted(CD)
        Cb = Cbasis if Cbasis is not None else jac_C_basis(d, CD)
        self.Cbasis = Cb
        K = len(Cb); self.K = K
        Af = to_float(A)
        AFij = [[fmul(Af, to_float(Fij[i][j])) for j in range(NV)] for i in range(NV)]
        RCH = sorted({chi_add(CHI_A, chi_add(c, c2)) for c in CD for c2 in CD})
        self.sectors = {}
        tot = 0
        for sgm in RCH:
            mons = mons_chi(2 * d, sgm)
            idx = {e: j for j, e in enumerate(mons)}
            gens = [g[2] for g in jac_gens(2 * d, [sgm])]
            JM = np.zeros((len(mons), len(gens)))
            for j, P in enumerate(gens):
                for e, v in P.items():
                    JM[idx[e], j] = float(v)
            if gens:
                U, s, _ = np.linalg.svd(JM, full_matrices=True)
                r = int(np.sum(s > 1e-9 * s[0]))
                Uc = U[:, r:]
            else:
                Uc = np.eye(len(mons))
            self.sectors[sgm] = (mons, idx, Uc, tot)
            tot += Uc.shape[1]
        self.dimR = tot
        T = np.zeros((tot, K, K))
        for k, (i, m) in enumerate(Cb):
            for l in range(k, K):
                (j, m2) = Cb[l]
                h = fshift(AFij[i][j], tuple(a + b for a, b in zip(m, m2)))
                if not h:
                    continue
                sgm = chi(next(iter(h)))
                mons, idx, Uc, o = self.sectors[sgm]
                v = np.zeros(len(mons))
                for e, val in h.items():
                    v[idx[e]] = val
                pv = Uc.T @ v
                T[o:o + Uc.shape[1], k, l] += pv
                if l != k:
                    T[o:o + Uc.shape[1], l, k] += pv
        self.T = T
        self.Tnorm = float(np.sqrt(np.sum(T ** 2)))
        basis_float = [fshift(to_float(Fi[i]), m) for (i, m) in Cb]
        self.basis = basis_float
        self.P, self.excluded_dim = _gauge_complement(basis_float, d, CD)
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self.name = 'hess'
        if verbose:
            print(f'[hess] d={d} K={K} dimR={tot} sectors={len(RCH)} gauge={self.excluded_dim} dimP={self.dimP} built {self.build_time:.1f}s', file=sys.stderr)

    def value(self, c):
        c = np.asarray(c, float); c = c / np.linalg.norm(c)
        G = np.einsum('rkl,l->rk', self.T, c)
        v = G @ c
        return float(v @ v)

    def value_and_grad(self, u):
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        c = a_full / na
        G = np.einsum('rkl,l->rk', self.T, c)
        v = G @ c
        f = float(v @ v)
        g = 4.0 * (G.T @ v)
        g_full = (g - float(g @ c) * c) / na
        return f, self.P.T @ g_full

    def bnorm(self, c):
        return self.Tnorm / max(self.K, 1)

    def search(self, n_restarts=10, seed=0, maxiter=500, gtol=1e-12, verbose=False, stop_rel=1e-13, init=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18})
            c = self.P @ res.x; c /= np.linalg.norm(c)
            f = self.value(c)
            bn = self.bnorm(c)
            rel = f / max(bn * bn, 1e-300)
            out.append({'value': f, 'rel': rel, 'a': c, 'iters': int(res.nit), 'restart': it, 'seed': seed})
            if verbose:
                print(f'[hess] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel:
                break
        out.sort(key=lambda c: c['rel'])
        return out

    def D_float(self, c):
        D = {}
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + c[k] * v
        return {e: v for e, v in D.items() if v != 0.0}


# ---------------------------------------------------------------- parallel restarts
_FAM = None


def _worker(args):
    seed, n, maxiter = args
    fam = _FAM
    t0 = time.time()
    cands = fam.search(n_restarts=n, seed=seed, maxiter=maxiter)
    rels = [c['rel'] for c in cands]
    best = cands[0]
    out = {'seed': seed, 'n': n, 'best_rel': best['rel'], 'best_value': best['value'], 'rels': rels,
           'iters_mean': float(np.mean([c['iters'] for c in cands])), 'seconds': round(time.time() - t0, 1),
           'best_a': best['a'].tolist()}
    out['low'] = [{'rel': c['rel'], 'a': c['a'].tolist(), 'seed': seed, 'restart': c['restart']} for c in cands if c['rel'] < 1e-8]
    return out


def rel_hist(rels):
    hist = {}
    for x in rels:
        b = 'lt1e-12' if x < 1e-12 else ('lt1e-10' if x < 1e-10 else ('lt1e-8' if x < 1e-8 else ('lt1e-6' if x < 1e-6 else ('lt1e-4' if x < 1e-4 else ('lt1e-2' if x < 1e-2 else 'ge1e-2')))))
        hist[b] = hist.get(b, 0) + 1
    return hist


def search_parallel(fam, n_total, workers, seed_base=0, per_task=None, maxiter=300, verbose=True):
    """Run n_total restarts over `workers` processes (fork; `fam` inherited).  Returns merged summary."""
    global _FAM
    _FAM = fam
    import multiprocessing as mp
    per_task = per_task or max(1, min(20, n_total // max(workers, 1)))
    tasks = []
    left = n_total; s = seed_base
    while left > 0:
        n = min(per_task, left)
        tasks.append((s, n, maxiter)); s += 1; left -= n
    t0 = time.time()
    results = []
    if workers <= 1:
        for t in tasks:
            results.append(_worker(t))
    else:
        ctx = mp.get_context('fork')
        with ctx.Pool(workers) as pool:
            for r in pool.imap_unordered(_worker, tasks):
                results.append(r)
                if verbose:
                    print(f'  task seed={r["seed"]} n={r["n"]} best={r["best_rel"]:.3e} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    rels = sorted(x for r in results for x in r['rels'])
    low = [c for r in results for c in r['low']]
    best = min(results, key=lambda r: r['best_rel'])
    quant = {k: rels[min(len(rels) - 1, int(f * len(rels)))] for k, f in (('min', 0), ('q10', 0.1), ('median', 0.5), ('q90', 0.9), ('max', 0.999))} if rels else {}
    return {'n_restarts': len(rels), 'workers': workers, 'best_rel': best['best_rel'], 'best_a': best['best_a'],
            'hist': rel_hist(rels), 'quantiles': quant, 'low': low, 'seconds': round(time.time() - t0, 1),
            'iters_mean': float(np.mean([r['iters_mean'] for r in results]))}


# ---------------------------------------------------------------- rational reconstruction + exact confirmation
def reconstruct_D(Dfloat, max_dens=(1, 6, 12, 60, 360, 2520, 10**4, 10**5, 10**6), tol=1e-7):
    """Rational candidates (exact polys) for a float D, normalized so its largest coefficient is 1."""
    if not Dfloat:
        return []
    j = max(Dfloat, key=lambda e: abs(Dfloat[e]))
    s = Dfloat[j]
    out = []
    for drop in (0.0, 1e-9, 1e-7):
        Dn = {e: v / s for e, v in Dfloat.items() if abs(v / s) > drop}
        for md in max_dens:
            cand = {}
            ok = True
            for e, v in Dn.items():
                fr = Fraction(v).limit_denominator(md)
                if abs(float(fr) - v) > tol:
                    ok = False; break
                if fr != 0:
                    cand[e] = fr
            if ok and cand and cand not in out:
                out.append(cand)
    return out


def confirm_exact(D, verbose=False):
    """Exact re-solve for a rational D: sector mod-p test (if eigen) then engine exact solve.  Returns (status, data)."""
    if not core.is_homogeneous(D):
        return 'inconsistent', {'note': 'inhomogeneous'}
    divF, _ = jacobian.divisible_by(F, D)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    cD = chi_of(D)
    if cD is not None:
        ok, info = sector_e1(D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails (sector)', 'e1': info}
        S = SectorSystem(D, verbose=verbose)
        r = S.solve(verbose=verbose)
        if not r['consistent']:
            return 'inconsistent', {'note': 'sector system inconsistent mod p', 'rank': r['rank'], 'nullity': r['nullity']}
    elif core.parity_of(D) is not None:
        ok, info = fixed_d.e1_solvable(F, A, q, D)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails (generic)', 'e1': info}
    if core.parity_of(D) is None:
        return 'unknown', {'note': 'mixed Y-parity D: engine exact solver needs pure parity; use D*tau(D)'}
    try:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
        if not out['consistent']:
            return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity'], 'note': 'engine compressed inconsistent'}
        if out.get('exact'):
            return 'solved', {'N': out['N'], 'Q': out['Q'], 'rank': out['rank'], 'nullity': out['nullity']}
    except Exception as ex:
        if verbose:
            print('[confirm] fast solver failed', ex, file=sys.stderr)
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
    return status, data


def save_witness(D, N, Q, note=''):
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    divF, _ = jacobian.divisible_by(F, D)
    w = {'d': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q),
         'verified': bool(ok), 'F_divides_D': bool(divF), 'note': note}
    with open(RES + '/witness-candidate.json', 'w') as f:
        json.dump(w, f)
    return ok and not divF, info


def process_candidates(low, fam, family, d, sector, verbose=True):
    """For every low-residual candidate: reconstruct D rationally, re-solve exactly, log."""
    results = []
    tried = set()
    for c in sorted(low, key=lambda c: c['rel']):
        a = np.asarray(c['a'], float)
        Df = fam.D_float(a)
        with open(CAND, 'a') as f:
            f.write(json.dumps({'family': family, 'd': d, 'sector': sector, 'rel': c['rel'], 'a': c['a'],
                                'Dfloat': {str(e): v for e, v in Df.items()}}) + '\n')
        for Dr in reconstruct_D(Df):
            key = core.to_str(Dr)
            if key in tried:
                continue
            tried.add(key)
            t0 = time.time()
            status, data = confirm_exact(Dr, verbose=verbose)
            rec = {'kind': 'confirm', 'family': family, 'd': d, 'sector': sector, 'rel': c['rel'], 'status': status,
                   'D': key[:3000], 'nterms': len(Dr), 'note': str(data.get('note', '')), 'seconds': round(time.time() - t0, 1)}
            wlog(**rec)
            results.append(rec)
            if status == 'solved':
                ok, info = save_witness(Dr, data['N'], data['Q'], note=f'{family} d={d} sector={sector}')
                wlog(kind='WITNESS', family=family, d=d, sector=sector, verified=ok, info=str(info)[:500])
                print('!!!!! WITNESS FOUND', key[:200], file=sys.stderr, flush=True)
                return results, True
    return results, False


# ---------------------------------------------------------------- mod-p colon dimensions (cross-check)
def colon_dims_modp(d, c):
    """dim S_d[c], dim J_d[c], dim (J:A)_d[c] computed mod p by ranks."""
    mons = mons_chi(d, c)
    n = len(mons)
    gens = [g[2] for g in jac_gens(d, [c])]
    dimJ = len(indep_subset_modp(gens, mons)) if gens else 0
    sigma = chi_add(c, CHI_A)
    tmons = mons_chi(d + 6, sigma)
    tidx = {e: j for j, e in enumerate(tmons)}
    tg = [g[2] for g in jac_gens(d + 6, [sigma])]
    cols = [shift(A, m, 1) for m in mons]
    def mat(polys):
        flat = [0] * (len(polys) * len(tmons))
        for r, P in enumerate(polys):
            for e, v in P.items():
                flat[r * len(tmons) + tidx[e]] = modp(v)
        return nmod_mat(len(polys), len(tmons), flat, p)
    rg = mat(tg).rank() if tg else 0
    ra = mat(cols + tg).rank() if (cols or tg) else 0
    colon = n - (ra - rg)
    return {'dimS': n, 'dimJ': dimJ, 'colon': colon, 'extra': colon - dimJ, 'dimR_d6': len(tmons) - rg}


def sector_list(spec, d):
    """Parse a sector spec: 'all' (one mixed family), 'each' (every character separately), 'even'/'odd'
    (characters of that Y-parity, separately), or comma-separated 4-digit strings."""
    chars = chars_of_degree(d)
    if spec == 'all':
        return [chars]
    if spec == 'each':
        return [[c] for c in chars]
    if spec in ('even', 'odd'):
        y = 0 if spec == 'even' else 1
        return [[c] for c in chars if c[3] == y]
    if spec == 'alleven' or spec == 'allodd':
        y = 0 if spec == 'alleven' else 1
        return [[c for c in chars if c[3] == y]]
    out = []
    for tok in spec.split(','):
        out.append([tuple(int(x) for x in tok)])
    return out


# ---------------------------------------------------------------- sector-restricted served-target space (q = 1)
_EPS_CACHE = {}


def eps_sector_basis():
    """Monomials of S_6[CHI_A] and the normal-form machinery of J_6[CHI_A] (mod p): returns (mons, idx, reducer rows)."""
    if 'b' in _EPS_CACHE:
        return _EPS_CACHE['b']
    mons = mons_chi(6, CHI_A)
    idx = {e: j for j, e in enumerate(mons)}
    gens = [g[2] for g in jac_gens(6, [CHI_A])]
    flat = [0] * (len(gens) * len(mons))
    for r, P in enumerate(gens):
        for e, v in P.items():
            flat[r * len(mons) + idx[e]] = modp(v)
    JM = nmod_mat(len(gens), len(mons), flat, p)
    JR, rk = JM.rref()
    Jrows = [[int(JR[i, j]) for j in range(len(mons))] for i in range(rk)]
    _EPS_CACHE['b'] = (mons, idx, Jrows)
    return _EPS_CACHE['b']


def served_sector(D, verbose=False):
    """Space of degree-6 targets A' in the CHI_A sector served by D (q=1), mod p, as row vectors over S_6[CHI_A]
    together with J_6[CHI_A]; returns dict(dim_mod_J, rows (list of int vectors incl. J rows), nrows, ncols, A_in)."""
    t0 = time.time()
    cD = chi_of(D)
    assert cD is not None, 'D must be a G_aff eigenvector'
    d = core.degree(D)
    Dp = core.reduce_mod(D, p)
    Fp_ = Fp
    Fi_ = Fip
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    FD = mul(Fp_, Dp, p)
    G = [scale(add(mul(Dp, Fi_[i], p), mul(Fp_, Di[i], p), p), p - 1, p) for i in range(NV)]
    F2 = scale(power(Fp_, 2, p), p - 1, p)
    D2 = mul(Dp, Dp, p)
    cols, polys = [], []
    for i in range(NV):
        ci = chi_add(CHI_A, chi_add(cD, CHI_Z[i]))
        for m in monomials(d + 1):
            if chi(m) != ci: continue
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); polys.append(col)
    cQ = chi_add(CHI_A, chi_add(cD, cD))
    if 2 * d - 6 >= 0:
        for m in monomials(2 * d - 6):
            if chi(m) == cQ:
                cols.append(('Q', -1, m)); polys.append(shift(F2, m, 1, p))
    amons, aidx, Jrows = eps_sector_basis()
    for m in amons:
        cols.append(('A', -1, m)); polys.append(shift(D2, m, 1, p))
    rows = {}
    for P in polys:
        for e in P:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nul = M.nullspace()
    a_idx = [j for j, c in enumerate(cols) if c[0] == 'A']
    vecs = []
    if nul:
        proj = nmod_mat(nul, len(amons), [int(X[j, k]) for k in range(nul) for j in a_idx], p)
        PR, rk = proj.rref()
        vecs = [[int(PR[i, j]) for j in range(len(amons))] for i in range(rk)]
    allrows = Jrows + vecs
    Mj = nmod_mat(len(Jrows), len(amons), [x for r in Jrows for x in r], p) if Jrows else None
    Ma = nmod_mat(len(allrows), len(amons), [x for r in allrows for x in r], p)
    rj = Mj.rank() if Mj is not None else 0
    ra = Ma.rank()
    Avec = [0] * len(amons)
    for e, v in Ap.items():
        Avec[aidx[e]] = v
    Mb = nmod_mat(len(allrows) + 1, len(amons), [x for r in allrows + [Avec] for x in r], p)
    A_in = (Mb.rank() == ra)
    return {'dim_mod_J': ra - rj, 'dim_raw': len(vecs), 'rows': vecs, 'nrows': R, 'ncols': C, 'A_in': A_in,
            'seconds': round(time.time() - t0, 1), 'd': d, 'chi': cstr(cD)}


def span_test(list_of_rows):
    """Given lists of served row vectors (mod p, over S_6[CHI_A]), return (dim span mod J, A in span)."""
    amons, aidx, Jrows = eps_sector_basis()
    allrows = list(Jrows)
    for rows in list_of_rows:
        allrows += rows
    Avec = [0] * len(amons)
    for e, v in Ap.items():
        Avec[aidx[e]] = v
    Ma = nmod_mat(len(allrows), len(amons), [x for r in allrows for x in r], p)
    Mb = nmod_mat(len(allrows) + 1, len(amons), [x for r in allrows + [Avec] for x in r], p)
    ra = Ma.rank()
    return ra - len(Jrows), Mb.rank() == ra
