"""Wave-2 library for klein-coupled-character (q = 1, A = X0^3 Y0^3).

Character-sector restricted joint (denominator + numerator) search in the Koszul
correction-map formulation (ch. 5.2), plus sector E1-kernel computations with
rational lifts, and the Hessian (mod-J) relaxation.

Conventions (kcc.py / kcc2.py, validated in wave 1 against accepted witnesses):
  sigma weights W = (1,2,4,0,0,0) mod 7; F has weight 0, A has weight a = 3.
  For a sigma-eigen denominator D of weight k and Y-parity pD:
     N_i : degree d+1, parity (pA - pD - PI_i) % 2, weight a + k + W_i
     p_ij: degree d-4, parity (pA + pD + PI_i + PI_j) % 2, weight a + k + W_i + W_j
     Q   : degree 2d-6, parity pA, weight a + 2k
  Koszul formulation (q = 1): with an E1 lift N0 (N0(F) = A D) every numerator is
  N = N0 + K(p), K(p)_i = sum_j p_ij F_j, p antisymmetric, and the certificate
  condition is   B_D(N0) + sum_{i<j} p_ij L_ij = F Q,
     B_D(N) = D div N - N(D),
     L_ij(m) = D[(d_i m) F_j - (d_j m) F_i] - m (F_j D_i - F_i D_j).
  For a family D(a) = sum_k a_k D_k with lifts N_k, N0(a) = sum a_k N_k and
     M(a) x = b(a),  M(a) = [ sum_k a_k Lmat_k | -F * Qmons ],  b(a) = -sum_{k,l} a_k a_l B_{D_k}(N_l).
"""
from __future__ import annotations
import os, sys, time, json, math, pickle
from fractions import Fraction
LANE = '/tmp/hodge/results/klein-coupled-character'
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, LANE)
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import minimize
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpz
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
import kcc, kcc2
from kcc import W, wt, sigma_weight, wmonomials, P0
from kcc2 import Ctx, e1_lift, dense_from_cols, rref_solve, PI, PAIRS

LOG = LANE + '/wave2-search.log.jsonl'
IID = 'klein-coupled-character'


def wlog(entry):
    entry.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('engine', 'w2-koszul-varpro')
    entry.setdefault('instance', IID)
    line = json.dumps(entry, default=str)
    with open(LOG, 'a') as f:
        f.write(line + '\n')


def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q


def fl(P):
    return {e: float(v) for e, v in P.items()}


def poly_norm(P):
    return math.sqrt(sum(float(v) ** 2 for v in P.values()))


# ----------------------------------------------------------------------------- sector E1 spaces
def sector_kernel_full(ctx, d, k, parD):
    """Canonical (rref) basis mod p of {D in S_d[parD,k] : A D in J_{d+6}} in D-monomial coordinates.

    Returns dict(Dm, kervecs (list of int lists), nl (dim (R/J)_{d+6} target sector), jd (dim J_d sector),
    jkeys (J_d generator keys (i,m)), jpolys)."""
    p = ctx.p
    Dm = wmonomials(d, parD, k)
    out = {'Dm': Dm, 'kervecs': [], 'nl': 0, 'jd': 0, 'jkeys': [], 'jpolys': []}
    if not Dm:
        return out
    tpar = (parD + ctx.pA) % 2
    tk = (k + ctx.a) % 7
    tgt = wmonomials(d + 6, tpar, tk)
    rows = {m: i for i, m in enumerate(tgt)}
    jkeys, jpolys = kcc2.jac_block_cols(ctx, d + 6, tpar, tk)
    if jpolys:
        MJ, _ = dense_from_cols(jpolys, None, p, rows)
        L, nl = MJ.transpose().nullspace()
    else:
        nl = len(rows)
        L = None
    AD, _ = dense_from_cols([shift(ctx.A, m, 1, p) for m in Dm], None, p, rows)
    if L is not None:
        Lt = nmod_mat(nl, len(rows), [int(L[r, c]) for c in range(nl) for r in range(len(rows))], p)
        Mult = Lt * AD
    else:
        Mult = AD
    K, nk = Mult.nullspace()
    kervecs = [[int(K[r, c]) for r in range(len(Dm))] for c in range(nk)]
    jdkeys, jdpolys = kcc2.jac_block_cols(ctx, d, parD, k)
    jd = 0
    if jdpolys:
        Jmat, _ = dense_from_cols(jdpolys, None, p, {m: i for i, m in enumerate(Dm)})
        jd = Jmat.rank()
    out.update({'kervecs': kervecs, 'nl': nl, 'jd': jd, 'jkeys': jdkeys, 'jpolys': jdpolys, 'nk': nk})
    return out


def crt_lift_vectors(getter, nprimes=12, verbose=False):
    """getter(p) -> list of int vectors (canonical mod p).  CRT + rational reconstruction.
    Returns list of Fraction vectors or None."""
    acc = None
    for idx, p in enumerate(core.PRIMES[:nprimes]):
        vecs = getter(p)
        if vecs is None:          # unlucky prime (pivot structure differs / E1 fails mod p): skip it
            continue
        if acc is None:
            acc = ([list(v) for v in vecs], p)
            continue
        xs, m = acc
        if len(xs) != len(vecs) or any(len(a) != len(b) for a, b in zip(xs, vecs)):
            if verbose:
                print('[crt] shape mismatch between primes; restarting', file=sys.stderr)
            acc = ([list(v) for v in vecs], p)
            continue
        xs = [[core.crt_pair(a, m, b, p)[0] for a, b in zip(v1, v2)] for v1, v2 in zip(xs, vecs)]
        acc = (xs, m * p)
        rat = [[core.ratrecon(a, acc[1]) for a in v] for v in acc[0]]
        if all(a is not None for v in rat for a in v):
            return rat
    return None


def lll_small_basis(int_rows):
    """LLL-reduce the lattice spanned by integer row vectors; returns reduced integer rows."""
    if not int_rows:
        return []
    M = fmpz_mat(len(int_rows), len(int_rows[0]), [int(x) for r in int_rows for x in r])
    R = M.lll()
    out = []
    for i in range(R.nrows()):
        row = [int(R[i, j]) for j in range(R.ncols())]
        if any(row):
            out.append(row)
    return out


def clear_denominators(vec):
    den = 1
    for v in vec:
        den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    return [int(Fraction(v) * den) for v in vec]


def sector_space(F, A, d, parD, k, nprimes=12, verbose=False, want_rational=True):
    """E1-admissible denominators of degree d, parity parD, sigma-weight k.

    Returns dict with dimensions and (if want_rational and the lift succeeds) small-height rational
    representatives of the kernel modulo J_d ('reps': list of exact polys) plus J generators."""
    t0 = time.time()
    ctx = Ctx(F, A, P0)
    S = sector_kernel_full(ctx, d, k, parD)
    Dm = S['Dm']
    info = {'d': d, 'parD': parD, 'k': k, 'nD': len(Dm), 'ring_dim_target': S['nl'], 'kernel_total': len(S['kervecs']),
            'J_d': S['jd'], 'nontrivial': len(S['kervecs']) - S['jd'], 'jkeys': S['jkeys'], 'p': P0}
    reps = []
    if want_rational and info['nontrivial'] > 0:
        def getter(p):
            c = Ctx(F, A, p)
            s = sector_kernel_full(c, d, k, parD)
            return s['kervecs']
        rat = crt_lift_vectors(getter, min(nprimes, 4), verbose)
        info['lift_ok'] = rat is not None
        if rat is None or (rat and max(max(abs(Fraction(x).numerator) for x in v) for v in rat) > 10**9):
            # canonical rref basis has large height: find short kernel vectors by LLL on the mod-m kernel lattice instead
            Dm2, good, nk2, jd2 = kernel_small_vectors(F, A, d, parD, k, nprimes=2, verbose=verbose)
            if good:
                rat = [[Fraction(x) for x in v] for v in good]
                info['lift_ok'] = True; info['lift_method'] = 'lll-lattice'
        if rat is not None:
            # lattice: J_d generator vectors (integral after x7) + cleared kernel vectors; LLL; pick reps independent mod J_d
            idx = {m: i for i, m in enumerate(Dm)}
            jrows = []
            for P in S['jpolys']:
                # exact: m * F_i over Q -> integer after clearing
                pass
            Fi = [pderiv(F, i) for i in range(NV)]
            jexact = [shift(Fi[i], m, 1) for (i, m) in S['jkeys']]
            for P in jexact:
                v = [Fraction(0)] * len(Dm)
                for e, c in P.items():
                    v[idx[e]] = Fraction(c)
                jrows.append(clear_denominators(v))
            krows = [clear_denominators(v) for v in rat]
            red = lll_small_basis(jrows + krows)
            # select vectors independent modulo J_d (mod p rank test)
            p = P0
            Jcols = [[x % p for x in r] for r in jrows]
            base = [r for r in Jcols]
            cur_rank = nmod_mat(len(base), len(Dm), [x for r in base for x in r], p).rank() if base else 0
            chosen = []
            # order candidates by norm (smallest first)
            red_sorted = sorted(red, key=lambda r: sum(x * x for x in r))
            for r in red_sorted:
                trial = base + [[x % p for x in r]]
                rk = nmod_mat(len(trial), len(Dm), [x for rr in trial for x in rr], p).rank()
                if rk > cur_rank:
                    base = trial; cur_rank = rk; chosen.append(r)
                if len(chosen) >= info['nontrivial']:
                    break
            for r in chosen:
                g = 0
                for x in r:
                    g = math.gcd(g, abs(x))
                g = g or 1
                reps.append({Dm[i]: Fraction(x // g) for i, x in enumerate(r) if x})
            info['rep_heights'] = [max(abs(int(c)) for c in P.values()) for P in reps]
            info['rep_terms'] = [len(P) for P in reps]
    info['reps'] = reps
    info['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(f"[sector] d={d} parD={parD} k={k}: nD={info['nD']} J_d={info['J_d']} kernel={info['kernel_total']} "
              f"nontrivial={info['nontrivial']} reps={len(reps)} heights={info.get('rep_heights')} [{info['seconds']}s]", file=sys.stderr)
    return info


# ----------------------------------------------------------------------------- E1 lifts
def e1_columns(ctx, d, pD, k, restrict=True):
    cols = []
    for i in range(NV):
        par = (ctx.pA - pD - PI[i]) % 2
        kN = (ctx.a + k + W[i]) % 7 if (restrict and k is not None) else None
        for m in wmonomials(d + 1, par, kN):
            cols.append((i, m))
    return cols


def rational_e1_lift(F, A, D, nprimes=12, restrict=True, verbose=False):
    """Exact rational E1 lift N (N(F) = A D) via canonical mod-p solves + CRT.  Returns list of 6 polys or None."""
    d = core.degree(D); pD = core.parity_of(D); k = sigma_weight(D) if restrict else None
    ctx0 = Ctx(F, A, P0)
    cols = e1_columns(ctx0, d, pD, k, restrict)
    piv_ref = [None]

    def getter(p):
        ctx = Ctx(F, A, p)
        ok, N, piv = e1_lift(ctx, D, weights=restrict)
        if not ok:
            return None
        if piv_ref[0] is None:
            piv_ref[0] = piv
        elif piv_ref[0] != piv:
            return None
        return [[int(N[i].get(m, 0)) for (i, m) in cols]]
    rat = crt_lift_vectors(getter, nprimes, verbose)
    if rat is None:
        return None
    N = [{} for _ in range(NV)]
    for (i, m), v in zip(cols, rat[0]):
        if v:
            N[i][m] = v
    # exact verification
    Fi = [pderiv(F, i) for i in range(NV)]
    acc = {}
    for i in range(NV):
        if N[i]:
            acc = add(acc, mul(N[i], Fi[i]))
    if sub(acc, mul(A, D)):
        if verbose:
            print('[lift] exact verification FAILED', file=sys.stderr)
        return None
    return N


def float_e1_lift(F, A, D, restrict=True):
    """Floating-point least-squares E1 lift (fallback).  Returns (N list of float polys, relative residual)."""
    d = core.degree(D); pD = core.parity_of(D); k = sigma_weight(D) if restrict else None
    ctx0 = Ctx(F, A, P0)
    cols = e1_columns(ctx0, d, pD, k, restrict)
    Fi = [fl(pderiv(F, i)) for i in range(NV)]
    colpolys = [shift(Fi[i], m, 1.0) for (i, m) in cols]
    rhs = fl(mul(A, D))
    rows = {}
    for c in colpolys:
        for e in c:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    M = np.zeros((len(rows), len(cols)))
    for j, c in enumerate(colpolys):
        for e, v in c.items():
            M[rows[e], j] = v
    b = np.zeros(len(rows))
    for e, v in rhs.items():
        b[rows[e]] = v
    x, *_ = sla.lstsq(M, b, cond=1e-12, check_finite=False)
    rel = np.linalg.norm(M @ x - b) / max(np.linalg.norm(b), 1e-300)
    N = [{} for _ in range(NV)]
    for (i, m), v in zip(cols, x):
        if abs(v) > 1e-14:
            N[i][m] = float(v)
    return N, float(rel)


def jacobian_family(F, A, d, parD, k, restrict=True):
    """Exact family of J(F)_d sector generators (m F_i) with canonical lifts A m e_i."""
    Fi = [pderiv(F, i) for i in range(NV)]
    fam = []
    for i in range(NV):
        if d - 5 < 0:
            break
        pm = (parD + PI[i]) % 2
        km = (k + W[i]) % 7 if (restrict and k is not None) else None
        for m in wmonomials(d - 5, pm, km):
            Dk = shift(Fi[i], m, 1)
            Nk = [{} for _ in range(NV)]
            Nk[i] = shift(A, m, 1)
            fam.append((Dk, Nk, ('J', i, m)))
    return fam


# ----------------------------------------------------------------------------- Koszul variable projection
class KoszulVarPro:
    """Joint search over D(a) = sum a_k D_k (all D_k E1-admissible with lifts N_k) in one sector.

    fam: list of (D_k exact poly, N_k list of 6 polys (exact Fractions or floats), tag).
    restrict: character-restrict the p_ij and Q unknowns (requires all D_k of the same sigma weight k).
    """

    def __init__(self, F, A, fam, restrict=True, k=None, exclude_F_multiples=True, cond=1e-9, verbose=False,
                 normalize=True, reg=0.0, reduce_F=True):
        t0 = time.time()
        self.F, self.A = F, A
        self.cond = cond
        self.reduce_F = reduce_F
        self._cache = None
        self.reg = reg      # >0: Tikhonov (relative lambda) via SVD instead of truncated-SVD lstsq
        self.fam = fam
        self.K = len(fam)
        if self.K == 0:
            raise ValueError('empty family')
        d = core.degree(fam[0][0]); self.d = d
        pD = core.parity_of(fam[0][0]); self.pD = pD
        for Dk, Nk, tag in fam:
            if core.degree(Dk) != d or core.parity_of(Dk) != pD:
                raise ValueError('family members must share degree and parity')
        pA = core.parity_of(A); self.pA = pA
        a = sigma_weight(A)
        if restrict and k is None:
            ks = {sigma_weight(Dk) for Dk, _, _ in fam}
            if len(ks) != 1 or None in ks:
                raise ValueError('restrict=True needs a common sigma weight')
            k = ks.pop()
        self.k = k; self.restrict = restrict
        Ff = fl(F)
        Fi = [fl(pderiv(F, i)) for i in range(NV)]
        # normalization of basis elements (unit coefficient norm) in floating point
        self.scales = []
        Dfl, Nfl = [], []
        for Dk, Nk, tag in fam:
            s = 1.0 / poly_norm(Dk) if normalize else 1.0
            self.scales.append(s)
            Dfl.append({e: float(v) * s for e, v in Dk.items()})
            Nfl.append([{e: float(v) * s for e, v in Ni.items()} for Ni in Nk])
        self.Dfl, self.Nfl = Dfl, Nfl
        # unknown columns
        pcols = []
        if d - 4 >= 0:
            for (i, j) in PAIRS:
                par = (pA + pD + PI[i] + PI[j]) % 2
                kp = (a + k + W[i] + W[j]) % 7 if restrict else None
                for m in wmonomials(d - 4, par, kp):
                    pcols.append((i, j, m))
        qcols = []
        if 2 * d - 6 >= 0:
            kQ = (a + 2 * k) % 7 if restrict else None
            qcols = wmonomials(2 * d - 6, pA, kQ)
        self.pcols, self.qcols = pcols, qcols
        self.nP, self.nQ = len(pcols), len(qcols)
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        # per-family p-blocks
        blocks = []
        for kk in range(self.K):
            D = Dfl[kk]
            Di = [pderiv(D, i) for i in range(NV)]
            DF = [mul(D, Fi[i]) for i in range(NV)]
            G = {}
            for (i, j) in PAIRS:
                G[(i, j)] = sub(mul(Fi[j], Di[i]), mul(Fi[i], Di[j]))
            ri, ci, vals = [], [], []
            for cidx, (i, j, m) in enumerate(pcols):
                col = scale(shift(G[(i, j)], m, 1.0), -1.0)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(DF[j], tuple(m2), float(m[i])))
                if m[j]:
                    m2 = list(m); m2[j] -= 1
                    col = sub(col, shift(DF[i], tuple(m2), float(m[j])))
                for e, v in col.items():
                    if v != 0.0:
                        ri.append(ridx(e)); ci.append(cidx); vals.append(v)
            blocks.append((ri, ci, vals))
        # Q block: -F * m
        qri, qci, qvals = [], [], []
        for cidx, m in enumerate(qcols):
            col = shift(Ff, m, -1.0)
            for e, v in col.items():
                qri.append(ridx(e)); qci.append(cidx); qvals.append(v)
        # rhs pair vectors v_kl = -B_{D_k}(N_l)
        pair = {}
        divN = [None] * self.K
        for l in range(self.K):
            acc = {}
            for i in range(NV):
                if Nfl[l][i]:
                    acc = add(acc, pderiv(Nfl[l][i], i))
            divN[l] = acc
        Dgrad = [[pderiv(Dfl[kk], i) for i in range(NV)] for kk in range(self.K)]
        for kk in range(self.K):
            for l in range(self.K):
                v = mul(Dfl[kk], divN[l])
                for i in range(NV):
                    if Nfl[l][i]:
                        v = sub(v, mul(Nfl[l][i], Dgrad[kk][i]))
                v = scale(v, -1.0)
                pair[(kk, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows)
        self.rows, self.row_index = R, rows
        K = self.K
        # stacked sparse blocks: Mbig = [M_0 | ... | M_{K-1}]  (R x K*nP);  M(a) = Mbig @ kron(a, I_nP)
        Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                            shape=(R, self.nP)) for ri, ci, vals in blocks]
        self.Mbig = sp.hstack(Mk, format='csr') if self.nP else sp.csr_matrix((R, 0))
        self.MQ = (sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                 shape=(R, self.nQ)).toarray() if self.nQ else np.zeros((R, 0)))
        # Bbig (R x K*K): column k*K + l = v_kl = -B_{D_k}(N_l);  b(a) = Bbig @ vec(a a^T)
        ri, ci, vals = [], [], []
        for kk in range(K):
            for l in range(K):
                for e, c in pair[(kk, l)].items():
                    ri.append(rows[e]); ci.append(kk * K + l); vals.append(c)
        self.Bbig = sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                  shape=(R, K * K))
        self.Bbig_full = self.Bbig
        if reduce_F and self.nQ:
            # eliminate Q: reduce every row monomial modulo F (leading term Y0^6); rows' = monomials not divisible by Y0^6
            Rmat, rrows = rem_F_matrix(rows, self.F)
            self.Mbig = (Rmat @ self.Mbig).tocsr()
            self.Bbig = (Rmat @ self.Bbig).tocsr()
            self.MQ = np.zeros((len(rrows), 0))
            self.nQ_eliminated = self.nQ
            self.nQ = 0
            self.rows_full = R
            R = len(rrows); self.rows, self.row_index = R, rrows
        self.nnz = int(self.Mbig.nnz + self.Bbig.nnz)
        self.ncols = self.nP + self.nQ
        self.P = np.eye(self.K); self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = self._complement_of_F_multiples()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[kvarpro] d={d} pD={pD} k={k} K={self.K} rows={R} cols={self.ncols} (p {self.nP}, Q {self.nQ}) '
                  f'nnz={self.nnz} excluded F-multiples dim={self.excluded_dim} dimP={self.dimP} built in {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------------
    def _complement_of_F_multiples(self):
        K = self.K
        dG = self.d - 6
        Gm = monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        if self.restrict:
            Gm = [m for m in Gm if wt(m) == self.k % 7]
        Ff = fl(self.F)
        Fcols = [shift(Ff, m, 1.0) for m in Gm]
        mons = {}
        for B in self.Dfl + Fcols:
            for e in B:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), K + len(Fcols)))
        for kk, B in enumerate(self.Dfl):
            for e, v in B.items():
                Mb[mons[e], kk] = v
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], K + j] = -v
        ns = sla.null_space(Mb)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return U[:, r:], r


    def _solve(self, M, B):
        """Least squares M X = B (B vector or matrix): truncated SVD (cond) or Tikhonov (reg > 0)."""
        if self.reg <= 0:
            X, *_ = sla.lstsq(M, B, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            return X
        U, s, Vt = sla.svd(M, full_matrices=False, check_finite=False, lapack_driver='gesdd')
        lam = self.reg * (s[0] if s.size else 1.0)
        filt = s / (s * s + lam * lam)
        UtB = U.T @ B
        if B.ndim == 1:
            return Vt.T @ (filt * UtB)
        return Vt.T @ (filt[:, None] * UtB)

    def assemble(self, a):
        K, nP = self.K, self.nP
        if nP:
            S = sp.kron(sp.csr_matrix(a.reshape(K, 1)), sp.identity(nP, format='csr'), format='csr')
            M = (self.Mbig @ S).toarray()
            if self.nQ:
                M = np.hstack([M, self.MQ])
        else:
            M = self.MQ
        b = self.Bbig @ np.outer(a, a).ravel()
        return M, b

    def value(self, a):
        a = np.asarray(a, dtype=float)
        na = np.linalg.norm(a)
        a = a / na
        M, b = self.assemble(a)
        x = self._solve(M, b)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r, b

    def value_and_grad(self, u):
        u = np.asarray(u, dtype=float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x = self._solve(M, b)
        r = M @ x - b
        f = float(r @ r)
        K = self.K
        # d/da_k of M(a) x = M_k x_P : columns of Mbig @ kron(I_K, x_P)
        if self.nP:
            xP = x[:self.nP]
            Mx = (self.Mbig @ sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(xP.reshape(self.nP, 1)), format='csr')).toarray()
        else:
            Mx = np.zeros((self.rows, K))
        # d b/da_k = sum_l a_l (v_kl + v_lk) : Bbig @ (e_k (x) a + a (x) e_k)
        Ia = sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(a.reshape(K, 1)), format='csr')
        aI = sp.kron(sp.csr_matrix(a.reshape(K, 1)), sp.identity(K, format='csr'), format='csr')
        Gb = (self.Bbig @ (Ia + aI)).toarray()
        g = 2.0 * (r @ (Mx - Gb))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def rel(self, a):
        """(f / |b_full|^2, f, |b_full|^2): residual relative to the UNreduced right-hand side -B_D(N0) (scale-free,
        comparable across sectors; in the F-reduced system the reduced rhs alone would give rel = 1 whenever p = 0)."""
        f, x, r, b = self.value(a)
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        bf = self.Bbig_full @ np.outer(a, a).ravel()
        bn = float(bf @ bf)
        return f / max(bn, 1e-300), f, bn

    def is_spurious(self, a, tol=1e-8):
        if self.d < 6:
            return False, 1.0
        a = np.asarray(a, float)
        D = {}
        for kk, B in enumerate(self.Dfl):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[kk] * v
        dG = self.d - 6
        Gm = monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Ff = fl(self.F)
        Fcols = [shift(Ff, m, 1.0) for m in Gm]
        mons = {}
        for C in Fcols + [D]:
            for e in C:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), len(Fcols)))
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], j] = v
        bv = np.zeros(len(mons))
        for e, v in D.items():
            bv[mons[e]] = v
        x, *_ = sla.lstsq(Mb, bv, check_finite=False)
        relr = float(np.linalg.norm(Mb @ x - bv) / max(np.linalg.norm(bv), 1e-300))
        return relr < tol, relr

    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, init=None, stop_rel=1e-13,
               callback=None):
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
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            relv, f, bn = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': relv, 'bnorm2': bn, 'a': a, 'spurious': spur, 'div_rel': srel,
                    'iters': int(res.nit), 'restart': it, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[kvarpro] restart {it}: f={f:.3e} rel={relv:.3e} iters={res.nit} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if callback is not None:
                callback(cand)
            if relv < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # ------------------------------------------------------------------ Gauss-Newton / LM (Kaufman varpro Jacobian)
    def _resid_and_jac(self, u, want_jac=True):
        """Residual r(u) = M(a) x - b(a) at a = P u/|P u|, x = lstsq, and the Kaufman Jacobian
        dr/du ~= P_perp (dM/da_k x - db/da_k) * (d a / d u)."""
        u = np.asarray(u, dtype=float)
        key = u.tobytes()
        if self._cache is not None and self._cache[0] == key:
            _, a, na, M, b, x, r = self._cache
        else:
            a_full = self.P @ u
            na = np.linalg.norm(a_full)
            a = a_full / na
            M, b = self.assemble(a)
            x = self._solve(M, b)
            r = M @ x - b
            self._cache = (key, a, na, M, b, x, r)
        self._nfev += 1
        if not want_jac:
            return r, None, a, x, b
        K = self.K
        if self.nP:
            xP = x[:self.nP]
            Mx = (self.Mbig @ sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(xP.reshape(self.nP, 1)), format='csr')).toarray()
        else:
            Mx = np.zeros((self.rows, K))
        Ia = sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(a.reshape(K, 1)), format='csr')
        aI = sp.kron(sp.csr_matrix(a.reshape(K, 1)), sp.identity(K, format='csr'), format='csr')
        Gb = (self.Bbig @ (Ia + aI)).toarray()
        V = Mx - Gb                                   # R x K : d r / d a (ignoring the x-dependence)
        # project onto the orthogonal complement of the column space of M (Kaufman)
        Y = self._solve(M, V)
        V = V - M @ Y
        # chain rule a = a_full/|a_full|, a_full = P u :  da/du = (I - a a^T)/na @ P
        Ja = (V - np.outer(V @ a, a)) / na
        J = Ja @ self.P
        return r, J, a, x, b

    def random_start(self, rng, sparse=0):
        """Random start in reduced coordinates; sparse>0 picks that many random family coordinates (sparse D)."""
        if sparse and sparse < self.K:
            a0 = np.zeros(self.K)
            idx = rng.choice(self.K, size=sparse, replace=False)
            a0[idx] = rng.standard_normal(sparse)
            u0 = self.P.T @ a0
            if np.linalg.norm(u0) < 1e-12:
                u0 = rng.standard_normal(self.dimP)
            return u0
        return rng.standard_normal(self.dimP)

    def search_lm(self, n_restarts=10, seed=0, max_nfev=300, verbose=False, stop_rel=1e-13, callback=None, init=None,
                  polish=True, sparse=0):
        """Random-restart Levenberg-Marquardt / TRF on the sphere with the Kaufman varpro Jacobian."""
        from scipy.optimize import least_squares
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = self.random_start(rng, sparse)
            u0 /= np.linalg.norm(u0)
            fun = lambda u: self._resid_and_jac(u, want_jac=False)[0]
            jac = lambda u: self._resid_and_jac(u, want_jac=True)[1]
            try:
                res = least_squares(fun, u0, jac=jac, method='trf', x_scale='jac', max_nfev=max_nfev,
                                    ftol=1e-15, xtol=1e-15, gtol=1e-15)
                uopt = res.x; nfev = int(res.nfev); status = int(res.status)
            except Exception as ex:
                uopt = u0; nfev = 0; status = -99
            a = self.P @ uopt; a /= np.linalg.norm(a)
            relv, f, bn = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': relv, 'bnorm2': bn, 'a': a, 'spurious': spur, 'div_rel': srel,
                    'iters': nfev, 'status': status, 'restart': it, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[kvarpro-lm] restart {it}: f={f:.3e} rel={relv:.3e} nfev={nfev} status={status} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if callback is not None:
                callback(cand)
            if relv < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out


    def search_joint(self, n_restarts=10, seed=0, max_nfev=300, verbose=False, stop_rel=1e-13, callback=None, init=None):
        """Random-restart LM/TRF over (u, x) jointly (no variable projection; no SVD cutoff)."""
        from scipy.optimize import least_squares
        rng = np.random.default_rng(seed)
        K, nP = self.K, self.nP
        out = []
        t0 = time.time()

        def unpack(z):
            u = z[:self.dimP]; x = z[self.dimP:]
            a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
            return u, x, a, na

        def fun(z):
            u, x, a, na = unpack(z)
            M, b = self.assemble(a)
            return M @ x - b

        def jac(z):
            u, x, a, na = unpack(z)
            M, b = self.assemble(a)
            if nP:
                xP = x[:nP]
                Mx = (self.Mbig @ sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(xP.reshape(nP, 1)), format='csr')).toarray()
            else:
                Mx = np.zeros((self.rows, K))
            Ia = sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(a.reshape(K, 1)), format='csr')
            aI = sp.kron(sp.csr_matrix(a.reshape(K, 1)), sp.identity(K, format='csr'), format='csr')
            Gb = (self.Bbig @ (Ia + aI)).toarray()
            V = Mx - Gb
            Ja = (V - np.outer(V @ a, a)) / na
            return np.hstack([Ja @ self.P, M])

        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            a0 = self.P @ u0; a0 /= np.linalg.norm(a0)
            M, b = self.assemble(a0)
            x0 = self._solve(M, b)
            z0 = np.concatenate([u0, x0])
            try:
                res = least_squares(fun, z0, jac=jac, method='trf', x_scale='jac', max_nfev=max_nfev,
                                    ftol=1e-15, xtol=1e-15, gtol=1e-15)
                z = res.x; nfev = int(res.nfev); status = int(res.status)
            except Exception as ex:
                z = z0; nfev = 0; status = -99
            u, x, a, na = unpack(z)
            relv, f, bn = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': relv, 'bnorm2': bn, 'a': a, 'spurious': spur, 'div_rel': srel,
                    'iters': nfev, 'status': status, 'restart': it, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[kvarpro-joint] restart {it}: f={f:.3e} rel={relv:.3e} nfev={nfev} status={status} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if callback is not None:
                callback(cand)
            if relv < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out


    # ------------------------------------------------------------------ two-stage: relaxation funnel, then full E2
    def ensure_relax(self, verbose=False):
        if getattr(self, '_relax', None) is None:
            self._relax = Relax(self, verbose=verbose)
        return self._relax

    def search_two_stage(self, n_restarts=50, seed=0, max_nfev=300, keep_frac=0.3, keep_rel=1e-6, verbose=False,
                         callback=None, relax_maxiter=600):
        """Stage 1: relaxation minima from n_restarts random starts; stage 2: full-E2 LM from the kept minima
        (all with relaxation rel < keep_rel, else the best keep_frac fraction).  Returns (E2 candidates, relax stats)."""
        rx = self.ensure_relax(verbose=verbose)
        t0 = time.time()
        r1 = rx.search(n_restarts, seed, maxiter=relax_maxiter)
        rels1 = np.array([c['rel'] for c in r1])
        nkeep = max(1, int(math.ceil(keep_frac * len(r1))))
        kept = [c for c in r1 if c['rel'] < keep_rel]
        if len(kept) < nkeep:
            kept = r1[:nkeep]
        out = []
        for j, c0 in enumerate(kept):
            cand = self.search_lm(n_restarts=1, init=c0['a'], max_nfev=max_nfev, seed=seed + j)[0]
            cand['relax_rel'] = c0['rel']; cand['restart'] = j
            cand['seconds'] = round(time.time() - t0, 1)
            out.append(cand)
            if callback is not None:
                callback(cand)
            if verbose:
                print(f'[two-stage] start {j}: relax rel={c0["rel"]:.2e} -> E2 rel={cand["rel"]:.2e} nfev={cand["iters"]}', file=sys.stderr)
            if cand['rel'] < 1e-13 and not cand['spurious']:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        stats = {'relax_restarts': len(r1), 'relax_best_rel': float(rels1.min()), 'relax_median_rel': float(np.median(rels1)),
                 'relax_n_lt_1e-6': int(np.sum(rels1 < 1e-6)), 'relax_n_lt_1e-10': int(np.sum(rels1 < 1e-10)),
                 'relax_nl': rx.nl, 'relax_dimP': self.dimP, 'kept': len(kept), 'seconds': round(time.time() - t0, 1)}
        return out, stats


    def search_clm(self, n_restarts=10, seed=0, max_nfev=300, mu=10.0, verbose=False, stop_rel=1e-13, callback=None, init=None):
        """LM on the combined residual [E2 varpro residual ; mu * relaxation residual L^T b(a)] (same zero set as E2;
        the relaxation part is cheap and shapes the landscape towards the mod-J solution set)."""
        from scipy.optimize import least_squares
        rx = self.ensure_relax()
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        scale_mu = mu * math.sqrt(max(1, self.rows) / max(1, rx.nl))

        def fun(u):
            r, _, a, x, b = self._resid_and_jac(u, want_jac=False)
            qv, _ = rx.q_and_grad(a)
            return np.concatenate([r, scale_mu * qv])

        def jac(u):
            r, J, a, x, b = self._resid_and_jac(u, want_jac=True)
            qv, dq = rx.q_and_grad(a)
            na = np.linalg.norm(self.P @ u)
            dqa = (dq - np.outer(dq @ a, a)) / na
            return np.vstack([J, scale_mu * (dqa @ self.P)])

        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            try:
                res = least_squares(fun, u0, jac=jac, method='trf', x_scale='jac', max_nfev=max_nfev,
                                    ftol=1e-15, xtol=1e-15, gtol=1e-15)
                uopt = res.x; nfev = int(res.nfev); status = int(res.status)
            except Exception:
                uopt = u0; nfev = 0; status = -99
            a = self.P @ uopt; a /= np.linalg.norm(a)
            relv, f, bn = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': relv, 'bnorm2': bn, 'a': a, 'spurious': spur, 'div_rel': srel, 'relax_rel': rx.rel(a)[0],
                    'iters': nfev, 'status': status, 'restart': it, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[kvarpro-clm] restart {it}: rel={relv:.3e} relax={cand["relax_rel"]:.2e} nfev={nfev} status={status} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if callback is not None:
                callback(cand)
            if relv < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    def solution_set_dim(self, a, tol=1e-7):
        """Numerical dimension of the local solution set at a (null space of the Kaufman Jacobian on the sphere)."""
        u = self.P.T @ np.asarray(a, float)
        r, J, a2, x, b = self._resid_and_jac(u)
        s = np.linalg.svd(J, compute_uv=False)
        if s.size == 0:
            return 0, s
        return int(np.sum(s < tol * s[0])), s

    # ------------------------------------------------------------------ exact side
    def D_float(self, a):
        a = np.asarray(a, float)
        D = {}
        for kk, B in enumerate(self.Dfl):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[kk] * v
        return {e: v for e, v in D.items() if v != 0.0}

    def reconstruct_D(self, a, max_den=1000, tol=1e-7):
        """Rational candidates for D(a) modulo F, via the canonical t-monic remainder chart rem_F."""
        Dr = rem_F_float(self.D_float(a), self.F)
        return reconstruct_poly(Dr, max_den=max_den, tol=tol)

    def coeffs_from_a(self, a):
        """Coefficients c_k on the *exact* basis D_k (undo normalization)."""
        a = np.asarray(a, float)
        return a * np.array(self.scales)

    def reconstruct(self, a, max_den=1000, tol=1e-7):
        c = self.coeffs_from_a(a)
        j = int(np.argmax(np.abs(c)))
        c = c / c[j]
        cands = []
        fr = [Fraction(float(v)).limit_denominator(max_den) for v in c]
        if max(abs(float(f) - v) for f, v in zip(fr, c)) < tol:
            cands.append(tuple(fr))
        for den in range(1, max_den + 1):
            v = c * den
            if np.max(np.abs(v - np.round(v))) < tol * den:
                cands.append(tuple(Fraction(int(round(x)), den) for x in v))
                break
        uniq = []
        for cc in cands:
            if cc not in uniq:
                uniq.append(cc)
        return uniq

    def D_exact(self, frac):
        D = {}
        for kk, (Dk, Nk, tag) in enumerate(self.fam):
            if frac[kk] == 0:
                continue
            D = add(D, scale(Dk, Fraction(frac[kk])))
        return D


def confirm_exact(F, A, D, verbose=False, restrict=True, with_engine=False):
    """Exact decision for a rational D: kcc2 Koszul test mod p (character restricted) then multi-prime lift.
    Returns (status, data) with status in {'inconsistent', 'solved', 'lift-failed', 'e1-fail'}."""
    t0 = time.time()
    if not D:
        return 'inconsistent', {'note': 'D = 0'}
    e = kcc2.full_test(F, A, D, 'w2-confirm', 'w2-confirm', weights=restrict, verbose=verbose, want_log=False)
    if not e.get('e1'):
        return 'e1-fail', {'e1': False, 'seconds_modp': round(time.time() - t0, 1)}
    data = {'rank': e.get('rank'), 'nullity': e.get('nullity'), 'rows': e.get('rows'), 'cols': e.get('cols'),
            'seconds_modp': round(time.time() - t0, 1)}
    if not e.get('consistent'):
        return 'inconsistent', data
    lift = kcc2.lift_certificate(F, A, D, weights=restrict, nprimes=12, verbose=verbose)
    if lift is None:
        return 'lift-failed', data
    data.update({'N': lift['N'], 'Q': lift['Q'], 'verify': lift['ok'], 'info': lift['info']})
    if lift['ok']:
        if with_engine:
            st, dd = fixed_d.solve_fixed_D(F, A, 1, D, verbose=verbose)
            data['engine_status'] = st
        return 'solved', data
    return 'lift-failed', data


# ----------------------------------------------------------------------------- Hessian (mod J) relaxation
class HessianRelaxation:
    """For a Jacobian-type family D = sum_i w_i F_i (w from jacobian_family generators), the necessary
    condition  [A * (w^T Hess(F) w)] = 0 in (R/J)_{2d}  (sector a+2k, parity pA).  Floating point:
    an orthonormal basis L of the left null space of the J_{2d} sector block gives f(w) = |L^T vec(A w^T H w)|^2.
    More generally for any family with lifts: f(a) = |L^T vec(B_{D(a)}(N0(a)))|^2 (same thing mod J)."""

    def __init__(self, F, A, fam, k, restrict=True, verbose=False):
        t0 = time.time()
        self.F, self.A, self.fam, self.k = F, A, fam, k
        d = core.degree(fam[0][0]); self.d = d
        pD = core.parity_of(fam[0][0]); pA = core.parity_of(A); a = sigma_weight(A)
        self.K = len(fam)
        Ff = fl(F)
        Fi = [fl(pderiv(F, i)) for i in range(NV)]
        tpar = pA
        tk = (a + 2 * k) % 7 if restrict else None
        tgt = wmonomials(2 * d, tpar, tk)
        rows = {m: i for i, m in enumerate(tgt)}
        # J block (floating) in degree 2d sector
        jcols = []
        for i in range(NV):
            pm = (tpar + PI[i]) % 2
            km = (tk + W[i]) % 7 if tk is not None else None
            for m in wmonomials(2 * d - 5, pm, km):
                jcols.append(shift(Fi[i], m, 1.0))
        J = np.zeros((len(rows), len(jcols)))
        for j, c in enumerate(jcols):
            for e, v in c.items():
                J[rows[e], j] = v
        # left null space (functionals vanishing on J): orthonormal basis of complement of col space
        U, s, Vt = np.linalg.svd(J, full_matrices=True)
        tol = max(J.shape) * np.finfo(float).eps * (s[0] if s.size else 1.0) * 1e3
        rank = int(np.sum(s > tol))
        self.L = U[:, rank:]           # rows x nl
        self.nl = self.L.shape[1]
        self.ring_dim_float = self.nl
        self.jrank = rank
        self.sv_gap = (float(s[rank - 1]), float(s[rank]) if rank < len(s) else 0.0) if rank > 0 else None
        # quadratic map: V[k][l] = L^T vec(-B_{D_k}(N_l)) ; f(a) = | sum a_k a_l V_kl |^2
        scales = [1.0 / poly_norm(Dk) for Dk, _, _ in fam]
        Dfl = [{e: float(v) * s_ for e, v in Dk.items()} for (Dk, _, _), s_ in zip(fam, scales)]
        Nfl = [[{e: float(v) * s_ for e, v in Ni.items()} for Ni in Nk] for (_, Nk, _), s_ in zip(fam, scales)]
        self.scales = scales
        divN = []
        for l in range(self.K):
            acc = {}
            for i in range(NV):
                if Nfl[l][i]:
                    acc = add(acc, pderiv(Nfl[l][i], i))
            divN.append(acc)
        Dgrad = [[pderiv(Dfl[kk], i) for i in range(NV)] for kk in range(self.K)]
        V = np.zeros((self.K, self.K, self.nl))
        for kk in range(self.K):
            for l in range(self.K):
                v = mul(Dfl[kk], divN[l])
                for i in range(NV):
                    if Nfl[l][i]:
                        v = sub(v, mul(Nfl[l][i], Dgrad[kk][i]))
                vec = np.zeros(len(rows))
                for e, c in v.items():
                    vec[rows[e]] = c
                V[kk, l, :] = self.L.T @ vec
        self.V = V
        self.Vs = V + V.transpose(1, 0, 2)   # symmetric part (for gradient)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[hess] d={d} k={k} K={self.K} rows={len(rows)} Jrank={rank} nl={self.nl} sv_gap={self.sv_gap} '
                  f'built in {self.build_time:.1f}s', file=sys.stderr)

    def value_and_grad(self, u):
        na = np.linalg.norm(u); a = u / na
        q = np.einsum('k,l,klm->m', a, a, self.V)
        f = float(q @ q)
        # d q_m / d a_k = sum_l a_l (V_klm + V_lkm)
        dq = np.einsum('l,klm->km', a, self.Vs)
        g = 2.0 * (dq @ q)
        g_full = (g - float(g @ a) * a) / na
        return f, g_full

    def search(self, n_restarts=1000, seed=0, maxiter=500, verbose=False):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            u0 = rng.standard_normal(self.K); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': 1e-14, 'ftol': 1e-20})
            a = res.x / np.linalg.norm(res.x)
            f, _ = self.value_and_grad(a)
            # relative scale: |V| applied to a (norm of the quadratic image before projection is unknown here);
            # use |a|=1 so f is already scale-free; also record the norm of the unprojected vector? keep f.
            out.append({'f': f, 'a': a, 'iters': int(res.nit)})
        out.sort(key=lambda c: c['f'])
        if verbose:
            print(f'[hess] {n_restarts} restarts best f={out[0]["f"]:.3e} ({time.time()-t0:.1f}s)', file=sys.stderr)
        return out


# ----------------------------------------------------------------------------- canonical chart modulo F
def rem_F_float(D, F):
    """Remainder of a floating-point poly D modulo F with respect to the leading monomial Y0^6 of F
    (F = Y0^6 + lower in Y0-degree).  Canonical representative of D mod F (same support on X)."""
    Ff = fl(F)
    t6 = (0, 0, 0, 6, 0, 0)
    assert abs(Ff.get(t6, 0.0) - 1.0) < 1e-12
    rest = {e: v for e, v in Ff.items() if e != t6}
    D = dict(D)
    changed = True
    while changed:
        changed = False
        for e in sorted(list(D.keys()), key=lambda e: -e[3]):
            if e[3] >= 6 and D.get(e, 0.0) != 0.0:
                c = D.pop(e)
                m = list(e); m[3] -= 6; m = tuple(m)
                for e2, v2 in rest.items():
                    e3 = tuple(x + y for x, y in zip(m, e2))
                    D[e3] = D.get(e3, 0.0) - c * v2
                changed = True
                break
    return {e: v for e, v in D.items() if abs(v) > 0.0}


def rem_F_exact(D, F):
    t6 = (0, 0, 0, 6, 0, 0)
    assert F.get(t6, 0) == 1
    rest = {e: v for e, v in F.items() if e != t6}
    D = dict(D)
    while True:
        big = [e for e in D if e[3] >= 6 and D[e] != 0]
        if not big:
            break
        e = max(big, key=lambda e: e[3])
        c = D.pop(e)
        m = list(e); m[3] -= 6; m = tuple(m)
        for e2, v2 in rest.items():
            e3 = tuple(x + y for x, y in zip(m, e2))
            D[e3] = D.get(e3, 0) - c * v2
    return {e: v for e, v in D.items() if v != 0}


def reconstruct_poly(Dr, max_den=1000, tol=1e-7, drop=1e-9):
    """Rational candidates (list of exact polys) for a floating poly, normalized so the largest |coeff| = 1."""
    if not Dr:
        return []
    items = sorted(Dr.items())
    vals = np.array([v for _, v in items])
    j = int(np.argmax(np.abs(vals)))
    vals = vals / vals[j]
    keep = np.abs(vals) > drop
    cands = []
    fr = [Fraction(float(v)).limit_denominator(max_den) if kp else Fraction(0) for v, kp in zip(vals, keep)]
    if max(abs(float(f) - v) for f, v in zip(fr, vals)) < tol:
        cands.append(fr)
    for den in range(1, max_den + 1):
        v = vals * den
        if np.max(np.abs(v - np.round(v))) < tol * den:
            cands.append([Fraction(int(round(x)), den) if kp else Fraction(0) for x, kp in zip(v, keep)])
            break
    out = []
    for c in cands:
        P = {e: f for (e, _), f in zip(items, c) if f != 0}
        if P and P not in out:
            out.append(P)
    return out


def rem_F_matrix(rows, F):
    """Sparse matrix R (reduced rows x rows) of the linear map 'remainder modulo F' (leading term Y0^6) on the
    span of the given row monomials; returns (R, reduced_row_index)."""
    t6 = (0, 0, 0, 6, 0, 0)
    assert F.get(t6, 0) == 1
    rest = {e: float(v) for e, v in F.items() if e != t6}
    memo = {}

    def red(e):
        if e[3] < 6:
            return {e: 1.0}
        if e in memo:
            return memo[e]
        m = list(e); m[3] -= 6; m = tuple(m)
        acc = {}
        for e2, v2 in rest.items():
            e3 = tuple(x + y for x, y in zip(m, e2))
            for e4, v4 in red(e3).items():
                acc[e4] = acc.get(e4, 0.0) - v2 * v4
        acc = {k: v for k, v in acc.items() if v != 0.0}
        memo[e] = acc
        return acc
    rrows = {}
    ri, ci, vals = [], [], []
    for e, j in rows.items():
        for e4, v in red(e).items():
            if e4 not in rrows:
                rrows[e4] = len(rrows)
            ri.append(rrows[e4]); ci.append(j); vals.append(v)
    R = sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                      shape=(len(rrows), len(rows)))
    return R, rrows


# ----------------------------------------------------------------------------- mod-J (Hessian) relaxation on a family
def jblock_float_sector(F, d, k, restrict=True):
    """Float J-block in degree 2d (Y-odd, weight 3+2k when restrict) over all sector monomials; returns (J, row index)."""
    a = 3; pA = 1
    tk = (a + 2 * k) % 7 if (restrict and k is not None) else None
    Fi = [fl(pderiv(F, i)) for i in range(NV)]
    tgt = wmonomials(2 * d, pA, tk)
    rows = {m: i for i, m in enumerate(tgt)}
    cols = []
    for i in range(NV):
        pm = (pA + PI[i]) % 2
        km = (tk + W[i]) % 7 if tk is not None else None
        for m in wmonomials(2 * d - 5, pm, km):
            cols.append(shift(Fi[i], m, 1.0))
    J = np.zeros((len(rows), len(cols)))
    for j, c in enumerate(cols):
        for e, v in c.items():
            J[rows[e], j] = v
    return J, rows


class Relax:
    """f(a) = |L^T b(a)|^2 where b(a) = -B_{D(a)}(N0(a)) (from a KoszulVarPro, reduced or not) and L = orthonormal
    functionals on the degree-2d sector vanishing on J(F).  Zero iff B_{D(a)}(N0(a)) in J, a necessary condition for E2."""

    def __init__(self, kv, verbose=False):
        t0 = time.time()
        self.kv = kv; self.K = kv.K; self.P = kv.P; self.dimP = kv.dimP
        J, frows = jblock_float_sector(kv.F, kv.d, kv.k, kv.restrict)
        U, s, Vt = np.linalg.svd(J, full_matrices=True)
        tol = max(J.shape) * np.finfo(float).eps * s[0] * 1e3
        rank = int(np.sum(s > tol))
        Lfull = U[:, rank:]
        self.nl = Lfull.shape[1]; self.jrank = rank
        self.sv_gap = (float(s[rank - 1]), float(s[rank]) if rank < len(s) else 0.0)
        sub_idx = np.array([frows[e] for e, j in sorted(kv.row_index.items(), key=lambda kv_: kv_[1])])
        self.L = Lfull[sub_idx, :]
        self.LtB = None
        if kv.Bbig.shape[1] * self.nl <= 3e7:
            self.LtB = (self.L.T @ kv.Bbig.toarray())          # nl x K^2
        self.build_time = time.time() - t0
        if verbose:
            print(f'[relax] d={kv.d} k={kv.k} K={self.K} nl={self.nl} jrank={rank} sv_gap={self.sv_gap} built in {self.build_time:.1f}s', file=sys.stderr)

    def q_and_grad(self, a):
        K = self.K
        v = np.outer(a, a).ravel()
        if self.LtB is not None:
            qv = self.LtB @ v
            G = self.LtB.reshape(-1, K, K)
            dq = np.einsum('mkl,l->mk', G, a) + np.einsum('mkl,k->ml', G, a)
        else:
            b = self.kv.Bbig @ v
            qv = self.L.T @ b
            Ia = sp.kron(sp.identity(K, format='csr'), sp.csr_matrix(a.reshape(K, 1)), format='csr')
            aI = sp.kron(sp.csr_matrix(a.reshape(K, 1)), sp.identity(K, format='csr'), format='csr')
            dq = self.L.T @ (self.kv.Bbig @ (Ia + aI)).toarray()
        return qv, dq

    def value_and_grad(self, u):
        a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
        qv, dq = self.q_and_grad(a)
        f = float(qv @ qv)
        g = 2.0 * (dq.T @ qv)
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def rel(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        b = self.kv.Bbig @ np.outer(a, a).ravel()
        qv = self.L.T @ b
        return float(qv @ qv) / max(float(b @ b), 1e-300), float(qv @ qv), float(b @ b)

    def search(self, n_restarts, seed, maxiter=600, init=None):
        rng = np.random.default_rng(seed)
        out = []
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': 1e-14, 'ftol': 1e-20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            relv, f, bn = self.rel(a)
            out.append({'f': f, 'rel': relv, 'bn': bn, 'a': a, 'iters': int(res.nit)})
        out.sort(key=lambda c: c['rel'])
        return out


def kernel_small_vectors(F, A, d, parD, k, nprimes=2, verbose=False, max_height=10**12):
    """Short integer vectors of the E1 kernel {D in S_d[parD,k] : A D in J} via LLL on the lattice
    {v in Z^nD : v mod m in ker_m}, m = product of nprimes primes (CRT of the canonical mod-p bases).
    Returns (Dm, list of integer vectors verified in the kernel mod an extra prime, nk, jd)."""
    ms = None; base = None; nk = jd = 0; Dm = None
    for p in core.PRIMES[:nprimes]:
        ctx = Ctx(F, A, p)
        S = sector_kernel_full(ctx, d, k, parD)
        Dm = S['Dm']; nk = len(S['kervecs']); jd = S['jd']
        if nk == 0:
            return Dm, [], 0, jd
        if base is None:
            base = [list(v) for v in S['kervecs']]; ms = p
        else:
            if len(S['kervecs']) != len(base):
                continue
            base = [[core.crt_pair(a, ms, b, p)[0] for a, b in zip(v1, v2)] for v1, v2 in zip(base, S['kervecs'])]
            ms *= p
    nD = len(Dm)
    rows = [[int(x) % ms for x in v] for v in base] + [[ms if i == j else 0 for i in range(nD)] for j in range(nD)]
    M = fmpz_mat(len(rows), nD, [x for r in rows for x in r])
    R = M.lll()
    cands = []
    for i in range(R.nrows()):
        v = [int(R[i, j]) for j in range(nD)]
        if any(v) and max(abs(x) for x in v) < max_height:
            cands.append(v)
    # verify mod an extra prime (not used in the lattice): functionals vanishing on J applied to A*D
    cands.sort(key=lambda v: sum(x * x for x in v))
    cands = cands[: nk + 10]
    pv = core.PRIMES[nprimes]
    ctx = Ctx(F, A, pv)
    tpar = (parD + ctx.pA) % 2; tk = (k + ctx.a) % 7
    tgt = wmonomials(d + 6, tpar, tk)
    rws = {m: i for i, m in enumerate(tgt)}
    jk, jp = kcc2.jac_block_cols(ctx, d + 6, tpar, tk)
    MJ, _ = dense_from_cols(jp, None, pv, rws)
    Lm, nl = MJ.transpose().nullspace()
    Lt = nmod_mat(nl, len(rws), [int(Lm[r, c]) for c in range(nl) for r in range(len(rws))], pv)
    good = []
    for v in cands:
        D = {m: c % pv for m, c in zip(Dm, v) if c % pv}
        AD = mul(ctx.A, D, pv)
        vec = [0] * len(rws)
        for e, c in AD.items():
            vec[rws[e]] = c
        if all(int(x) == 0 for x in (Lt * nmod_mat(len(rws), 1, vec, pv)).entries()):
            good.append(v)
    good.sort(key=lambda v: sum(x * x for x in v))
    if verbose:
        print(f'[kernel-lll] d={d} parD={parD} k={k}: nk={nk} jd={jd} lattice candidates {len(cands)} verified {len(good)} '
              f'heights {[max(abs(x) for x in v) for v in good[:8]]}', file=sys.stderr)
    return Dm, good, nk, jd


def float_kernel_reps(F, A, d, parD, k, tol=1e-9, verbose=False):
    """Floating-point basis of the E1 kernel {D in S_d[parD,k] : A D in J} orthogonal to J_d (nontrivial directions),
    for sectors whose rational representatives have astronomically large height.  Returns (reps as float polys,
    lifts as float polys, info).  Uses the functionals L (left null space of the degree-(d+6) J-block) in float."""
    t0 = time.time()
    a = 3; pA = 1
    Dm = wmonomials(d, parD, k)
    if not Dm:
        return [], [], {'nD': 0}
    tpar = (parD + pA) % 2; tk = (k + a) % 7
    tgt = wmonomials(d + 6, tpar, tk)
    rows = {m: i for i, m in enumerate(tgt)}
    Fi = [fl(pderiv(F, i)) for i in range(NV)]
    jcols = []
    for i in range(NV):
        pm = (tpar + PI[i]) % 2; km = (tk + W[i]) % 7
        for m in wmonomials(d + 1, pm, km):
            jcols.append(shift(Fi[i], m, 1.0))
    J = np.zeros((len(rows), len(jcols)))
    for j, c in enumerate(jcols):
        for e, v in c.items():
            J[rows[e], j] = v
    U, s, Vt = np.linalg.svd(J, full_matrices=True)
    cut = max(J.shape) * np.finfo(float).eps * s[0] * 1e3
    rank = int(np.sum(s > cut))
    Lf = U[:, rank:]                       # functionals vanishing on J_{d+6}
    Af = fl(A)
    AD = np.zeros((len(rows), len(Dm)))
    for j, m in enumerate(Dm):
        for e, v in shift(Af, m, 1.0).items():
            AD[rows[e], j] = v
    Mult = Lf.T @ AD                        # nl x nD
    U2, s2, Vt2 = np.linalg.svd(Mult, full_matrices=True)
    cut2 = max(Mult.shape) * np.finfo(float).eps * (s2[0] if s2.size else 1.0) * 1e3
    r2 = int(np.sum(s2 > cut2))
    Kb = Vt2[r2:, :].T                      # nD x nk orthonormal kernel basis
    nk = Kb.shape[1]
    # J_d part in D coordinates
    idx = {m: i for i, m in enumerate(Dm)}
    jd_cols = []
    for i in range(NV):
        pm = (parD + PI[i]) % 2; km = (k + W[i]) % 7
        for m in wmonomials(d - 5, pm, km):
            jd_cols.append(shift(Fi[i], m, 1.0))
    Jd = np.zeros((len(Dm), len(jd_cols)))
    for j, c in enumerate(jd_cols):
        for e, v in c.items():
            Jd[idx[e], j] = v
    if jd_cols:
        Q, Rr = np.linalg.qr(Jd)
        sj = np.linalg.svd(Jd, compute_uv=False)
        rj = int(np.sum(sj > max(Jd.shape) * np.finfo(float).eps * sj[0] * 1e3))
        Uj, _, _ = np.linalg.svd(Jd, full_matrices=False)
        Uj = Uj[:, :rj]
        Kperp = Kb - Uj @ (Uj.T @ Kb)       # component orthogonal to J_d
    else:
        rj = 0; Kperp = Kb
    U3, s3, Vt3 = np.linalg.svd(Kperp, full_matrices=False)
    cut3 = (s3[0] if s3.size else 1.0) * 1e-8
    r3 = int(np.sum(s3 > cut3))
    reps_vec = U3[:, :r3]                   # nD x nontrivial, orthonormal, orthogonal to J_d
    reps = []
    for j in range(r3):
        v = reps_vec[:, j]
        reps.append({Dm[i]: float(v[i]) for i in range(len(Dm)) if abs(v[i]) > 1e-15})
    lifts = []; lift_res = []
    for D in reps:
        N, relr = float_e1_lift(F, A, D, restrict=True)
        lifts.append(N); lift_res.append(relr)
    info = {'nD': len(Dm), 'nl': Lf.shape[1], 'kernel_total': nk, 'J_d': rj, 'nontrivial': r3, 'lift_res': lift_res,
            'sv_kernel_gap': [float(s2[r2 - 1]) if r2 > 0 else None, float(s2[r2]) if r2 < len(s2) else None],
            'seconds': round(time.time() - t0, 1)}
    if verbose:
        print(f'[float-kernel] d={d} parD={parD} k={k}: nD={len(Dm)} nl={Lf.shape[1]} kernel={nk} J_d={rj} nontrivial={r3} '
              f'lift residuals {lift_res} [{info["seconds"]}s]', file=sys.stderr)
    return reps, lifts, info
