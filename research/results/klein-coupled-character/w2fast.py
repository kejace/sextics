"""Fast sector-restricted Koszul variable-projection engine for klein-coupled-character (q = 1).

Complements w2lib.py (exact sector kernels / rational reps / exact confirmation).  Same mathematics
(correction-map form, ch. 5.2):

  sigma weights W = (1,2,4,0,0,0) mod 7; F character 0; A = X0^3 Y0^3 character aA = 3, Y-odd.
  Sector (d, parD, chi):  D in S_d[parD, chi].  Unknown sectors (averaging, ch. 1):
     N_i  deg d+1,  parity pA+parD+PI_i,        char aA+chi+W_i
     p_ij deg d-4,  parity pA+parD+PI_i+PI_j,   char aA+chi+W_i+W_j
     Q    deg 2d-6, parity pA,                  char aA+2chi
     rows deg 2d,   parity pA,                  char aA+2chi
  For monomials D = z^mu, N_i = z^m:   B_D(N) = D div N - N(D) = (m_i - mu_i) z^{mu+m-e_i}
     L_ij(m) = (m_i - mu_i) z^{mu+m-e_i} F_j - (m_j - mu_j) z^{mu+m-e_j} F_i.
  System:  T(c) p - F Q = b(c, n),  b = -B_D(N0) = sum (mu_i - m_i) c_mu n_(i,m) z^{mu+m-e_i}.
  T(c) = sum_mu c_mu T_mu is stored as one sparse matrix Tall ((rows*nP) x nD); b as Vb (rows x nD*nN).
  The constant F*Q block is projected out once (QR), so each evaluation is one dense lstsq of size
  rows x nP (gelsd, fixed relative cutoff).

Families (Bc: nD x K orthonormal D-coefficient basis, Bn: nN x K lifts with E1 exact):
  'full': float null space of the sector E1 system  (= (J:A)_{d,sector}, lifts for free)
  'w'   : J_d sector (D = sum w_i F_i, lifts A w, minimal-norm w)
  Hessian relaxation: family 'w' with use_koszul=False (N = A w exactly, Q free).
Spurious F-multiples (D in F*S_{d-6}) are projected out of the parameter space.

Validation: validate() rebuilds the functional for accepted q=1 witnesses of other instances
(unrestricted sectors) and checks f ~ 0 at the witness denominator.
"""
from __future__ import annotations
import sys, time, json, os, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d
from hodge.core import NV, monomials, pderiv, mul, add, sub, scale, shift

LANE = 'klein-coupled-character'
LANE_DIR = '/tmp/hodge/results/klein-coupled-character'
BATCH_LOG = os.path.join(LANE_DIR, 'wave2-search.log.jsonl')
W = (1, 2, 4, 0, 0, 0)
PI = (0, 0, 0, 1, 1, 1)
PAIRS = [(i, j) for i in range(NV) for j in range(i + 1, NV)]


def wt(e):
    return (e[0] + 2 * e[1] + 4 * e[2]) % 7


def sector_monomials(deg, par, chi):
    if deg < 0:
        return []
    if chi is None:
        return list(monomials(deg, par))
    return [e for e in monomials(deg, par) if wt(e) == chi % 7]


def enc(E):
    E = np.asarray(E, dtype=np.int64).reshape(-1, 6)
    return ((((E[:, 0] * 64 + E[:, 1]) * 64 + E[:, 2]) * 64 + E[:, 3]) * 64 + E[:, 4]) * 64 + E[:, 5]


class Index:
    def __init__(self, mons):
        self.mons = list(mons)
        self.n = len(self.mons)
        if self.n:
            keys = enc(np.array(self.mons, dtype=np.int64).reshape(-1, 6))
        else:
            keys = np.zeros(0, np.int64)
        order = np.argsort(keys)
        self.keys = keys[order]
        self.order = order

    def lookup(self, E, strict=True):
        k = enc(E)
        if self.n == 0:
            if strict and len(k):
                raise KeyError('empty index')
            return np.zeros(0, np.int64), np.zeros(len(k), bool)
        pos = np.searchsorted(self.keys, k)
        pos = np.minimum(pos, self.n - 1)
        ok = self.keys[pos] == k
        if strict and not np.all(ok):
            raise KeyError(f'{int((~ok).sum())} monomials not in index')
        return self.order[pos], ok


def poly_to_vec(P, index):
    v = np.zeros(index.n)
    if P:
        E = np.array(list(P.keys()), dtype=np.int64)
        idx, _ = index.lookup(E)
        v[idx] = [float(x) for x in P.values()]
    return v


def batch_log(entry, path=BATCH_LOG):
    entry.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('engine', 'w2fast-koszul-varpro')
    entry.setdefault('instance', LANE)
    with open(path, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')


def load():
    F, A, q, rec = core.load_instance(LANE)
    return F, A, q


def char_of(P):
    ws = {wt(e) for e in P}
    return ws.pop() if len(ws) == 1 else None


# ==============================================================================================
class Sector:
    """Monomial bookkeeping + matrices for one (d, parD, chi) sector.  chi=None: no character restriction."""

    def __init__(self, F, A, d, parD, chi, verbose=False):
        t0 = time.time()
        self.F, self.A, self.d, self.parD, self.chi = F, A, d, parD, chi
        self.pA = core.parity_of(A)
        self.aA = char_of(A) if chi is not None else None
        pA, aA = self.pA, self.aA
        if chi is not None and (char_of(F) != 0 or aA is None):
            raise ValueError('character restriction needs char(F)=0 and A an eigenvector')
        self.Fi = [pderiv(F, i) for i in range(NV)]
        cc = (lambda x: (x % 7)) if chi is not None else (lambda x: None)
        self.Dm = sector_monomials(d, parD, chi)
        self.nD = len(self.Dm)
        self.Ncols, self.Nblocks = [], []
        for i in range(NV):
            ms = sector_monomials(d + 1, (pA + parD + PI[i]) % 2, cc(aA + chi + W[i]) if chi is not None else None)
            self.Nblocks.append((len(self.Ncols), len(ms)))
            self.Ncols += [(i, m) for m in ms]
        self.nN = len(self.Ncols)
        self.Pcols, self.Pblocks = [], {}
        if d >= 4:
            for (i, j) in PAIRS:
                ms = sector_monomials(d - 4, (pA + parD + PI[i] + PI[j]) % 2,
                                      cc(aA + chi + W[i] + W[j]) if chi is not None else None)
                self.Pblocks[(i, j)] = (len(self.Pcols), len(ms))
                self.Pcols += [(i, j, m) for m in ms]
        self.nP = len(self.Pcols)
        qchi = cc(aA + 2 * chi) if chi is not None else None
        self.Qm = sector_monomials(2 * d - 6, pA, qchi) if 2 * d - 6 >= 0 else []
        self.nQ = len(self.Qm)
        self.rows = sector_monomials(2 * d, pA, qchi)
        self.ridx = Index(self.rows)
        self.nrows = len(self.rows)
        self.Didx = Index(self.Dm)
        self.verbose = verbose
        self._built = False
        self.t_init = time.time() - t0
        if verbose:
            print(f'[sector d={d} par={parD} chi={chi}] nD={self.nD} nN={self.nN} nP={self.nP} nQ={self.nQ} '
                  f'rows={self.nrows} ({self.t_init:.1f}s)', file=sys.stderr)

    def dims(self):
        return {'d': self.d, 'parD': self.parD, 'chi': self.chi, 'nD': self.nD, 'nN': self.nN, 'nP': self.nP,
                'nQ': self.nQ, 'rows': self.nrows}

    # ------------------------------------------------------------------ E1 space
    def e1_matrix(self):
        d, parD, chi = self.d, self.parD, self.chi
        tchi = (self.aA + chi) % 7 if chi is not None else None
        rows6 = sector_monomials(d + 6, (parD + self.pA) % 2, tchi)
        idx6 = Index(rows6)
        E = np.zeros((len(rows6), self.nD + self.nN))
        for k, mu in enumerate(self.Dm):
            E[:, k] = poly_to_vec(shift(self.A, mu, 1), idx6)
        for k, (i, m) in enumerate(self.Ncols):
            E[:, self.nD + k] = -poly_to_vec(shift(self.Fi[i], m, 1), idx6)
        return E, idx6

    def e1_space(self, rcond=1e-10, verbose=None):
        """Orthonormal basis Bc (nD x K) of (J:A)_{d,sector} and lifts Bn (nN x K) with A D = sum F_i N_i."""
        verbose = self.verbose if verbose is None else verbose
        t0 = time.time()
        if self.nD == 0:
            return np.zeros((0, 0)), np.zeros((self.nN, 0)), {'K': 0}
        E, idx6 = self.e1_matrix()
        ns = sla.null_space(E, rcond=rcond)
        ZD, ZN = ns[:self.nD], ns[self.nD:]
        if ZD.size == 0:
            return np.zeros((self.nD, 0)), np.zeros((self.nN, 0)), {'K': 0, 'nullity_e1': int(ns.shape[1])}
        U, s, Vt = np.linalg.svd(ZD, full_matrices=False)
        r = int(np.sum(s > 1e-7))   # null vectors are unit norm: genuine D-directions have s = O(1)
        Bc = U[:, :r]
        Bn = (ZN @ Vt[:r].T) / s[:r]
        chk = float(np.abs(E @ np.vstack([Bc, Bn])).max()) if r else 0.0
        info = {'K': r, 'nullity_e1': int(ns.shape[1]), 'syz': int(ns.shape[1] - r), 'e1_check': chk,
                'e1_rows': int(E.shape[0]), 'e1_cols': int(E.shape[1]), 'seconds': round(time.time() - t0, 1),
                'sv_gap': [float(s[r - 1]), float(s[r]) if r < len(s) else 0.0] if r else None}
        if verbose:
            print(f'[e1 d={self.d} par={self.parD} chi={self.chi}] dim(J:A)={r} syz={info["syz"]} check={chk:.1e} '
                  f'({info["seconds"]}s)', file=sys.stderr)
        return Bc, Bn, info

    def w_generators(self):
        """Columns m F_i of J_d in the sector (D-coeff matrix nD x g) with lifts A m e_i (nN x g)."""
        gens_c, gens_n, keys = [], [], []
        for i in range(NV):
            if self.d - 5 < 0:
                break
            mchi = (self.chi + W[i]) % 7 if self.chi is not None else None
            ms = sector_monomials(self.d - 5, (self.parD + PI[i]) % 2, mchi)
            start, cnt = self.Nblocks[i]
            blk = Index([m for (_, m) in self.Ncols[start:start + cnt]])
            for m in ms:
                c = poly_to_vec(shift(self.Fi[i], m, 1), self.Didx)
                n = np.zeros(self.nN)
                Am = shift(self.A, m, 1)
                E = np.array(list(Am.keys()), dtype=np.int64)
                pos, _ = blk.lookup(E)
                n[start + pos] = [float(v) for v in Am.values()]
                gens_c.append(c); gens_n.append(n); keys.append((i, m))
        if not gens_c:
            return np.zeros((self.nD, 0)), np.zeros((self.nN, 0)), []
        return np.array(gens_c).T, np.array(gens_n).T, keys

    def w_space(self, verbose=None):
        verbose = self.verbose if verbose is None else verbose
        Jc, Jn, keys = self.w_generators()
        if Jc.shape[1] == 0:
            return np.zeros((self.nD, 0)), np.zeros((self.nN, 0)), {'K': 0, 'gens': 0}
        U, s, Vt = np.linalg.svd(Jc, full_matrices=False)
        r = int(np.sum(s > 1e-8 * max(s[0], 1e-300)))
        Bc = U[:, :r]
        Bn = (Jn @ Vt[:r].T) / s[:r]
        info = {'K': r, 'gens': int(Jc.shape[1]), 'syz': int(Jc.shape[1] - r)}
        if verbose:
            print(f'[w d={self.d} par={self.parD} chi={self.chi}] dim J_d={r} gens={Jc.shape[1]}', file=sys.stderr)
        return Bc, Bn, info

    def f_multiples(self):
        if self.d < 6:
            return np.zeros((self.nD, 0))
        ms = sector_monomials(self.d - 6, self.parD, self.chi)
        if not ms:
            return np.zeros((self.nD, 0))
        return np.array([poly_to_vec(shift(self.F, m, 1), self.Didx) for m in ms]).T

    def family_from_polys(self, pairs):
        """Family from explicit (D, N0) pairs (exact or float polys, N0 an exact E1 lift).  Returns Bc, Bn."""
        Bc = np.array([poly_to_vec(D, self.Didx) for D, _ in pairs]).T
        Bn = np.zeros((self.nN, len(pairs)))
        blks = []
        for i in range(NV):
            start, cnt = self.Nblocks[i]
            blks.append((start, Index([m for (_, m) in self.Ncols[start:start + cnt]])))
        for k, (D, N0) in enumerate(pairs):
            for i in range(NV):
                if N0[i]:
                    start, blk = blks[i]
                    E = np.array(list(N0[i].keys()), dtype=np.int64)
                    pos, _ = blk.lookup(E)
                    Bn[start + pos, k] = [float(v) for v in N0[i].values()]
        return Bc, Bn

    # ------------------------------------------------------------------ matrices
    def build(self, verbose=None):
        verbose = self.verbose if verbose is None else verbose
        if self._built:
            return
        t0 = time.time()
        Dm = np.array(self.Dm, dtype=np.int64).reshape(-1, 6)
        nD, nP, nN, nQ, nrows = self.nD, self.nP, self.nN, self.nQ, self.nrows
        Fterms = [(np.array(list(Fi.keys()), dtype=np.int64).reshape(-1, 6),
                   np.array([float(v) for v in Fi.values()])) for Fi in self.Fi]
        R_, C_, V_ = [], [], []
        for (i, j), (start, cnt) in self.Pblocks.items():
            if cnt == 0 or nD == 0:
                continue
            Pm = np.array([m for (_, _, m) in self.Pcols[start:start + cnt]], dtype=np.int64).reshape(-1, 6)
            MU = Dm[:, None, :]; MM = Pm[None, :, :]
            for (a, b, sgn) in ((i, j, 1.0), (j, i, -1.0)):
                coef = sgn * (MM[:, :, a] - MU[:, :, a]).astype(float)
                base = (MU + MM).copy(); base[:, :, a] -= 1
                valid = (coef != 0) & (base[:, :, a] >= 0)
                if not valid.any():
                    continue
                dsel, psel = np.nonzero(valid)
                bexp = base[dsel, psel]
                cval = coef[dsel, psel]
                Ft, Fv = Fterms[b]
                for t in range(len(Fv)):
                    ridx, _ = self.ridx.lookup(bexp + Ft[t])
                    R_.append(ridx.astype(np.int64) * nP + (start + psel))
                    C_.append(dsel.astype(np.int64))
                    V_.append(cval * Fv[t])
        if R_:
            R = np.concatenate(R_); C = np.concatenate(C_); V = np.concatenate(V_)
            self.Tall = sp.csr_matrix((V, (R, C)), shape=(nrows * nP, nD))
        else:
            self.Tall = sp.csr_matrix((nrows * nP, nD))
        self.Tall.sum_duplicates()
        R_, C_, V_ = [], [], []
        for i in range(NV):
            start, cnt = self.Nblocks[i]
            if cnt == 0 or nD == 0:
                continue
            Nm = np.array([m for (_, m) in self.Ncols[start:start + cnt]], dtype=np.int64).reshape(-1, 6)
            MU = Dm[:, None, :]; MM = Nm[None, :, :]
            coef = (MU[:, :, i] - MM[:, :, i]).astype(float)
            base = (MU + MM).copy(); base[:, :, i] -= 1
            valid = (coef != 0) & (base[:, :, i] >= 0)
            dsel, nsel = np.nonzero(valid)
            if len(dsel) == 0:
                continue
            ridx, _ = self.ridx.lookup(base[dsel, nsel])
            R_.append(ridx.astype(np.int64)); C_.append(dsel.astype(np.int64) * nN + (start + nsel)); V_.append(coef[dsel, nsel])
        if R_:
            R = np.concatenate(R_); C = np.concatenate(C_); V = np.concatenate(V_)
            self.Vb = sp.csr_matrix((V, (R, C)), shape=(nrows, nD * nN))
        else:
            self.Vb = sp.csr_matrix((nrows, max(nD * nN, 0)))
        self.Vb.sum_duplicates()
        MFQ = np.zeros((nrows, nQ))
        for k, m in enumerate(self.Qm):
            MFQ[:, k] = poly_to_vec(shift(self.F, m, 1), self.ridx)
        self.MFQ = MFQ
        if nQ:
            Qf, Rf = np.linalg.qr(MFQ)
            self.Qf = Qf
        else:
            self.Qf = np.zeros((nrows, 0))
        self._built = True
        self.t_build = time.time() - t0
        if verbose:
            print(f'[build d={self.d} par={self.parD} chi={self.chi}] Tall nnz={self.Tall.nnz} Vb nnz={self.Vb.nnz} '
                  f'MFQ {MFQ.shape} ({self.t_build:.1f}s)', file=sys.stderr)

    def D_from_c(self, c, tol=0.0):
        return {self.Dm[k]: float(c[k]) for k in range(self.nD) if abs(c[k]) > tol}


# ==============================================================================================
class KVP:
    """Koszul variable-projection functional over the family (Bc, Bn) of a Sector."""

    def __init__(self, sec: Sector, Bc, Bn, use_koszul=True, exclude_F_multiples=True, cond=1e-9):
        sec.build()
        self.sec, self.Bc, self.Bn = sec, np.asarray(Bc, float), np.asarray(Bn, float)
        self.K = self.Bc.shape[1]
        self.use_koszul = use_koszul and sec.nP > 0
        self.cond = cond
        self.P = np.eye(self.K); self.excluded_dim = 0
        self.Fc = sec.f_multiples()
        if exclude_F_multiples and self.K and self.Fc.shape[1]:
            # coordinates of the F-multiples in the family (least squares; keep those inside the family)
            S, *_ = sla.lstsq(self.Bc, self.Fc, check_finite=False)
            resid = np.linalg.norm(self.Bc @ S - self.Fc, axis=0) / np.maximum(np.linalg.norm(self.Fc, axis=0), 1e-300)
            S = S[:, resid < 1e-6]
            if S.shape[1]:
                U, s, _ = np.linalg.svd(S, full_matrices=True)
                r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
                self.P = U[:, r:]; self.excluded_dim = r
        self.dimP = self.P.shape[1]
        self._nfev = 0

    def solve(self, a):
        sec = self.sec
        c = self.Bc @ a; n = self.Bn @ a
        b = sec.Vb @ np.outer(c, n).ravel()
        Qf = sec.Qf
        bp = b - Qf @ (Qf.T @ b) if Qf.shape[1] else b
        if self.use_koszul:
            T = (sec.Tall @ c).reshape(sec.nrows, sec.nP)
            Tp = T - Qf @ (Qf.T @ T) if Qf.shape[1] else T
            xP, _, rank, sv = sla.lstsq(Tp, bp, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            r = Tp @ xP - bp
        else:
            xP = np.zeros(0); r = -bp; rank = 0
        self._nfev += 1
        return float(r @ r), c, n, b, xP, r, int(rank)

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        f, c, n, b, xP, r, rank = self.solve(a)
        return f, float(b @ b), xP, r, rank

    def value_and_grad(self, u):
        sec = self.sec
        u = np.asarray(u, float)
        af = self.P @ u; na = np.linalg.norm(af); a = af / na
        f, c, n, b, xP, r, rank = self.solve(a)
        G = (sec.Vb.T @ r).reshape(sec.nD, sec.nN)
        g = -(self.Bc.T @ (G @ n)) - (self.Bn.T @ (G.T @ c))
        if self.use_koszul:
            gT = sec.Tall.T @ np.outer(r, xP).ravel()
            g = g + self.Bc.T @ gT
        g = 2.0 * g
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def is_spurious(self, a, tol=1e-8):
        a = np.asarray(a, float)
        c = self.Bc @ (a / np.linalg.norm(a))
        if self.Fc.shape[1] == 0:
            return False, 1.0
        x, *_ = sla.lstsq(self.Fc, c, check_finite=False)
        rel = float(np.linalg.norm(self.Fc @ x - c) / max(np.linalg.norm(c), 1e-300))
        return rel < tol, rel

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-10, verbose=False, init=None, stop_rel=1e-12,
               keep_a_below=1e-6, callback=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            t1 = time.time()
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-16, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, bb, xP, r, rank = self.value(a)
            rel = f / max(bb, 1e-300)
            spur, srel = self.is_spurious(a)
            cand = {'f': f, 'rel': rel, 'b2': bb, 'spurious': bool(spur), 'div_rel': srel, 'iters': int(res.nit),
                    'nfev': int(res.nfev), 'restart': it, 'seed': seed, 'rankT': int(rank),
                    'seconds': round(time.time() - t1, 2), 'msg': str(res.message)[:50]}
            if rel < keep_a_below or it == 0:
                cand['a'] = a
            out.append(cand)
            if callback:
                callback(cand)
            if verbose:
                print(f'[kvp] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} nfev={res.nfev} spur={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    def polish(self, a, maxiter=3000, gtol=1e-14):
        u0 = self.P.T @ np.asarray(a, float)
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 30})
        a = self.P @ res.x; a /= np.linalg.norm(a)
        f, bb, xP, r, rank = self.value(a)
        return a, f, f / max(bb, 1e-300), int(res.nit)


# ==============================================================================================
def reconstruct_D(sec: Sector, c, max_dens=(10, 100, 1000, 10**4, 10**5, 10**6)):
    """Rational candidates for D from its float coefficient vector (normalized by the largest entry)."""
    c = np.asarray(c, float)
    j = int(np.argmax(np.abs(c)))
    c = c / c[j]
    cands = []
    for md in max_dens:
        fr = [Fraction(float(v)).limit_denominator(md) for v in c]
        err = max(abs(float(f) - v) for f, v in zip(fr, c))
        if err < 1e-9:
            D = {sec.Dm[k]: f for k, f in enumerate(fr) if f != 0}
            if D and D not in cands:
                cands.append(D)
    return cands


def normalize_e1(F, A, D, N):
    """Given a certificate numerator N (E = F^2 Q), return N' = N - (h/6) Z with N'(F) = A D exactly
    (h = (N(F) - A D)/F), which leaves E unchanged.  Exact arithmetic."""
    from hodge import jacobian
    Fi = [pderiv(F, i) for i in range(NV)]
    NF = {}
    for i in range(NV):
        if N[i]:
            NF = add(NF, mul(N[i], Fi[i]))
    diff = sub(NF, mul(A, D))
    if not diff:
        return [dict(n) for n in N]
    ok, h = jacobian.divisible_by(F, diff)
    if not ok:
        raise ValueError('N(F) - A D not divisible by F')
    Np = []
    for i in range(NV):
        e = [0] * 6; e[i] = 1
        Np.append(sub(N[i], scale(shift(h, tuple(e), 1), Fraction(1, 6))))
    return Np


def validate(iids=('cone-ladder-linear', 'two-reflection-cubic-surface-d2', 'f3-block-x0-quartic-strict',
                   'hesse-radial-hesse-cofactor-d4-rational'), verbose=True):
    """Rebuild the functional (unrestricted sector) for accepted q=1 witnesses; f must be ~ 0 at the witness."""
    out = {}
    for iid in iids:
        F, A, q, rec = core.load_instance(iid)
        if q != 1:
            continue
        d, D, N, Q = core.load_witness(iid)
        t0 = time.time()
        Np = normalize_e1(F, A, D, N)
        sec = Sector(F, A, d, core.parity_of(D), None)
        Bc, Bn = sec.family_from_polys([(D, Np)])
        # add a random second direction to make the family non-trivial
        kvp = KVP(sec, Bc, Bn, exclude_F_multiples=False)
        f, bb, xP, r, rank = kvp.value(np.array([1.0]))
        out[iid] = {'d': d, 'f': f, 'b2': bb, 'rel': f / max(bb, 1e-300), 'rankT': rank, 'rows': sec.nrows,
                    'nP': sec.nP, 'nQ': sec.nQ, 'seconds': round(time.time() - t0, 1)}
        if verbose:
            print(iid, out[iid], file=sys.stderr)
    return out


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--validate', action='store_true')
    args = ap.parse_args()
    if args.validate:
        print(json.dumps(validate(), indent=1))
