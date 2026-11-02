"""(U,V)-tensor-frame Hessian relaxation for Jacobian-type denominators of F = (f(U)+f(V))/2.

R_F/J = R_f (x) R_f (bigraded).  For cofactor cubics c_i(U,V) (i = 0,1,2) the two parity sectors are
    D_+ = T + swap(T),  D_- = T - swap(T),   T = sum_i c_i f_i(U),
and both have the SAME Hessian form  Hess_F(v,v) ∝ H + swap(H),  H = sum_ij c_i c_j f_ij(U).
The relaxation  A * Hess in J(F)  is therefore one condition on c for both sectors, invariant under the
Euler directions c_i -> c_i + U_i h (h any quadric; these are exactly F-multiples and the spurious g*h family).
Everything here is exact (Q or mod p) and small: normal forms are tensor products of ternary normal forms.
"""
import sys, time, itertools, json, os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np, scipy.sparse as sp
from fractions import Fraction as Fr
import w2lib as L
from hodge import core

def tmons(k):
    return [(a, b, k - a - b) for a in range(k + 1) for b in range(k + 1 - a)]

def tmul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2]); out[e] = out.get(e, 0) + c1 * c2
    return {e: v for e, v in out.items() if v}

def tderiv(P, i):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; out[tuple(e2)] = c * e[i]
    return out

class TernaryRing:
    """R_f = C[z0,z1,z2]/J(f) over Q (Fractions), degrees 0..12 (socle 12; R_k = 0 for k > 12)."""
    def __init__(self, f3, maxdeg=12):
        self.f = f3; self.fi = [tderiv(f3, i) for i in range(3)]
        self.maxdeg = maxdeg
        self.std = {}; self.nf = {}   # nf[(monomial)] -> dict std-index -> Fraction
        for k in range(maxdeg + 1):
            mons = tmons(k); idx = {m: i for i, m in enumerate(mons)}
            gens = []
            if k >= 5:
                for m in tmons(k - 5):
                    for Pi in self.fi:
                        gens.append(tmul({m: Fr(1)}, Pi))
            # RREF over Q
            rows = [[g.get(m, Fr(0)) for m in mons] for g in gens]
            piv = []; R = []
            for r in rows:
                r = list(r)
                for (pc, pr) in zip(piv, R):
                    if r[pc]:
                        c = r[pc]; r = [x - c * y for x, y in zip(r, pr)]
                j = next((j for j, x in enumerate(r) if x), None)
                if j is None: continue
                c = r[j]; r = [x / c for x in r]
                # eliminate j from existing rows
                for t, (pc, pr) in enumerate(zip(piv, R)):
                    if pr[j]:
                        cc = pr[j]; R[t] = [x - cc * y for x, y in zip(pr, r)]
                piv.append(j); R.append(r)
            pivset = set(piv)
            std = [m for j, m in enumerate(mons) if j not in pivset]
            self.std[k] = std
            sidx = {m: i for i, m in enumerate(std)}
            # normal form of each monomial: standard -> itself; pivot -> -(rest of its RREF row)
            for j, m in enumerate(mons):
                if j in pivset:
                    r = R[piv.index(j)]
                    self.nf[m] = {sidx[mons[t]]: -r[t] for t in range(len(mons)) if t != j and r[t]}
                else:
                    self.nf[m] = {sidx[m]: Fr(1)}
        self.dims = [len(self.std[k]) for k in range(maxdeg + 1)]

    def NF(self, m):
        k = sum(m)
        if k > self.maxdeg: return {}
        return self.nf[m]

class TensorBasis:
    """Basis of (R_f (x) R_f)_n: triples (a, i, j), a = U-degree, i in std[a], j in std[n-a]."""
    def __init__(self, R, n):
        self.R, self.n = R, n
        self.index = {}; self.labels = []
        for a in range(n + 1):
            b = n - a
            if a > R.maxdeg or b > R.maxdeg: continue
            for i in range(R.dims[a]):
                for j in range(R.dims[b]):
                    self.index[(a, i, j)] = len(self.labels); self.labels.append((a, i, j))
        self.dim = len(self.labels)
    def nf_uv(self, e):
        """normal form of the UV monomial e (6-tuple) as dict index -> Fraction."""
        a = e[0] + e[1] + e[2]; b = e[3] + e[4] + e[5]
        if a > self.R.maxdeg or b > self.R.maxdeg: return {}
        nu = self.R.NF((e[0], e[1], e[2])); nv = self.R.NF((e[3], e[4], e[5]))
        out = {}
        for i, x in nu.items():
            for j, y in nv.items():
                out[self.index[(a, i, j)]] = x * y
        return out

def uv_poly_to_tensor(P, TB):
    out = {}
    for e, c in P.items():
        for k, v in TB.nf_uv(e).items():
            out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}

class TensorRelaxation:
    """Quadrics of the relaxation in cofactor coordinates c (K = 3 * #cofactor UV monomials), vectorised.
    modp=None -> float64 tensor Wt; modp=p -> additionally exact int64 quadric coefficients mod p (self.Wmod)."""
    def __init__(self, d, verbose=True, modp=None, mixed=False):
        t0 = time.time()
        self.mixed = mixed
        F, A, q = L.load()
        f6 = L.f_ternary()
        f3 = {(e[0], e[1], e[2]): Fr(c) for e, c in f6.items()}
        self.R = TernaryRing(f3)
        self.d = d; n = 2 * d
        self.TB = TensorBasis(self.R, n)
        self.Auv = L.xy_to_uv(A)
        fij3 = [[tderiv(self.R.fi[i], j) for j in range(3)] for i in range(3)]
        cm = core.monomials(d - 5, None)
        self.cmons = cm
        self.basis = [(i, m) for i in range(3) for m in cm]
        if mixed:
            self.basis = [(i, m, 'U') for i in range(3) for m in cm] + [(i, m, 'V') for i in range(3) for m in cm]
        K = len(self.basis); self.K = K
        # --- degree-n UV monomial index and NF matrix (sparse)
        mons_n = core.monomials(n, None)
        En = np.array(mons_n, dtype=np.int64); keys_n = En @ (32 ** np.arange(6, dtype=np.int64))
        order = np.argsort(keys_n); skeys = keys_n[order]
        def lookup(E):
            k = E @ (32 ** np.arange(6, dtype=np.int64)); pos = np.searchsorted(skeys, k)
            assert (skeys[pos] == k).all(); return order[pos]
        rows, cols, vals, valsp = [], [], [], []
        p = modp
        for r, e in enumerate(mons_n):
            a = e[0] + e[1] + e[2]; b = n - a
            if a > self.R.maxdeg or b > self.R.maxdeg: continue
            nu = self.R.NF((e[0], e[1], e[2])); nv = self.R.NF((e[3], e[4], e[5]))
            for i, x in nu.items():
                for j, y in nv.items():
                    rows.append(r); cols.append(self.TB.index[(a, i, j)]); v = x * y; vals.append(float(v))
                    if p: valsp.append((v.numerator * pow(v.denominator, -1, p)) % p)
        NF = sp.csr_matrix((vals, (rows, cols)), shape=(len(mons_n), self.TB.dim))
        NFp = sp.csr_matrix((np.array(valsp, dtype=np.int64), (rows, cols)), shape=(len(mons_n), self.TB.dim)) if p else None
        t1 = time.time()
        # --- pair polynomials  P_kl = A * (T + swap T),  T = m_k m_l f_{ik il}(U)
        Aexp = np.array(list(self.Auv.keys()), dtype=np.int64); Acoef = np.array([float(v) for v in self.Auv.values()])
        Acoefp = np.array([(Fr(v).numerator * pow(Fr(v).denominator, -1, p)) % p for v in self.Auv.values()], dtype=np.int64) if p else None
        fexp = {}; fcoef = {}; fcoefp = {}
        for i in range(3):
            for j in range(3):
                items = list(fij3[i][j].items())
                fexp[(i, j)] = np.array([(e[0], e[1], e[2], 0, 0, 0) for e, _ in items], dtype=np.int64).reshape(-1, 6)
                fcoef[(i, j)] = np.array([float(c) for _, c in items])
                if p: fcoefp[(i, j)] = np.array([(c.numerator * pow(c.denominator, -1, p)) % p for _, c in items], dtype=np.int64)
        pairs = [(k, l) for k in range(K) for l in range(k, K)]
        self.pairs = pairs
        tr, tc, tv, tvp = [], [], [], []
        swap_idx = [3, 4, 5, 0, 1, 2]
        for pr, (k, l) in enumerate(pairs):
            if mixed:
                ik, mk, bk = self.basis[k]; il, ml, bl = self.basis[l]
                if bk != bl: continue                          # block-diagonal Hessian: no U/V cross terms
            else:
                ik, mk = self.basis[k]; il, ml = self.basis[l]; bk = 'U'
            fe = fexp[(ik, il)]
            if fe.shape[0] == 0: continue
            mm = np.array(mk, dtype=np.int64) + np.array(ml, dtype=np.int64)
            T = fe + mm[None, :]                              # (nf, 6)   f_ij(U) m_k m_l
            if mixed:
                TT = T if bk == 'U' else T[:, swap_idx] if False else None
                # block U: f_ij(U) m_k m_l ; block V: f_ij(V) m_k m_l  (f_ij(V) = swap of f_ij(U) exponents; m's stay)
                if bk == 'U':
                    TT = T
                else:
                    TT = fe[:, swap_idx] + mm[None, :]
                cc = fcoef[(ik, il)]
            else:
                Ts = T[:, swap_idx]
                TT = np.concatenate([T, Ts], axis=0); cc = np.concatenate([fcoef[(ik, il)], fcoef[(ik, il)]])
            E = (TT[:, None, :] + Aexp[None, :, :]).reshape(-1, 6)
            co = (cc[:, None] * Acoef[None, :]).reshape(-1)
            idx = lookup(E)
            tr.extend([pr] * len(idx)); tc.extend(idx.tolist()); tv.extend(co.tolist())
            if p:
                ccp = fcoefp[(ik, il)] if mixed else np.concatenate([fcoefp[(ik, il)], fcoefp[(ik, il)]])
                cop = ((ccp[:, None] * Acoefp[None, :]) % p).reshape(-1)
                tvp.extend(cop.tolist())
        Terms = sp.csr_matrix((tv, (tr, tc)), shape=(len(pairs), len(mons_n)))
        W = (Terms @ NF).toarray()                            # pairs x tensor dim
        self.Wpairs = W
        if p:
            Termsp = sp.csr_matrix((np.array(tvp, dtype=np.int64), (tr, tc)), shape=(len(pairs), len(mons_n)))
            self.Wmod = L.spmod(Termsp, NFp, p)               # exact quadric coefficients mod p (pairs x tensor dim)
            self.p = p
        t2 = time.time()
        # independent coordinates of the swap-antisymmetric part: blocks (a,b) with a < b fully, block (d,d) with i < j
        nzc = []
        for t, (a, i, j) in enumerate(self.TB.labels):
            b = n - a
            if mixed or a < b or (a == b and i < j):
                nzc.append(t)
        nzc = np.array(nzc)
        self.coord_rows = nzc; self.nquot = len(nzc)
        Wt = np.zeros((self.nquot, K, K))
        Wn = W[:, nzc]
        # sanity: the symmetric part must vanish identically (A is swap-antisymmetric)
        symchk = 0.0 if mixed else max((abs(W[:, t]).max() for t, (a, i, j) in enumerate(self.TB.labels) if a == n - a and i == j), default=0.0)
        self.sym_part_max = float(symchk)
        for pr, (k, l) in enumerate(pairs):
            if k == l: Wt[:, k, k] = Wn[pr]
            else: Wt[:, k, l] = Wn[pr]; Wt[:, l, k] = Wn[pr]     # (k,l) and (l,k) both occur in the quadratic form
        self.Wt = Wt
        self.W2 = np.ascontiguousarray(Wt.reshape(self.nquot * K, K))   # for BLAS gemv
        E = []
        self.quadrics_uv = core.monomials(d - 6, None)
        if mixed:
            pos = {(i, m, b): k for k, (i, m, b) in enumerate(self.basis)}
            for blk, off in (('U', 0), ('V', 3)):
                for h in self.quadrics_uv:
                    v = np.zeros(K)
                    for i in range(3):
                        m = list(h); m[off + i] += 1; v[pos[(i, tuple(m), blk)]] = 1.0   # c_i = U_i h  or  ct_i = V_i h
                    E.append(v)
        else:
            pos = {(i, m): k for k, (i, m) in enumerate(self.basis)}
            for h in self.quadrics_uv:
                v = np.zeros(K)
                for i in range(3):
                    m = list(h); m[i] += 1; v[pos[(i, tuple(m))]] = 1.0
                E.append(v)
        self.Euler = np.array(E).T
        self.P = L.complement(self.Euler, K)
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        if verbose:
            print(f'[TensorRelax] d={d} K={K} tensor dim {self.TB.dim} nonzero coords {self.nquot} Euler dirs {self.Euler.shape[1]} dimP {self.dimP} '
                  f'R_f dims {self.R.dims} sym-part max {self.sym_part_max:.1e} NF {t1-t0:.1f}s pairs {t2-t1:.1f}s total {self.build_time:.1f}s', file=sys.stderr, flush=True)

    def r_and_J(self, c):
        Wa = (self.W2 @ c).reshape(self.nquot, self.K)
        return Wa @ c, 2 * Wa

    def value(self, u):
        c = self.P @ u
        r, J = self.r_and_J(c)
        return float(r @ r), r, J @ self.P

    def D_xy(self, c, sector=0):
        """denominator in the XY frame for cofactor vector c and sector (0: +swap, 1: -swap)."""
        T = {}
        for k, (i, m) in enumerate(self.basis):
            if abs(c[k]) < 1e-15: continue
            for e, v in self.R.fi[i].items():
                ee = (e[0] + m[0], e[1] + m[1], e[2] + m[2], m[3], m[4], m[5]); T[ee] = T.get(ee, 0.0) + c[k] * float(v)
        D = dict(T)
        for e, v in T.items():
            es = (e[3], e[4], e[5], e[0], e[1], e[2]); D[es] = D.get(es, 0.0) + (v if sector == 0 else -v)
        return L.uv_to_xy({e: Fr(v).limit_denominator(10**12) for e, v in D.items() if abs(v) > 1e-15})

    def sector_subspaces(self):
        """Euler directions by quadric swap-parity: symmetric h (12) = F-multiples of D_+ / spurious g*h of D_-;
        antisymmetric h (9) = spurious g*h of D_+ / F-multiples of D_-."""
        sym, anti = [], []
        for t, h in enumerate(self.quadrics_uv):
            hs = (h[3], h[4], h[5], h[0], h[1], h[2])
            if hs == h: sym.append(self.Euler[:, t])
            elif h < hs:
                ts = self.quadrics_uv.index(hs)
                sym.append(self.Euler[:, t] + self.Euler[:, ts]); anti.append(self.Euler[:, t] - self.Euler[:, ts])
        return np.array(sym).T, np.array(anti).T
