"""Wave-2 SEARCH lens library for lane cycle-q2 (engine v2 on sys.path; engine files untouched).

Pieces
------
Inst                      instance + diagonal group + character keys (tuple of weights on the generators)
sector_kernel(I,d,key)    exact rational basis (multi-prime CRT + rational reconstruction) of
                          K_d^key = {D in S_d^{Y-even, key} : A*D in J(F)}  (= E1-admissible denominators),
                          plus dim J_d^key and the non-Jacobian quotient dimension
wspace_basis(I,d,key)     the Jacobian part parametrised by w: D = sum_i w_i F_i, w_i sector-restricted
SectorVarPro              variable-projection functional restricted to the character sector
                            mode='family'  : D = sum a_k B_k, unknowns N (sector), Q (sector)      [full E2]
                            mode='wspace'  : D = sum_k a_k m_k F_{i_k}, N = A w/q + Koszul(c), unknowns c, Q
                                             equation divided by F:  D div N - N(D) = F^2 Q          [full E2]
                            mode='hess'    : same without Koszul corrections (Hessian relaxation)
run_search(...)           multiprocess random restarts, incremental JSONL logging, candidate handling
reconstruct_and_confirm   rational reconstruction with increasing denominators + exact sector re-solve
"""
import os, sys, time, json, math, hashlib
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from flint import nmod_mat
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, reduce_mod

LANE = 'cycle-q2'
LANEDIR = '/tmp/hodge/results/cycle-q2'
LOG = LANEDIR + '/wave2-search.log.jsonl'
CANDDIR = LANEDIR + '/wave2-cands'
KERNDIR = LANEDIR + '/wave2-kernels'
PRIMES = core.PRIMES


def log_line(rec, path=LOG):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('lane', LANE)
    rec.setdefault('engine', 3)
    with open(path, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')


def keystr(key):
    return '-'.join(str(k) for k in key)


# ----------------------------------------------------------------------------------------------
class Inst:
    def __init__(self, iid=LANE):
        F, A, q, meta = core.load_instance(iid)
        self.iid = iid
        self.F, self.A, self.q = F, A, q
        self.G = symmetry.diagonal_group(F)
        self.n = self.G.order
        self.gens = [tuple(int(x) for x in g) for g in self.G.generators]
        self.elements = [tuple(int(x) for x in g) for g in self.G.elements]
        self.ng = len(self.gens)
        self.Fi = [pderiv(F, i) for i in range(NV)]
        self.keyF = self.key(next(iter(F)))
        self.keyA = self.key(next(iter(A)))
        self.keyFi = [self.key(next(iter(self.Fi[i]))) for i in range(NV)]
        self.parFi = [core.parity_of(self.Fi[i]) for i in range(NV)]
        self.pA = core.parity_of(A)
        self.pF = core.parity_of(F)
        # generators must separate the characters of the group
        # the generator keys must determine the full character (i.e. the generators generate G):
        # test on all degree-6 monomials (every character of G occurs there)
        fullmap = {}
        for e in monomials(6):
            k = self.key(e)
            fc = tuple(sum(a * b for a, b in zip(e, w)) % self.n for w in self.elements)
            fullmap.setdefault(k, set()).add(fc)
        self.gens_separate = all(len(v) == 1 for v in fullmap.values())
        self.n_characters = len({fc for v in fullmap.values() for fc in v})
        self._mcache = {}

    def key(self, e):
        return tuple(sum(a * b for a, b in zip(e, g)) % self.n for g in self.gens)

    def kadd(self, *ks):
        out = [0] * self.ng
        for k in ks:
            for j in range(self.ng):
                out[j] = (out[j] + k[j]) % self.n
        return tuple(out)

    def kneg(self, k):
        return tuple((-x) % self.n for x in k)

    def kmul(self, c, k):
        return tuple((c * x) % self.n for x in k)

    def gkey(self, i):
        """key contribution of the variable weight w_i (used by the derivative d/dz_i)."""
        return tuple(g[i] % self.n for g in self.gens)

    def sector_monomials(self, deg, parity, key):
        """Monomials of degree deg, Y-parity parity, character key (vectorised)."""
        if deg < 0:
            return []
        ck = (deg, parity)
        if ck not in self._mcache:
            ms = monomials(deg, parity)
            if ms:
                E = np.array(ms, dtype=np.int64)
                W = np.array(self.gens, dtype=np.int64).T
                T = (E @ W) % self.n
                groups = {}
                for m, row in zip(ms, T):
                    groups.setdefault(tuple(int(x) for x in row), []).append(m)
            else:
                groups = {}
            self._mcache[ck] = groups
        return self._mcache[ck].get(tuple(key), [])

    def keys_of_degree(self, deg, parity=0):
        self.sector_monomials(deg, parity, ())  # populate
        return sorted(self._mcache[(deg, parity)].keys())

    # characters of the certificate unknowns for a denominator of character keyD
    def keyN(self, i, keyD):
        return self.kadd(self.keyA, keyD, self.kneg(self.keyF), self.gkey(i))

    def keyQ(self, keyD):
        return self.kadd(self.keyA, self.kmul(2, keyD), self.kneg(self.kmul(self.q + 1, self.keyF)))

    def keyC(self, a, b, keyD):
        return self.kadd(self.keyA, keyD, self.kneg(self.kmul(2, self.keyF)), self.gkey(a), self.gkey(b))

    def keyW(self, i, keyD):
        return self.kadd(keyD, self.kneg(self.keyF), self.gkey(i))


# ----------------------------------------------------------------------------------------------
def _nmod_from_cols(cols, p, rows=None):
    if rows is None:
        rows = {}
        for col in cols:
            for e in col:
                rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p), rows


def _rref_rows(Pr, rk, K):
    pivots, vals = [], []
    for i in range(rk):
        row = [int(Pr[i, j]) for j in range(K)]
        piv = next(j for j in range(K) if row[j])
        pivots.append(piv); vals.append(row)
    return pivots, vals


def sector_kernel(I, d, key, parity=0, nprimes=2, exact_try=True, verbose=False):
    """Basis of K_d^key (see module docstring).  Dimensions from two primes; the basis is exact rational
    (CRT + rational reconstruction, verified at a third prime) when the heights are small, otherwise a float
    orthonormal basis (SVD) -- the rational rref bases typically have ~170-digit entries in this lane.
    Returns dict with keys d, key, nmons, kernel_dim, jac_dim, nonjac_dim, ring_dim, basis (list of polys),
    exact (bool), mons, pivots, primes, shape, seconds."""
    t0 = time.time()
    Dm = I.sector_monomials(d, parity, key)
    K = len(Dm)
    out = {'d': d, 'key': list(key), 'nmons': K, 'kernel_dim': 0, 'jac_dim': 0, 'nonjac_dim': 0,
           'ring_dim': K, 'basis': [], 'exact': True, 'mons': Dm, 'pivots': [], 'primes': [], 'shape': [0, 0]}
    if K == 0:
        out['seconds'] = round(time.time() - t0, 2)
        return out
    tkey = I.kadd(I.keyA, key); tpar = (I.pA + parity) % 2
    gens = []
    for i in range(NV):
        mkey = I.kadd(tkey, I.kneg(I.keyFi[i])); mpar = (tpar - I.parFi[i]) % 2
        for m in I.sector_monomials(d + 7, mpar, mkey):
            gens.append((i, m))
    jg = []
    for i in range(NV):
        mkey = I.kadd(key, I.kneg(I.keyFi[i])); mpar = (parity - I.parFi[i]) % 2
        for m in I.sector_monomials(d - 5, mpar, mkey):
            jg.append((i, m))
    Dindex = {m: j for j, m in enumerate(Dm)}
    cols_exact = [shift(I.A, m, 1) for m in Dm] + [shift(I.Fi[i], m, 1) for (i, m) in gens]
    rows = {}
    for col in cols_exact:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols_exact)

    def compute(p):
        flat = [0] * (R * C)
        for j, col in enumerate(cols_exact):
            for e, v in col.items():
                flat[rows[e] * C + j] = (v.numerator * pow(v.denominator, -1, p)) % p
        M = nmod_mat(R, C, flat, p)
        X, nullity = M.nullspace()
        if nullity == 0:
            pivots, vals = [], []
        else:
            flat2 = [int(X[j, k]) for k in range(nullity) for j in range(K)]
            proj = nmod_mat(nullity, K, flat2, p)
            Pr, rk = proj.rref()
            pivots, vals = _rref_rows(Pr, rk, K)
        if jg:
            jcols = [shift(reduce_mod(I.Fi[i], p), m, 1, p) for (i, m) in jg]
            JM, _ = _nmod_from_cols(jcols, p, rows=dict(Dindex))
            jrank = JM.rank()
        else:
            jrank = 0
        return pivots, vals, jrank

    ps = PRIMES[:nprimes + 1]
    per = {p: compute(p) for p in ps[:nprimes]}
    dims = [len(per[p][0]) for p in ps[:nprimes]]
    jr = [per[p][2] for p in ps[:nprimes]]
    agree = (len(set(dims)) == 1 and len(set(jr)) == 1 and all(per[p][0] == per[ps[0]][0] for p in ps[:nprimes]))
    rk = min(dims); jrank = min(jr)
    basis = []; exact = False
    if exact_try and rk > 0 and agree:
        try:
            basis_rows = []
            ok = True
            for r in range(rk):
                vals_frac = []
                for j in range(K):
                    a, m = per[ps[0]][1][r][j], ps[0]
                    for p in ps[1:nprimes]:
                        a, m = core.crt_pair(a, m, per[p][1][r][j], p)
                    fr = core.ratrecon(a, m)
                    if fr is None:
                        ok = False; break
                    vals_frac.append(fr)
                if not ok:
                    break
                basis_rows.append(vals_frac)
            if ok:
                pv = ps[nprimes]
                per[pv] = compute(pv)
                ok = per[pv][0] == per[ps[0]][0]
                for r in range(rk):
                    if not ok:
                        break
                    for j in range(K):
                        fr = basis_rows[r][j]
                        if (fr.numerator * pow(fr.denominator, -1, pv)) % pv != per[pv][1][r][j] % pv:
                            ok = False; break
            if ok:
                basis = [{Dm[j]: v for j, v in enumerate(row) if v != 0} for row in basis_rows]
                exact = True
        except Exception:
            exact = False
    if rk > 0 and not exact:
        Mf = np.zeros((R, C))
        for j, col in enumerate(cols_exact):
            for e, v in col.items():
                Mf[rows[e], j] = float(v)
        Nf = sla.null_space(Mf)
        if Nf.shape[1]:
            U, s, _ = np.linalg.svd(Nf[:K, :], full_matrices=False)
            r = int(np.sum(s > 1e-9 * max(1.0, s[0])))
            if r != rk:
                out['warning'] = f'float kernel rank {r} != modular {rk}'
            basis = [{Dm[j]: float(U[j, i]) for j in range(K) if abs(U[j, i]) > 1e-14} for i in range(min(r, rk) if r >= rk else r)]
    out.update({'kernel_dim': rk, 'jac_dim': jrank, 'nonjac_dim': rk - jrank, 'ring_dim': K - jrank,
                'basis': basis, 'exact': exact, 'pivots': per[ps[0]][0], 'primes': list(per.keys()),
                'n_jgens': len(gens), 'shape': [R, C], 'agree': agree, 'seconds': round(time.time() - t0, 2)})
    if verbose:
        print(f'[kernel] d={d} key={keystr(key)} mons={K} ker={rk} J={jrank} nonJ={rk - jrank} exact={exact} '
              f'shape={out["shape"]} ({out["seconds"]}s)', file=sys.stderr)
    return out


def kernel_path(d, key, parity=0):
    return f'{KERNDIR}/d{d}-p{parity}-{keystr(key)}.json'


def save_kernel(kd, parity=0):
    rec = dict(kd)
    if kd['exact']:
        rec['basis'] = [{core.to_str({m: Fraction(1)}): str(v) for m, v in B.items()} for B in kd['basis']]
    else:
        rec['basis'] = [{core.to_str({m: Fraction(1)}): float(v) for m, v in B.items()} for B in kd['basis']]
    rec['mons'] = [list(m) for m in kd['mons']]
    with open(kernel_path(kd['d'], tuple(kd['key']), parity), 'w') as f:
        json.dump(rec, f)


def load_kernel(I, d, key, parity=0, compute=True, **kw):
    path = kernel_path(d, key, parity)
    if os.path.exists(path):
        rec = json.load(open(path))
        basis = []
        for B in rec['basis']:
            P = {}
            for ms, v in B.items():
                e = next(iter(core.parse(ms)))
                P[e] = Fraction(v) if rec['exact'] else float(v)
            basis.append(P)
        rec['basis'] = basis
        rec['mons'] = [tuple(m) for m in rec['mons']]
        rec['key'] = tuple(rec['key'])
        return rec
    if not compute:
        return None
    kd = sector_kernel(I, d, key, parity=parity, **kw)
    save_kernel(kd, parity)
    return kd


def wspace_basis(I, d, keyD, parity=0):
    """[(i, m)] with D_k = m F_i in sector keyD (deg m = d-5)."""
    out = []
    for i in range(NV):
        mpar = (parity - I.parFi[i]) % 2
        for m in I.sector_monomials(d - 5, mpar, I.keyW(i, keyD)):
            out.append((i, m))
    return out


def cone_subbasis(kd, variables):
    """Sub-family of span(basis) supported on the given variable subset (float orthonormal basis)."""
    vs = set(variables)
    basis = kd['basis']; mons = kd['mons']
    if not basis:
        return []
    bad = [m for m in mons if any(m[i] for i in range(NV) if i not in vs)]
    if not bad:
        return list(basis)
    Bm = np.zeros((len(bad), len(basis)))
    for r, m in enumerate(bad):
        for c, B in enumerate(basis):
            Bm[r, c] = float(B.get(m, 0))
    ns = sla.null_space(Bm)
    out = []
    for i in range(ns.shape[1]):
        D = {}
        for c, B in enumerate(basis):
            if abs(ns[c, i]) > 1e-14:
                for e, v in B.items():
                    D[e] = D.get(e, 0.0) + ns[c, i] * float(v)
        D = {e: v for e, v in D.items() if abs(v) > 1e-13}
        if D:
            out.append(D)
    return out


# ----------------------------------------------------------------------------------------------
def _to_float_cols(cols, rows):
    ri, ci, vals = [], [], []
    for j, col in enumerate(cols):
        for e, v in col.items():
            ri.append(rows[e]); ci.append(j); vals.append(float(v))
    return ri, ci, vals


class SectorVarPro:
    """Variable projection restricted to a character sector (see module docstring)."""

    def __init__(self, I, keyD, mode='family', basis=None, wbasis=None, parity=0, exclude=True,
                 cond=1e-9, verbose=False, ls='auto'):
        t0 = time.time()
        self.I = I; self.mode = mode; self.keyD = tuple(keyD); self.parity = parity; self.cond = cond
        F, A, q = I.F, I.A, I.q
        self.ls = ls
        if mode == 'family':
            assert basis, 'family mode needs a basis'
            degs = {core.degree(B) for B in basis}
            assert len(degs) == 1, 'family members must share one degree'
            d = degs.pop()
            self.Dk = [dict(B) for B in basis]
            self.labels = [f'B{k}' for k in range(len(basis))]
        else:
            assert wbasis, 'wspace mode needs a w-basis'
            d = core.degree(wbasis[0][1]) + 5
            self.wbasis = list(wbasis)
            self.Dk = [shift(I.Fi[i], m, 1) for (i, m) in wbasis]
            self.labels = [f'{core.to_str({m: Fraction(1)})}*F_{i}' for (i, m) in wbasis]
        self.d = d; self.K = len(self.Dk)
        pD = parity
        self.pars = fixed_d.numerator_parities({next(iter(self.Dk[0])): Fraction(1)}, A)
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        self.degN, self.degQ = degN, degQ
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        blocks = []  # per k: (ri, ci, vals) for the linear-in-a columns
        if mode == 'family':
            Ncols = []
            for i in range(NV):
                for m in I.sector_monomials(degN, self.pars[i], I.keyN(i, self.keyD)):
                    Ncols.append((i, m))
            self.cols = Ncols
            Fi = I.Fi
            for Bk in self.Dk:
                FBk = mul(F, Bk)
                Bki = [pderiv(Bk, i) for i in range(NV)]
                G = [scale(add(scale(mul(Bk, Fi[i]), q), mul(F, Bki[i])), -1) for i in range(NV)]
                ri, ci, vals = [], [], []
                for j, (i, m) in enumerate(Ncols):
                    col = shift(G[i], m, 1)
                    if m[i]:
                        m2 = list(m); m2[i] -= 1
                        col = add(col, shift(FBk, tuple(m2), m[i]))
                    for e, v in col.items():
                        ri.append(ridx(e)); ci.append(j); vals.append(float(v))
                blocks.append((ri, ci, vals))
            Qcols = I.sector_monomials(degQ, I.pA, I.keyQ(self.keyD)) if degQ >= 0 else []
            Fq1 = scale(power(F, q + 1), -1)
            self.Qcols = Qcols
            qcols = [shift(Fq1, m, 1) for m in Qcols]
            # rhs pairs
            AB = [mul(A, Bk) for Bk in self.Dk]
            pair = {}
            for k in range(self.K):
                for l in range(k, self.K):
                    v = scale(mul(AB[k], self.Dk[l]), -1)
                    pair[(k, l)] = v
                    for e in v:
                        ridx(e)
        else:
            # equation divided by F: D div N - N(D) - F^2 Q = 0, N = A w / q + K_c
            Ccols = []
            if mode == 'wspace':
                for a in range(NV):
                    for b in range(a + 1, NV):
                        cpar = (self.pars[a] - I.parFi[b]) % 2
                        for m in I.sector_monomials(d + 6 * q - 10, cpar, I.keyC(a, b, self.keyD)):
                            Ccols.append((a, b, m))
            self.cols = Ccols
            Fi = I.Fi
            Dki = [[pderiv(Dk, i) for i in range(NV)] for Dk in self.Dk]
            # U_{ab,m} = (d_a m) F_b - (d_b m) F_a   (column-only part)
            Ucol = []
            for (a, b, m) in Ccols:
                mp = {m: Fraction(1)}
                U = sub(mul(pderiv(mp, a), Fi[b]), mul(pderiv(mp, b), Fi[a]))
                Ucol.append(U)
            for k, Dk in enumerate(self.Dk):
                Vpair = {}
                ri, ci, vals = [], [], []
                for j, (a, b, m) in enumerate(Ccols):
                    if (a, b) not in Vpair:
                        Vpair[(a, b)] = sub(mul(Fi[b], Dki[k][a]), mul(Fi[a], Dki[k][b]))
                    col = sub(mul(Dk, Ucol[j]), shift(Vpair[(a, b)], m, 1))
                    for e, v in col.items():
                        ri.append(ridx(e)); ci.append(j); vals.append(float(v))
                blocks.append((ri, ci, vals))
            Qcols = I.sector_monomials(degQ, I.pA, I.keyQ(self.keyD)) if degQ >= 0 else []
            self.Qcols = Qcols
            F2 = scale(power(F, 2), -1)
            qcols = [shift(F2, m, 1) for m in Qcols]
            # quadratic part: term(k,l) = (1/q)[D_k d_{i_l}(A m_l) - A m_l d_{i_l} D_k]
            pair = {}
            Am = [shift(A, m, 1) for (i, m) in self.wbasis]
            for k in range(self.K):
                for l in range(self.K):
                    il, ml = self.wbasis[l]
                    t = sub(mul(self.Dk[k], pderiv(Am[l], il)), mul(Am[l], Dki[k][il]))
                    t = scale(t, Fraction(-1, q))
                    kk, ll = min(k, l), max(k, l)
                    pair[(kk, ll)] = add(pair.get((kk, ll), {}), scale(t, Fraction(1, 2) if k != l else Fraction(1)))
            for v in pair.values():
                for e in v:
                    ridx(e)
            # pair[(k,l)] for k<l now holds (term(k,l)+term(l,k))/2 * (-1) ; for k==l holds -term(k,k)
        qri, qci, qvals = [], [], []
        for j, col in enumerate(qcols):
            for e, v in col.items():
                qri.append(ridx(e)); qci.append(j); qvals.append(float(v))
        R = len(rows); self.rows = R; self.row_index = rows
        self.nX = len(self.cols); self.nQ = len(self.Qcols)
        self.Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nX)) for ri, ci, vals in blocks]
        self.MQ = sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nQ)).toarray() if self.nQ else np.zeros((R, 0))
        self.Vk = []
        for k in range(self.K):
            ri, ci, vals = [], [], []
            for l in range(self.K):
                v = pair.get((min(k, l), max(k, l)), {})
                for e, c in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(float(c))
            self.Vk.append(sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, self.K)))
        self.ncols = self.nX + self.nQ
        self.P = np.eye(self.K); self.excluded_dim = 0
        if exclude:
            self.P, self.excluded_dim = self._complement_of_spurious()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[svp] mode={mode} d={d} key={keystr(self.keyD)} K={self.K} rows={R} cols={self.ncols} '
                  f'(X {self.nX}, Q {self.nQ}) excluded={self.excluded_dim} dimP={self.dimP} '
                  f'built in {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------------
    def _complement_of_spurious(self):
        """Complement of {a : D(a) in F*S_{d-6}} (includes a with D(a) = 0 for redundant w)."""
        K = self.K
        dG = self.d - 6
        Fcols = []
        if dG >= 0:
            Gm = self.I.sector_monomials(dG, (self.parity - self.I.pF) % 2, self.I.kadd(self.keyD, self.I.kneg(self.I.keyF)))
            Fcols = [shift(self.I.F, m, 1) for m in Gm]
        mons = {}
        for B in self.Dk + Fcols:
            for e in B:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((max(len(mons), 1), K + len(Fcols)))
        for k, B in enumerate(self.Dk):
            for e, v in B.items():
                Mb[mons[e], k] = float(v)
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], K + j] = -float(v)
        ns = sla.null_space(Mb)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        if r == 0:
            return np.eye(K), 0
        return U[:, r:], r

    # ------------------------------------------------------------------
    def assemble(self, a):
        M = self.Mk[0] * a[0]
        for k in range(1, self.K):
            M = M + self.Mk[k] * a[k]
        M = M.toarray()
        if self.nQ:
            M = np.hstack([M, self.MQ])
        b = np.zeros(self.rows)
        for k in range(self.K):
            b += a[k] * (self.Vk[k] @ a)
        return M, b

    def _lstsq(self, M, b, exact=False):
        if exact or self.ls == 'gelsd':
            x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
            return x
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        return x

    def value(self, a, exact=True):
        a = np.asarray(a, dtype=float)
        na = np.linalg.norm(a)
        a = a / na
        M, b = self.assemble(a)
        x = self._lstsq(M, b, exact=exact)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r, float(b @ b)

    def value_and_grad(self, u):
        u = np.asarray(u, dtype=float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x = self._lstsq(M, b)
        r = M @ x - b
        f = float(r @ r)
        xX = x[:self.nX]
        g = np.zeros(self.K)
        for k in range(self.K):
            dM = self.Mk[k] @ xX
            db = 2.0 * (self.Vk[k] @ a)
            g[k] = 2.0 * float(r @ (dM - db))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def is_spurious(self, a, tol=1e-8):
        """D(a) numerically divisible by F (or zero)?"""
        a = np.asarray(a, float)
        D = {}
        for k, B in enumerate(self.Dk):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * float(v)
        nD = math.sqrt(sum(v * v for v in D.values()))
        if nD < 1e-10:
            return True, 0.0
        if self.d < 6:
            return False, 1.0
        dG = self.d - 6
        Gm = self.I.sector_monomials(dG, (self.parity - self.I.pF) % 2, self.I.kadd(self.keyD, self.I.kneg(self.I.keyF)))
        Fcols = [shift(self.I.F, m, 1) for m in Gm]
        if not Fcols:
            return False, 1.0
        mons = {}
        for C in Fcols + [D]:
            for e in C:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), len(Fcols)))
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                Mb[mons[e], j] = float(v)
        bv = np.zeros(len(mons))
        for e, v in D.items():
            bv[mons[e]] = v
        x, *_ = sla.lstsq(Mb, bv, check_finite=False)
        rel = float(np.linalg.norm(Mb @ x - bv) / max(np.linalg.norm(bv), 1e-300))
        return rel < tol, rel

    def D_float(self, a):
        D = {}
        for k, B in enumerate(self.Dk):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * float(v)
        return D

    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.Dk):
            if frac[k] == 0:
                continue
            D = add(D, scale(B, Fraction(frac[k])))
        return D

    # ------------------------------------------------------------------
    def one_restart(self, rng, maxiter=300, gtol=1e-10, init=None, polish=True):
        if init is not None:
            u0 = self.P.T @ np.asarray(init, float)
        else:
            u0 = rng.standard_normal(self.dimP)
        u0 /= np.linalg.norm(u0)
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-16})
        nit = int(res.nit)
        u = res.x
        a = self.P @ u; a /= np.linalg.norm(a)
        f, x, r, bb = self.value(a, exact=True)
        rel = f / max(bb, 1e-300)
        if polish and rel < 1e-6:
            # tighten
            res2 = minimize(self.value_and_grad, u, jac=True, method='L-BFGS-B',
                            options={'maxiter': 2 * maxiter, 'gtol': 1e-14, 'ftol': 1e-18})
            a2 = self.P @ res2.x; a2 /= np.linalg.norm(a2)
            f2, x2, r2, bb2 = self.value(a2, exact=True)
            rel2 = f2 / max(bb2, 1e-300)
            if rel2 < rel:
                a, f, x, rel, bb = a2, f2, x2, rel2, bb2
            nit += int(res2.nit)
        spur, srel = self.is_spurious(a)
        return {'value': f, 'rel': rel, 'a': a, 'spurious': bool(spur), 'div_rel': srel, 'iters': nit,
                'bnorm2': bb}


# ----------------------------------------------------------------------------------------------
def reconstruct_candidates(a, max_dens=(1, 2, 6, 12, 60, 120, 1000, 10**4, 10**5, 10**6, 10**8), tol=1e-6):
    """Rational candidates for a (normalised: largest |a_k| = 1), increasing denominators."""
    a = np.asarray(a, float)
    j = int(np.argmax(np.abs(a)))
    a = a / a[j]
    cands = []
    for md in max_dens:
        fr = tuple(Fraction(float(v)).limit_denominator(md) for v in a)
        err = max(abs(float(f) - v) for f, v in zip(fr, a))
        if fr not in [c for c, _ in cands]:
            cands.append((fr, err))
    return cands


def exact_confirm(I, D, verbose=False, use_symmetric=True):
    """Exact re-solve for a rational D.  Returns (status, data)."""
    if not D:
        return 'inconsistent', {'note': 'D = 0'}
    divF, _ = jacobian.divisible_by(I.F, D)
    if divF:
        return 'spurious', {'note': 'F divides D'}
    ok, info = fixed_d.e1_solvable(I.F, I.A, I.q, D)
    if not ok:
        return 'e1-fails', {'e1': str(info)[:200]}
    try:
        if use_symmetric:
            status, data = symmetry.solve_fixed_D_symmetric(I.F, I.A, I.q, D, verbose=verbose, group=I.G)
        else:
            status, data = fixed_d.solve_fixed_D(I.F, I.A, I.q, D, verbose=verbose)
    except Exception as ex:
        return 'error', {'note': repr(ex)[:300]}
    if status == 'solved' or (isinstance(data, dict) and data.get('N') is not None and status not in ('inconsistent',)):
        N, Q = data['N'], data['Q']
        ok2 = core.verify_witness(I.F, I.A, I.q, D, N, Q)
        return ('solved' if ok2 else 'verify-failed'), {'N': N, 'Q': Q, 'rank': data.get('rank'), 'nullity': data.get('nullity')}
    return status, {k: v for k, v in (data.items() if isinstance(data, dict) else []) if k in ('rank', 'nullity', 'note')}


def handle_candidate(svp, cand, tag, job, extra=None, confirm=True):
    """Save a low-residual candidate and try rational reconstruction + exact re-solve."""
    I = svp.I
    a = np.asarray(cand['a'], float)
    Df = svp.D_float(a / a[int(np.argmax(np.abs(a)))])
    h = hashlib.sha1(json.dumps([round(float(v), 8) for v in a]).encode()).hexdigest()[:10]
    path = f'{CANDDIR}/cand-{tag}-{h}.json'
    rec = {'tag': tag, 'job': job, 'mode': svp.mode, 'd': svp.d, 'key': list(svp.keyD), 'rel': cand['rel'],
           'value': cand['value'], 'spurious': cand['spurious'], 'div_rel': cand['div_rel'],
           'a': [float(v) for v in a], 'labels': svp.labels,
           'D_float': {core.to_str({e: Fraction(1)}): v for e, v in Df.items()}, 'attempts': []}
    if extra:
        rec.update(extra)
    if confirm and not cand['spurious']:
        tried = set()
        for fr, err in reconstruct_candidates(a):
            D = svp.D_exact(fr)
            Ds = core.to_str(D) if D else '0'
            if Ds in tried:
                continue
            tried.add(Ds)
            t0 = time.time()
            status, data = exact_confirm(I, D)
            att = {'max_err': err, 'status': status, 'seconds': round(time.time() - t0, 1),
                   'D': Ds[:4000], 'rank': data.get('rank'), 'nullity': data.get('nullity'), 'note': data.get('note')}
            rec['attempts'].append(att)
            log_line({'step': 'confirm', 'job': job, 'tag': tag, 'mode': svp.mode, 'd': svp.d, 'key': list(svp.keyD),
                      'rel': cand['rel'], 'status': status, 'max_err': err, 'seconds': att['seconds'],
                      'rank': data.get('rank'), 'nullity': data.get('nullity'), 'D': Ds[:2000]})
            if status == 'solved':
                wit = {'d': svp.d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']],
                       'Q': core.to_str(data['Q']) if data['Q'] else '0'}
                with open(f'{LANEDIR}/witness-candidate.json', 'w') as f:
                    json.dump(wit, f)
                rec['WITNESS'] = True
                log_line({'step': 'WITNESS', 'job': job, 'tag': tag, 'd': svp.d, 'D': core.to_str(D)})
                break
            if status in ('spurious',):
                break
    with open(path, 'w') as f:
        json.dump(rec, f)
    return rec


# ----------------------------------------------------------------------------------------------
_SVP = None


def _worker(args):
    seed, n, maxiter = args
    rng = np.random.default_rng(seed)
    out = []
    for it in range(n):
        c = _SVP.one_restart(rng, maxiter=maxiter)
        c['seed'] = seed; c['restart'] = it
        out.append(c)
    return out


def hist_decades(rels):
    """Counts of log10(rel) by decade from >= 0 down to < -16."""
    h = {}
    for r in rels:
        if r <= 0:
            k = '-inf'
        else:
            k = str(int(math.floor(math.log10(r))))
        h[k] = h.get(k, 0) + 1
    return h


def run_search(svp, tag, job, restarts, cpus=1, seed0=0, chunk=None, maxiter=300, cand_tol=1e-10,
               save_tol=1e-8, extra=None, verbose=True):
    """Random-restart search with incremental logging.  Returns summary dict."""
    global _SVP
    _SVP = svp
    import multiprocessing as mp
    chunk = chunk or max(1, min(64, restarts // max(1, cpus)))
    per = max(1, chunk // max(1, cpus))
    done = 0
    allrels = []
    best = None
    t0 = time.time()
    seed = seed0
    while done < restarts:
        n_this = min(chunk, restarts - done)
        tasks = []
        rem = n_this
        while rem > 0:
            k = min(per, rem)
            tasks.append((seed, k, maxiter)); seed += 1; rem -= k
        if cpus > 1 and len(tasks) > 1:
            ctx = mp.get_context('fork')
            with ctx.Pool(min(cpus, len(tasks))) as pool:
                results = pool.map(_worker, tasks)
        else:
            results = [_worker(t) for t in tasks]
        cands = [c for res in results for c in res]
        done += len(cands)
        rels = [c['rel'] for c in cands]
        allrels.extend(rels)
        nonspur = [c for c in cands if not c['spurious']]
        cb = min(nonspur, key=lambda c: c['rel']) if nonspur else None
        if cb is not None and (best is None or cb['rel'] < best['rel']):
            best = cb
        rec = {'step': 'search', 'job': job, 'tag': tag, 'mode': svp.mode, 'd': svp.d, 'key': list(svp.keyD),
               'family_dim': svp.K, 'dimP': svp.dimP, 'excluded': svp.excluded_dim, 'rows': svp.rows,
               'cols': svp.ncols, 'nX': svp.nX, 'nQ': svp.nQ, 'restarts': len(cands), 'restarts_total': done,
               'best_rel': (cb['rel'] if cb else None), 'best_rel_sofar': (best['rel'] if best else None),
               'n_spurious': len(cands) - len(nonspur),
               'rel_quantiles': [float(np.quantile(rels, qq)) for qq in (0.0, 0.1, 0.5, 0.9, 1.0)] if rels else None,
               'hist': hist_decades([c['rel'] for c in nonspur]), 'iters_mean': float(np.mean([c['iters'] for c in cands])),
               'seconds': round(time.time() - t0, 1), 'seed_range': [seed0, seed - 1]}
        if extra:
            rec.update(extra)
        log_line(rec)
        if verbose:
            print(f'[search] {tag} {done}/{restarts} best={rec["best_rel"]:.3e} hist={rec["hist"]} ({rec["seconds"]}s)',
                  file=sys.stderr, flush=True)
        for c in nonspur:
            if c['rel'] < save_tol:
                handle_candidate(svp, c, tag, job, extra=extra, confirm=(c['rel'] < cand_tol))
    return {'best': best, 'rels': allrels, 'seconds': time.time() - t0}


# ----------------------------------------------------------------------------------------------
class HessRelax(SectorVarPro):
    """Hessian relaxation (NECESSARY condition for a Jacobian-type denominator D = sum_i w_i F_i):
         quad(w) := D div(A w)/q - (A/q) w.grad(D)  must lie in J(F)_{2d+6} (sector, parity of A),
       because the full E2 reads quad(w) + T_D(c) - F^2 Q = 0 with T_D(c), F^2 Q in J(F).
       f(w) = || Pi_{J-perp} quad(w) ||^2 ; the number of independent conditions is dim (S/J)_{2d+6}^{sector}
       (zero for 2d+6 > 24, i.e. d >= 10; tiny for d = 8, 9)."""

    def __init__(self, I, keyD, wbasis, parity=0, exclude=True, verbose=False):
        t0 = time.time()
        self.I = I; self.mode = 'hessJ'; self.keyD = tuple(keyD); self.parity = parity; self.cond = 1e-9
        self.ls = 'auto'
        q = I.q
        self.wbasis = list(wbasis); self.K = len(wbasis)
        d = core.degree(wbasis[0][1]) + 5; self.d = d
        self.Dk = [shift(I.Fi[i], m, 1) for (i, m) in wbasis]
        self.labels = [f'{core.to_str({m: Fraction(1)})}*F_{i}' for (i, m) in wbasis]
        self.pars = fixed_d.numerator_parities({next(iter(self.Dk[0])): Fraction(1)}, I.A)
        Dki = [[pderiv(Dk, i) for i in range(NV)] for Dk in self.Dk]
        Am = [shift(I.A, m, 1) for (i, m) in wbasis]
        pair = {}
        for k in range(self.K):
            for l in range(self.K):
                il, ml = wbasis[l]
                t = sub(mul(self.Dk[k], pderiv(Am[l], il)), mul(Am[l], Dki[k][il]))
                t = scale(t, Fraction(1, q))
                kk, ll = min(k, l), max(k, l)
                pair[(kk, ll)] = add(pair.get((kk, ll), {}), scale(t, Fraction(1, 2) if k != l else Fraction(1)))
        degE = 2 * d + 6
        keyE = I.kadd(I.keyA, I.kmul(2, self.keyD), I.kneg(I.keyF))
        jcols = []
        for i in range(NV):
            mkey = I.kadd(keyE, I.kneg(I.keyFi[i])); mpar = (I.pA - I.parFi[i]) % 2
            for m in I.sector_monomials(degE - 5, mpar, mkey):
                jcols.append(shift(I.Fi[i], m, 1))
        rows = {}
        for v in pair.values():
            for e in v:
                rows.setdefault(e, len(rows))
        for c in jcols:
            for e in c:
                rows.setdefault(e, len(rows))
        R = len(rows); self.rows = R; self.row_index = rows
        Jm = np.zeros((R, len(jcols)))
        for j, c in enumerate(jcols):
            for e, v in c.items():
                Jm[rows[e], j] = float(v)
        # orthonormal basis of the orthogonal complement of col(Jm)
        if len(jcols):
            U, s, _ = np.linalg.svd(Jm, full_matrices=True)
            rk = int(np.sum(s > 1e-9 * max(1.0, s[0])))
            C = U[:, rk:]
        else:
            rk = 0; C = np.eye(R)
        self.jrank = rk; self.ncond = C.shape[1]; self.nJ = len(jcols)
        self.Cmat = C
        # Wk = C^T Vk  (ncond x K)
        self.Wk = []
        for k in range(self.K):
            V = np.zeros((R, self.K))
            for l in range(self.K):
                v = pair.get((min(k, l), max(k, l)), {})
                for e, c in v.items():
                    V[rows[e], l] = float(c)
            self.Wk.append(C.T @ V)
        self.nX = 0; self.nQ = 0; self.ncols = self.nJ
        self.P = np.eye(self.K); self.excluded_dim = 0
        if exclude:
            self.P, self.excluded_dim = self._complement_of_spurious()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[hessJ] d={d} key={keystr(self.keyD)} K={self.K} rows={R} J-gens={self.nJ} rankJ={rk} '
                  f'conditions={self.ncond} excluded={self.excluded_dim} dimP={self.dimP} built in {self.build_time:.1f}s',
                  file=sys.stderr)

    def _b(self, a):
        b = np.zeros(self.ncond)
        for k in range(self.K):
            b += a[k] * (self.Wk[k] @ a)
        return b

    def value(self, a, exact=True):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        b = self._b(a)
        # normaliser: squared norm of quad(w) itself is not available after projection; use |a|^4 = 1
        self._nfev += 1
        return float(b @ b), None, None, 1.0

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
        b = self._b(a)
        f = float(b @ b)
        g = np.zeros(self.K)
        for k in range(self.K):
            g[k] = 4.0 * float(b @ (self.Wk[k] @ a))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full
