#!/usr/bin/env python3
"""Wave-2 joint search library for lane reflection-coupled-residual (q=1, A = X0^3 Y0^3).

Klein group chi(e) = (e5 mod 2, (e3+e4) mod 2); chi(F) = (0,0); chi(A) = (0,1).
Certificate (Euler normalised, fixed D):  E1: N(F) = A D,  E2: B_D(N) := D div N - N(D) in (F).
E2-mod-J obstruction  ob(D) = [B_D(N)] in R_{2d} = S_{2d}/J_{2d}  is independent of the E1 solution N.
For D = sum g_i F_i (Jacobian type) with N = A g:   ob(g) = -[A * g^T Hess(F) g]   ("Hessian relaxation").
For D = D0 + sum g_i F_i, N = N0 + A g:
   ob(g) = [B(D0,N0)] + [D0 (g(A) + A div g) - A g(D0) - N0^T H g] - [A g^T H g].

Subcommands (run inside the lane dir, engine on PYTHONPATH):
  prep    --deg DEG --sec ab           cokernel basis of J_DEG^{sec} (float SVD, exact rank mod p) -> coker_DEG_ab.npz
  e1      --d D --sec ab               E1-admissible space of degree D in sector ab (exact dims mod p + float basis)
  relax   --d D --sec ab|e|o [--affine kernel.json:idx] --restarts N --seed S   LM search on the relaxation
  msolve  --d D --sec ab               exact decision of the relaxation mod p with msolve
  varpro  --d D --sec ab --family J|E1|K8J --restarts N --seed S   full-system variable projection
  confirm --D EXPR                     exact pipeline for one denominator
"""
from __future__ import annotations
import argparse, json, math, os, sys, time, hashlib, subprocess
from fractions import Fraction
import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares, minimize
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat

LANE = 'reflection-coupled-residual'
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'wave2-search.log.jsonl')
P = core.PRIMES[0]
SECS = [(0, 0), (0, 1), (1, 0), (1, 1)]
# control-instance override (positive controls on other (F, A) with the same Klein structure):
#   HODGE_RCR_INSTANCE=<catalogue id>  or  HODGE_RCR_INSTANCE=custom:<name> with HODGE_RCR_F / HODGE_RCR_A expressions
INST_ID = os.environ.get('HODGE_RCR_INSTANCE', LANE)
CACHE_PREFIX = '' if INST_ID == LANE else INST_ID.replace(':', '_') + '_'


def chi(e):
    return (e[5] & 1, (e[3] + e[4]) & 1)


def chi_var(i):
    return (0, 0) if i < 3 else ((1, 0) if i == 5 else (0, 1))


def sadd(a, b):
    return ((a[0] + b[0]) & 1, (a[1] + b[1]) & 1)


def sec_str(s):
    return f'{s[0]}{s[1]}'


def parse_sec(s):
    if s in ('e', 'even'):
        return 'e'
    if s in ('o', 'odd'):
        return 'o'
    if s == 'all':
        return 'all'
    return (int(s[0]), int(s[1]))


def sectors_of(spec):
    """'e' -> Y-even Klein sectors, 'o' -> Y-odd, 'all' -> all four, else the single sector."""
    if spec == 'all':
        return list(SECS)
    if spec == 'e':
        return [(0, 0), (1, 1)]
    if spec == 'o':
        return [(0, 1), (1, 0)]
    return [spec]


def sector_monomials(deg, sec):
    return [e for e in monomials(deg) if chi(e) == sec]


def log_event(rec):
    rec = dict(rec)
    rec.setdefault('lane', LANE)
    if INST_ID != LANE:
        rec.setdefault('control_instance', INST_ID)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', fixed_d.ENGINE_VERSION)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')


_INST = None


def instance():
    global _INST
    if _INST is None:
        if INST_ID.startswith('custom:'):
            F = core.parse(os.environ['HODGE_RCR_F']); A = core.parse(os.environ['HODGE_RCR_A']); q = 1
        else:
            F, A, q, rec = core.load_instance(INST_ID)
        assert q == 1 and chi(next(iter(A))) == (0, 1) and len({chi(e) for e in F}) == 1 and chi(next(iter(F))) == (0, 0)
        Fi = [pderiv(F, i) for i in range(NV)]
        Hs = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
        _INST = (F, A, q, Fi, Hs)
    return _INST


# ----------------------------------------------------------------------------
# Jacobian generators / cokernel bases
# ----------------------------------------------------------------------------

def gens_in_sector(deg, sec, p=None):
    """(i, m, m*F_i) for all generators of J_deg lying in Klein sector sec."""
    F, A, q, Fi, Hs = instance()
    Fip = [core.reduce_mod(f, p) for f in Fi] if p else Fi
    out = []
    if deg < 5:
        return out
    for i in range(NV):
        need = sadd(sec, chi_var(i))
        for m in sector_monomials(deg - 5, need):
            out.append((i, m, shift(Fip[i], m, 1, p)))
    return out


def coker_path(deg, sec):
    return os.path.join(HERE, f'{CACHE_PREFIX}coker_{deg}_{sec_str(sec)}.npz')


def prep_coker(deg, sec, force=False, verbose=True):
    """Float orthonormal basis B (rows) of the left null space of the generator matrix of J_deg^sec,
    i.e. coordinates on R_deg^sec = S_deg^sec / J_deg^sec.  Rank taken exactly mod p."""
    path = coker_path(deg, sec)
    if os.path.exists(path) and not force:
        return load_coker(deg, sec)
    t0 = time.time()
    mons = sector_monomials(deg, sec)
    idx = {e: i for i, e in enumerate(mons)}
    gens = gens_in_sector(deg, sec)
    R, C = len(mons), len(gens)
    # exact rank mod p
    gp = gens_in_sector(deg, sec, p=P)
    flat = [0] * (R * C)
    for j, (_, _, g) in enumerate(gp):
        for e, v in g.items():
            flat[idx[e] * C + j] = v
    rank = nmod_mat(R, C, flat, P).rank() if C else 0
    del flat
    t1 = time.time()
    G = np.zeros((R, C))
    for j, (_, _, g) in enumerate(gens):
        for e, v in g.items():
            G[idx[e], j] = float(v)
    if C:
        # left null space of G via the Gram matrix (R x R): eigenvectors of the R-rank smallest eigenvalues
        Gram = G @ G.T
        w, U = sla.eigh(Gram, check_finite=False, driver='evd')
        B = np.ascontiguousarray(U[:, :R - rank].T)
        # one refinement step: orthogonalise against the column space more precisely (Gram squares the conditioning)
        gap = float(w[R - rank] / max(w[R - rank - 1], 1e-300)) if R - rank - 1 >= 0 and R - rank < R else float('inf')
        quality = float(np.linalg.norm(B @ G) / np.linalg.norm(G)) if R - rank else 0.0
    else:
        B = np.eye(R); gap = float('inf'); quality = 0.0
    tmp = f'{path}.{os.getpid()}.tmp.npz'
    np.savez_compressed(tmp, B=B, mons=np.array(mons, dtype=np.int16), rank=rank)
    os.replace(tmp, path)
    info = {'event': 'coker', 'deg': deg, 'sec': sec_str(sec), 'rows': R, 'gens': C, 'rank_modp': rank,
            'dimR': R - rank, 'eig_gap': gap, 'quality_BG': quality, 't_rank': round(t1 - t0, 1), 't_total': round(time.time() - t0, 1)}
    log_event(info)
    if verbose:
        print('[coker]', info, flush=True)
    return B, mons, rank


def load_coker(deg, sec):
    z = np.load(coker_path(deg, sec))
    mons = [tuple(int(x) for x in row) for row in z['mons']]
    return z['B'], mons, int(z['rank'])


# ----------------------------------------------------------------------------
# E1-admissible denominators (exact dims mod p, float basis)
# ----------------------------------------------------------------------------

def e1_path(d, sec):
    return os.path.join(HERE, f'{CACHE_PREFIX}e1_{d}_{sec_str(sec)}')


def e1_space(d, sec, force=False, verbose=True):
    """E1-admissible denominators of degree d in Klein sector sec: {D : A D in J_{d+6}}.
    Exact (mod p): kernel_dim (incl. J_d), dim J_d, non-Jacobian classes (reduced basis mod J_d).
    Float: orthonormal basis of the full admissible subspace in the monomial coordinates of S_d^sec."""
    base = e1_path(d, sec)
    if os.path.exists(base + '.json') and os.path.exists(base + '.npz') and not force:
        return json.load(open(base + '.json'))
    F, A, q, Fi, Hs = instance()
    t0 = time.time()
    Dm = sector_monomials(d, sec)
    K = len(Dm)
    tsec = sadd(sec, chi(next(iter(A))))
    Ap = core.reduce_mod(A, P)
    gens = gens_in_sector(d + 6, tsec, p=P)
    colsp = [shift(Ap, m, 1, P) for m in Dm] + [g[2] for g in gens]
    tm = sector_monomials(d + 6, tsec)
    tidx = {e: i for i, e in enumerate(tm)}
    R, C = len(tm), len(colsp)
    flat = [0] * (R * C)
    for j, col in enumerate(colsp):
        for e, v in col.items():
            flat[tidx[e] * C + j] = v
    M = nmod_mat(R, C, flat, P)
    del flat
    X, nullity = M.nullspace()
    del M
    vecs = []
    if nullity:
        proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], P)
        del X
        Pr, rk = proj.rref()
        for i in range(rk):
            vecs.append({Dm[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])})
    kernel_dim = len(vecs)
    # reduce modulo J_d^sec
    gd = gens_in_sector(d, sec, p=P)
    didx = {e: i for i, e in enumerate(Dm)}
    if gd:
        flat = [0] * (len(gd) * K)
        for r, (_, _, g) in enumerate(gd):
            for e, v in g.items():
                flat[r * K + didx[e]] = v
        Jr, jrank = nmod_mat(len(gd), K, flat, P).rref()
        pivots, rows = [], []
        for i in range(jrank):
            j = 0
            while j < K and int(Jr[i, j]) == 0:
                j += 1
            pivots.append(j)
            rows.append({Dm[k]: int(Jr[i, k]) for k in range(j, K) if int(Jr[i, k])})
    else:
        jrank, pivots, rows = 0, [], []

    def nf(Pdict):
        Pd = dict(Pdict)
        for piv, row in zip(pivots, rows):
            c = Pd.get(Dm[piv], 0)
            if c:
                Pd = sub(Pd, scale(row, c, P), P)
        return Pd
    reduced = [v for v in (nf(v) for v in vecs) if v]
    basis = []
    if reduced:
        flat = [0] * (len(reduced) * K)
        for r, v in enumerate(reduced):
            for e, c in v.items():
                flat[r * K + didx[e]] = c
        Qr, qrk = nmod_mat(len(reduced), K, flat, P).rref()
        for i in range(qrk):
            basis.append({Dm[j]: int(Qr[i, j]) for j in range(K) if int(Qr[i, j])})
    t1 = time.time()
    # float basis via the cokernel of J_{d+6}^{tsec}
    B, tmons, trank = prep_coker(d + 6, tsec, verbose=verbose)
    assert tmons == tm
    V = np.zeros((R, K))
    for k, m in enumerate(Dm):
        for e, v in shift(A, m, 1).items():
            V[tidx[e], k] = float(v)
    BV = B @ V
    _, s, Vt = sla.svd(BV, full_matrices=True, check_finite=False)
    fdim = kernel_dim
    Kb = np.ascontiguousarray(Vt[K - fdim:, :]) if fdim else np.zeros((0, K))
    sgap = float(s[K - fdim - 1] / s[K - fdim]) if (0 < K - fdim <= len(s) - 1 and fdim and s[K - fdim] > 0) else float('inf')
    sval_min_nonzero = float(s[K - fdim - 1]) if K - fdim - 1 >= 0 and K - fdim - 1 < len(s) else None
    vecs_arr = np.zeros((len(vecs), K), dtype=np.int64)
    for r, v in enumerate(vecs):
        for e, c in v.items():
            vecs_arr[r, didx[e]] = c
    tmp = f'{base}.{os.getpid()}.tmp.npz'
    np.savez_compressed(tmp, Kb=Kb, mons=np.array(Dm, dtype=np.int16), vecs_modp=vecs_arr, p=P)
    os.replace(tmp, base + '.npz')
    def sym(v):
        return {core.to_str({e: Fraction(c if c <= P // 2 else c - P)}) : None for e, c in v.items()}
    info = {'event': 'e1space', 'd': d, 'sec': sec_str(sec), 'dimS': K, 'dimJ': jrank, 'dimR': K - jrank,
            'kernel_dim': kernel_dim, 'nonJ': len(basis), 'target_rows': R, 'target_cols': C,
            'float_gap': sgap, 'sval_min_nonzero': sval_min_nonzero,
            'basis_modp': [{'|'.join(f'{core.VARS[i]}^{k}' for i, k in enumerate(e) if k) or '1':
                            (c if c <= P // 2 else c - P) for e, c in v.items()} for v in basis],
            't_modp': round(t1 - t0, 1), 't_total': round(time.time() - t0, 1), 'p': P}
    json.dump(info, open(base + '.json', 'w'))
    log_event({k: v for k, v in info.items() if k != 'basis_modp'} | {'nonJ_basis_terms': [len(v) for v in basis]})
    if verbose:
        print('[e1space]', {k: v for k, v in info.items() if k != 'basis_modp'}, flush=True)
    return info


def e1_member_modp(D):
    """Exact mod-p test: is D (pure Klein sector) in the E1-admissible space?  Uses the saved rref kernel basis."""
    secs = {chi(e) for e in D}
    d = core.degree(D)
    if len(secs) != 1:
        # mixed Klein character: test each component separately (E1 is linear and sector-graded)
        comps = {}
        for e, c in D.items():
            comps.setdefault(chi(e), {})[e] = c
        return all(e1_member_modp(c) for c in comps.values())
    sec = secs.pop()
    z = np.load(e1_path(d, sec) + '.npz')
    mons = [tuple(int(x) for x in row) for row in z['mons']]
    idx = {e: i for i, e in enumerate(mons)}
    V = z['vecs_modp']
    Dp = core.reduce_mod(D, P)
    row = [0] * len(mons)
    for e, c in Dp.items():
        row[idx[e]] = c
    if V.shape[0] == 0:
        return not any(row)
    flat = [int(x) for x in V.ravel()] + row
    M = nmod_mat(V.shape[0] + 1, len(mons), flat, P)
    return M.rank() == V.shape[0]


def load_e1_basis(d, sec):
    z = np.load(e1_path(d, sec) + '.npz')
    mons = [tuple(int(x) for x in row) for row in z['mons']]
    return z['Kb'], mons


# ----------------------------------------------------------------------------
# g-bases and relaxation tensors
# ----------------------------------------------------------------------------

def g_basis(d, sec):
    """Basis (i, m) of {g : sum g_i F_i in S_d^sec}: deg m = d-5, chi(m) = sec + chi(z_i)."""
    out = []
    if d < 5:
        return out
    for i in range(NV):
        need = sadd(sec, chi_var(i))
        for m in sector_monomials(d - 5, need):
            out.append((i, m))
    return out


def vec_in(poly, idx, n):
    v = np.zeros(n)
    for e, c in poly.items():
        v[idx[e]] += float(c)
    return v


class Relaxation:
    """ob(g) = c0 + L g + g^T T g in float cokernel coordinates; D(g) = D0 + sum g_i F_i."""

    def __init__(self, d, dsecs, D0=None, verbose=True):
        F, A, q, Fi, Hs = instance()
        self.d = d
        self.dsecs = list(dsecs)
        self.D0 = D0
        gb = []
        for s in self.dsecs:
            gb += g_basis(d, s)
        self.gb = gb
        K = len(gb); self.K = K
        # output sectors: chi(A) + sec_k + sec_l ; collect needed cokernels at degree 2d
        outsecs = sorted({sadd((0, 1), sadd(sk, sl)) for sk in self.dsecs for sl in self.dsecs}
                         | ({sadd((0, 1), s) for s in self.dsecs} if D0 is not None else set()))
        self.outsecs = outsecs
        self.cok = {}
        off = 0
        self.offsets = {}
        for s in outsecs:
            B, mons, rank = prep_coker(2 * d, s, verbose=verbose)
            idx = {e: i for i, e in enumerate(mons)}
            self.cok[s] = (B, mons, idx)
            self.offsets[s] = off
            off += B.shape[0]
        self.m = off
        t0 = time.time()
        AH = [[mul(A, Hs[i][j]) for j in range(NV)] for i in range(NV)]
        T = np.zeros((K, K, off))
        for k, (i, mk) in enumerate(gb):
            for l in range(k, K):
                j, ml = gb[l]
                if not AH[i][j]:
                    continue
                poly = shift(AH[i][j], tuple(a + b for a, b in zip(mk, ml)), 1)
                s = chi(next(iter(poly)))
                B, mons, idx = self.cok[s]
                v = -(B @ vec_in(poly, idx, len(mons)))   # ordered pair (k,l): coefficient of g_k g_l in -[A g^T H g]
                o = self.offsets[s]
                T[k, l, o:o + len(v)] += v
                if k != l:
                    T[l, k, o:o + len(v)] += v
        self.T = T
        self.c0 = np.zeros(off)
        self.L = np.zeros((K, off))
        if D0 is not None:
            self._affine_part(D0)
        self.scale = float(np.sqrt(np.mean(np.sum(T ** 2, axis=2)))) or 1.0
        # spurious subspace Z' = {g : D(g) in F*S_{d-6}} (contains the syzygies); ob is invariant under Z'
        self.P, self.dimZ = self._complement_of_spurious()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        if verbose:
            print(f'[relax] d={d} secs={[sec_str(s) for s in self.dsecs]} K={K} m={off} affine={D0 is not None} '
                  f'spurious dim={self.dimZ} search dim={self.dimP} built in {self.build_time:.1f}s', flush=True)

    def _complement_of_spurious(self):
        F, A, q, Fi, Hs = instance()
        K = self.K
        dmons = {}
        cols = []
        for (i, m) in self.gb:
            cols.append(shift(Fi[i], m, 1))
        Fcols = []
        if self.d >= 6:
            for s in self.dsecs:
                Fcols += [shift(F, h, 1) for h in sector_monomials(self.d - 6, s)]
        for Cc in cols + Fcols:
            for e in Cc:
                dmons.setdefault(e, len(dmons))
        Mb = np.zeros((len(dmons), K + len(Fcols)))
        for k, Cc in enumerate(cols):
            for e, v in Cc.items():
                Mb[dmons[e], k] = float(v)
        for j, Cc in enumerate(Fcols):
            for e, v in Cc.items():
                Mb[dmons[e], K + j] = -float(v)
        ns = sla.null_space(Mb)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return np.ascontiguousarray(U[:, r:]), r

    def _affine_part(self, D0):
        F, A, q, Fi, Hs = instance()
        d = self.d
        # E1 particular solution N0 (float least squares, sector restricted)
        sD = chi(next(iter(D0)))
        Ncols = []
        for i in range(NV):
            s = sadd(sadd((0, 1), sD), chi_var(i))
            Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
        tsec = sadd((0, 1), sD)
        tm = sector_monomials(d + 6, tsec)
        tidx = {e: i for i, e in enumerate(tm)}
        M = np.zeros((len(tm), len(Ncols)))
        for j, (i, m) in enumerate(Ncols):
            for e, v in shift(Fi[i], m, 1).items():
                M[tidx[e], j] = float(v)
        b = vec_in(mul(A, D0), tidx, len(tm))
        x, *_ = sla.lstsq(M, b, cond=1e-12, lapack_driver='gelsd', check_finite=False)
        self.e1_residual = float(np.linalg.norm(M @ x - b) / np.linalg.norm(b))
        N0 = [dict() for _ in range(NV)]
        for j, (i, m) in enumerate(Ncols):
            if abs(x[j]) > 1e-14:
                N0[i][m] = float(x[j])
        D0f = {e: float(v) for e, v in D0.items()}
        Af = {e: float(v) for e, v in A.items()}
        Fif = [{e: float(v) for e, v in f.items()} for f in Fi]
        Hsf = [[{e: float(v) for e, v in h.items()} for h in row] for row in Hs]

        def cls(poly):
            out = np.zeros(self.m)
            bysec = {}
            for e, c in poly.items():
                bysec.setdefault(chi(e), {})[e] = c
            for s, pp in bysec.items():
                if s not in self.cok:
                    raise RuntimeError(f'unexpected output sector {s}')
                B, mons, idx = self.cok[s]
                v = B @ vec_in(pp, idx, len(mons))
                o = self.offsets[s]
                out[o:o + len(v)] += v
            return out
        divN = {}
        for i in range(NV):
            divN = add(divN, pderiv(N0[i], i))
        ND = {}
        for i in range(NV):
            ND = add(ND, mul(N0[i], pderiv(D0f, i)))
        self.c0 = cls(sub(mul(D0f, divN), ND))
        for k, (i, mk) in enumerate(self.gb):
            g = {mk: 1.0}
            # D0 (g(A) + A div g) - A g d_i D0 - m_k sum_j N0_j F_ij
            term = mul(D0f, add(mul(g, pderiv(Af, i)), mul(Af, pderiv(g, i))))
            term = sub(term, mul(mul(Af, g), pderiv(D0f, i)))
            acc = {}
            for j in range(NV):
                if N0[j] and Hsf[i][j]:
                    acc = add(acc, mul(N0[j], Hsf[i][j]))
            term = sub(term, mul(g, acc))
            self.L[k] = cls(term)

    # residual in normalised form
    def ob(self, g):
        return self.c0 + g @ self.L + np.einsum('k,kl->l', g, np.einsum('klm,l->km', self.T, g))

    def ob_and_jac(self, g):
        Tg = np.einsum('klm,l->km', self.T, g)        # K x m
        ob = self.c0 + g @ self.L + g @ Tg
        J = self.L + 2.0 * Tg                           # K x m  (d ob / d g_k)
        return ob, J

    def residual_u(self, u):
        g = self.P @ u
        r, J = self.residual_fun(g)
        return r, J @ self.P

    def residual_fun(self, g):
        if self.D0 is None:
            n2 = float(g @ g)
            ob, J = self.ob_and_jac(g)
            r = ob / n2
            Jr = J / n2 - 2.0 * np.outer(g, ob) / (n2 * n2)
            # extra residual keeps g on the unit sphere (the objective is scale invariant)
            r = np.concatenate([r / self.scale, [n2 - 1.0]])
            Jf = np.vstack([Jr.T / self.scale, 2.0 * g[None, :]])
            return r, Jf
        ob, J = self.ob_and_jac(g)
        return ob / self.scale, J.T / self.scale

    def search(self, restarts, seed, maxfev=2000, verbose=False, stop_rel=1e-13):
        rng = np.random.default_rng(seed)
        results = []
        t0 = time.time()
        for it in range(restarts):
            u0 = rng.standard_normal(self.dimP)
            if self.D0 is None:
                u0 /= np.linalg.norm(u0)
            method = 'lm' if self.m + 1 >= self.dimP else 'trf'
            res = least_squares(lambda u: self.residual_u(u)[0], u0,
                                jac=lambda u: self.residual_u(u)[1], method=method,
                                xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxfev)
            g = self.P @ res.x
            if self.D0 is None:
                g = g / np.linalg.norm(g)
            r = float(np.linalg.norm(self.ob(g) / (g @ g if self.D0 is None else 1.0))) / self.scale
            results.append({'rel': r, 'g': g, 'nfev': int(res.nfev), 'restart': it})
            if verbose and (it < 5 or it % 50 == 0):
                print(f'[relax] restart {it}: rel={r:.3e} nfev={res.nfev} ({time.time()-t0:.1f}s)', flush=True)
            if r < stop_rel:
                break
        results.sort(key=lambda c: c['rel'])
        return results

    def D_of(self, g, exact=False):
        F, A, q, Fi, Hs = instance()
        D = dict(self.D0) if (self.D0 is not None and exact) else ({e: float(v) for e, v in self.D0.items()} if self.D0 is not None else {})
        for k, (i, m) in enumerate(self.gb):
            if g[k]:
                D = add(D, shift(Fi[i], m, g[k]))
        return D


def reconstruct_vector(g, max_dens=(1, 2, 3, 4, 6, 12, 24, 60, 120, 360, 720, 2520, 5040, 10000, 100000, 1000000), tol=1e-7, normalize=True):
    """Rational candidates for a float vector (optionally normalised by its largest entry)."""
    g = np.asarray(g, float)
    if normalize:
        j = int(np.argmax(np.abs(g)))
        if g[j] == 0:
            return []
        g = g / g[j]
    cands = []
    for md in max_dens:
        fr = [Fraction(float(v)).limit_denominator(md) for v in g]
        err = max(abs(float(f) - v) for f, v in zip(fr, g))
        if err < tol:
            t = tuple(fr)
            if t not in cands:
                cands.append(t)
    return cands


def exact_confirm(D, note='', symmetric=True, verbose=True):
    """Full exact pipeline for one denominator: divisibility, E1, fixed-D solve (sector restricted), verify."""
    F, A, q, Fi, Hs = instance()
    t0 = time.time()
    rec = {'event': 'confirm', 'note': note, 'd': core.degree(D), 'D': core.to_str(D)}
    if not D:
        rec['status'] = 'zero'; log_event(rec); return rec
    divF, _ = jacobian.divisible_by(F, D)
    rec['F_divides_D'] = bool(divF)
    if divF:
        rec['status'] = 'spurious (F | D)'; log_event(rec); return rec
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    rec['e1'] = bool(ok)
    if not ok:
        rec['status'] = 'e1 fails'; rec['seconds'] = round(time.time() - t0, 1); log_event(rec); return rec
    if symmetric and core.parity_of(D) is not None and len({chi(e) for e in D}) == 1:
        status, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, verbose=verbose, compressed=True)
    else:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
        status = 'inconsistent' if not out['consistent'] else ('solved' if out.get('exact') else 'lift-failed')
        data = out
    rec['status'] = status
    rec.update({k: data.get(k) for k in ('rank', 'nullity', 'rows', 'cols')})
    rec['seconds'] = round(time.time() - t0, 1)
    if status == 'solved':
        N, Q = data['N'], data['Q']
        ok2, info2 = core.verify_witness(F, A, q, D, N, Q)
        rec['verified'] = bool(ok2)
        if ok2:
            w = {'instance': LANE, 'd': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in N],
                 'Q': core.to_str(Q), 'note': note}
            h = hashlib.sha1(w['D'].encode()).hexdigest()[:10]
            path = os.path.join(HERE, f'witness-candidate-{h}.json')
            json.dump(w, open(path, 'w'), indent=1)
            if INST_ID == LANE:
                json.dump(w, open(os.path.join(HERE, 'witness-candidate.json'), 'w'), indent=1)
                print('!!!!!!!! WITNESS FOUND', path, flush=True)
            else:
                print(f'[control {INST_ID}] exact witness found: {path}', flush=True)
            rec['witness_file'] = path
    log_event(rec)
    if verbose:
        print('[confirm]', {k: v for k, v in rec.items() if k != 'D'}, flush=True)
    return rec


def exact_relaxation_check(D, g_exact_polys=None):
    """Exact mod-p test of the relaxation for D = sum g_i F_i: A g^T H g in J_{2d}."""
    F, A, q, Fi, Hs = instance()
    acc = {}
    for i in range(NV):
        for j in range(NV):
            if g_exact_polys[i] and g_exact_polys[j] and Hs[i][j]:
                acc = add(acc, mul(mul(g_exact_polys[i], g_exact_polys[j]), Hs[i][j]))
    if not acc:
        return True, {'note': 'zero'}
    return jacobian.in_jacobian(F, mul(A, acc))


# ----------------------------------------------------------------------------
# sector-restricted variable projection (full system)
# ----------------------------------------------------------------------------

class SectorVarPro:
    """Variable projection over D(a) = sum a_k M_k with all M_k in one Klein sector (unknowns restricted)."""

    def __init__(self, basis, verbose=True, cond=1e-10):
        F, A, q, Fi, Hs = instance()
        self.F, self.A, self.q = F, A, q
        self.basis = basis
        K = len(basis); self.K = K
        degs = {core.degree(B) for B in basis}
        assert len(degs) == 1
        d = degs.pop(); self.d = d
        secs = {chi(e) for B in basis for e in B}
        assert len(secs) == 1, f'basis spans several sectors {secs}'
        sD = secs.pop(); self.sD = sD
        degN, degQ = d + 1, 2 * d - 6
        Ncols = []
        for i in range(NV):
            s = sadd(sadd((0, 1), sD), chi_var(i))
            Ncols += [(i, m) for m in sector_monomials(degN, s)]
        Qcols = sector_monomials(degQ, (0, 1)) if degQ >= 0 else []
        self.Ncols, self.Qcols = Ncols, Qcols
        self.nN, self.nQ = len(Ncols), len(Qcols)
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        Fq1 = scale(power(F, 2), -1)
        blocks = []
        for Bk in basis:
            FBk = mul(F, Bk)
            G = [scale(add(mul(Bk, Fi[i]), mul(F, pderiv(Bk, i))), -1) for i in range(NV)]
            ri, ci, vals = [], [], []
            for j, (i, m) in enumerate(Ncols):
                col = shift(G[i], m, 1)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FBk, tuple(m2), m[i]))
                for e, v in col.items():
                    ri.append(ridx(e)); ci.append(j); vals.append(float(v))
            blocks.append((ri, ci, vals))
        qri, qci, qvals = [], [], []
        for j, m in enumerate(Qcols):
            for e, v in shift(Fq1, m, 1).items():
                qri.append(ridx(e)); qci.append(j); qvals.append(float(v))
        AB = [mul(A, Bk) for Bk in basis]
        pair = {}
        for k in range(K):
            for l in range(k, K):
                v = scale(mul(AB[k], basis[l]), -1)
                pair[(k, l)] = v
                for e in v:
                    ridx(e)
        R = len(rows); self.rows = R
        import scipy.sparse as sp
        self.Mk = [sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                 shape=(R, self.nN)) for ri, ci, vals in blocks]
        self.MQ = sp.csr_matrix((np.array(qvals), (np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64))),
                                shape=(R, self.nQ)).toarray() if self.nQ else np.zeros((R, 0))
        self.Vk = []
        for k in range(K):
            ri, ci, vals = [], [], []
            for l in range(K):
                v = pair[(min(k, l), max(k, l))]
                for e, c in v.items():
                    ri.append(rows[e]); ci.append(l); vals.append(float(c))
            self.Vk.append(sp.csr_matrix((np.array(vals), (np.array(ri, dtype=np.int64), np.array(ci, dtype=np.int64))),
                                         shape=(R, K)))
        self.ncols = self.nN + self.nQ
        self.cond = cond
        # exclude F-multiples
        self.P = np.eye(K); self.excluded = 0
        if d >= 6:
            Gm = sector_monomials(d - 6, sD)
            Fcols = [shift(F, m, 1) for m in Gm]
            mons = {}
            for B in basis + Fcols:
                for e in B:
                    mons.setdefault(e, len(mons))
            Mb = np.zeros((len(mons), K + len(Fcols)))
            for k, B in enumerate(basis):
                for e, v in B.items():
                    Mb[mons[e], k] = float(v)
            for j, C in enumerate(Fcols):
                for e, v in C.items():
                    Mb[mons[e], K + j] = -float(v)
            ns = sla.null_space(Mb)
            S = ns[:K, :] if ns.size else np.zeros((K, 0))
            if S.shape[1]:
                U, s, _ = np.linalg.svd(S, full_matrices=True)
                r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
                self.P, self.excluded = U[:, r:], r
        self.dimP = self.P.shape[1]
        self._nfev = 0
        if verbose:
            print(f'[svarpro] d={d} sec={sec_str(sD)} K={K} rows={R} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'excluded F-multiples={self.excluded}', flush=True)

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

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r, float(b @ b)

    def value_and_grad(self, u):
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        f = float(r @ r)
        xN = x[:self.nN]
        g = np.zeros(self.K)
        for k in range(self.K):
            g[k] = 2.0 * float(r @ (self.Mk[k] @ xN - 2.0 * (self.Vk[k] @ a)))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def search(self, n_restarts, seed, maxiter=400, verbose=False, stop_rel=1e-13, init=None):
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
                           options={'maxiter': maxiter, 'gtol': 1e-12, 'ftol': 1e-16})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r, bb = self.value(a)
            rel = f / max(bb, 1e-300)
            out.append({'value': f, 'rel': rel, 'a': a, 'iters': int(res.nit), 'restart': it})
            if verbose and (it < 5 or it % 20 == 0):
                print(f'[svarpro] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)', flush=True)
            if rel < stop_rel:
                break
        out.sort(key=lambda c: c['rel'])
        return out

    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.basis):
            if frac[k]:
                D = add(D, scale(B, Fraction(frac[k])))
        return D


def family_basis(d, sec, family):
    """Exact (rational) basis polynomials of the requested family in one sector."""
    F, A, q, Fi, Hs = instance()
    if family == 'J':
        return [shift(Fi[i], m, 1) for (i, m) in g_basis(d, sec)]
    if family == 'K8J':
        assert d == 8 and sec == (0, 0)
        K8 = core.parse(json.load(open(os.path.join(HERE, 'kernel_d8_X03Y03.json')))[0])
        return [K8] + [shift(Fi[i], m, 1) for (i, m) in g_basis(d, sec)]
    if family == 'E1':
        # float orthonormal basis of the E1-admissible subspace (not rational); used only numerically
        Kb, mons = load_e1_basis(d, sec)
        return [{m: float(Kb[r, j]) for j, m in enumerate(mons) if abs(Kb[r, j]) > 1e-15} for r in range(Kb.shape[0])]
    if family == 'mon':
        return [{m: Fraction(1)} for m in sector_monomials(d, sec)]
    raise ValueError(family)


# ----------------------------------------------------------------------------
# msolve export of the relaxation quadrics (exact mod p)
# ----------------------------------------------------------------------------

def relaxation_quadrics_modp(d, dsecs):
    """Exact quadrics mod p: coordinates of [A g^T H g] in standard monomials of R_{2d}^{out sectors}."""
    F, A, q, Fi, Hs = instance()
    gb = []
    for s in dsecs:
        gb += g_basis(d, s)
    K = len(gb)
    Ap = core.reduce_mod(A, P)
    Hp = [[core.reduce_mod(h, P) for h in row] for row in Hs]
    outsecs = sorted({sadd((0, 1), sadd(sk, sl)) for sk in dsecs for sl in dsecs})
    reducers = {}
    for s in outsecs:
        mons = sector_monomials(2 * d, s)
        gens = gens_in_sector(2 * d, s, p=P)
        idx = {e: i for i, e in enumerate(mons)}
        flat = [0] * (len(gens) * len(mons))
        for r, (_, _, g) in enumerate(gens):
            for e, v in g.items():
                flat[r * len(mons) + idx[e]] = v
        Jr, rank = nmod_mat(len(gens), len(mons), flat, P).rref()
        pivots, rows = [], []
        for i in range(rank):
            j = 0
            while j < len(mons) and int(Jr[i, j]) == 0:
                j += 1
            pivots.append(j)
            rows.append({mons[k]: int(Jr[i, k]) for k in range(j, len(mons)) if int(Jr[i, k])})
        standard = [m for j, m in enumerate(mons) if j not in set(pivots)]
        reducers[s] = (mons, pivots, rows, standard)

    def nf(poly, s):
        mons, pivots, rows, standard = reducers[s]
        Pd = dict(poly)
        for piv, row in zip(pivots, rows):
            c = Pd.get(mons[piv], 0)
            if c:
                Pd = sub(Pd, scale(row, c, P), P)
        return Pd
    # quadrics: for each standard monomial coordinate, a quadratic form in g
    coords = {}  # (sector, standard monomial) -> dict {(k,l): coeff}
    for k, (i, mk) in enumerate(gb):
        for l in range(k, K):
            j, ml = gb[l]
            if not Hp[i][j]:
                continue
            poly = shift(mul(Ap, Hp[i][j], P), tuple(a + b for a, b in zip(mk, ml)), 1, P)
            s = chi(next(iter(poly)))
            red = nf(poly, s)
            for e, c in red.items():
                coords.setdefault((s, e), {})[(k, l)] = c
    return gb, coords


def spurious_pivots_modp(d, dsecs, gb):
    """Pivot coordinates (rref mod p) of Z' = {g : D(g) in F*S_{d-6}}: z*h/6 (h monomial) and Koszul syzygies."""
    F, A, q, Fi, Hs = instance()
    gidx = {g: k for k, g in enumerate(gb)}
    K = len(gb)
    vecs = []
    if d >= 6:
        for s in dsecs:
            for h in sector_monomials(d - 6, s):
                v = {}
                for i in range(NV):
                    m = list(h); m[i] += 1; m = tuple(m)
                    if (i, m) in gidx:
                        v[gidx[(i, m)]] = 1
                if len(v) == NV:
                    vecs.append(v)
    if d >= 10:
        Fip = [core.reduce_mod(f, P) for f in Fi]
        for i in range(NV):
            for j in range(i + 1, NV):
                # syzygy: F_j e_i - F_i e_j, times monomials h of degree d-10
                for h in monomials(d - 10):
                    v = {}
                    ok = True
                    for (slot, poly, sgn) in ((i, Fip[j], 1), (j, Fip[i], P - 1)):
                        for e, c in shift(poly, h, 1, P).items():
                            key = (slot, e)
                            if key not in gidx:
                                ok = False; break
                            v[gidx[key]] = (v.get(gidx[key], 0) + sgn * c) % P
                        if not ok:
                            break
                    if ok and v:
                        vecs.append(v)
    if not vecs:
        return []
    flat = [0] * (len(vecs) * K)
    for r, v in enumerate(vecs):
        for k, c in v.items():
            flat[r * K + k] = c % P
    Zr, rank = nmod_mat(len(vecs), K, flat, P).rref()
    piv = []
    for i in range(rank):
        j = 0
        while j < K and int(Zr[i, j]) == 0:
            j += 1
        piv.append(j)
    return piv


def write_msolve(d, dsecs, path):
    gb, coords = relaxation_quadrics_modp(d, dsecs)
    K = len(gb)
    names = [f'g{k}' for k in range(K)]
    lines = [','.join(names), str(P)]
    polys = []
    piv = spurious_pivots_modp(d, dsecs, gb)
    for j in piv:
        polys.append(names[j])
    for (s, e), q in coords.items():
        terms = []
        for (k, l), c in q.items():
            terms.append(f'{c}*{names[k]}*{names[l]}' if k != l else f'{c}*{names[k]}^2')
        polys.append('+'.join(terms))
    with open(path, 'w') as f:
        f.write(lines[0] + '\n' + lines[1] + '\n' + ',\n'.join(polys) + '\n')
    return K, len(polys)


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def cmd_prep(a):
    for deg in a.deg:
        for sec in sectors_of(parse_sec(a.sec)):
            prep_coker(deg, sec, force=a.force)


def cmd_e1(a):
    for d in a.d:
        for sec in sectors_of(parse_sec(a.sec)):
            e1_space(d, sec, force=a.force)


def cmd_relax(a):
    d = a.d
    secspec = parse_sec(a.sec)
    dsecs = sectors_of(secspec)
    D0 = None
    if a.affine:
        path, idx = a.affine.split(':')
        D0 = core.parse(json.load(open(os.path.join(HERE, path)))[int(idx)])
        dsecs = [chi(next(iter(D0)))]
    rel = Relaxation(d, dsecs, D0=D0)
    t0 = time.time()
    res = rel.search(a.restarts, a.seed, verbose=True, maxfev=a.maxfev)
    rels = np.array([c['rel'] for c in res])
    hist = {f'<1e-{k}': int(np.sum(rels < 10.0 ** (-k))) for k in (2, 4, 6, 8, 10, 12)}
    rec = {'event': 'relax', 'd': d, 'secs': [sec_str(s) for s in dsecs], 'affine': a.affine, 'K': rel.K,
           'dimZ': rel.dimZ, 'dimP': rel.dimP, 'm': rel.m, 'restarts': len(res), 'seed': a.seed, 'best_rel': float(rels.min()),
           'median_rel': float(np.median(rels)), 'hist': hist, 'seconds': round(time.time() - t0, 1),
           'e1_residual': getattr(rel, 'e1_residual', None)}
    print('[relax] summary', rec, flush=True)
    log_event(rec)
    best = res[:a.keep]
    json.dump([{'rel': c['rel'], 'g': [float(x) for x in c['g']]} for c in best],
              open(os.path.join(HERE, f'relax_d{d}_{a.sec}_{"aff" if a.affine else "hom"}_s{a.seed}.json'), 'w'))
    # candidates below threshold: reconstruct + exact
    for c in res:
        if c['rel'] > a.thresh:
            break
        g = c['g']
        for fr in reconstruct_vector(g, normalize=(D0 is None)):
            gpolys = [dict() for _ in range(NV)]
            for k, (i, m) in enumerate(rel.gb):
                if fr[k]:
                    gpolys[i][m] = Fraction(fr[k])
            if D0 is None:
                ok, info = exact_relaxation_check(None, gpolys)
                print(f'[relax] exact relaxation check of reconstruction: {ok}', flush=True)
                log_event({'event': 'relax-reconstruct', 'd': d, 'rel': c['rel'], 'exact_relaxation': bool(ok),
                           'g': {f'{i}:{m}': str(fr[k]) for k, (i, m) in enumerate(rel.gb) if fr[k]}})
                if not ok:
                    continue
            D = rel.D_of([Fraction(x) for x in fr], exact=True)
            exact_confirm(D, note=f'relax d={d} secs={a.sec} rel={c["rel"]:.2e}')


def cmd_msolve(a):
    d = a.d
    dsecs = sectors_of(parse_sec(a.sec))
    path = os.path.join(HERE, f'relax_d{d}_{a.sec}.ms')
    K, n = write_msolve(d, dsecs, path)
    print(f'[msolve] d={d} secs={a.sec}: {n} polys (quadrics + spurious pivot kills) in {K} unknowns -> {path}', flush=True)
    log_event({'event': 'msolve-export', 'd': d, 'secs': a.sec, 'K': K, 'nquadrics': n, 'path': path})
    if a.run:
        # dehomogenise with a random affine chart: add sum r_k g_k - 1
        import random
        random.seed(a.seed)
        txt = open(path).read().rstrip('\n')
        names = txt.split('\n')[0].split(',')
        chart = '+'.join(f'{random.randrange(1, P)}*{nm}' for nm in names) + f'+{P-1}'
        path2 = path.replace('.ms', '_chart.ms')
        open(path2, 'w').write(txt + ',\n' + chart + '\n')
        t0 = time.time()
        out = path2.replace('.ms', '.out')
        r = subprocess.run(['msolve', '-v', '1', '-t', str(a.threads), '-f', path2, '-o', out], capture_output=True, text=True)
        res = open(out).read() if os.path.exists(out) else ''
        empty = res.strip().startswith('[-1]')
        rec = {'event': 'msolve', 'd': d, 'secs': a.sec, 'K': K, 'nquadrics': n, 'empty_projective': bool(empty),
               'out_head': res[:300], 'stderr_tail': r.stderr[-500:], 'seconds': round(time.time() - t0, 1)}
        print('[msolve]', rec, flush=True)
        log_event(rec)


def cmd_varpro(a):
    d = a.d
    sec = parse_sec(a.sec)
    basis = family_basis(d, sec, a.family)
    if a.family == 'E1':
        bb = [{e: Fraction(v).limit_denominator(10 ** 6) for e, v in B.items()} for B in basis]
        basis = bb
    fam = SectorVarPro(basis)
    t0 = time.time()
    res = fam.search(a.restarts, a.seed, maxiter=a.maxiter, verbose=True)
    rels = np.array([c['rel'] for c in res])
    hist = {f'<1e-{k}': int(np.sum(rels < 10.0 ** (-k))) for k in (2, 4, 6, 8, 10, 12)}
    rec = {'event': 'varpro', 'd': d, 'sec': a.sec, 'family': a.family, 'K': fam.K, 'rows': fam.rows,
           'cols': fam.ncols, 'restarts': len(res), 'seed': a.seed, 'best_rel': float(rels.min()),
           'median_rel': float(np.median(rels)), 'hist': hist, 'seconds': round(time.time() - t0, 1)}
    print('[varpro] summary', rec, flush=True)
    log_event(rec)
    json.dump([{'rel': c['rel'], 'a': [float(x) for x in c['a']]} for c in res[:a.keep]],
              open(os.path.join(HERE, f'varpro_d{d}_{a.sec}_{a.family}_s{a.seed}.json'), 'w'))
    for c in res:
        if c['rel'] > a.thresh:
            break
        for fr in reconstruct_vector(c['a']):
            D = fam.D_exact(fr)
            exact_confirm(D, note=f'varpro d={d} sec={a.sec} fam={a.family} rel={c["rel"]:.2e}')


def cmd_screen(a):
    """E1-screen a candidate file (lines: tag<TAB>expr) against the saved E1 spaces; confirm survivors exactly."""
    F, A, q, Fi, Hs = instance()
    for line in open(a.file):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        tag, expr = line.split('\t', 1) if '\t' in line else ('cand', line)
        try:
            D = core.parse(expr)
        except Exception as ex:
            log_event({'event': 'screen', 'tag': tag, 'error': str(ex)}); continue
        d = core.degree(D)
        if not D or not core.is_homogeneous(D) or core.parity_of(D) is None:
            log_event({'event': 'screen', 'tag': tag, 'd': d, 'note': 'not homogeneous / mixed Y-parity'}); continue
        t0 = time.time()
        try:
            ok = e1_member_modp(D)
        except FileNotFoundError:
            ok2, info = fixed_d.e1_solvable(F, A, q, D)
            ok = bool(ok2)
        rec = {'event': 'screen', 'tag': tag, 'd': d, 'terms': len(D), 'e1': bool(ok), 'seconds': round(time.time() - t0, 1),
               'D': core.to_str(D)[:200]}
        print('[screen]', rec, flush=True)
        log_event(rec)
        if ok:
            exact_confirm(D, note=f'screen {tag}')


def cmd_confirm(a):
    D = core.parse(a.D)
    exact_confirm(D, note=a.note, symmetric=not a.full)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('prep'); p.add_argument('--deg', type=int, nargs='+', required=True); p.add_argument('--sec', required=True); p.add_argument('--force', action='store_true')
    p = sp.add_parser('e1'); p.add_argument('--d', type=int, nargs='+', required=True); p.add_argument('--sec', default='all'); p.add_argument('--force', action='store_true')
    p = sp.add_parser('relax'); p.add_argument('--d', type=int, required=True); p.add_argument('--sec', required=True)
    p.add_argument('--affine', default=None); p.add_argument('--restarts', type=int, default=100); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--maxfev', type=int, default=2000); p.add_argument('--thresh', type=float, default=1e-9); p.add_argument('--keep', type=int, default=20)
    p = sp.add_parser('msolve'); p.add_argument('--d', type=int, required=True); p.add_argument('--sec', required=True); p.add_argument('--run', action='store_true'); p.add_argument('--threads', type=int, default=8); p.add_argument('--seed', type=int, default=1)
    p = sp.add_parser('varpro'); p.add_argument('--d', type=int, required=True); p.add_argument('--sec', required=True); p.add_argument('--family', default='J')
    p.add_argument('--restarts', type=int, default=20); p.add_argument('--seed', type=int, default=0); p.add_argument('--maxiter', type=int, default=400)
    p.add_argument('--thresh', type=float, default=1e-9); p.add_argument('--keep', type=int, default=20)
    p = sp.add_parser('screen'); p.add_argument('--file', required=True)
    p = sp.add_parser('confirm'); p.add_argument('--D', required=True); p.add_argument('--note', default=''); p.add_argument('--full', action='store_true')
    a = ap.parse_args(argv)
    {'prep': cmd_prep, 'e1': cmd_e1, 'relax': cmd_relax, 'msolve': cmd_msolve, 'varpro': cmd_varpro, 'confirm': cmd_confirm, 'screen': cmd_screen}[a.cmd](a)




# ----------------------------------------------------------------------------
# Generalised E1-family relaxation: pairs (D_k, N_k) with N_k an E1 solution for D_k;
# ob(c) = [B(D(c), N(c))],  B(D, N) = D div N - N(D), bilinear  ->  tensor T[k,l] = Pi B(D_k, N_l)
# ----------------------------------------------------------------------------

def _pack(E):
    """int64 keys of exponent rows (base 64)."""
    E = np.asarray(E, dtype=np.int64)
    return E @ (64 ** np.arange(6, dtype=np.int64))


def jacobian_pairs(d, sec):
    """Exact pairs (m F_i, A m e_i) spanning J_d^sec."""
    F, A, q, Fi, Hs = instance()
    pairs = []
    for (i, m) in g_basis(d, sec):
        N = [dict() for _ in range(NV)]
        N[i] = shift(A, m, 1)
        pairs.append((shift(Fi[i], m, 1), N, f'J:{i}:{m}'))
    return pairs


def e1_particular_solutions(Ds, sec, verbose=False):
    """Float least-squares E1 solutions N (sector restricted) for a list of float D's of one sector."""
    F, A, q, Fi, Hs = instance()
    d = core.degree(Ds[0])
    Ncols = []
    for i in range(NV):
        s = sadd(sadd((0, 1), sec), chi_var(i))
        Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
    tsec = sadd((0, 1), sec)
    tm = sector_monomials(d + 6, tsec)
    tidx = {e: i for i, e in enumerate(tm)}
    M = np.zeros((len(tm), len(Ncols)))
    for j, (i, m) in enumerate(Ncols):
        for e, v in shift(Fi[i], m, 1).items():
            M[tidx[e], j] = float(v)
    Bm = np.zeros((len(tm), len(Ds)))
    for c, D in enumerate(Ds):
        for e, v in mul(A, D).items():
            Bm[tidx[e], c] += float(v)
    X, *_ = sla.lstsq(M, Bm, cond=1e-12, lapack_driver='gelsd', check_finite=False)
    resid = np.linalg.norm(M @ X - Bm, axis=0) / np.maximum(np.linalg.norm(Bm, axis=0), 1e-300)
    Ns = []
    for c in range(len(Ds)):
        N = [dict() for _ in range(NV)]
        for j, (i, m) in enumerate(Ncols):
            if abs(X[j, c]) > 1e-14:
                N[i][m] = float(X[j, c])
        Ns.append(N)
    if verbose:
        print(f'[e1sol] {len(Ds)} systems {M.shape}: max rel residual {resid.max():.2e}', flush=True)
    return Ns, resid


def kernel_pairs(d, sec, verbose=True):
    """Float pairs for a complement of J_d^sec inside the E1-admissible space (the non-Jacobian classes)."""
    F, A, q, Fi, Hs = instance()
    info = e1_space(d, sec, verbose=verbose)
    Kb, mons = load_e1_basis(d, sec)
    if info['nonJ'] == 0 or Kb.shape[0] == 0:
        return []
    idx = {e: i for i, e in enumerate(mons)}
    # J_d^sec in the same coordinates
    gens = gens_in_sector(d, sec)
    G = np.zeros((len(mons), len(gens)))
    for j, (_, _, g) in enumerate(gens):
        for e, v in g.items():
            G[idx[e], j] = float(v)
    # project the E1 basis rows onto the orthogonal complement of colspace(G), take the top nonJ directions
    if len(gens):
        Qg, _ = np.linalg.qr(G)
        rk = info['dimJ']
        Qg = Qg[:, :rk]
        R = Kb - (Kb @ Qg) @ Qg.T
    else:
        R = Kb
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    nonJ = info['nonJ']
    W = Vt[:nonJ, :]            # nonJ x |mons|, orthonormal, spans E1 space mod J (numerically)
    Ds = [{m: float(W[r, j]) for j, m in enumerate(mons) if abs(W[r, j]) > 1e-14} for r in range(nonJ)]
    Ns, resid = e1_particular_solutions(Ds, sec, verbose=verbose)
    if verbose:
        print(f'[kernel_pairs] d={d} sec={sec_str(sec)} nonJ={nonJ} singular values of complement: '
              f'{s[:nonJ].round(3)} next {s[nonJ] if len(s) > nonJ else 0:.2e}; E1 residuals max {resid.max():.1e}', flush=True)
    return [(Ds[r], Ns[r], f'ker{r}') for r in range(nonJ)]


class FamilyRelaxation:
    def __init__(self, d, pairs, verbose=True):
        F, A, q, Fi, Hs = instance()
        self.d = d
        self.pairs = pairs
        K = len(pairs); self.K = K
        self.tags = [p[2] for p in pairs]
        dsecs = sorted({chi(e) for p in pairs for e in p[0]})
        self.dsecs = dsecs
        outsecs = sorted({sadd((0, 1), sadd(sk, sl)) for sk in dsecs for sl in dsecs})
        self.cok = {}
        self.offsets = {}
        off = 0
        for s in outsecs:
            B, mons, rank = prep_coker(2 * d, s, verbose=verbose)
            self.cok[s] = (B, mons, {e: i for i, e in enumerate(mons)})
            self.offsets[s] = off
            off += B.shape[0]
        self.m = off
        t0 = time.time()
        # global row index over the union of output-sector monomials; stacked cokernel matrix
        allmons = []
        for s in outsecs:
            allmons += self.cok[s][1]
        keys = _pack(allmons)
        order = np.argsort(keys)
        skeys = keys[order]
        Bstack = np.zeros((off, len(allmons)))
        col0 = 0
        for s in outsecs:
            B, mons, _ = self.cok[s]
            o = self.offsets[s]
            Bstack[o:o + B.shape[0], col0:col0 + len(mons)] = B
            col0 += len(mons)
        # N-monomial columns: union over slots of all monomials of degree d+1 of the needed sectors
        nsecs = sorted({sadd(sadd((0, 1), s), chi_var(i)) for s in dsecs for i in range(NV)})
        Ncols = []
        for i in range(NV):
            for s in sorted({sadd(sadd((0, 1), sd), chi_var(i)) for sd in dsecs}):
                Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
        nidx = {c: j for j, c in enumerate(Ncols)}
        nN = len(Ncols)
        Em = [np.array([m for (i, m) in Ncols if i == ii], dtype=np.int64) for ii in range(NV)]
        Ecol = [np.array([nidx[(ii, m)] for (i, m) in Ncols if i == ii], dtype=np.int64) for ii in range(NV)]
        # N matrix (nN x K)
        Nmat = np.zeros((nN, K))
        for l, (Dl, Nl, tag) in enumerate(pairs):
            for i in range(NV):
                for m, v in Nl[i].items():
                    Nmat[nidx[(i, m)], l] += float(v)
        import scipy.sparse as sp
        T = np.zeros((K, K, off))
        for k, (Dk, Nk, tag) in enumerate(pairs):
            Ed = np.array(list(Dk.keys()), dtype=np.int64)
            cd = np.array([float(v) for v in Dk.values()])
            rows_all, cols_all, vals_all = [], [], []
            for i in range(NV):
                if Em[i].shape[0] == 0:
                    continue
                # exponent of (t, m): Ed[t] + Em[m] - e_i ; coefficient cd[t]*(m_i - Ed[t,i])
                E = Ed[:, None, :] + Em[i][None, :, :]
                E[:, :, i] -= 1
                coef = cd[:, None] * (Em[i][None, :, i] - Ed[:, None, i])
                valid = (E[:, :, i] >= 0) & (coef != 0)
                kk = _pack(E[valid])
                pos = np.searchsorted(skeys, kk)
                pos = np.clip(pos, 0, len(skeys) - 1)
                ok = skeys[pos] == kk
                if not np.all(ok):
                    raise RuntimeError('monomial outside the output sectors')
                rows_all.append(order[pos])
                cols_all.append(np.broadcast_to(Ecol[i][None, :], valid.shape)[valid])
                vals_all.append(coef[valid])
            Mk = sp.csr_matrix((np.concatenate(vals_all), (np.concatenate(rows_all), np.concatenate(cols_all))),
                               shape=(len(allmons), nN))
            PM = (Mk.T @ Bstack.T).T          # off x nN
            T[k] = (PM @ Nmat).T               # K x off : row l = Pi B(D_k, N_l)
        self.T = 0.5 * (T + T.transpose(1, 0, 2))
        self.c0 = np.zeros(off)
        self.L = np.zeros((K, off))
        self.D0 = None
        self.scale = float(np.sqrt(np.mean(np.sum(self.T ** 2, axis=2)))) or 1.0
        # spurious subspace {c : D(c) in F*S}
        self.P, self.dimZ = self._complement()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        if verbose:
            print(f'[famrelax] d={d} secs={[sec_str(s) for s in dsecs]} K={K} m={off} nN={nN} spurious={self.dimZ} '
                  f'search dim={self.dimP} built in {self.build_time:.1f}s', flush=True)

    def _complement(self):
        F, A, q, Fi, Hs = instance()
        K = self.K
        dmons = {}
        cols = [p[0] for p in self.pairs]
        Fcols = []
        if self.d >= 6:
            for s in self.dsecs:
                Fcols += [shift(F, h, 1) for h in sector_monomials(self.d - 6, s)]
        for Cc in cols + Fcols:
            for e in Cc:
                dmons.setdefault(e, len(dmons))
        Mb = np.zeros((len(dmons), K + len(Fcols)))
        for k, Cc in enumerate(cols):
            for e, v in Cc.items():
                Mb[dmons[e], k] = float(v)
        for j, Cc in enumerate(Fcols):
            for e, v in Cc.items():
                Mb[dmons[e], K + j] = -float(v)
        ns = sla.null_space(Mb, rcond=1e-10)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return np.ascontiguousarray(U[:, r:]), r

    ob = Relaxation.ob
    ob_and_jac = Relaxation.ob_and_jac
    residual_fun = Relaxation.residual_fun
    residual_u = Relaxation.residual_u
    search = Relaxation.search

    def D_of(self, c):
        D = {}
        for k, (Dk, Nk, tag) in enumerate(self.pairs):
            if c[k]:
                D = add(D, scale({e: float(v) for e, v in Dk.items()}, float(c[k])))
        return D


def reconstruct_poly(Df, max_dens=(1, 2, 3, 4, 6, 7, 12, 14, 21, 28, 42, 84, 168, 294, 588, 1176, 2352, 4116, 8232, 10000, 100000, 1000000), tol=1e-7):
    """Rational reconstruction of a float polynomial (normalised so that its largest coefficient is 1)."""
    if not Df:
        return []
    e0 = max(Df, key=lambda e: abs(Df[e]))
    s = Df[e0]
    vals = {e: v / s for e, v in Df.items() if abs(v / s) > 1e-10}
    out = []
    for md in max_dens:
        fr = {e: Fraction(v).limit_denominator(md) for e, v in vals.items()}
        err = max(abs(float(fr[e]) - vals[e]) for e in vals)
        if err < tol:
            D = {e: f for e, f in fr.items() if f}
            if D not in out:
                out.append(D)
    return out


def cmd_famrelax(a):
    d = a.d
    sec = parse_sec(a.sec)
    secs = sectors_of(sec)
    pairs = []
    for s in secs:
        if 'J' in a.family:
            pairs += jacobian_pairs(d, s)
        if 'K' in a.family:
            pairs += kernel_pairs(d, s)
    rel = FamilyRelaxation(d, pairs)
    t0 = time.time()
    res = rel.search(a.restarts, a.seed, verbose=True, maxfev=a.maxfev)
    rels = np.array([c['rel'] for c in res])
    hist = {f'<1e-{k}': int(np.sum(rels < 10.0 ** (-k))) for k in (2, 4, 6, 8, 10, 12)}
    rec = {'event': 'famrelax', 'd': d, 'secs': [sec_str(s) for s in secs], 'family': a.family, 'K': rel.K,
           'dimZ': rel.dimZ, 'dimP': rel.dimP, 'm': rel.m, 'restarts': len(res), 'seed': a.seed,
           'best_rel': float(rels.min()), 'median_rel': float(np.median(rels)), 'hist': hist,
           'seconds': round(time.time() - t0, 1), 'build_seconds': round(rel.build_time, 1)}
    print('[famrelax] summary', rec, flush=True)
    log_event(rec)
    json.dump([{'rel': c['rel'], 'g': [float(x) for x in c['g']]} for c in res[:a.keep]],
              open(os.path.join(HERE, f'famrelax_d{d}_{a.sec}_{a.family}_s{a.seed}.json'), 'w'))
    for c in res:
        if c['rel'] > a.thresh:
            break
        Df = rel.D_of(c['g'])
        for D in reconstruct_poly(Df):
            exact_confirm(D, note=f'famrelax d={d} secs={a.sec} fam={a.family} rel={c["rel"]:.2e}')



# ----------------------------------------------------------------------------
# exact (mod p) quadrics of the family relaxation for msolve
# ----------------------------------------------------------------------------

def _parse_key(key):
    e = [0] * 6
    if key != '1':
        for part in key.split('|'):
            v, k = part.split('^')
            e[core.VARS.index(v)] = int(k)
    return tuple(e)


def family_pairs_modp(d, sec, family):
    """Pairs (D_k mod p, N_k mod p, tag): Jacobian pairs and/or the non-Jacobian E1 classes (exact mod p)."""
    F, A, q, Fi, Hs = instance()
    Ap = core.reduce_mod(A, P)
    Fip = [core.reduce_mod(f, P) for f in Fi]
    pairs = []
    if 'J' in family:
        for (i, m) in g_basis(d, sec):
            N = [dict() for _ in range(NV)]
            N[i] = shift(Ap, m, 1, P)
            pairs.append((shift(Fip[i], m, 1, P), N, f'J:{i}:{m}'))
    if 'K' in family:
        info = e1_space(d, sec, verbose=False)
        classes = [{_parse_key(k): v % P for k, v in b.items()} for b in info['basis_modp']]
        if classes:
            Ncols = []
            for i in range(NV):
                s = sadd(sadd((0, 1), sec), chi_var(i))
                Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
            tsec = sadd((0, 1), sec)
            tm = sector_monomials(d + 6, tsec)
            tidx = {e: i for i, e in enumerate(tm)}
            R, C = len(tm), len(Ncols)
            nc = len(classes)
            flat = [0] * (R * (C + nc))
            for j, (i, m) in enumerate(Ncols):
                for e, v in shift(Fip[i], m, 1, P).items():
                    flat[tidx[e] * (C + nc) + j] = v
            for c, Dk in enumerate(classes):
                for e, v in mul(Ap, Dk, P).items():
                    flat[tidx[e] * (C + nc) + C + c] = v
            Mr, rk = nmod_mat(R, C + nc, flat, P).rref()
            pivots = []
            for r in range(rk):
                j = 0
                while j < C + nc and int(Mr[r, j]) == 0:
                    j += 1
                pivots.append(j)
            assert all(j < C for j in pivots), 'E1 inconsistent for a kernel class (mod p)'
            # particular solutions: x_piv = -Mr[r, C+c]  (system  M x = b  written as [M | -b]?  we stored +b, so x = -coef)
            for c, Dk in enumerate(classes):
                N = [dict() for _ in range(NV)]
                for r, j in enumerate(pivots):
                    v = int(Mr[r, C + c])
                    if v:
                        i, m = Ncols[j]
                        N[i][m] = (-v) % P
                # verify
                acc = {}
                for i in range(NV):
                    acc = add(acc, mul(N[i], Fip[i], P), P)
                assert sub(acc, mul(Ap, Dk, P), P) == {}, 'E1 particular solution check failed'
                pairs.append((Dk, N, f'ker{c}'))
    return pairs


def family_quadrics_modp(d, pairs):
    """Exact quadrics mod p: for ob(c) = sum_{k,l} c_k c_l NF[B(D_k, N_l)] over the output sectors.
    Returns (coords dict {(sector, std monomial): {(k,l): coeff}} , reducers)."""
    import scipy.sparse as sp
    K = len(pairs)
    dsecs = sorted({chi(e) for p in pairs for e in p[0]})
    outsecs = sorted({sadd((0, 1), sadd(sk, sl)) for sk in dsecs for sl in dsecs})
    # exact normal-form maps NF_s (dimR_s x |S_2d^s|)
    NFs, stds, allmons = {}, {}, []
    for s in outsecs:
        mons = sector_monomials(2 * d, s)
        gens = gens_in_sector(2 * d, s, p=P)
        idx = {e: i for i, e in enumerate(mons)}
        flat = [0] * (len(gens) * len(mons))
        for r, (_, _, g) in enumerate(gens):
            for e, v in g.items():
                flat[r * len(mons) + idx[e]] = v
        Jr, rank = nmod_mat(len(gens), len(mons), flat, P).rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j < len(mons) and int(Jr[i, j]) == 0:
                j += 1
            pivots.append(j)
        pset = set(pivots)
        std = [j for j in range(len(mons)) if j not in pset]
        sidx = {j: r for r, j in enumerate(std)}
        NF = np.zeros((len(std), len(mons)), dtype=np.int64)
        for r, j in enumerate(std):
            NF[r, j] = 1
        for i, j in enumerate(pivots):
            for jj in std:
                v = int(Jr[i, jj])
                if v:
                    NF[sidx[jj], j] = (-v) % P
        NFs[s] = NF
        stds[s] = [mons[j] for j in std]
        allmons += mons
    keys = _pack(allmons)
    order = np.argsort(keys)
    skeys = keys[order]
    # stacked NF over the union of monomials
    dimR = sum(NFs[s].shape[0] for s in outsecs)
    NFstack = np.zeros((dimR, len(allmons)), dtype=np.int64)
    r0 = 0; c0 = 0
    rowsec = []
    for s in outsecs:
        NF = NFs[s]
        NFstack[r0:r0 + NF.shape[0], c0:c0 + NF.shape[1]] = NF
        rowsec += [(s, m) for m in stds[s]]
        r0 += NF.shape[0]; c0 += NF.shape[1]
    NF_lo = NFstack & 0xFFFF
    NF_hi = NFstack >> 16
    # N columns
    Ncols = []
    for i in range(NV):
        for s in sorted({sadd(sadd((0, 1), sd), chi_var(i)) for sd in dsecs}):
            Ncols += [(i, m) for m in sector_monomials(d + 1, s)]
    nidx = {c: j for j, c in enumerate(Ncols)}
    nN = len(Ncols)
    Em = [np.array([m for (i, m) in Ncols if i == ii], dtype=np.int64).reshape(-1, 6) for ii in range(NV)]
    Ecol = [np.array([nidx[(ii, m)] for (i, m) in Ncols if i == ii], dtype=np.int64) for ii in range(NV)]
    Nmat = np.zeros((nN, K), dtype=np.int64)
    for l, (Dl, Nl, tag) in enumerate(pairs):
        for i in range(NV):
            for m, v in Nl[i].items():
                Nmat[nidx[(i, m)], l] = (Nmat[nidx[(i, m)], l] + v) % P
    N_lo = Nmat & 0xFFFF
    N_hi = Nmat >> 16
    T = np.zeros((K, K, dimR), dtype=np.int64)
    for k, (Dk, Nk, tag) in enumerate(pairs):
        Ed = np.array(list(Dk.keys()), dtype=np.int64).reshape(-1, 6)
        cd = np.array([int(v) % P for v in Dk.values()], dtype=np.int64)
        rows_all, cols_all, vals_all = [], [], []
        for i in range(NV):
            if Em[i].shape[0] == 0:
                continue
            E = Ed[:, None, :] + Em[i][None, :, :]
            E[:, :, i] -= 1
            coef = (cd[:, None] * ((Em[i][None, :, i] - Ed[:, None, i]) % P)) % P
            valid = (E[:, :, i] >= 0) & (coef != 0)
            kk = _pack(E[valid])
            pos = np.clip(np.searchsorted(skeys, kk), 0, len(skeys) - 1)
            assert np.all(skeys[pos] == kk), 'monomial outside the output sectors'
            rows_all.append(order[pos]); cols_all.append(np.broadcast_to(Ecol[i][None, :], valid.shape)[valid]); vals_all.append(coef[valid])
        Mk = sp.csr_matrix((np.concatenate(vals_all), (np.concatenate(rows_all), np.concatenate(cols_all))),
                           shape=(len(allmons), nN), dtype=np.int64)
        Mk.sum_duplicates(); Mk.data %= P
        MkT = Mk.T.tocsr()
        PM = ((MkT @ NF_lo.T) % P + (((MkT @ NF_hi.T) % P) << 16) % P) % P      # nN x dimR
        PM_lo = PM & 0xFFFF; PM_hi = PM >> 16
        Tk = ((N_lo.T @ PM_lo) % P + (((N_lo.T @ PM_hi) % P) << 16) % P + (((N_hi.T @ PM_lo) % P) << 16) % P
              + (((N_hi.T @ PM_hi) % P) << 32) % P) % P                        # K x dimR : row l = NF B(D_k, N_l)
        T[k] = Tk
    coords = {}
    for r, (s, m) in enumerate(rowsec):
        qd = {}
        for k in range(K):
            v = int(T[k, k, r])
            if v:
                qd[(k, k)] = v
            for l in range(k + 1, K):
                v = (int(T[k, l, r]) + int(T[l, k, r])) % P
                if v:
                    qd[(k, l)] = v
        if qd:
            coords[(s, m)] = qd
    return coords


def spurious_pivots_family_modp(d, pairs):
    """Pivots (rref mod p) of {c : D(c) in F*S_{d-6}} for a general pair family."""
    F, A, q, Fi, Hs = instance()
    Fp = core.reduce_mod(F, P)
    K = len(pairs)
    dsecs = sorted({chi(e) for p in pairs for e in p[0]})
    Fcols = []
    if d >= 6:
        for s in dsecs:
            Fcols += [shift(Fp, h, 1, P) for h in sector_monomials(d - 6, s)]
    cols = [p[0] for p in pairs] + Fcols
    mons = {}
    for Cc in cols:
        for e in Cc:
            mons.setdefault(e, len(mons))
    R, C = len(mons), len(cols)
    flat = [0] * (R * C)
    for j, Cc in enumerate(cols):
        for e, v in Cc.items():
            flat[mons[e] * C + j] = v % P
    X, nullity = nmod_mat(R, C, flat, P).nullspace()
    if nullity == 0:
        return []
    proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], P)
    Pr, rk = proj.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < K and int(Pr[i, j]) == 0:
            j += 1
        if j < K:
            piv.append(j)
    return piv


def write_msolve_family(d, sec, family, path):
    pairs = []
    for s in sectors_of(sec):
        pairs += family_pairs_modp(d, s, family)
    K = len(pairs)
    coords = family_quadrics_modp(d, pairs)
    piv = spurious_pivots_family_modp(d, pairs)
    names = [f'g{k}' for k in range(K)]
    polys = [names[j] for j in piv]
    for (s, e), qd in coords.items():
        terms = [f'{c}*{names[k]}*{names[l]}' if k != l else f'{c}*{names[k]}^2' for (k, l), c in qd.items()]
        polys.append('+'.join(terms))
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n' + str(P) + '\n' + ',\n'.join(polys) + '\n')
    return K, len(polys), len(piv), [p[2] for p in pairs]


def run_msolve_file(path, threads=8, seed=1, timeout=None):
    import random
    random.seed(seed)
    txt = open(path).read().rstrip('\n')
    names = txt.split('\n')[0].split(',')
    chart = '+'.join(f'{random.randrange(1, P)}*{nm}' for nm in names) + f'+{P-1}'
    path2 = path.replace('.ms', '_chart.ms')
    open(path2, 'w').write(txt + ',\n' + chart + '\n')
    out = path2.replace('.ms', '.out')
    if os.path.exists(out):
        os.remove(out)
    t0 = time.time()
    try:
        r = subprocess.run(['msolve', '-v', '1', '-t', str(threads), '-f', path2, '-o', out], capture_output=True, text=True, timeout=timeout)
        err = r.stderr[-400:]
    except subprocess.TimeoutExpired:
        err = 'TIMEOUT'
    res = open(out).read() if os.path.exists(out) else ''
    return {'empty_projective': bool(res.strip().startswith('[-1]')), 'out_head': res[:200], 'stderr_tail': err,
            'seconds': round(time.time() - t0, 1), 'timed_out': err == 'TIMEOUT'}


def cmd_msolvef(a):
    path = os.path.join(HERE, f'famrelax_d{a.d}_{a.sec}_{a.family}.ms')
    t0 = time.time()
    K, n, npiv, tags = write_msolve_family(a.d, parse_sec(a.sec), a.family, path)
    rec = {'event': 'msolve-family-export', 'd': a.d, 'secs': a.sec, 'family': a.family, 'K': K, 'npolys': n,
           'spurious_pivots': npiv, 'kernel_classes': sum(1 for t in tags if t.startswith('ker')), 'seconds': round(time.time() - t0, 1)}
    print('[msolvef]', rec, flush=True)
    log_event(rec)
    if a.run:
        res = run_msolve_file(path, threads=a.threads, seed=a.seed, timeout=a.timeout)
        rec2 = {'event': 'msolve-family', 'd': a.d, 'secs': a.sec, 'family': a.family, 'K': K, 'npolys': n, **res}
        print('[msolvef]', rec2, flush=True)
        log_event(rec2)


_old_main = main


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == 'msolvef':
        ap = argparse.ArgumentParser()
        ap.add_argument('--d', type=int, required=True); ap.add_argument('--sec', required=True)
        ap.add_argument('--family', default='JK'); ap.add_argument('--run', action='store_true')
        ap.add_argument('--threads', type=int, default=8); ap.add_argument('--seed', type=int, default=1)
        ap.add_argument('--timeout', type=float, default=None)
        return cmd_msolvef(ap.parse_args(argv[1:]))
    if argv and argv[0] == 'famrelax':
        ap = argparse.ArgumentParser()
        ap.add_argument('--d', type=int, required=True); ap.add_argument('--sec', required=True)
        ap.add_argument('--family', default='JK'); ap.add_argument('--restarts', type=int, default=200)
        ap.add_argument('--seed', type=int, default=0); ap.add_argument('--maxfev', type=int, default=2000)
        ap.add_argument('--thresh', type=float, default=1e-9); ap.add_argument('--keep', type=int, default=20)
        return cmd_famrelax(ap.parse_args(argv[1:]))
    return _old_main(argv)


if __name__ == '__main__':
    main()
