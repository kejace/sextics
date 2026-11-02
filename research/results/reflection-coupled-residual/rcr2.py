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
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', fixed_d.ENGINE_VERSION)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')


_INST = None


def instance():
    global _INST
    if _INST is None:
        F, A, q, rec = core.load_instance(LANE)
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
    return os.path.join(HERE, f'coker_{deg}_{sec_str(sec)}.npz')


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
    return os.path.join(HERE, f'e1_{d}_{sec_str(sec)}')


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
    tsec = sadd(sec, chi((3, 0, 0, 3, 0, 0)))
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
    tmp = f'{base}.{os.getpid()}.tmp.npz'
    np.savez_compressed(tmp, Kb=Kb, mons=np.array(Dm, dtype=np.int16))
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
        self.build_time = time.time() - t0
        if verbose:
            print(f'[relax] d={d} secs={[sec_str(s) for s in self.dsecs]} K={K} m={off} affine={D0 is not None} '
                  f'built in {self.build_time:.1f}s', flush=True)

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
            g0 = rng.standard_normal(self.K)
            if self.D0 is None:
                g0 /= np.linalg.norm(g0)
            method = 'lm' if self.m + 1 >= self.K else 'trf'
            res = least_squares(lambda g: self.residual_fun(g)[0], g0,
                                jac=lambda g: self.residual_fun(g)[1], method=method,
                                xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxfev)
            g = res.x
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
            json.dump(w, open(os.path.join(HERE, 'witness-candidate.json'), 'w'), indent=1)
            rec['witness_file'] = path
            print('!!!!!!!! WITNESS FOUND', path, flush=True)
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


def write_msolve(d, dsecs, path):
    gb, coords = relaxation_quadrics_modp(d, dsecs)
    K = len(gb)
    names = [f'g{k}' for k in range(K)]
    lines = [','.join(names), str(P)]
    polys = []
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
           'm': rel.m, 'restarts': len(res), 'seed': a.seed, 'best_rel': float(rels.min()),
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
    print(f'[msolve] d={d} secs={a.sec}: {n} quadrics in {K} unknowns -> {path}', flush=True)
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
        empty = res.strip().startswith('[1]') or res.strip() == '[1]:' or '[1]' == res.strip()[:3]
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
    p = sp.add_parser('confirm'); p.add_argument('--D', required=True); p.add_argument('--note', default=''); p.add_argument('--full', action='store_true')
    a = ap.parse_args(argv)
    {'prep': cmd_prep, 'e1': cmd_e1, 'relax': cmd_relax, 'msolve': cmd_msolve, 'varpro': cmd_varpro, 'confirm': cmd_confirm}[a.cmd](a)


if __name__ == '__main__':
    main()
