"""Wave-2 library for the symmetric-sextic-residual lane.

Character-reduced (S3xS3 or a subgroup H, sign character chi) *float* E-systems for a linear family of
chi-semi-invariant denominators, used for a variable-projection (varpro) joint search over the full
E1-admissible family K_chi(d, pD) = { D of type chi, degree d, Y-parity pD : A*D in J(F) }.

Conventions follow charsolve.CharSystem (same orbit folding, same column meaning): the unknown attached to
column ('N', s, m) is the chi-twisted symmetrised field (1/|H|) sum_g chi(g) g.(m d_s), the Q-column ('Q', m)
is (1/|H|) sum_g g.m; one equation per H-orbit of degree-(2d+6) Y-odd monomials.  The varpro objective is
f(a) = min_x |M(a) x - b(a)|^2 on the unit sphere of the admissible subspace (E1 kernel, F-multiples removed).

Only an exact re-solve (charsolve.solve_char -> core.verify_witness) counts as evidence.
"""
from __future__ import annotations
import sys, os, time, json, math, itertools
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from flint import nmod_mat
from hodge import core, fixed_d, jacobian
from charsolve import Group, orbit_reps, live_reps, pN_of, act, CharSystem, solve_char, e1_char, J_char_basis

LANE = 'symmetric-sextic-residual'
RES = f'/tmp/hodge/results/{LANE}'
WLOG = f'{RES}/wave2-search.log.jsonl'
P0 = core.PRIMES[0]
ENGINE = 'w2lib-v1'


def wlog(rec, path=None):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', ENGINE)
    rec.setdefault('lane', LANE)
    with open(path or WLOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec


def load():
    F, A, q, rec = core.load_instance(LANE)
    return F, A, q


# ---------------------------------------------------------------- symmetric-function helpers
def var(i):
    e = [0] * 6; e[i] = 1
    return {tuple(e): Fraction(1)}

X = [var(0), var(1), var(2)]
Y = [var(3), var(4), var(5)]

def e_k(vs, k):
    out = {}
    for c in itertools.combinations(vs, k):
        t = {(0,) * 6: Fraction(1)}
        for v in c:
            t = core.mul(t, v)
        out = core.add(out, t)
    return out

def p_k(vs, k):
    out = {}
    for v in vs:
        out = core.add(out, core.power(v, k))
    return out

def delta(vs):
    return core.mul(core.mul(core.sub(vs[0], vs[1]), core.sub(vs[0], vs[2])), core.sub(vs[1], vs[2]))

GROUPS = {
    'triv': lambda: Group([[0, 1, 2], [3, 4, 5]], [1, 1]),
    'sgnX': lambda: Group([[0, 1, 2], [3, 4, 5]], [-1, 1]),
    'sgnY': lambda: Group([[0, 1, 2], [3, 4, 5]], [1, -1]),
    'sgnsgn': lambda: Group([[0, 1, 2], [3, 4, 5]], [-1, -1]),
    'S22_++': lambda: Group([[1, 2], [4, 5]], [1, 1]),
    'S22_-+': lambda: Group([[1, 2], [4, 5]], [-1, 1]),
    'S22_+-': lambda: Group([[1, 2], [4, 5]], [1, -1]),
    'S22_--': lambda: Group([[1, 2], [4, 5]], [-1, -1]),
    'S3S2_++': lambda: Group([[0, 1, 2], [4, 5]], [1, 1]),
    'S3S2_-+': lambda: Group([[0, 1, 2], [4, 5]], [-1, 1]),
    'S3S2_+-': lambda: Group([[0, 1, 2], [4, 5]], [1, -1]),
    'S3S2_--': lambda: Group([[0, 1, 2], [4, 5]], [-1, -1]),
    'S2S3_++': lambda: Group([[1, 2], [3, 4, 5]], [1, 1]),
    'S2S3_-+': lambda: Group([[1, 2], [3, 4, 5]], [-1, 1]),
    'S2S3_+-': lambda: Group([[1, 2], [3, 4, 5]], [1, -1]),
    'S2S3_--': lambda: Group([[1, 2], [3, 4, 5]], [-1, -1]),
}


def fpoly(P):
    """exact dict poly -> (E int64 n x 6, C float64 n)"""
    if not P:
        return np.zeros((0, 6), dtype=np.int64), np.zeros(0)
    items = list(P.items())
    return np.array([e for e, _ in items], dtype=np.int64), np.array([float(c) for _, c in items])


def fmul(a, b):
    """product of two float array-polys WITHOUT merging duplicates (folds merge later)."""
    (E1, C1), (E2, C2) = a, b
    if E1.shape[0] == 0 or E2.shape[0] == 0:
        return np.zeros((0, 6), dtype=np.int64), np.zeros(0)
    E = (E1[:, None, :] + E2[None, :, :]).reshape(-1, 6)
    C = (C1[:, None] * C2[None, :]).ravel()
    return E, C


def fscale(a, c):
    return a[0], a[1] * c


def fadd(a, b):
    return np.concatenate([a[0], b[0]]), np.concatenate([a[1], b[1]])


def fderiv(a, i):
    E, C = a
    keep = E[:, i] > 0
    E2 = E[keep].copy(); C2 = C[keep] * E2[:, i]
    E2[:, i] -= 1
    return E2, C2


def coords_in_sos(P, reps):
    """coordinates of a chi-type poly in the signed-orbit-sum basis of `reps` (coefficient at each rep)."""
    return np.array([float(P.get(r, 0)) for r in reps])


# ---------------------------------------------------------------- E1 kernel of an arbitrary family (mod p and float)
def e1_family_modp(F, A, q, basis, G, pD, p=P0):
    """Kernel of D -> [A D] in R/J restricted to span(basis) (chi-type polys of degree d, parity pD), mod p.
    Returns (dimK, kernel vectors as rows over the basis (ints mod p), nN, nrows)."""
    from charsolve import _tfold_factory, to_arr
    d = core.degree(basis[0]); pA = core.parity_of(A); degN = d + 6 * q - 5
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    tfold, nrows = _tfold_factory(G, d + 6 * q, (pA + pD) % 2, p)
    cols = []
    for s in G.slots:
        FE, FC = to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            cols.append(tfold(FE + np.array(m, dtype=np.int64), FC))
    nN = len(cols)
    for B in basis:
        E_, C_ = to_arr(core.mul(Ap, core.reduce_mod(B, p), p), p)
        cols.append(tfold(E_, C_))
    K = len(basis)
    M = np.stack(cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    if nullity == 0:
        return 0, [], nN, nrows
    proj = nmod_mat(K, nullity, [int(Xn[nN + j, k]) for j in range(K) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    return rk, [[int(Tr[i, j]) for j in range(K)] for i in range(rk)], nN, nrows


def e1_family_float(F, A, q, basis, G, pD, expect=None, tol=None):
    """Orthonormal basis (K x dimK) of the E1 kernel inside span(basis), by SVD of the twisted E1 system.
    `expect` (from e1_family_modp) fixes the kernel dimension; otherwise a gap heuristic is used."""
    d = core.degree(basis[0]); pA = core.parity_of(A); degN = d + 6 * q - 5
    rowmons = np.array(core.monomials(d + 6 * q, (pA + pD) % 2), dtype=np.int64)
    codes = G.canon(rowmons); signs = G.tsign(rowmons)
    live = np.unique(codes[signs != 0]); R = len(live)
    Fi = [fpoly(core.pderiv(F, i)) for i in range(6)]
    cols = []   # list of (row idx array, val array)
    for s in G.slots:
        sb = G.stab_blocks(s)
        FE, FC = Fi[s]
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            E = FE + np.array(m, dtype=np.int64)
            sg = G.tsign(E); keep = sg != 0
            idx = np.searchsorted(live, G.canon(E[keep]))
            cols.append((idx, q * FC[keep] * sg[keep]))
    nN = len(cols)
    Af = fpoly(A)
    for B in basis:
        E, C = fmul(Af, fpoly(B))
        sg = G.tsign(E); keep = sg != 0
        idx = np.searchsorted(live, G.canon(E[keep]))
        cols.append((idx, C[keep] * sg[keep]))
    K = len(basis)
    ri = np.concatenate([c[0] for c in cols]); ci = np.concatenate([np.full(len(c[0]), j) for j, c in enumerate(cols)])
    vals = np.concatenate([c[1] for c in cols])
    M = sp.coo_matrix((vals, (ri, ci)), shape=(R, nN + K)).toarray()
    # nullspace restricted to D-coordinates: use the SVD of the full matrix
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    if expect is not None:
        # the kernel of the FULL matrix has dimension nullity_full; its projection to D-coords has rank `expect`.
        pass
    smax = s[0] if len(s) else 1.0
    thr = tol if tol is not None else 1e-9 * smax
    rank = int(np.sum(s > thr))
    null = Vt[rank:, :].T            # (nN+K) x nullity
    Dpart = null[nN:, :]             # K x nullity
    # orthonormal basis of the column space of Dpart
    Ud, sd, _ = np.linalg.svd(Dpart, full_matrices=False)
    if expect is not None:
        r2 = expect
        gap = (sd[r2 - 1] if r2 > 0 else None, sd[r2] if r2 < len(sd) else None)
    else:
        r2 = int(np.sum(sd > 1e-9 * max(1.0, sd[0] if len(sd) else 1.0)))
        gap = None
    return Ud[:, :r2], {'rows': R, 'nN': nN, 'K': K, 'rank_full': rank, 'nullity_full': null.shape[1],
                        'dimK': r2, 'sv_gap': gap, 'sv_small': s[max(rank - 3, 0):rank + 3].tolist()}


# ---------------------------------------------------------------- character-reduced float family
class CharFamily:
    """Float E-system for D(a) = sum_k a_k B_k, all B_k of character chi (group G), degree d, parity pD."""

    def __init__(self, F, A, q, G, basis, reps=None, verbose=False, cond=1e-9, label=''):
        t0 = time.time()
        assert q == 1
        self.F, self.A, self.q, self.G = F, A, q, G
        self.cond = cond; self.label = label
        self.basis = basis; K = len(basis); self.K = K
        self.reps = reps
        d = core.degree(basis[0]); self.d = d
        pD = core.parity_of(basis[0]); self.pD = pD
        pA = core.parity_of(A); self.pA = pA
        degN, degQ, degE = d + 6 * q - 5, 2 * d - 6, 2 * d + 6 * q
        self.degN, self.degQ = degN, degQ
        rowmons = np.array(core.monomials(degE, pA), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); R = len(self.orbit_codes); self.rows = R
        Ff = fpoly(F); Fi = [fpoly(core.pderiv(F, i)) for i in range(6)]
        Bf = [fpoly(B) for B in basis]
        # N columns
        Ncols = []
        for s in G.slots:
            sb = G.stab_blocks(s)
            for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                Ncols.append((s, m))
        self.Ncols = Ncols; nN = len(Ncols); self.nN = nN
        # per slot: rep arrays
        slot_reps = {}
        for j, (s, m) in enumerate(Ncols):
            slot_reps.setdefault(s, []).append((j, m))
        tri = []  # (row, col = j*K + k ... ) we build Mstack of shape (R*nN, K): entry (r, j) of M_k at row r*nN + j, col k
        for k, B in enumerate(Bf):
            FB = fmul(Ff, B)
            for s, lst in slot_reps.items():
                Gs = fadd(fscale(fmul(B, Fi[s]), -q), fscale(fmul(Ff, fderiv(B, s)), -1.0))   # -(q B F_s + F dB/ds)
                J = np.array([j for j, _ in lst]); Mm = np.array([m for _, m in lst], dtype=np.int64)
                # part 1: shift(Gs, m) for all m
                E1 = (Gs[0][None, :, :] + Mm[:, None, :]).reshape(-1, 6)
                C1 = np.broadcast_to(Gs[1][None, :], (len(lst), Gs[1].shape[0])).ravel()
                r1 = np.searchsorted(self.orbit_codes, G.canon(E1))
                j1 = np.repeat(J, Gs[1].shape[0])
                tri.append((r1 * nN + j1, np.full(len(r1), k), C1))
                # part 2: m[s] * shift(FB, m - e_s) for m with m[s] > 0
                sel = Mm[:, s] > 0
                if sel.any():
                    M2 = Mm[sel].copy(); M2[:, s] -= 1; J2 = J[sel]; mult = Mm[sel][:, s].astype(float)
                    E2 = (FB[0][None, :, :] + M2[:, None, :]).reshape(-1, 6)
                    C2 = (mult[:, None] * FB[1][None, :]).ravel()
                    r2 = np.searchsorted(self.orbit_codes, G.canon(E2))
                    j2 = np.repeat(J2, FB[1].shape[0])
                    tri.append((r2 * nN + j2, np.full(len(r2), k), C2))
        ri = np.concatenate([t[0] for t in tri]); ci = np.concatenate([t[1] for t in tri]); vals = np.concatenate([t[2] for t in tri])
        self.Mstack = sp.coo_matrix((vals, (ri, ci)), shape=(R * nN, K)).tocsr()
        self.MstackT = self.Mstack.T.tocsr()
        # Q columns (constant): -F^2 * Sym(m)
        Qcols = orbit_reps(degQ, pA, G.canon) if degQ >= 0 else []
        self.Qcols = Qcols; nQ = len(Qcols); self.nQ = nQ
        MQ = np.zeros((R, nQ))
        if nQ:
            F2 = fmul(Ff, Ff)
            Qm = np.array(Qcols, dtype=np.int64)
            E = (F2[0][None, :, :] + Qm[:, None, :]).reshape(-1, 6)
            C = -np.broadcast_to(F2[1][None, :], (nQ, F2[1].shape[0])).ravel()
            r = np.searchsorted(self.orbit_codes, G.canon(E)); c = np.repeat(np.arange(nQ), F2[1].shape[0])
            np.add.at(MQ, (r, c), C)
        self.MQ = MQ
        # rhs pair vectors: v_kl = -A B_k B_l   -> Vstack (R, K*K)
        Af = fpoly(A)
        AB = [fmul(Af, B) for B in Bf]
        tri = []
        for k in range(K):
            for l in range(K):
                E, C = fmul(AB[k], Bf[l])
                r = np.searchsorted(self.orbit_codes, G.canon(E))
                tri.append((r, np.full(len(r), k * K + l), -C))
        ri = np.concatenate([t[0] for t in tri]); ci = np.concatenate([t[1] for t in tri]); vals = np.concatenate([t[2] for t in tri])
        self.Vstack = sp.coo_matrix((vals, (ri, ci)), shape=(R, K * K)).tocsr()
        self.VstackT = self.Vstack.T.tocsr()
        self.ncols = nN + nQ
        self.P = np.eye(K); self.dimP = K
        self._nfev = 0
        self.build_time = time.time() - t0
        if verbose:
            print(f'[charfamily {label}] d={d} pD={pD} K={K} rows={R} cols={self.ncols} (N {nN}, Q {nQ}) '
                  f'nnz(M)={self.Mstack.nnz} built in {self.build_time:.1f}s', file=sys.stderr)

    # ------------------------------------------------------------ parameter subspace
    def set_subspace(self, P):
        """restrict a = P u (P: K x dimP with orthonormal columns)."""
        self.P = np.asarray(P, float); self.dimP = self.P.shape[1]

    def f_multiple_vectors(self):
        """vectors (in basis coordinates) spanning {a : D(a) = F * G, G of type chi, degree d-6}; needs reps."""
        if self.d < 6 or self.reps is None:
            return np.zeros((self.K, 0))
        dG = self.d - 6
        Gl = live_reps(dG, (self.pD - core.parity_of(self.F)) % 2, self.G)
        vecs = []
        for _, Gp in Gl:
            prod = core.mul(self.F, Gp)
            vecs.append(coords_in_sos(prod, self.reps))
        return np.array(vecs).T if vecs else np.zeros((self.K, 0))

    # ------------------------------------------------------------ objective
    def assemble(self, a):
        M = (self.Mstack @ a).reshape(self.rows, self.nN)
        if self.nQ:
            M = np.hstack([M, self.MQ])
        b = self.Vstack @ np.outer(a, a).ravel()
        return M, b

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r

    def rel(self, a):
        f, x, r = self.value(a)
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, b = self.assemble(a)
        return f / max(float(b @ b), 1e-300), f, x

    def value_and_grad(self, u):
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        f = float(r @ r)
        xN = x[:self.nN]
        gM = self.MstackT @ np.outer(r, xN).ravel()          # K: r . (M_k xN)
        W = (self.VstackT @ r).reshape(self.K, self.K)       # W[k,l] = r . v_kl
        gb = 2.0 * (W @ a)                                   # r . db/da_k
        g = 2.0 * (gM - gb)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def is_spurious(self, a, tol=1e-8):
        """numerically: is D(a) divisible by F?"""
        if self.d < 6:
            return False, 1.0
        a = np.asarray(a, float)
        D = {}
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * float(v)
        dG = self.d - 6
        Gm = core.monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Fcols = [core.shift(self.F, m, 1) for m in Gm]
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

    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, inits=(), stop_rel=1e-13,
               callback=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if it < len(inits):
                u0 = self.P.T @ np.asarray(inits[it], float)
                if np.linalg.norm(u0) < 1e-12:
                    u0 = rng.standard_normal(self.dimP)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            rel, f, x = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': rel, 'a': a, 'spurious': spur, 'div_rel': srel, 'iters': int(res.nit),
                    'restart': it, 'seed': seed, 'nfev': int(res.nfev)}
            out.append(cand)
            if callback:
                callback(cand)
            if verbose:
                print(f'[search {self.label}] restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # ------------------------------------------------------------ exact confirmation
    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.basis):
            if frac[k]:
                D = core.add(D, core.scale(B, Fraction(frac[k])))
        return D

    def reconstruct(self, a, max_dens=(12, 100, 1000, 10000, 100000), tol_scale=1e-7):
        """rational candidates for a (normalised: largest |a_k| = 1), increasing denominators."""
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a)))
        a = a / a[j]
        cands = []
        for md in max_dens:
            fr = tuple(Fraction(float(v)).limit_denominator(md) for v in a)
            err = max(abs(float(f) - v) for f, v in zip(fr, a))
            if err < tol_scale * max(1.0, md ** 0.5) or md == max_dens[-1]:
                if fr not in cands:
                    cands.append(fr)
        # common-denominator scaling
        for den in range(1, 2001):
            v = a * den
            if np.max(np.abs(v - np.round(v))) < 1e-6 * den:
                fr = tuple(Fraction(int(round(x)), den) for x in v)
                if fr not in cands:
                    cands.insert(0, fr)
                break
        return cands

    def confirm(self, frac, verbose=False):
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        if self.d >= 6:
            divF, _ = jacobian.divisible_by(self.F, D)
            if divF:
                return 'spurious', {'note': 'F divides D', 'D': D}
        ok, info = e1_char(self.F, self.A, self.q, D, self.G)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails', 'e1': info, 'D': D}
        status, data = solve_char(self.F, self.A, self.q, D, self.G, verbose=verbose)
        data['D'] = D
        return status, data


def build_sector_family(F, A, q, chi, d, pD, p=P0, verbose=False, exclude_F=True, label=None):
    """CharFamily on the signed-orbit-sum basis of character chi, restricted (via P) to the E1 kernel minus
    F-multiples.  Returns (fam, info) or (None, info) if the admissible subspace is trivial."""
    t0 = time.time()
    G = GROUPS[chi]()
    L = live_reps(d, pD, G)
    reps = [r for r, _ in L]; basis = [P for _, P in L]
    info = {'chi': chi, 'd': d, 'pD': pD, 'nD': len(basis)}
    if not basis:
        info.update(dimK=0, dimP=0); return None, info
    dimK, kvecs, nN1, rows1 = e1_family_modp(F, A, q, basis, G, pD, p)
    info.update(dimK=dimK, e1_rows=rows1, e1_nN=nN1)
    # J part (mod p): dimension of J^chi_d inside the family
    Jp = J_char_basis(core.reduce_mod(F, p), d, pD, G, p)
    if Jp:
        Jm = np.array([[Pj.get(r, 0) % p for r in reps] for Pj in Jp], dtype=np.int64)
        dimJ = nmod_mat(Jm.shape[0], Jm.shape[1], Jm.ravel().tolist(), p).rank()
    else:
        dimJ = 0
    info['dimJ'] = dimJ; info['ker_beyond_J'] = dimK - dimJ
    if dimK == 0:
        info['dimP'] = 0; return None, info
    Kb, finfo = e1_family_float(F, A, q, basis, G, pD, expect=dimK)
    info['sv_gap'] = finfo['sv_gap']
    fam = CharFamily(F, A, q, G, basis, reps=reps, verbose=verbose, label=label or f'{chi}/d{d}/p{pD}')
    Pm = Kb
    nF = 0
    if exclude_F and d >= 6:
        Fv = fam.f_multiple_vectors()
        if Fv.shape[1]:
            S = Kb.T @ Fv
            resid = np.linalg.norm(Fv - Kb @ S) / max(np.linalg.norm(Fv), 1e-300)
            info['F_in_kernel_resid'] = float(resid)
            U, s, _ = np.linalg.svd(S, full_matrices=True)
            nF = int(np.sum(s > 1e-9 * max(1.0, s[0])))
            Pm = Kb @ U[:, nF:]
    info['dimF'] = nF
    fam.set_subspace(Pm)
    info['dimP'] = fam.dimP
    info['rows'] = fam.rows; info['nN'] = fam.nN; info['nQ'] = fam.nQ
    info['build_s'] = round(time.time() - t0, 1)
    fam.info = info
    return fam, info


def witness_json(D, N, Q):
    return {'d': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}


def try_confirm(fam, cand, logrec, verbose=True):
    """rational reconstruction + exact re-solve of a candidate; writes witness-candidate.json on success."""
    results = []
    for fr in fam.reconstruct(cand['a']):
        status, data = fam.confirm(fr, verbose=verbose)
        rec = {'status': status, 'den': max(f.denominator for f in fr), 'nz': sum(1 for f in fr if f)}
        if status == 'solved':
            D, N, Q = data['D'], data['N'], data['Q']
            ok, info = core.verify_witness(fam.F, fam.A, fam.q, D, N, Q)
            rec['verified'] = bool(ok)
            if ok:
                w = witness_json(D, N, Q)
                w.update({'instance': LANE, 'family': fam.label, 'rel': cand['rel'], 'engine': ENGINE})
                path = f'{RES}/witness-candidate.json'
                with open(path, 'w') as f:
                    json.dump(w, f, indent=1)
                wlog(dict(logrec, event='WITNESS', status='solved', D=w['D'][:2000], d=w['d'], path=path))
                print('WITNESS FOUND', path, flush=True)
                results.append(rec)
                return True, results
        results.append(rec)
    return False, results


# ---------------------------------------------------------------- generic helpers (added after validation)
def generic_f_multiples(F, basis):
    """K x nF matrix of coefficient vectors a with D(a) in F * R_{d-6} (monomial-space nullspace)."""
    d = core.degree(basis[0]); pD = core.parity_of(basis[0]); K = len(basis)
    if d < 6:
        return np.zeros((K, 0))
    Gm = core.monomials(d - 6, (pD - core.parity_of(F)) % 2)
    Fcols = [core.shift(F, m, 1) for m in Gm]
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
    if S.shape[1] == 0:
        return np.zeros((K, 0))
    U, s, _ = np.linalg.svd(S, full_matrices=False)
    r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
    return U[:, :r]


def build_family(F, A, q, G, basis, label, p=P0, verbose=False, exclude_F=True):
    """CharFamily on an arbitrary chi-type basis, restricted to the E1 kernel minus F-multiples."""
    t0 = time.time()
    d = core.degree(basis[0]); pD = core.parity_of(basis[0])
    info = {'label': label, 'd': d, 'pD': pD, 'nD': len(basis)}
    dimK, kvecs, nN1, rows1 = e1_family_modp(F, A, q, basis, G, pD, p)
    info.update(dimK=dimK, e1_rows=rows1, e1_nN=nN1)
    if dimK == 0:
        info['dimP'] = 0; return None, info
    Kb, finfo = e1_family_float(F, A, q, basis, G, pD, expect=dimK)
    info['sv_gap'] = finfo['sv_gap']
    fam = CharFamily(F, A, q, G, basis, reps=None, verbose=verbose, label=label)
    Pm = Kb; nF = 0
    if exclude_F and d >= 6:
        Fv = generic_f_multiples(F, basis)
        if Fv.shape[1]:
            S = Kb.T @ Fv
            info['F_in_kernel_resid'] = float(np.linalg.norm(Fv - Kb @ S) / max(np.linalg.norm(Fv), 1e-300))
            U, s, _ = np.linalg.svd(S, full_matrices=True)
            nF = int(np.sum(s > 1e-9 * max(1.0, s[0])))
            Pm = Kb @ U[:, nF:]
    info['dimF'] = nF
    fam.set_subspace(Pm)
    info.update(dimP=fam.dimP, rows=fam.rows, nN=fam.nN, nQ=fam.nQ, build_s=round(time.time() - t0, 1))
    fam.info = info
    return fam, info


def restricted_family(polys, G, p=P0):
    """chi-symmetrise the given exact polys (common degree/parity) and extract a basis of their span."""
    syms = []
    for P in polys:
        S = G.sym_chi(P)
        if S:
            syms.append(S)
    if not syms:
        return []
    mons = {}
    for S in syms:
        for e in S:
            mons.setdefault(e, len(mons))
    rows = []
    for S in syms:
        Sp = core.reduce_mod(S, p)
        rows.append([Sp.get(e, 0) for e in mons])
    M = nmod_mat(len(rows), len(mons), [v for r in rows for v in r], p)
    # greedy independent subset via rref of the transpose (pivot columns = independent polys)
    Mt = M.transpose()
    Mr, rk = Mt.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < Mt.ncols() and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    return [syms[j] for j in piv]


def solve_single(F, A, q, D, G, p=P0, lift=False):
    """fixed-D character-reduced solve (mod p); D may be an exact or a mod-p (int coefficient) poly."""
    t0 = time.time()
    ok, info = e1_char(F, A, q, D, G, p)
    if not ok:
        return {'e1': False, 'consistent': False, 'rank': info['rank'], 'nullity': info['cols'] - info['rank'],
                'rows': info['rows'], 'cols': info['cols'], 'seconds': round(time.time() - t0, 2)}
    S = CharSystem(F, A, q, D, G, p=p)
    r = S.solve()
    out = {'e1': True, 'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
           'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 2)}
    if r['consistent'] and lift:
        status, data = solve_char(F, A, q, D, G)
        out['lift_status'] = status
        if status == 'solved':
            out['N'], out['Q'] = data['N'], data['Q']
    return out


def kernel_classes(F, A, q, chi, d, pD, p=P0):
    """E1-kernel classes modulo J^chi_d (mod p; signed-int coefficient polys) via charsolve.e1_kernel_char."""
    from charsolve import e1_kernel_char
    G = GROUPS[chi]()
    info, polys = e1_kernel_char(F, A, q, d, pD, G, p)
    return G, info, polys
