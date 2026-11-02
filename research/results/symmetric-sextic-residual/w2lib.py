"""w2lib -- wave-2 search library for the symmetric-sextic-residual lane.

Character-twisted (H = S3(X) x S3(Y), chi = sign character per block) variable-projection search over exact
rational linear families of chi-semi-invariant denominators  D(a) = sum_k a_k M_k  (E1 kernel classes +
J(F)_d^chi), with the numerator field N chi-equivariant and Q H-invariant (lossless for such D, see
charsolve.py).  Rows are folded over H-orbits, so the dense least-squares problem is ~36x smaller than
hodge.varpro.VarProFamily.  Candidates with small residual are rationally reconstructed and re-solved
EXACTLY with charsolve.solve_char (reduced, CRT-lifted) and core.verify_witness.

Engine (/home/kejace/devel/sextics-work/hodge) is not modified.
"""
import sys, os, time, json, math, itertools, pickle
sys.path.insert(0, '/home/kejace/devel/sextics-work')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.linalg as sla, scipy.sparse as sp
from scipy.optimize import minimize
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from flint import nmod_mat
import charsolve as cs
from charsolve import Group, orbit_reps, live_reps, to_arr, act, pN_of

LANE = HERE
IID = 'symmetric-sextic-residual'
WLOG = f'{LANE}/wave2-search.log.jsonl'
PRIMES = core.PRIMES

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def wlog(rec):
    rec = dict(rec); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); rec.setdefault('engine', 2)
    with open(WLOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec

CHARS = {'triv': (1, 1), 'sgnX': (-1, 1), 'sgnY': (1, -1), 'sgnsgn': (-1, -1)}
BLOCKS = {'S3S3': [[0, 1, 2], [3, 4, 5]], 'S2S2': [[1, 2], [4, 5]], 'S3S2': [[0, 1, 2], [4, 5]], 'S2S3': [[1, 2], [3, 4, 5]]}

def group(chi='triv', blocks='S3S3'):
    return Group(BLOCKS[blocks], CHARS[chi])

# ----------------------------------------------------------------------------------------------------------
# exact rational families: kernel classes (CRT + ratrecon) and an independent basis of J_d^chi
# ----------------------------------------------------------------------------------------------------------
def _poly_float(P):
    return {e: float(v) for e, v in P.items()}

def _kernel_vecs_modp(F, A, q, d, pD, G, p):
    """Returns (Dreps, Dbasis(signed orbit sums), Kvecs (rref rows mod p in Dbasis coords), Jvecs)."""
    pA = core.parity_of(A); degN = d + 6 * q - 5
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    tfold, nrows = cs._tfold_factory(G, d + 6 * q, (pA + pD) % 2, p)
    cols = []
    for s in G.slots:
        FE, FC = to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            cols.append(tfold(FE + np.array(m, dtype=np.int64), FC))
    nN = len(cols)
    Dl = live_reps(d, pD, G)
    Dreps = [r for r, _ in Dl]; Dbasis = [P for _, P in Dl]
    for P in Dbasis:
        E_, C_ = to_arr(core.mul(Ap, core.reduce_mod(P, p), p), p)
        cols.append(tfold(E_, C_))
    nD = len(Dbasis)
    Jpolys = cs.J_char_basis(Fp, d, pD, G, p)
    Jvecs = [[Pj.get(r, 0) % p for r in Dreps] for Pj in Jpolys]
    if nD == 0:
        return Dreps, Dbasis, [], Jvecs, {'rows': nrows, 'nN': nN, 'nD': 0}
    M = np.stack(cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    if nullity == 0:
        return Dreps, Dbasis, [], Jvecs, {'rows': nrows, 'nN': nN, 'nD': nD}
    proj = nmod_mat(nD, nullity, [int(Xn[nN + j, k]) for j in range(nD) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    Kvecs = [[int(Tr[i, j]) for j in range(nD)] for i in range(rk)]
    return Dreps, Dbasis, Kvecs, Jvecs, {'rows': nrows, 'nN': nN, 'nD': nD}

def _rref_pivots(vecs, n, p):
    if not vecs:
        return [], []
    M = nmod_mat(len(vecs), n, [v % p for row in vecs for v in row], p)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < n and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    rows = [[int(Mr[i, j]) for j in range(n)] for i in range(rk)]
    return rows, piv

def _crt_lift(rows_by_prime, primes):
    """rows_by_prime: list (per prime) of equal-shape int matrices (lists). Returns Fraction rows or None."""
    acc = [list(r) for r in rows_by_prime[0]]; mod = primes[0]
    for rows, p in zip(rows_by_prime[1:], primes[1:]):
        acc = [[core.crt_pair(a, mod, b, p)[0] for a, b in zip(ra, rb)] for ra, rb in zip(acc, rows)]
        mod *= p
    out = []
    for r in acc:
        rr = [core.ratrecon(a, mod) for a in r]
        if any(v is None for v in rr):
            return None
        out.append(rr)
    return out

def _primelist(n=24):
    import sympy
    ps = []; p = 2 ** 31
    while len(ps) < n:
        p = sympy.prevprime(p); ps.append(p)
    return ps
PRIMELIST = _primelist(24)

def family(F, A, q, d, pD, chi='triv', blocks='S3S3', max_primes=23, verbose=False, cache=True):
    """Exact rational chi-sector family at degree d, parity pD:
       returns dict(kernel=[polys], jac=[polys], dims..., G).  Kernel polys = E1 classes beyond J (rref canonical,
       CRT over as many primes as needed + ratrecon, verified by the twisted E1 test at a fresh prime),
       jac polys = independent basis of J_d^chi (exact, from J_char_basis)."""
    fn = f'{LANE}/w2fam_{blocks}_{chi}_d{d}_p{pD}.pkl'
    if cache and os.path.exists(fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)
    t0 = time.time()
    G = group(chi, blocks)
    primes = PRIMELIST[:max_primes + 1]
    Krows_by_p = []; used = []; piv0 = None
    Dreps = Dbasis = None; info = None
    kernel = None; nker = 0; verified = False
    for idx, p in enumerate(primes[:-1]):
        Dreps, Dbasis, Kvecs, Jvecs, info = _kernel_vecs_modp(F, A, q, d, pD, G, p)
        nD = len(Dreps)
        Jr, Jpiv = _rref_pivots(Jvecs, nD, p)
        Mr, piv = _rref_pivots(Jr + Kvecs, nD, p)
        extra = [row for row, pv in zip(Mr, piv) if pv not in Jpiv]
        sig = (tuple(Jpiv), tuple(pv for pv in piv if pv not in Jpiv))
        if piv0 is not None and sig != piv0:
            if verbose: print(f'[family] pivot mismatch at prime {p}; restarting CRT', file=sys.stderr)
            Krows_by_p = []; used = []
        piv0 = sig
        Krows_by_p.append(extra); used.append(p)
        nker = len(extra)
        if nker == 0:
            kernel = []; verified = True; break
        Kfrac = _crt_lift(Krows_by_p, used)
        if Kfrac is None:
            continue
        polys = []
        for row in Kfrac:
            P = {}
            for j, v in enumerate(row):
                if v:
                    P = core.add(P, core.scale(Dbasis[j], Fraction(v)))
            polys.append(P)
        pv = primes[idx + 1]
        if all(cs.e1_char(F, A, q, P, G, p=pv)[0] for P in polys):
            kernel = polys; verified = True
            if verbose: print(f'[family] kernel ratrecon ok with {len(used)} primes, verified mod {pv}', file=sys.stderr)
            break
        elif verbose:
            print(f'[family] ratrecon candidate failed verification with {len(used)} primes; continuing', file=sys.stderr)
    if kernel is None:
        if verbose: print('[family] ratrecon failed; kernel left EMPTY (ratrecon_ok=False)', file=sys.stderr)
        kernel = []
    # exact J basis: independent subset of J_char_basis (exact Fractions) by mod-p pivots
    Jexact = cs.J_char_basis(F, d, pD, G, None)
    p = primes[0]
    Jv = [[_modp(Pj.get(r, 0), p) for r in Dreps] for Pj in Jexact]
    jac = []; cur = []; rk = 0
    for Pj, v in zip(Jexact, Jv):
        test = cur + [v]
        rows, piv = _rref_pivots(test, len(Dreps), p)
        if len(rows) > rk:
            rk = len(rows); cur = test; jac.append(Pj)
    out = {'d': d, 'pD': pD, 'chi': chi, 'blocks': blocks, 'kernel': kernel, 'jac': jac, 'ker': nker,
           'dimJ': len(jac), 'nD': len(Dreps), 'primes_used': used, 'ratrecon_ok': verified,
           'kernel_verified_modp': verified if nker else None,
           'rows_e1': info['rows'], 'nN_e1': info['nN'], 'seconds': round(time.time() - t0, 1)}
    if cache and verified:
        with open(fn, 'wb') as f:
            pickle.dump(out, f)
    return out

def _modp(v, p):
    v = Fraction(v)
    return (v.numerator * pow(v.denominator, -1, p)) % p


def pow2_scale(P):
    """power-of-2 factor s (Fraction) with |s|*||P||_2 ~ 1 (keeps rationality)."""
    nrm = math.sqrt(sum(float(v) ** 2 for v in P.values())) if P else 1.0
    if nrm == 0: return Fraction(1)
    e = int(round(math.log2(nrm)))
    return Fraction(1, 2 ** e) if e >= 0 else Fraction(2 ** (-e))

def scale_pow2(P):
    return core.scale(P, pow2_scale(P))

def e1_lift_float(F, A, q, D, G):
    """Float chi-equivariant N with N(F) ~ A D / q (twisted reduced E1 least squares).  Returns (N list of float polys, relres)."""
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A); degN = d + 6 * q - 5
    Ff = _poly_float(F); Df = _poly_float(D); Af = _poly_float(A)
    Fi = [core.pderiv(Ff, i) for i in range(6)]
    rowmons = np.array(core.monomials(d + 6 * q, (pA + pD) % 2), dtype=np.int64)
    codes = G.canon(rowmons); signs = G.tsign(rowmons)
    live = np.unique(codes[signs != 0]); n = len(live)
    def tfold(P):
        if not P: return np.zeros(n)
        items = list(P.items()); E = np.array([e for e, _ in items], dtype=np.int64); C = np.array([float(c) for _, c in items])
        s = G.tsign(E); keep = s != 0
        idx = np.searchsorted(live, G.canon(E[keep]))
        return np.bincount(idx, weights=C[keep] * s[keep], minlength=n)
    cols = []; info = []
    for s in G.slots:
        sb = G.stab_blocks(s)
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            cols.append(tfold(core.shift(core.scale(Fi[s], q), m, 1))); info.append((s, m))
    M = np.stack(cols, axis=1); rhs = tfold(core.mul(Af, Df))
    x, *_ = sla.lstsq(M, rhs, cond=1e-12, lapack_driver='gelsd', check_finite=False)
    relres = float(np.linalg.norm(M @ x - rhs) / max(np.linalg.norm(rhs), 1e-300))
    N = [{} for _ in range(6)]
    for (s, m), v in zip(info, x):
        if abs(v) < 1e-300: continue
        mono = {m: v / G.order}
        for perm, c in G.perms:
            N[perm[s]] = core.add(N[perm[s]], core.scale(act(mono, perm), c))
    return N, relres

def reconstruct_poly(Dfloat, max_den=1000, tol=1e-7):
    """Rational reconstruction of a float polynomial normalised by its largest |coefficient|."""
    if not Dfloat: return None
    mx = max(abs(v) for v in Dfloat.values())
    P = {}
    for e, v in Dfloat.items():
        x = v / mx
        if abs(x) < tol: continue
        fr = Fraction(x).limit_denominator(max_den)
        if abs(float(fr) - x) > tol: return None
        P[e] = fr
    return P

# ----------------------------------------------------------------------------------------------------------
# folded (H-reduced) float variable projection
# ----------------------------------------------------------------------------------------------------------
class CharVarPro:
    """D(a) = sum a_k M_k (all M_k chi-semi-invariant, same degree/parity); N chi-equivariant; Q invariant.
    Rows = H-orbits of degree 2d+6q monomials of parity pA (plain fold).  float64."""
    def __init__(self, F, A, q, basis, G, exclude_F_multiples=True, verbose=False, cond=1e-9, drop_dead=True):
        t0 = time.time()
        assert q == 1
        self.F, self.A, self.q, self.G = F, A, q, G
        self.cond = cond
        self.basis = [scale_pow2(B) for B in basis]
        K = len(basis); self.K = K
        degs = {core.degree(B) for B in basis}; assert len(degs) == 1
        d = degs.pop(); self.d = d
        pars = {core.parity_of(B) for B in basis}; assert len(pars) == 1 and None not in pars
        pD = pars.pop(); self.pD = pD
        pA = core.parity_of(A); self.pA = pA
        degN = d + 6 * q - 5; degQ = 2 * d - 6; degE = 2 * d + 6 * q
        self.degN, self.degQ = degN, degQ
        for B in self.basis:
            assert G.is_semi_invariant(B), 'family member not chi-semi-invariant'
        Ff = _poly_float(F); Af = _poly_float(A)
        Fi = [core.pderiv(Ff, i) for i in range(6)]
        rowmons = np.array(core.monomials(degE, pA), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); R = len(self.orbit_codes); self.rows = R
        def fold(E, C):
            if E.shape[0] == 0:
                return np.zeros(R)
            idx = np.searchsorted(self.orbit_codes, G.canon(E))
            return np.bincount(idx, weights=C, minlength=R)
        def arr(P):
            if not P:
                return np.zeros((0, 6), dtype=np.int64), np.zeros(0)
            items = list(P.items())
            return np.array([e for e, _ in items], dtype=np.int64), np.array([float(c) for _, c in items])
        # N columns
        Ncols = []
        for s in G.slots:
            sb = G.stab_blocks(s)
            for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                Ncols.append((s, m))
        self.Ncols = Ncols
        blocks = []
        for Bk in self.basis:
            Bf = _poly_float(Bk)
            FB = core.mul(Ff, Bf); FBE, FBC = arr(FB)
            Bi = [core.pderiv(Bf, i) for i in range(6)]
            cols = np.zeros((R, len(Ncols)))
            for s in G.slots:
                G_s = core.scale(core.add(core.scale(core.mul(Bf, Fi[s]), q), core.mul(Ff, Bi[s])), -1)
                GE, GC = arr(G_s)
                for j, (ss, m) in enumerate(Ncols):
                    if ss != s: continue
                    mv = np.array(m, dtype=np.int64)
                    E1 = GE + mv; C1 = GC
                    if m[s]:
                        m2 = mv.copy(); m2[s] -= 1
                        E1 = np.concatenate([E1, FBE + m2]); C1 = np.concatenate([C1, FBC * m[s]])
                    cols[:, j] = fold(E1, C1)
            blocks.append(cols)
        # dead N columns (zero in every block)
        if drop_dead:
            alive = np.zeros(len(Ncols), dtype=bool)
            for cols in blocks:
                alive |= (np.abs(cols).sum(axis=0) > 0)
            self.Ncols = [c for c, a in zip(Ncols, alive) if a]
            blocks = [cols[:, alive] for cols in blocks]
        self.nN = len(self.Ncols)
        self.Mk = [sp.csr_matrix(cols) for cols in blocks]
        del blocks
        # Q columns
        Qcols = []
        if degQ >= 0:
            F2E, F2C = arr(core.scale(core.mul(Ff, Ff), -1))
            Qreps = orbit_reps(degQ, pA, G.canon)
            MQ = np.zeros((R, len(Qreps)))
            for j, m in enumerate(Qreps):
                MQ[:, j] = fold(F2E + np.array(m, dtype=np.int64), F2C)
            Qcols = Qreps
        else:
            MQ = np.zeros((R, 0))
        self.Qcols = Qcols; self.nQ = len(Qcols); self.MQ = MQ
        # rhs pair vectors v_kl = -A M_k M_l folded
        AB = [core.mul(Af, _poly_float(Bk)) for Bk in self.basis]
        self.V = np.zeros((K, K, R))
        for k in range(K):
            for l in range(k, K):
                E_, C_ = arr(core.mul(AB[k], _poly_float(self.basis[l])))
                v = -fold(E_, C_)
                self.V[k, l] = v; self.V[l, k] = v
        self.ncols = self.nN + self.nQ
        self.P = np.eye(K); self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = self._complement_of_F_multiples()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[cvp] d={d} pD={pD} K={K} rows={R} cols={self.ncols} (N {self.nN}, Q {self.nQ}) '
                  f'excl F-mult dim={self.excluded_dim} dimP={self.dimP} built {self.build_time:.1f}s', file=sys.stderr)

    def _complement_of_F_multiples(self):
        K = self.K; dG = self.d - 6
        Gm = core.monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Fcols = [core.shift(_poly_float(self.F), m, 1) for m in Gm]
        mons = {}
        for B in self.basis + Fcols:
            for e in B: mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), K + len(Fcols)))
        for k, B in enumerate(self.basis):
            for e, v in B.items(): Mb[mons[e], k] = float(v)
        for j, C in enumerate(Fcols):
            for e, v in C.items(): Mb[mons[e], K + j] = -float(v)
        ns = sla.null_space(Mb)
        S = ns[:K, :] if ns.size else np.zeros((K, 0))
        if S.shape[1] == 0:
            return np.eye(K), 0
        U, s, _ = np.linalg.svd(S, full_matrices=True)
        r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
        return U[:, r:], r

    def assemble(self, a):
        M = self.Mk[0] * a[0]
        for k in range(1, self.K):
            M = M + self.Mk[k] * a[k]
        M = M.toarray()
        if self.nQ:
            M = np.hstack([M, self.MQ])
        b = np.einsum('k,kl,lr->r', a, np.ones((self.K, self.K)), np.zeros((self.K, self.rows))) if False else np.tensordot(a, np.tensordot(a, self.V, axes=(0, 0)), axes=(0, 0))
        return M, b

    def value(self, a):
        a = np.asarray(a, dtype=float); na = np.linalg.norm(a); a = a / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r

    def value_and_grad(self, u):
        u = np.asarray(u, dtype=float)
        a_full = self.P @ u; na = np.linalg.norm(a_full); a = a_full / na
        M, b = self.assemble(a)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver='gelsd', check_finite=False)
        r = M @ x - b
        f = float(r @ r)
        xN = x[:self.nN]
        Va = np.tensordot(a, self.V, axes=(0, 0))  # K x R : row k = sum_l a_l v_kl
        g = np.zeros(self.K)
        for k in range(self.K):
            dM = self.Mk[k] @ xN
            db = 2.0 * Va[k]
            g[k] = 2.0 * float(r @ (dM - db))
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def bnorm(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, b = self.assemble(a)
        return float(np.linalg.norm(b))

    def is_spurious(self, a, tol=1e-8):
        if self.d < 6:
            return False, 1.0
        a = np.asarray(a, float)
        D = {}
        for k, B in enumerate(self.basis):
            for e, v in B.items():
                D[e] = D.get(e, 0.0) + a[k] * float(v)
        dG = self.d - 6
        Gm = core.monomials(dG, (self.pD - core.parity_of(self.F)) % 2)
        Fcols = [core.shift(_poly_float(self.F), m, 1) for m in Gm]
        mons = {}
        for C in Fcols + [D]:
            for e in C: mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), len(Fcols)))
        for j, C in enumerate(Fcols):
            for e, v in C.items(): Mb[mons[e], j] = float(v)
        bv = np.zeros(len(mons))
        for e, v in D.items(): bv[mons[e]] = v
        x, *_ = sla.lstsq(Mb, bv, check_finite=False)
        rel = float(np.linalg.norm(Mb @ x - bv) / max(np.linalg.norm(bv), 1e-300))
        return rel < tol, rel

    def search(self, n_restarts=10, seed=0, maxiter=300, gtol=1e-10, verbose=False, init=None, stop_rel=1e-12,
               on_candidate=None):
        rng = np.random.default_rng(seed)
        out = []; t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-16})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f, x, r = self.value(a)
            bn = self.bnorm(a)
            rel = f / max(bn * bn, 1e-300)
            spur, srel = self.is_spurious(a)
            cand = {'value': f, 'rel': rel, 'a': a, 'spurious': spur, 'div_rel': srel,
                    'iters': int(res.nit), 'restart': it, 'seed': seed, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if on_candidate: on_candidate(cand)
            if verbose:
                print(f'[cvp] seed {seed} restart {it}: f={f:.3e} rel={rel:.3e} iters={res.nit} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # exact confirmation --------------------------------------------------------------------------------
    def reconstruct(self, a, max_den=200, tol=1e-6):
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a))); a = a / a[j]
        cands = []
        fr = [Fraction(float(v)).limit_denominator(max_den) for v in a]
        if max(abs(float(f) - v) for f, v in zip(fr, a)) < tol:
            cands.append(tuple(fr))
        for den in range(1, max_den + 1):
            v = a * den
            if np.max(np.abs(v - np.round(v))) < tol * den:
                cands.append(tuple(Fraction(int(round(x)), den) for x in v)); break
        uniq = []
        for c in cands:
            if c not in uniq: uniq.append(c)
        return uniq

    def D_exact(self, frac):
        D = {}
        for k, B in enumerate(self.basis):
            if frac[k] == 0: continue
            D = core.add(D, core.scale(B, Fraction(frac[k])))
        return D

    def confirm(self, frac, verbose=False):
        """Exact (reduced, CRT-lifted) re-solve for the rational candidate."""
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        ok, info = cs.e1_char(self.F, self.A, self.q, D, self.G)
        if not ok:
            return 'inconsistent', {'note': 'e1 fails (twisted)', 'e1': info, 'D': D}
        status, data = cs.solve_char(self.F, self.A, self.q, D, self.G, verbose=verbose)
        data['D'] = D
        if status == 'solved':
            okv, vinfo = core.verify_witness(self.F, self.A, self.q, D, data['N'], data['Q'])
            data['verified'] = okv
            divF, _ = jacobian.divisible_by(self.F, D)
            data['F_divides_D'] = divF
            if not okv or divF:
                status = 'spurious' if divF else 'lift-bad'
        return status, data


# ----------------------------------------------------------------------------------------------------------
# Koszul-reduced w-space variable projection (D = sum a_k D_k, N = sum a_k N0_k + K(c)),  E = F * B
# ----------------------------------------------------------------------------------------------------------
PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]

def e1_lift_char(F, A, q, D, G, max_primes=5, verbose=False):
    """Exact chi-equivariant N with N(F) = A D (twisted reduced E1 system, CRT + ratrecon).  Returns list of 6 polys or None."""
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A); degN = d + 6 * q - 5
    acc = None; mod = None; piv0 = None; cols_info = None
    for p in PRIMES[:max_primes]:
        Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        tfold, nrows = cs._tfold_factory(G, d + 6 * q, (pA + pD) % 2, p)
        cols = []; info = []
        for s in G.slots:
            FE, FC = to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
            for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                cols.append(tfold(FE + np.array(m, dtype=np.int64), FC)); info.append((s, m))
        rE, rC = to_arr(core.mul(Ap, Dp, p), p)
        rhs = tfold(rE, rC)
        M = np.concatenate([np.stack(cols, axis=1), rhs[:, None]], axis=1)
        Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
        Mr, rank = Mn.rref()
        C = len(cols); piv = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        if piv and piv[-1] == C:
            return None  # inconsistent
        x = [0] * C
        for i, j in enumerate(piv): x[j] = int(Mr[i, C])
        if piv0 is None or piv0 != piv:
            acc, mod, piv0, cols_info = x, p, piv, info
        else:
            acc = [core.crt_pair(a, mod, b, p)[0] for a, b in zip(acc, x)]; mod *= p
        rat = [core.ratrecon(a, mod) for a in acc]
        if all(v is not None for v in rat):
            N = [{} for _ in range(6)]
            for (s, m), v in zip(cols_info, rat):
                if not v: continue
                mono = {m: Fraction(v) / G.order}
                for perm, c in G.perms:
                    N[perm[s]] = core.add(N[perm[s]], core.scale(act(mono, perm), c))
            # verify N(F) = q^{-1} A D  (E1 normalisation of the engine: q N(F) = A D)
            lhs = {}
            for i in range(6):
                lhs = core.add(lhs, core.mul(N[i], core.pderiv(F, i)))
            lhs = core.scale(lhs, q)
            if core.sub(lhs, core.mul(A, D)):
                if verbose: print('[lift] verification failed; continuing CRT', file=sys.stderr)
                continue
            return N
    return None

def jac_family_with_lifts(F, A, q, d, pD, G):
    """(D_k, N0_k) for a spanning set of J_d^chi: D_k = Sym_chi(m F_s), N0_k = A * Sym_chi(m e_s)."""
    out = []
    if d < 5: return out
    Fi = [core.pderiv(F, i) for i in range(6)]
    for s in G.slots:
        sb = G.stab_blocks(s)
        for m in orbit_reps(d - 5, (pD + (1 if s >= 3 else 0)) % 2, lambda E: G.canon(E, sb)):
            Dk = G.sym_chi(core.shift(Fi[s], m, 1))
            if not Dk: continue
            w = [{} for _ in range(6)]
            for perm, c in G.perms:
                w[perm[s]] = core.add(w[perm[s]], core.scale(act({m: Fraction(1)}, perm), c))
            Nk = [core.mul(A, wi) for wi in w]
            out.append((Dk, Nk))
    return out

class KoszulVarPro:
    """Parameters a over (D_k, N0_k) pairs with N0_k(F) = A D_k (q=1).  N = N0(a) + K(c),
    K(c) = sum_{i<j} c_ij (F_j d_i - F_i d_j), deg c_ij = d-4, chi-equivariant; Q invariant of degree 2d-6.
    Residual B(a,c,Q) = D div N - N(D) - F Q  (degree 2d, parity pA), rows folded over H-orbits.
    koszul=False gives the 'Hessian relaxation' (N = N0(a), only Q free)."""
    def __init__(self, F, A, q, pairs, G, koszul=True, exclude_F_multiples=True, verbose=False, cond=1e-9):
        t0 = time.time(); assert q == 1
        self.F, self.A, self.q, self.G = F, A, q, G; self.cond = cond
        self.basis = []; self.N0 = []
        for Dk, Nk in pairs:
            s = pow2_scale(Dk)
            self.basis.append(core.scale(Dk, s)); self.N0.append([core.scale(x, s) if x else {} for x in Nk])
        K = len(pairs); self.K = K
        d = core.degree(self.basis[0]); self.d = d
        pD = core.parity_of(self.basis[0]); self.pD = pD; pA = core.parity_of(A); self.pA = pA
        degB = 2 * d; degQ = 2 * d - 6; degc = d - 4
        Ff = _poly_float(F); Fi = [core.pderiv(Ff, i) for i in range(6)]
        rowmons = np.array(core.monomials(degB, pA), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); R = len(self.orbit_codes); self.rows = R
        def fold(E, C):
            if E.shape[0] == 0: return np.zeros(R)
            idx = np.searchsorted(self.orbit_codes, G.canon(E))
            return np.bincount(idx, weights=C, minlength=R)
        def arr(P):
            if not P: return np.zeros((0, 6), dtype=np.int64), np.zeros(0)
            items = list(P.items())
            return np.array([e for e, _ in items], dtype=np.int64), np.array([float(c) for _, c in items])
        def foldpoly(P):
            E, C = arr(P); return fold(E, C)
        Bf = [_poly_float(B) for B in self.basis]
        Bi = [[core.pderiv(B, i) for i in range(6)] for B in Bf]
        N0f = [[_poly_float(x) for x in Nk] for Nk in self.N0]
        # pair vectors: v_kl = fold( D_k div N0_l - N0_l(D_k) ), symmetrised over (k,l)
        self.V = np.zeros((K, K, R))
        def Bterm(k, l):
            P = {}
            for i in range(6):
                P = core.add(P, core.mul(Bf[k], core.pderiv(N0f[l][i], i)))
                P = core.sub(P, core.mul(N0f[l][i], Bi[k][i]))
            return P
        for k in range(K):
            for l in range(K):
                self.V[k, l] += foldpoly(Bterm(k, l))
        self.V = -(self.V + np.transpose(self.V, (1, 0, 2))) / 2.0   # b(a) = -sum a_k a_l (..)
        # Koszul columns, per basis k: col(pair (i,j), m) = D_k (m_i F_j - m_j F_i) - m (F_j D_k,i - F_i D_k,j)
        self.ccols = []
        blocks = []
        if koszul and degc >= 0:
            # orbit reps of (pair, m): pairs under H: take reps (0,1),(3,4),(0,3) [+ (1,2)? no: (0,1)~(1,2)] for S3S3;
            # for general blocks enumerate all pairs and dedupe by column equality is expensive; use all pairs and drop dead/duplicate later
            pair_reps = PAIRS
            for (i, j) in pair_reps:
                # c_ij parity: term D c (m_i F_j) has parity pA = pD + pc + par(F_j) + ... ; F_j parity = (j>=3), m_i removes nothing -> pc = pA - pD - par(F_i) - par(F_j)... derive: B parity pA; D c_i F_j: pD + pc + (j>=3) - ... monomial parity of d_i m = par(m) - (i>=3)
                pc = (pA - pD - (1 if j >= 3 else 0) + (1 if i >= 3 else 0)) % 2
                for m in orbit_reps(degc, pc, lambda E: G.canon(E, [[x for x in b if x not in (i, j)] for b in G.blocks])):
                    self.ccols.append((i, j, m))
            for k in range(K):
                cols = np.zeros((R, len(self.ccols)))
                for jj, (i, j, m) in enumerate(self.ccols):
                    mono = {m: 1.0}
                    mi = core.pderiv(mono, i); mj = core.pderiv(mono, j)
                    P = core.sub(core.mul(Bf[k], core.sub(core.mul(mi, Fi[j]), core.mul(mj, Fi[i]))),
                                 core.mul(mono, core.sub(core.mul(Fi[j], Bi[k][i]), core.mul(Fi[i], Bi[k][j]))))
                    cols[:, jj] = foldpoly(P)
                blocks.append(cols)
            alive = np.zeros(len(self.ccols), dtype=bool)
            for cols in blocks: alive |= (np.abs(cols).sum(axis=0) > 0)
            blocks = [cols[:, alive] for cols in blocks]
            self.ccols = [c for c, a in zip(self.ccols, alive) if a]
            # drop duplicate columns (same up to sign across all blocks): hash stacked columns
            if blocks:
                St = np.vstack(blocks)
                keep = []; seen = {}
                for jj in range(St.shape[1]):
                    col = St[:, jj]; nz = np.flatnonzero(col)
                    if nz.size == 0: continue
                    sgn = np.sign(col[nz[0]]); key = tuple(np.round(col * sgn, 6)[nz[:40]]) + (int(nz.size),)
                    if key in seen: continue
                    seen[key] = jj; keep.append(jj)
                blocks = [cols[:, keep] for cols in blocks]
                self.ccols = [self.ccols[jj] for jj in keep]
        self.nC = len(self.ccols)
        self.Mk = [sp.csr_matrix(cols) for cols in blocks] if blocks else [sp.csr_matrix((R, 0)) for _ in range(K)]
        # Q columns: -F m
        Qreps = orbit_reps(degQ, pA, G.canon) if degQ >= 0 else []
        MQ = np.zeros((R, len(Qreps)))
        FE, FC = arr(core.scale(Ff, -1))
        for jj, m in enumerate(Qreps):
            MQ[:, jj] = fold(FE + np.array(m, dtype=np.int64), FC)
        self.Qcols = Qreps; self.nQ = len(Qreps); self.MQ = MQ
        self.ncols = self.nC + self.nQ
        self.P = np.eye(K); self.excluded_dim = 0
        if exclude_F_multiples and d >= 6:
            self.P, self.excluded_dim = CharVarPro._complement_of_F_multiples(self)
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0; self._nfev = 0
        if verbose:
            print(f'[kvp] d={d} pD={pD} K={K} rows={R} cols={self.ncols} (c {self.nC}, Q {self.nQ}) koszul={koszul} '
                  f'excl F-mult dim={self.excluded_dim} dimP={self.dimP} built {self.build_time:.1f}s', file=sys.stderr)

    nN = property(lambda self: self.nC)
    assemble = CharVarPro.assemble
    value = CharVarPro.value
    value_and_grad = CharVarPro.value_and_grad
    bnorm = CharVarPro.bnorm
    is_spurious = CharVarPro.is_spurious
    search = CharVarPro.search
    reconstruct = CharVarPro.reconstruct
    D_exact = CharVarPro.D_exact
    confirm = CharVarPro.confirm


# ----------------------------------------------------------------------------------------------------------
# rational reconstruction with increasing denominators + exact confirmation
# ----------------------------------------------------------------------------------------------------------
def try_reconstruct_confirm(vp, a, dens=(10, 100, 1000, 10000, 100000), tol_scale=50.0, verbose=False, log_tag=''):
    """For each max denominator: continued-fraction reconstruction of a (normalised by its largest entry);
    if the reconstruction error is below tol, confirm exactly.  Returns (status, data, frac) of first solved or last tried."""
    a = np.asarray(a, float); j = int(np.argmax(np.abs(a))); a = a / a[j]
    tried = set(); last = ('none', {}, None)
    for den in dens:
        fr = tuple(Fraction(float(v)).limit_denominator(den) for v in a)
        err = max(abs(float(f) - v) for f, v in zip(fr, a))
        if err > tol_scale / den ** 2 and den < dens[-1]:
            continue
        if fr in tried: continue
        tried.add(fr)
        status, data = vp.confirm(fr, verbose=verbose)
        last = (status, data, fr)
        if verbose: print(f'[recon{log_tag}] den<={den} err={err:.2e} -> {status}', file=sys.stderr)
        if status == 'solved':
            return status, data, fr
    # route (ii): reconstruct the monomial coefficients of D(a) directly
    Df = {}
    for k, B in enumerate(vp.basis):
        for e, v in B.items():
            Df[e] = Df.get(e, 0.0) + a[k] * float(v)
    for den in dens:
        P = reconstruct_poly(Df, max_den=den, tol=max(1e-9, tol_scale / den ** 2))
        if P is None: continue
        key = core.to_str(P)
        if key in tried: continue
        tried.add(key)
        if not vp.G.is_semi_invariant(P): continue
        ok, info = cs.e1_char(vp.F, vp.A, vp.q, P, vp.G)
        if verbose: print(f'[recon{log_tag}] poly-route den<={den}: e1={ok}', file=sys.stderr)
        if not ok: continue
        status, data = cs.solve_char(vp.F, vp.A, vp.q, P, vp.G, verbose=verbose)
        data['D'] = P; last = (status, data, None)
        if status == 'solved':
            okv, _ = core.verify_witness(vp.F, vp.A, vp.q, P, data['N'], data['Q'])
            divF, _ = jacobian.divisible_by(vp.F, P)
            data['verified'] = okv; data['F_divides_D'] = divF
            if okv and not divF:
                return status, data, None
    return last

def save_witness_candidate(D, N, Q, path=f'{LANE}/witness-candidate.json', extra=None):
    rec = {'d': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(x) if x else '0' for x in N], 'Q': core.to_str(Q) if Q else '0'}
    if extra: rec.update(extra)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    return rec
