"""Wave-2 search library for quadric-contact-residual (q = 1, A = X0^3 Y0^3).

Everything here is lane-local (the shared engine in ~/devel/sextics-work is not modified).

Contents
--------
* G_aff characters chi(e) = (e0+e3 mod 2, e2+e5 mod 2, e1-e4 mod 6, e3+e4+e5 mod 2)  (fqlx's bookkeeping;
  the 4th component is the Y-parity).  F has character 0, A has CHI_A = (0,0,0,1).
* e1_family(d, c, vars_mask): EXACT (over Q, python-flint fmpz nullspace) basis of
      (J:A)_d[c] = { D in S_d[c] : A D in J(F)_{d+6} }   (optionally restricted to a variable subset)
  split into the Jacobian part J_d[c] (sparse generators m*F_i) and the non-Jacobian "extras".
* KoszulVarPro: numerical variable-projection search over a linear family D(a) = sum a_k D_k of
  E1-admissible denominators with E1-lifts N_k (N_k(F) = A D_k).  Every E1 solution is
  N = sum a_k N_k + K with K a Koszul field K = sum_{r<s} Q_rs (F_s e_r - F_r e_s), and the
  certificate identity reduces (q = 1) to   B_D(N) = F Q,  B_D(N) = D div N - N(D).
  Unknowns x = (Koszul coefficients Q_rs, quotient Q) enter linearly; the family parameters a enter
  quadratically through D(a) and sum a_k N_k.  f(a) = min_x ||M(a) x - b(a)||^2 is minimized on the
  sphere (Golub-Pereyra gradient).  Character sectors restrict all unknown monomials losslessly when
  D is a G_aff-eigenvector; mode 'parity' / 'mixed' allow mixed characters.
  Spurious zeros (D(a) = 0 directions, D(a) in F*S_{d-6}) are projected out of the parameter space.
* reconstruction helpers: monomial-coordinate rational reconstruction with increasing denominators,
  exact re-solve through the engine (symmetry.solve_fixed_D_symmetric / fixed_d.solve_fixed_D_fast),
  core.verify_witness, coprimality check (gcd with F via univariate specialization).
"""
from __future__ import annotations
import sys, os, json, time, math, itertools, random
os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['MKL_NUM_THREADS'] = '1'
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.optimize import minimize
from hodge import core, fixed_d, jacobian, symmetry
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
import flint
from flint import nmod_mat, fmpz_mat, fmpq_mat

IID = 'quadric-contact-residual'
RES = '/tmp/hodge/results/quadric-contact-residual'
LOG = os.environ.get('HODGE_W2LOG', RES + '/wave2-search.log.jsonl')   # remote jobs set a per-job file (rsync-safe)
P0 = core.PRIMES[0]

F, A, q, REC = core.load_instance(IID)
assert q == 1
Fi = [pderiv(F, i) for i in range(NV)]
Ff = {e: float(v) for e, v in F.items()}
Fif = [{e: float(v) for e, v in P.items()} for P in Fi]
Af = {e: float(v) for e, v in A.items()}
UNIT = [tuple(1 if j == i else 0 for j in range(NV)) for i in range(NV)]


def log_line(rec, path=LOG):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', 2)
    rec.setdefault('host', os.uname().nodename)
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec


# ----------------------------------------------------------------------------- characters
# SECTOR_MODE 'Gaff': 4-component character of G_aff = <sigma0, sigma2, h1, tau_Y> (last component = Y-parity).
# SECTOR_MODE 'H'   : 3-component character of the index-2 subgroup H = <sigma0, sigma2, h1> (no parity):
#                     needed for subspaces that are not tau_Y-stable (e.g. cones with vertex p+ = [1:0:0:1:0:0]).
SECTOR_MODE = 'Gaff'
_MODS = (2, 2, 6, 2)


def set_sector_mode(mode):
    global SECTOR_MODE
    assert mode in ('Gaff', 'H')
    SECTOR_MODE = mode


def chi(e):
    full = ((e[0] + e[3]) % 2, (e[2] + e[5]) % 2, (e[1] - e[4]) % 6, (e[3] + e[4] + e[5]) % 2)
    return full if SECTOR_MODE == 'Gaff' else full[:3]


def chi_add(a, b):
    return tuple((x + y) % m for x, y, m in zip(a, b, _MODS))


def chi_neg(a):
    return tuple((-x) % m for x, m in zip(a, _MODS))


def chi_of(P):
    cs = {chi(e) for e in P}
    return cs.pop() if len(cs) == 1 else None


class _ChiConst:
    """Mode-aware character constants (recomputed on access so set_sector_mode() takes effect)."""
    @property
    def A(self): return chi_of(A)
    @property
    def Z(self): return [chi(u) for u in UNIT]
    @property
    def FI(self): return [chi_neg(c) for c in self.Z]
    @property
    def zero(self): return tuple(0 for _ in chi(UNIT[0]))
_CC = _ChiConst()


class _Lazy:
    def __init__(self, fn): self.fn = fn
    def __getitem__(self, i): return self.fn()[i]
    def __iter__(self): return iter(self.fn())
    def __len__(self): return len(self.fn())
    def __eq__(self, other): return self.fn() == other
    def __hash__(self): return hash(self.fn())
    def __repr__(self): return repr(self.fn())
    def __add__(self, other): return self.fn() + other


def _chiA(): return _CC.A
CHI_Z = _Lazy(lambda: _CC.Z)
CHI_FI = _Lazy(lambda: _CC.FI)
CHI_A = _Lazy(lambda: _CC.A)
CHI0 = (0, 0, 0, 0)


def chars_in_degree(d):
    return sorted({chi(e) for e in monomials(d)})


def chi_str(c):
    return ''.join(str(x) for x in c)


def parse_chi(s):
    s = s.strip().strip('()').replace(',', '').replace(' ', '')
    return tuple(int(x) for x in s)


def mons_sector(d, c=None, vars_mask=None, parity=None):
    out = []
    for e in monomials(d):
        if c is not None and chi(e) != c:
            continue
        if parity is not None and core.yparity(e) != parity:
            continue
        if vars_mask is not None and any(e[j] and not vars_mask[j] for j in range(NV)):
            continue
        out.append(e)
    return out


def jgens(deg, c=None, parity=None):
    """Jacobian generators m*F_i of degree `deg` with character c (or Y-parity)."""
    out = []
    if deg < 5:
        return out
    for i in range(NV):
        for m in monomials(deg - 5):
            if c is not None and chi_add(chi(m), CHI_FI[i]) != c:
                continue
            if parity is not None and (core.yparity(m) + (1 if i >= 3 else 0)) % 2 != parity:
                continue
            out.append((i, m))
    return out


def fpoly(P):
    out = {}
    for e, v in P.items():
        try:
            out[e] = float(v)
        except OverflowError:
            out[e] = float(Fraction(v).limit_denominator(10**9)) if abs(Fraction(v)) < 1e300 else (1e300 if v > 0 else -1e300)
    return out


def to_int_poly(P):
    """Clear denominators of a Fraction poly -> integer poly (primitive)."""
    den = 1
    for v in P.values():
        den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    g = 0
    out = {}
    for e, v in P.items():
        n = int(Fraction(v) * den)
        out[e] = n
        g = math.gcd(g, abs(n))
    if g > 1:
        out = {e: v // g for e, v in out.items()}
    return out


# ----------------------------------------------------------------------------- exact E1 families
def e1_family(d, c, vars_mask=None, verbose=False, parity=None):
    """Exact basis over Q of {D in span(mons) : A D in J_{d+6}} where mons = monomials of degree d with
    character c (or Y-parity `parity` when c is None), optionally supported on a variable subset.

    Returns dict with keys: mons, kernel (list of integer polys, rref basis), dimJ (dim of J_d cap span),
    extras (kernel elements independent modulo J, integer polys), gens (list of (i,m) Jacobian
    generators of degree d in the sector, only when vars_mask is None), dims, seconds."""
    t0 = time.time()
    mons = mons_sector(d, c, vars_mask, parity)
    K = len(mons)
    if c is not None:
        c6 = chi_add(c, CHI_A)
        gens6 = jgens(d + 6, c6)
    else:
        par6 = None if parity is None else (parity + 1) % 2
        gens6 = jgens(d + 6, parity=par6)
    cols = [shift(A, m, 1) for m in mons] + [shift(Fi[i], m, 1) for i, m in gens6]
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    # integer matrix (scale by 7: F has coefficients with denominator 7)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = int(Fraction(v) * 7)
    M = fmpz_mat(R, C, flat)
    t1 = time.time()
    X, nul = M.nullspace()
    t2 = time.time()
    kernel = []
    if nul:
        # project kernel vectors onto the D coordinates, rref over Q, then LLL-reduce the integer basis of the
        # row space (the rref representatives have astronomically large entries at d >= 10; LLL makes them small)
        proj = fmpq_mat(nul, K, [int(X[j, k]) for k in range(nul) for j in range(K)])
        Pr, rk = proj.rref()
        int_rows = []
        for i in range(rk):
            P = {}
            for j in range(K):
                v = Pr[i, j]
                if v != 0:
                    P[mons[j]] = Fraction(int(v.p), int(v.q))
            Pi = to_int_poly(P)
            int_rows.append([Pi.get(m, 0) for m in mons])
        for row in int_rows:
            P = {mons[j]: row[j] for j in range(K) if row[j]}
            if P:
                kernel.append(P)
    # Jacobian part of degree d inside span(mons)
    gens = []
    if vars_mask is None:
        gens = jgens(d, c, parity)
        gpolys = [shift(Fi[i], m, 1) for i, m in gens]
    else:
        gens_all = jgens(d, c, parity)
        gpolys_all = [shift(Fi[i], m, 1) for i, m in gens_all]
        # combinations supported on mons: nullspace of the "outside" rows
        allm = {}
        for P in gpolys_all:
            for e in P:
                allm.setdefault(e, len(allm))
        inside = set(mons)
        outside = [e for e in allm if e not in inside]
        oidx = {e: i for i, e in enumerate(outside)}
        if outside and gpolys_all:
            flat = [0] * (len(outside) * len(gpolys_all))
            for j, P in enumerate(gpolys_all):
                for e, v in P.items():
                    if e in oidx:
                        flat[oidx[e] * len(gpolys_all) + j] = int(Fraction(v) * 7)
            Mo = fmpz_mat(len(outside), len(gpolys_all), flat)
            Y, ny = Mo.nullspace()
            gpolys = []
            for k in range(ny):
                P = {}
                for j, G in enumerate(gpolys_all):
                    cj = int(Y[j, k])
                    if cj:
                        P = add(P, scale(G, Fraction(cj)))
                if P:
                    gpolys.append(P)
        else:
            gpolys = gpolys_all
    midx = {m: i for i, m in enumerate(mons)}
    p = P0

    def vec_mod(P):
        v = [0] * K
        for e, cv in P.items():
            v[midx[e]] = int(Fraction(cv).numerator * pow(Fraction(cv).denominator, -1, p)) % p
        return v
    dimJ = 0
    Jrows = [vec_mod(P) for P in gpolys]
    if Jrows:
        JM = nmod_mat(len(Jrows), K, [x for v in Jrows for x in v], p)
        dimJ = JM.rank()
    # extras: kernel vectors independent modulo J (mod p rank test)
    extras = []
    cur = list(Jrows)
    currank = dimJ
    for P in kernel:
        if len(extras) >= len(kernel) - dimJ:
            break
        test = cur + [vec_mod(P)]
        TM = nmod_mat(len(test), K, [x for v in test for x in v], p)
        rk = TM.rank()
        if rk > currank:
            cur = test; currank = rk; extras.append(P)
    out = {'d': d, 'chi': None if c is None else chi_str(c), 'parity': parity, 'vars_mask': vars_mask,
           'dimS': K, 'rows': R, 'cols': C, 'nullity': nul, 'dim_kernel': len(kernel), 'dimJ': dimJ,
           'n_extra': len(extras), 'dim_gauge': len(mons_sector(d - 6, c, vars_mask, parity)) if d >= 6 else 0,
           'seconds': round(time.time() - t0, 1), 't_null': round(t2 - t1, 1)}
    if verbose:
        print('[e1_family]', {k: v for k, v in out.items() if k != 'vars_mask'}, file=sys.stderr, flush=True)
    return {'info': out, 'mons': mons, 'kernel': kernel, 'extras': extras, 'gens': gens, 'gpolys': gpolys}


def e1_family_span(d, polys, c, verbose=False):
    """Exact kernel of D -> NF(A D) on the span of the given degree-d polynomials (all of character c, or
    c=None for no sector restriction).  Returns (kernel polys as Fraction dicts, info)."""
    t0 = time.time()
    K = len(polys)
    if K == 0:
        return [], {'d': d, 'K': 0, 'dim_kernel': 0}
    c6 = None if c is None else chi_add(c, CHI_A)
    gens6 = jgens(d + 6, c6)
    cols = [mul(A, P) for P in polys] + [shift(Fi[i], m, 1) for i, m in gens6]
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    den = 1
    for col in cols[:K]:
        for v in col.values():
            den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    den = den * 7 // math.gcd(den, 7)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = int(Fraction(v) * den)
    M = fmpz_mat(R, C, flat)
    X, nul = M.nullspace()
    kernel = []
    if nul:
        proj = fmpq_mat(nul, K, [int(X[j, k]) for k in range(nul) for j in range(K)])
        Pr, rk = proj.rref()
        int_rows = []
        for i in range(rk):
            row = [Fraction(int(Pr[i, j].p), int(Pr[i, j].q)) for j in range(K)]
            den = 1
            for v in row:
                den = den * v.denominator // math.gcd(den, v.denominator)
            int_rows.append([int(v * den) for v in row])
        for row in int_rows:
            P = {}
            for j in range(K):
                if row[j]:
                    P = add(P, scale(polys[j], Fraction(row[j])))
            if P:
                kernel.append(to_int_poly(P))
    info = {'d': d, 'K': K, 'rows': R, 'cols': C, 'nullity': nul, 'dim_kernel': len(kernel), 'seconds': round(time.time() - t0, 1)}
    if verbose:
        print('[e1_family_span]', info, file=sys.stderr, flush=True)
    return kernel, info


def numeric_kernel(d, c=None, parity=None, vars_mask=None, polys=None, rcond=1e-10):
    """Float orthonormal basis of {D in span : A D in J_{d+6}} where span = sector monomials (vars_mask) or the
    given polys.  Returns (basis polys (float, unit norm), monomial list / polys, singular-value gap info)."""
    if polys is None:
        mons = mons_sector(d, c, vars_mask, parity)
        base = [{m: 1.0} for m in mons]
    else:
        base = [fpoly(P) for P in polys]
    K = len(base)
    if K == 0:
        return [], base, {}
    c6 = None if c is None else chi_add(c, CHI_A)
    par6 = None if parity is None else (parity + 1) % 2
    gens6 = jgens(d + 6, c6, par6)
    cols = [mul(Af, B) for B in base] + [shift(Fif[i], m, 1.0) for i, m in gens6]
    rows = {}
    for col in cols:
        for e in col:
            rows.setdefault(e, len(rows))
    M = np.zeros((len(rows), len(cols)))
    for j, col in enumerate(cols):
        for e, v in col.items():
            M[rows[e], j] = v
    N = sla.null_space(M, rcond=rcond)
    Pk = N[:K, :]
    if Pk.size == 0:
        return [], base, {'nullity': int(N.shape[1]), 'rank': 0}
    U, s, Vt = np.linalg.svd(Pk, full_matrices=False)
    tol = max(s[0], 1e-300) * 1e-8
    r = int(np.sum(s > tol))
    gap = float(s[r - 1] / s[r]) if r < len(s) and r > 0 and s[r] > 0 else float('inf')
    out = []
    for i in range(r):
        P = {}
        for j, B in enumerate(base):
            cj = U[j, i]
            if abs(cj) > 1e-14:
                for e, v in B.items():
                    P[e] = P.get(e, 0.0) + cj * v
        out.append(P)
    return out, base, {'nullity': int(N.shape[1]), 'rank': r, 'gap': gap, 'smin_kept': float(s[r - 1]) if r else None,
                       'smax_dropped': float(s[r]) if r < len(s) else 0.0}


def numeric_extras(kernel_polys, gens_polys):
    """Orthonormal basis of the part of span(kernel) orthogonal to span(gens) (float): the non-Jacobian directions."""
    mons = {}
    for P in kernel_polys + gens_polys:
        for e in P:
            mons.setdefault(e, len(mons))
    Kmat = np.zeros((len(mons), len(kernel_polys)))
    for j, P in enumerate(kernel_polys):
        for e, v in P.items():
            Kmat[mons[e], j] = v
    if gens_polys:
        G = np.zeros((len(mons), len(gens_polys)))
        for j, P in enumerate(gens_polys):
            for e, v in P.items():
                G[mons[e], j] = v
        Qg, Rg = np.linalg.qr(G, mode='reduced')
        keep = np.abs(np.diag(Rg)) > 1e-10 * max(1.0, np.abs(Rg).max())
        Qg = Qg[:, keep]
        Kmat = Kmat - Qg @ (Qg.T @ Kmat)
    U, s, Vt = np.linalg.svd(Kmat, full_matrices=False)
    r = int(np.sum(s > 1e-8 * max(s[0], 1e-300))) if s.size else 0
    inv = {i: e for e, i in mons.items()}
    out = []
    for i in range(r):
        out.append({inv[j]: float(U[j, i]) for j in range(len(mons)) if abs(U[j, i]) > 1e-14})
    return out


# ----------------------------------------------------------------------------- numerical E1 lifts
class LiftSolver:
    """Least-squares lifts N (deg d+1) with N(F) = T for T of degree d+6 in a sector (float)."""

    def __init__(self, deg6, c=None, parity=None):
        self.gens = jgens(deg6, c, parity)
        self.polys = [shift(Fif[i], m, 1.0) for i, m in self.gens]
        rows = {}
        for P in self.polys:
            for e in P:
                rows.setdefault(e, len(rows))
        self.rows = rows
        M = np.zeros((len(rows), len(self.polys)))
        for j, P in enumerate(self.polys):
            for e, v in P.items():
                M[rows[e], j] = v
        self.M = M
        # pseudo-inverse pieces via SVD once
        U, s, Vt = np.linalg.svd(M, full_matrices=False)
        tol = s[0] * 1e-12
        r = int(np.sum(s > tol))
        self.U, self.s, self.Vt, self.rank = U[:, :r], s[:r], Vt[:r, :], r

    def lift(self, T):
        b = np.zeros(len(self.rows))
        for e, v in T.items():
            if e not in self.rows:
                raise ValueError('target monomial outside the generator span')
            b[self.rows[e]] = float(v)
        x = self.Vt.T @ ((self.U.T @ b) / self.s)
        res = float(np.linalg.norm(self.M @ x - b) / max(np.linalg.norm(b), 1e-300))
        N = [{} for _ in range(NV)]
        for (i, m), v in zip(self.gens, x):
            if v != 0.0:
                N[i][m] = N[i].get(m, 0.0) + float(v)
        return N, res


# ----------------------------------------------------------------------------- <s0, r1> symmetrization
PERM_ID = (0, 1, 2, 3, 4, 5); PERM_S0 = (3, 1, 2, 0, 4, 5); PERM_R1 = (0, 4, 2, 3, 1, 5); PERM_S0R1 = (3, 4, 2, 0, 1, 5)
GRP_SR = [PERM_ID, PERM_S0, PERM_R1, PERM_S0R1]


def perm_poly(P, pi):
    return {tuple(e[pi[j]] for j in range(NV)): c for e, c in P.items()}


def perm_field(N, g):
    return [perm_poly(N[g[i]], g) for i in range(NV)]


def symmetrize_member(Dk, Nk, grp=GRP_SR):
    """Average (D, N) over the coordinate exchanges (F and A are invariant, so lifts stay lifts)."""
    Ds = {}
    Ns = [{} for _ in range(NV)]
    for g in grp:
        Ds = add(Ds, perm_poly(Dk, g))
        gN = perm_field(Nk, g)
        Ns = [add(Ns[i], gN[i]) for i in range(NV)]
    s = 1.0 / len(grp)
    return scale(Ds, s), [scale(P, s) for P in Ns]


# ----------------------------------------------------------------------------- family construction
def build_family(d, c=None, parity=None, vars_mask=None, mode='gens', verbose=False, fam=None, sym=False):
    """Members (D_k, N_k) spanning (J:A)_d[c] (or the masked subspace).

    mode 'gens': sparse Jacobian generators m*F_i (lift A m e_i exact) + extras (numerical lifts).
    mode 'basis': the exact rref kernel basis (dense; numerical lifts for all)."""
    t0 = time.time()
    fam = fam or e1_family(d, c, vars_mask, verbose=verbose, parity=parity)
    members = []   # (D poly float, N list of 6 float polys, tag)
    need_lift = []
    kinfo = {}
    if mode == 'gens' and vars_mask is None:
        gpolys = []
        for i, m in fam['gens']:
            Dk = shift(Fif[i], m, 1.0)
            Nk = [{} for _ in range(NV)]
            Nk[i] = shift(Af, m, 1.0)
            members.append((Dk, Nk, f'J:{i}:{m}'))
            gpolys.append(Dk)
        if fam['info']['n_extra'] > 0:
            kb, _, kinfo = numeric_kernel(d, c, parity, vars_mask)
            need_lift = numeric_extras(kb, gpolys)
            if len(need_lift) != fam['info']['n_extra']:
                print(f'[build_family] WARNING numeric extras {len(need_lift)} != exact {fam["info"]["n_extra"]} (kernel info {kinfo})', file=sys.stderr, flush=True)
    else:
        need_lift, _, kinfo = numeric_kernel(d, c, parity, vars_mask)
        if len(need_lift) != fam['info']['dim_kernel']:
            print(f'[build_family] WARNING numeric kernel {len(need_lift)} != exact {fam["info"]["dim_kernel"]} (kernel info {kinfo})', file=sys.stderr, flush=True)
    lift_res = []
    if need_lift:
        deg6 = d + 6
        c6 = None if c is None else chi_add(c, CHI_A)
        par6 = None if parity is None else (parity + 1) % 2
        LS = LiftSolver(deg6, c6, par6)
        for k, E in enumerate(need_lift):
            T = mul(Af, E)
            Nk, res = LS.lift(T)
            lift_res.append(res)
            members.append((E, Nk, f'E:{k}'))
    if sym:
        # average over <s0, r1>; dependent / vanishing members are removed by the null-space projection later
        sm = []
        seen = set()
        for Dk, Nk, tag in members:
            Ds, Ns = symmetrize_member(Dk, Nk)
            Ds = {e: v for e, v in Ds.items() if abs(v) > 1e-14}
            if not Ds:
                continue
            key = tuple(sorted(Ds))
            if key in seen:
                continue
            seen.add(key)
            sm.append((Ds, Ns, 'sym:' + tag))
        members = sm
    # normalize every member to unit coefficient norm (D_k and its lift N_k by the same factor): the exact
    # kernel extras come as integer polynomials with huge entries, which otherwise wreck the SVD tolerances
    # of parameter_complement (observed at d = 10: 108 of 128 directions wrongly excluded) and the float lstsq.
    normed = []
    for Dk, Nk, tag in members:
        vals = list(Dk.values())
        if all(isinstance(v, (int, Fraction)) for v in vals):
            nrm2 = sum(Fraction(v) ** 2 for v in vals)
            nrm = Fraction(math.isqrt(nrm2.numerator), math.isqrt(nrm2.denominator)) if nrm2 else Fraction(0)
            if nrm == 0:
                continue
            Dk = {e: float(Fraction(v) / nrm) for e, v in Dk.items()}
            Nk = [{e: float(Fraction(v) / nrm) for e, v in P.items()} for P in Nk]
        else:
            nrm = math.sqrt(sum(float(v) ** 2 for v in vals))
            if nrm == 0.0:
                continue
            Dk = scale(Dk, 1.0 / nrm); Nk = [scale(P, 1.0 / nrm) for P in Nk]
        normed.append((Dk, Nk, tag))
    members = normed
    info = dict(fam['info'])
    info.update({'mode': mode, 'sym': sym, 'n_members': len(members), 'n_lifted': len(need_lift), 'numeric_kernel': kinfo,
                 'max_lift_res': max(lift_res) if lift_res else 0.0, 't_family': round(time.time() - t0, 1)})
    if verbose:
        print('[build_family]', {k: v for k, v in info.items() if k != 'vars_mask'}, file=sys.stderr, flush=True)
    return members, info, fam


def family_matrix(members):
    fm = {}
    for Dk, _, _ in members:
        for e in Dk:
            fm.setdefault(e, len(fm))
    Fam = np.zeros((len(fm), len(members)))
    for k, (Dk, _, _) in enumerate(members):
        for e, v in Dk.items():
            Fam[fm[e], k] = v
    return Fam, fm


def parameter_complement(Fam, fmons, d, c=None, parity=None):
    """Orthonormal basis of {a : a in rowspace(Fam) and Fam a orthogonal to F*S_{d-6}}; returns (P, n_excluded)."""
    K = Fam.shape[1]
    U, s, Vt = np.linalg.svd(Fam, full_matrices=True)
    tol = max(1.0, s[0]) * 1e-10 if s.size else 0.0
    r = int(np.sum(s > tol))
    Null = Vt[r:, :]              # rows span the null space of Fam
    constraints = [Null]
    dG = d - 6
    if dG >= 0:
        if c is not None:
            Gm = [m for m in monomials(dG) if chi(m) == c]
        elif parity is not None:
            Gm = monomials(dG, parity)
        else:
            Gm = monomials(dG)
        Fcols = [shift(Ff, m, 1.0) for m in Gm]
        fm = fmons
        G = np.zeros((len(fm), len(Fcols)))
        missing = 0
        for j, C in enumerate(Fcols):
            for e, v in C.items():
                if e in fm:
                    G[fm[e], j] = v
                else:
                    missing += 1
        if Fcols and missing == 0:
            # F-multiples that lie INSIDE the family span: joint null space of [Fam | -G] -> s-components
            NS = sla.null_space(np.hstack([Fam, -G]), rcond=1e-10)
            Sg = NS[K:, :]
            if Sg.size and np.linalg.norm(Sg) > 1e-12:
                Gin = G @ Sg                           # in-family F-multiples (monomial coordinates)
                Qg, Rg = np.linalg.qr(Gin, mode='reduced')
                keep = np.abs(np.diag(Rg)) > 1e-10 * max(1.0, np.abs(Rg).max())
                Qg = Qg[:, keep]
                if Qg.shape[1]:
                    constraints.append(Qg.T @ Fam)     # monomial-orthogonality of D(a) to those multiples
        elif Fcols:
            # some F-multiples have monomials outside the family's support (masked families): only the
            # multiples fully inside the span matter; handle through the joint null space on the union support
            allm = dict(fm)
            for C in Fcols:
                for e in C:
                    allm.setdefault(e, len(allm))
            FamU = np.zeros((len(allm), K)); FamU[:len(fm), :] = Fam
            GU = np.zeros((len(allm), len(Fcols)))
            for j, C in enumerate(Fcols):
                for e, v in C.items():
                    GU[allm[e], j] = v
            NS = sla.null_space(np.hstack([FamU, -GU]), rcond=1e-10)
            Sg = NS[K:, :]
            if Sg.size and np.linalg.norm(Sg) > 1e-12:
                Gin = (GU @ Sg)[:len(fm), :]
                Qg, Rg = np.linalg.qr(Gin, mode='reduced')
                keep = np.abs(np.diag(Rg)) > 1e-10 * max(1.0, np.abs(Rg).max())
                Qg = Qg[:, keep]
                if Qg.shape[1]:
                    constraints.append(Qg.T @ Fam)
    Cm = np.vstack([c for c in constraints if c.size]) if any(c.size for c in constraints) else np.zeros((0, K))
    if Cm.shape[0] == 0:
        return np.eye(K), 0
    P = sla.null_space(Cm, rcond=1e-10)
    return P, K - P.shape[1]


# ----------------------------------------------------------------------------- Koszul variable projection
class KoszulVarPro:
    def __init__(self, members, d, c=None, parity=None, verbose=False, cond=1e-9, with_Q=True, driver='gelsy'):
        t0 = time.time()
        self.cond = cond
        self.driver = driver
        self.d = d; self.c = c; self.parity = parity
        self.members = members
        K = len(members); self.K = K
        degK = d - 4           # Koszul coefficients Q_rs
        degQ = 2 * d - 6       # quotient
        pD = 0 if parity is None else parity
        if c is not None and len(c) == 4:
            pD = c[3]
        # unknown columns
        Kcols = []
        for r in range(NV):
            for s in range(r + 1, NV):
                if degK < 0:
                    continue
                if c is not None:
                    req = chi_add(CHI_A, chi_add(c, chi_add(CHI_Z[r], CHI_Z[s])))
                    ms = [m for m in monomials(degK) if chi(m) == req]
                elif parity is not None:
                    req = (1 + pD + (1 if r >= 3 else 0) + (1 if s >= 3 else 0)) % 2
                    ms = monomials(degK, req)
                else:
                    ms = monomials(degK)
                Kcols += [(r, s, m) for m in ms]
        if degQ >= 0 and with_Q:
            if c is not None:
                reqQ = chi_add(CHI_A, chi_add(c, c))
                Qcols = [m for m in monomials(degQ) if chi(m) == reqQ]
            elif parity is not None:
                Qcols = monomials(degQ, 1)
            else:
                Qcols = monomials(degQ)
        else:
            Qcols = []
        self.Kcols, self.Qcols = Kcols, Qcols
        nK, nQ = len(Kcols), len(Qcols)
        self.nK, self.nQ = nK, nQ
        rows = {}

        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        # blocks: for member k, column (r,s,m): D_k (d_r m F_s - d_s m F_r) - m (F_s d_r D_k - F_r d_s D_k)
        bri, bci, bval = [], [], []     # stacked: flat row index = row*nK + col ; column = k
        for k, (Dk, Nk, tag) in enumerate(members):
            Dki = [pderiv(Dk, i) for i in range(NV)]
            DF = [mul(Dk, Fif[s]) for s in range(NV)]
            cache_Q = {}
            for j, (r, s, m) in enumerate(Kcols):
                col = {}
                if m[r]:
                    m2 = list(m); m2[r] -= 1
                    col = add(col, shift(DF[s], tuple(m2), float(m[r])))
                if m[s]:
                    m2 = list(m); m2[s] -= 1
                    col = sub(col, shift(DF[r], tuple(m2), float(m[s])))
                if (r, s) not in cache_Q:
                    cache_Q[(r, s)] = sub(mul(Fif[s], Dki[r]), mul(Fif[r], Dki[s]))
                col = sub(col, shift(cache_Q[(r, s)], m, 1.0))
                for e, v in col.items():
                    if v != 0.0:
                        bri.append(ridx(e) * nK + j); bci.append(k); bval.append(v)
        # Q block: -F m
        qri, qci, qval = [], [], []
        for j, m in enumerate(Qcols):
            for e, v in shift(Ff, m, -1.0).items():
                qri.append(ridx(e)); qci.append(j); qval.append(v)
        # pair vectors: v_kl = -B_{D_k}(N_l) = -(D_k div N_l - N_l(D_k)), symmetrized
        pair = {}
        for l, (Dl, Nl, _) in enumerate(members):
            divN = {}
            for i in range(NV):
                if Nl[i]:
                    divN = add(divN, pderiv(Nl[i], i))
            for k, (Dk, Nk, _) in enumerate(members):
                Dki = None
                Bv = mul(Dk, divN)
                for i in range(NV):
                    if Nl[i]:
                        Bv = sub(Bv, mul(Nl[i], pderiv(Dk, i)))
                if Bv:
                    pair[(k, l)] = Bv
                    for e in Bv:
                        ridx(e)
        R = len(rows); self.rows = R; self.row_index = rows
        self.Mstack = sp.csr_matrix((np.array(bval), (np.array(bri, dtype=np.int64), np.array(bci, dtype=np.int64))),
                                    shape=(R * nK, K)) if nK else None
        self.MQ = np.zeros((R, nQ))
        if nQ:
            self.MQ[np.array(qri, dtype=np.int64), np.array(qci, dtype=np.int64)] = np.array(qval)
        # orthonormal basis of range(MQ) (multiplication by F is injective, so MQ has full column rank);
        # the Q-unknowns are eliminated by projecting onto the orthogonal complement.
        if nQ:
            Qb, Rq = np.linalg.qr(self.MQ, mode='reduced')
            keep = np.abs(np.diag(Rq)) > 1e-10 * max(1.0, np.abs(Rq).max())
            self.Qb = Qb[:, keep]
        else:
            self.Qb = np.zeros((R, 0))
        vri, vci, vval = [], [], []
        for (k, l), Bv in pair.items():
            for e, v in Bv.items():
                # symmetrized: contributes -v/2 to (k,l) and (l,k)
                vri.append(rows[e] * K + k); vci.append(l); vval.append(-0.5 * v)
                vri.append(rows[e] * K + l); vci.append(k); vval.append(-0.5 * v)
        self.Vstack = sp.csr_matrix((np.array(vval), (np.array(vri, dtype=np.int64), np.array(vci, dtype=np.int64))),
                                    shape=(R * K, K))
        self.ncols = nK + nQ
        # family matrix in monomial coordinates (for null directions / gauge / reconstruction)
        self.Fam, self.fmons = family_matrix(members)
        self.P, self.excluded = self._parameter_complement()
        self.dimP = self.P.shape[1]
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[koszul-varpro] d={d} chi={None if c is None else chi_str(c)} K={K} rows={R} '
                  f'cols={self.ncols} (K {nK}, Q {nQ}) Mstack nnz={0 if self.Mstack is None else self.Mstack.nnz} '
                  f'Vstack nnz={self.Vstack.nnz} dimP={self.dimP} excluded={self.excluded} built {self.build_time:.1f}s',
                  file=sys.stderr, flush=True)

    # --------------------------------------------------------------
    def _parameter_complement(self):
        return parameter_complement(self.Fam, self.fmons, self.d, self.c, self.parity)

    # --------------------------------------------------------------
    def assemble(self, a):
        """Projected K-block P_Q M_K(a), projected rhs P_Q b(a), and W (db/da = 2W)."""
        M = (self.Mstack @ a).reshape(self.rows, self.nK) if self.nK else np.zeros((self.rows, 0))
        W = (self.Vstack @ a).reshape(self.rows, self.K)
        b = W @ a
        if self.Qb.shape[1]:
            M = M - self.Qb @ (self.Qb.T @ M)
            b = b - self.Qb @ (self.Qb.T @ b)
        return M, b, W

    def _solve(self, M, b):
        if M.shape[1] == 0:
            return np.zeros(0)
        x, *_ = sla.lstsq(M, b, cond=self.cond, lapack_driver=self.driver, check_finite=False)
        return x

    def value(self, a):
        """f(a/|a|) = min_x ||P_Q (M_K(a) x - b(a))||^2 ; returns (f, xK, r, b_projected)."""
        a = np.asarray(a, float)
        na = np.linalg.norm(a)
        a = a / na
        M, b, W = self.assemble(a)
        x = self._solve(M, b)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r, b

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M, b, W = self.assemble(a)
        x = self._solve(M, b)
        r = M @ x - b
        f = float(r @ r)
        g = np.zeros(self.K)
        if self.nK:
            xK = x[:self.nK]
            g += 2.0 * (self.Mstack.T @ np.outer(r, xK).ravel())
        g -= 4.0 * (W.T @ r)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def rel(self, a):
        """Relative residual f/|b|^2 with b the unprojected quadratic data -B_D(N)."""
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        f, x, r, bp = self.value(a)
        W = (self.Vstack @ a).reshape(self.rows, self.K)
        b = W @ a
        bn = float(b @ b)
        return f / max(bn, 1e-300), f, bn, x

    def D_of(self, a):
        a = np.asarray(a, float)
        D = {}
        for k, (Dk, _, _) in enumerate(self.members):
            if a[k] == 0.0:
                continue
            for e, v in Dk.items():
                D[e] = D.get(e, 0.0) + a[k] * v
        return D

    def is_spurious(self, a, tol=1e-8):
        """Numerically: is D(a) divisible by F?  (relative residual of the LS division)."""
        if self.d < 6:
            return False, 1.0
        D = self.D_of(a)
        dG = self.d - 6
        Gm = monomials(dG, core.parity_of({e: 1 for e in D}) if core.parity_of({e: 1 for e in D}) is not None else None)
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
        rel = float(np.linalg.norm(Mb @ x - bv) / max(np.linalg.norm(bv), 1e-300))
        return rel < tol, rel

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
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            rel, f, bn, x = self.rel(a)
            spur, srel = self.is_spurious(a)
            cand = {'rel': rel, 'f': f, 'bnorm2': bn, 'a': a, 'spurious': spur, 'div_rel': srel,
                    'iters': int(res.nit), 'restart': it, 'seed': seed, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[search seed={seed}] restart {it}: rel={rel:.3e} f={f:.3e} iters={res.nit} spurious={spur} '
                      f'({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
            if rel < stop_rel and not spur:
                break
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out


# ----------------------------------------------------------------------------- reconstruction / exact re-solve
def reconstruct_monomial(D, max_den):
    """Rational reconstruction of a float poly after normalizing the largest coefficient to 1."""
    if not D:
        return None
    e0 = max(D, key=lambda e: abs(D[e]))
    scale0 = D[e0]
    out = {}
    err = 0.0
    for e, v in D.items():
        x = v / scale0
        fr = Fraction(x).limit_denominator(max_den)
        err = max(err, abs(float(fr) - x))
        if fr != 0:
            out[e] = fr
    return out, err


def coprime_with_F(D, tries=20, seed=1):
    """Heuristic exact coprimality test: specialize all variables but one at random integer points
    (several times); if some univariate gcd is 1 then gcd(D, F) = 1 over Q[Z] (F is monic in each Z_j)."""
    import sympy as spy
    x = spy.symbols('x')
    rng = random.Random(seed)
    for j in range(NV):
        for _ in range(tries):
            pt = [rng.randint(-5, 5) for _ in range(NV)]
            def uni(P):
                poly = 0
                for e, cv in P.items():
                    t = spy.Rational(Fraction(cv).numerator, Fraction(cv).denominator)
                    for i in range(NV):
                        if e[i]:
                            t *= (x if i == j else spy.Integer(pt[i])) ** e[i]
                    poly += t
                return spy.Poly(poly, x)
            fF, fD = uni(F), uni(D)
            if fD.is_zero:
                continue
            g = spy.gcd(fF, fD)
            if g.degree() == 0:
                return True, {'j': j, 'point': pt}
    return False, {}


def exact_resolve(D, verbose=False, prefer_symmetric=True):
    """Exact fixed-D solve through the engine.  Returns (status, data)."""
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    if not ok:
        return 'e1-fail', {'e1': info, 'seconds': round(time.time() - t0, 1)}
    status, data = None, None
    if prefer_symmetric and chi_of(D) is not None:
        try:
            status, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, verbose=verbose)
        except Exception as ex:
            if verbose:
                print('[exact_resolve] symmetric solver failed:', ex, file=sys.stderr)
            status = None
    if status is None:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
        if not out['consistent']:
            status, data = 'inconsistent', out
        elif out.get('exact'):
            status, data = 'solved', out
        else:
            status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
    data = dict(data or {})
    data['seconds'] = round(time.time() - t0, 1)
    if status == 'solved':
        ok, vinfo = core.verify_witness(F, A, q, D, data['N'], data['Q'])
        data['verified'] = bool(ok)
        data['verify_info'] = str(vinfo)[:200]
        cp, cinfo = coprime_with_F(D)
        data['coprime'] = bool(cp)
        data['coprime_info'] = cinfo
    return status, data


def save_witness_candidate(D, N, Q, path=RES + '/witness-candidate.json', extra=None):
    rec = {'instanceId': IID, 'd': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(P) for P in N],
           'Q': core.to_str(Q), 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    if extra:
        rec.update(extra)
    with open(path, 'w') as f:
        json.dump(rec, f)
    return path


def try_reconstruct_and_confirm(vp, a, dens=(1, 2, 6, 12, 60, 420, 10**3, 10**4, 10**5, 10**6), verbose=False, tag=''):
    """For a low-residual parameter vector a: reconstruct D(a) rationally with increasing denominators
    and exactly re-solve.  Returns list of attempt records; writes witness-candidate.json on success."""
    D = vp.D_of(a)
    attempts = []
    tried = set()
    for md in dens:
        rec_, err = reconstruct_monomial(D, md)
        key = core.to_str(rec_)
        if key in tried:
            continue
        tried.add(key)
        Dq = {e: Fraction(v) for e, v in rec_.items()}
        if core.parity_of(Dq) is None and chi_of(Dq) is None:
            pass
        t0 = time.time()
        status, data = exact_resolve(Dq, verbose=verbose)
        att = {'max_den': md, 'recon_err': err, 'status': status, 'seconds': round(time.time() - t0, 1),
               'D_len': len(Dq), 'tag': tag}
        for k in ('rank', 'nullity', 'rows', 'cols', 'verified', 'coprime'):
            if k in data:
                att[k] = data[k]
        attempts.append(att)
        if verbose:
            print('[reconstruct]', att, file=sys.stderr, flush=True)
        if status == 'solved' and data.get('verified'):
            path = save_witness_candidate(Dq, data['N'], data['Q'], extra={'coprime': data.get('coprime'), 'tag': tag})
            att['witness_path'] = path
            print('WITNESS FOUND (exact re-solve verified) ->', path, file=sys.stderr, flush=True)
            break
    return attempts
