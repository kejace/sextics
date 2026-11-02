"""Wave-2 exclusion/rigor library for lane reflection-coupled-residual.

Klein group (Z/2)^2 = {1, r, rho, tau}: character of a monomial m is
    chi(m) = (m5 mod 2, (m3+m4) mod 2)            (Y2-parity, (Y0+Y1)-parity)
F has chi = (0,0); A = X0^3 Y0^3 has chi = (0,1).  chi(F_i) = chi(z_i).

All linear algebra is dense python-flint nmod_mat mod p.  Every receipt records the
prime and the good-reduction check (rank of J_n^chi mod p == dim S_n^chi - Fermat count),
which is what makes "kernel dim mod p >= kernel dim over Q" rigorous.
"""
import sys, time, json, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv, mul, add, sub, scale
from flint import nmod_mat
import numpy as np

IID = 'reflection-coupled-residual'
WORK = f'/tmp/hodge/results/{IID}'
OUT = f'{WORK}/wave2-excl'
F, A, q, rec = core.load_instance(IID)
if os.environ.get('W2LIB_F'):          # positive-control override (e.g. the Fermat sextic); A, q unchanged
    F = core.parse(os.environ['W2LIB_F'])
    OUT = f'{WORK}/wave2-excl/control-' + os.environ.get('W2LIB_TAG', 'alt')
    os.makedirs(OUT, exist_ok=True)
PRIMES = core.PRIMES
CHARS = [(0, 0), (0, 1), (1, 0), (1, 1)]
CHI_Z = [(0, 0), (0, 0), (0, 0), (0, 1), (0, 1), (1, 0)]
CHI_A = (0, 1)


def chi(m):
    return (m[5] & 1, (m[3] + m[4]) & 1)


def cadd(a, b):
    return ((a[0] + b[0]) & 1, (a[1] + b[1]) & 1)


def mons_chi(n, c):
    """Monomials of degree n with character c (c=None: all)."""
    return [m for m in monomials(n) if c is None or chi(m) == c]


def fermat_count(n, c):
    """dim (S/J)_n^c for the Fermat sextic (all exponents <= 4) = expected dim for any smooth sextic."""
    cnt = 0
    for e in itertools.product(range(5), repeat=NV):
        if sum(e) == n and (c is None or chi(e) == c):
            cnt += 1
    return cnt


def receipt(job, **kw):
    """Append one JSON line to OUT/<job>.jsonl and to the lane-level wave2 log; print it."""
    kw.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    kw.setdefault('job', job)
    line = json.dumps(kw, default=str)
    os.makedirs(OUT, exist_ok=True)
    with open(f'{OUT}/{job}{os.environ.get("W2LIB_SUFFIX", "")}.jsonl', 'a') as f:
        f.write(line + '\n')
    print('RECEIPT', line, flush=True)


def _mat(rows_polys, idx, C, p):
    R = len(rows_polys)
    flat = [0] * (R * C)
    for r, poly in enumerate(rows_polys):
        base = r * C
        for e, v in poly.items():
            flat[base + idx[e]] = v
    return nmod_mat(R, C, flat, p)


class CharReducer:
    """J(F)_n restricted to character c (mod p): rref, pivots, class coordinates.

    cls(poly) -> list of `nul` ints: coordinates of [poly] in (S/J)_n^c (zero iff poly in J).
    nf(poly)  -> canonical normal form (pivot monomials eliminated) as a dict poly.
    """

    def __init__(self, n, c, p, need_W=True, need_nf=False, verbose=True):
        t0 = time.time()
        self.n, self.c, self.p = n, c, p
        Fp = core.reduce_mod(F, p)
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        self.Fi = Fi
        self.mons = mons_chi(n, c)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        C = len(self.mons)
        gens = []
        if n >= 5:
            for i in range(NV):
                ci = CHI_Z[i]
                for m in monomials(n - 5):
                    if cadd(chi(m), ci) == c:
                        gens.append(shift(Fi[i], m, 1, p))
        self.ngens = len(gens)
        self.fermat = fermat_count(n, c)
        if gens:
            M = _mat(gens, self.idx, C, p)
            if need_W:
                W, nul = M.nullspace()
                self.nul = nul
                self.rank = C - nul
                self.W = nmod_mat(C, nul, [int(W[i, k]) for i in range(C) for k in range(nul)], p)
            else:
                self.rank = M.rank(); self.nul = C - self.rank; self.W = None
            if need_nf:
                Mr, rk = M.rref()
                assert rk == self.rank
                piv = []
                for i in range(rk):
                    j = 0
                    while int(Mr[i, j]) == 0:
                        j += 1
                    piv.append(j)
                self.pivots = piv
                self.rref_rows = [{k: int(Mr[i, k]) for k in range(C) if int(Mr[i, k])} for i in range(rk)]
                self.standard = [m for j, m in enumerate(self.mons) if j not in set(piv)]
        else:
            self.rank, self.nul = 0, C
            self.W = nmod_mat(C, C, [1 if i == j else 0 for i in range(C) for j in range(C)], p) if need_W else None
            self.pivots, self.rref_rows, self.standard = [], [], list(self.mons)
        self.good = (self.rank == C - self.fermat)
        self.dimS = C
        self.t = round(time.time() - t0, 1)
        if verbose:
            print(f'[reducer n={n} chi={c} p={p}] dim S={C} gens={self.ngens} rank J={self.rank} '
                  f'nul={self.nul} fermat={self.fermat} good_reduction={self.good} ({self.t}s)', flush=True)

    def row(self, poly):
        r = [0] * self.dimS
        for e, v in poly.items():
            if chi(e) != self.c or sum(e) != self.n:
                raise ValueError(f'monomial {e} not in degree {self.n} character {self.c}')
            r[self.idx[e]] = v % self.p
        return r

    def cls_many(self, polys):
        """nmod_mat (len(polys) x nul) of class coordinates."""
        R = len(polys)
        flat = []
        for P in polys:
            flat.extend(self.row(P))
        M = nmod_mat(R, self.dimS, flat, self.p)
        return M * self.W

    def cls(self, poly):
        v = self.cls_many([poly])
        return [int(v[0, k]) for k in range(self.nul)]

    def nf(self, poly):
        """Canonical normal form (dict) using the rref rows; requires need_nf=True."""
        v = self.row(poly)
        p = self.p
        for j, r in zip(self.pivots, self.rref_rows):
            if v[j]:
                cf = v[j]
                for k, a in r.items():
                    v[k] = (v[k] - cf * a) % p
        return {self.mons[k]: v[k] for k in range(self.dimS) if v[k]}


def e1_kernel(d, c, p, red_target=None, verbose=True):
    """Kernel of [A]: (S/J)_d^c -> (S/J)_{d+6}^{c+chi(A)} mod p.

    Returns dict with dims, ranks, good-reduction flags, kernel dim (in S_d^c) and modulo J_d^c,
    and canonical representatives (normal forms w.r.t. J_d^c rref) of the kernel classes."""
    t0 = time.time()
    Ap = core.reduce_mod(A, p)
    ct = cadd(c, CHI_A)
    RT = red_target or CharReducer(d + 6, ct, p, verbose=verbose)
    Dm = mons_chi(d, c)
    K = RT.cls_many([shift(Ap, m, 1, p) for m in Dm])     # |Dm| x nul
    # kernel: vectors x (over Dm) with x^T K = 0  -> nullspace of K^T
    X, nullity = K.transpose().nullspace()
    kers = []
    for k in range(nullity):
        kers.append({Dm[j]: int(X[j, k]) for j in range(len(Dm)) if int(X[j, k])})
    out = {'d': d, 'chi': c, 'p': p, 'dim_Sd': len(Dm), 'target_chi': ct, 'dim_S_target': RT.dimS,
           'rankJ_target': RT.rank, 'fermat_target': RT.fermat, 'good_target': RT.good,
           'dim_RJ_target': RT.nul, 'rank_A': len(Dm) - nullity, 'e1_kernel_dim': nullity}
    if d >= 5:
        Rd = CharReducer(d, c, p, need_W=False, need_nf=True, verbose=verbose)
        out.update({'rankJ_d': Rd.rank, 'fermat_d': Rd.fermat, 'good_d': Rd.good, 'dim_RJ_d': Rd.nul})
        nfs = [Rd.nf(v) for v in kers]
        nfs = [v for v in nfs if v]
        reps = []
        if nfs:
            M = _mat(nfs, Rd.idx, Rd.dimS, p)
            Mr, rk = M.rref()
            for i in range(rk):
                reps.append({Rd.mons[j]: int(Mr[i, j]) for j in range(Rd.dimS) if int(Mr[i, j])})
        out['kernel_mod_J'] = len(reps)
        out['reps'] = reps
    else:
        out.update({'rankJ_d': 0, 'good_d': True, 'dim_RJ_d': len(Dm), 'kernel_mod_J': nullity, 'reps': kers})
    out['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(f'[e1 d={d} chi={c} p={p}] dim S_d^c={len(Dm)} rankA={out["rank_A"]} ker={nullity} '
              f'ker mod J={out["kernel_mod_J"]} good=({out["good_d"]},{RT.good}) ({out["seconds"]}s)', flush=True)
    return out


def sym_str(P, p):
    return core.to_str({e: (v if v <= p // 2 else v - p) for e, v in P.items()})


def ratrecon_polys(polys_by_prime):
    """polys_by_prime: list of (p, poly dict mod p) with the same canonical form.  CRT + ratrecon."""
    keys = set()
    for _, P in polys_by_prime:
        keys |= set(P)
    out = {}
    ok = True
    for e in sorted(keys):
        x, m = polys_by_prime[0][1].get(e, 0), polys_by_prime[0][0]
        for p, P in polys_by_prime[1:]:
            x, m = core.crt_pair(x, m, P.get(e, 0), p)
        r = core.ratrecon(x, m)
        if r is None:
            ok = False
        out[e] = r
    return ok, out


# ---------------------------------------------------------------- Hessian relaxation

def hessian(p):
    Fp = core.reduce_mod(F, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    H = {}
    for i in range(NV):
        for j in range(i, NV):
            H[(i, j)] = pderiv(Fi[i], j, p)
    return Fi, H


def gbasis(d, c=None):
    """Basis (i, m) of g-space at degree d: g_i in S_{d-5}; character of g_i d_i is chi(m)+chi(z_i)."""
    out = []
    for i in range(NV):
        for m in monomials(d - 5):
            if c is None or cadd(chi(m), CHI_Z[i]) == c:
                out.append((i, m))
    return out


def gchar(i, m):
    return cadd(chi(m), CHI_Z[i])


def trivial_g(d, c=None):
    """Vectors (as dict {(i,m): coeff}) spanning T = {z G : G in S_{d-6}} (character c of z G = chi(G))."""
    out = []
    if d < 6:
        return out
    for G in monomials(d - 6):
        if c is not None and chi(G) != c:
            continue
        v = {}
        for i in range(NV):
            m = list(G); m[i] += 1
            v[(i, tuple(m))] = 1
        out.append(v)
    return out


def bil_poly(ga, gb, H, Ap, p):
    """A * sum_{i,j} ga_i gb_j F_ij  (mod p); ga, gb: list of 6 polys."""
    acc = {}
    for i in range(NV):
        if not ga[i]:
            continue
        for j in range(NV):
            if not gb[j]:
                continue
            Hij = H[(min(i, j), max(i, j))]
            if not Hij:
                continue
            acc = add(acc, mul(mul(ga[i], gb[j], p), Hij, p), p)
    return mul(Ap, acc, p)
