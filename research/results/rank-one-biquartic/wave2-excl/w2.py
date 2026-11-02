"""Wave-2 exclusion library for lane rank-one-biquartic (q=1, A = X0^3 Y0^3).

Independent code path from the wave-1 scripts (lib.py/e1kernel.py) and from the dossier's
rank_one_compute.py: polynomial arithmetic with python-flint nmod_mpoly, linear algebra with nmod_mat.

Conventions: variables X0,X1,X2,Y0,Y1,Y2 (indices 0..5); Y-parity of a monomial = (e3+e4+e5) % 2.
Certificate (q=1): D (deg d), N_i (deg d+1), Q (deg 2d-6):
    E1:  N(F) = A D            (sum_i N_i F_i = A D)
    E2:  B_D(N) := D div N - N(D) = F C,  C = Q  (deg 2d-6)
Koszul corrections: K = sum_{i<j} c_ij (F_j d_i - F_i d_j), deg c_ij = d-4, B_D(K) = T_D(c) in J(F).
Hessian relaxation (necessary): for D in the E1 kernel with lift n, B_D(n) in J(F)_{2d}.
"""
import sys, json, time, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from flint import nmod_mpoly_ctx, nmod_mpoly, nmod_mat
from fractions import Fraction

NV = 6
NAMES = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')
IID = 'rank-one-biquartic'
WORK = '/tmp/hodge/results/rank-one-biquartic'
LOG = WORK + '/wave2-excl.log.jsonl'
PACKET = WORK + '/wave2-excl'


def log(entry):
    entry.setdefault('ts', time.strftime('%Y-%m-%d %H:%M:%S'))
    entry.setdefault('lens', 'excl')
    entry.setdefault('engine', 'w2.py (nmod_mpoly/nmod_mat), independent of hodge.fixed_d')
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(entry, default=str) + '\n')


def yparity(e):
    return (e[3] + e[4] + e[5]) % 2


def monomials(deg, par=None, pred=None):
    """Exponent tuples of total degree deg with given Y-parity (None = both) and optional predicate."""
    out = []
    if deg < 0:
        return out
    for c in itertools.combinations(range(deg + NV - 1), NV - 1):
        e = []
        prev = -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + NV - 2 - prev)
        e = tuple(e)
        if par is not None and yparity(e) != par:
            continue
        if pred is not None and not pred(e):
            continue
        out.append(e)
    return out


class Ring:
    def __init__(self, p):
        self.p = p
        self.ctx = nmod_mpoly_ctx.get(NAMES, ordering='degrevlex', modulus=p)
        self.gens = self.ctx.gens()

    def from_dict(self, P):
        """engine dict poly (Fraction/int coefficients) -> nmod_mpoly mod p."""
        p = self.p
        d = {}
        for e, c in P.items():
            if isinstance(c, Fraction):
                v = (c.numerator % p) * pow(c.denominator % p, -1, p) % p
            else:
                v = int(c) % p
            if v:
                d[tuple(e)] = v
        return self.ctx.from_dict(d)

    def parse(self, s):
        return self.from_dict(core.parse(s))

    def zero(self):
        return self.ctx.from_dict({})

    def mono(self, e, c=1):
        return self.ctx.from_dict({tuple(e): c % self.p})

    def deriv(self, P, i):
        return P.derivative(i)

    def div(self, vec):
        """divergence of a vector field (list of 6 polys)."""
        s = self.zero()
        for i in range(NV):
            s += vec[i].derivative(i)
        return s

    def apply(self, vec, P):
        """vec(P) = sum_i vec_i dP/dz_i."""
        s = self.zero()
        for i in range(NV):
            s += vec[i] * P.derivative(i)
        return s

    def B(self, D, N):
        """B_D(N) = D div N - N(D)."""
        return D * self.div(N) - self.apply(N, D)

    def parity(self, P):
        pars = {yparity(e) for e in P.to_dict()}
        if len(pars) == 1:
            return pars.pop()
        return None if pars else 0


class Index:
    """Monomial index for one degree (and parity / predicate)."""
    def __init__(self, deg, par=None, pred=None):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par, pred)
        self.idx = {e: i for i, e in enumerate(self.mons)}
        self.n = len(self.mons)

    def vec(self, P):
        v = [0] * self.n
        for e, c in P.to_dict().items():
            v[self.idx[e]] = c
        return v


def build_mat(cols, rowidx, p, rhs=None):
    """nmod_mat whose columns are the given polys (vectorized on rowidx). rhs appended as last column."""
    R = rowidx.n
    C = len(cols) + (1 if rhs is not None else 0)
    flat = [0] * (R * C)
    for j, P in enumerate(cols):
        for e, c in P.to_dict().items():
            flat[rowidx.idx[e] * C + j] = c
    if rhs is not None:
        for e, c in rhs.to_dict().items():
            flat[rowidx.idx[e] * C + C - 1] = c
    return nmod_mat(R, C, flat, p)


def pivots_of_rref(Mr, rank):
    C = Mr.ncols()
    piv = []
    for i in range(rank):
        j = 0
        while j < C and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    return piv


# --------------------------------------------------------------------------------------
# instance
# --------------------------------------------------------------------------------------

def load(p):
    R = Ring(p)
    F, A, q, rec = core.load_instance(IID)
    Fp, Ap = R.from_dict(F), R.from_dict(A)
    Fi = [Fp.derivative(i) for i in range(NV)]
    return R, Fp, Ap, Fi, F, A


def num_parities(pD, pA=1):
    """parity of N_i given parity of D (pure parity D)."""
    return [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]


def koszul_parity(pD, i, j, pA=1):
    """parity of c_ij: c_ij F_j is the i-th component, so par(c_ij) = par(N_i) - par(F_j)."""
    pN = num_parities(pD, pA)
    pj = 0 if j < 3 else 1
    return (pN[i] - pj) % 2


# --------------------------------------------------------------------------------------
# E1 family (kernel of multiplication by A, with lifts)
# --------------------------------------------------------------------------------------

def not_div_X06(e):
    return e[0] < 6


def e1_family(R, Fp, Ap, Fi, d, pD, dpred=None, verbose=True):
    """All D of degree d and Y-parity pD (None = mixed) with D monomials restricted by dpred, such that
    A D in J(F); returns basis [(D, n)] with n(F) = A D (n = list of 6 polys), plus rank data.

    Columns: D monomials, then n_i monomials (degree d+1, parity par(N_i)); rows: degree d+6, parity 1+pD.
    Equation: A D - sum_i n_i F_i = 0."""
    t0 = time.time()
    p = R.p
    prow = None if pD is None else (1 + pD) % 2
    rows = Index(d + 6, prow)
    Dmons = monomials(d, pD, dpred)
    cols = [Ap * R.mono(e) for e in Dmons]
    ncols_n = []
    for i in range(NV):
        pi = None if pD is None else num_parities(pD)[i]
        ms = monomials(d + 1, pi)
        ncols_n.append(ms)
        for e in ms:
            cols.append(-(Fi[i] * R.mono(e)))
    M = build_mat(cols, rows, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    nD = len(Dmons)
    C = M.ncols()
    # rows of X^T restricted; rref with D-columns first (they are first already)
    if nullity == 0:
        return dict(basis=[], dimK=0, nD=nD, rows=rows.n, cols=C, nullity=0, t=round(time.time() - t0, 1),
                    Dmons=Dmons)
    XT = nmod_mat(nullity, C, [int(X[j, k]) for k in range(nullity) for j in range(C)], p)
    Xr, rk = XT.rref()
    piv = pivots_of_rref(Xr, rk)
    basis = []
    for r, j in enumerate(piv):
        if j >= nD:
            break
        Dd = {Dmons[k]: int(Xr[r, k]) for k in range(nD) if int(Xr[r, k])}
        D = R.ctx.from_dict(Dd)
        n = []
        off = nD
        for i in range(NV):
            ms = ncols_n[i]
            nd = {ms[k]: int(Xr[r, off + k]) for k in range(len(ms)) if int(Xr[r, off + k])}
            n.append(R.ctx.from_dict(nd))
            off += len(ms)
        basis.append((D, n))
    # J_d (parity pD, restricted to dpred-span): rank of Jacobian generators projected... we report dim of J_d cap span
    info = dict(basis=basis, dimK=len(basis), nD=nD, rows=rows.n, cols=C, nullity=nullity,
                nJcols=C - nD, t_build=round(t1 - t0, 1), t_null=round(t2 - t1, 1), t=round(time.time() - t0, 1),
                Dmons=Dmons)
    if verbose:
        print(f'[e1_family d={d} pD={pD}] rows={rows.n} cols={C} (D {nD}) nullity={nullity} dimK={len(basis)} '
              f'({info["t"]}s)', file=sys.stderr, flush=True)
    return info


def jacobian_dim(R, Fi, d, pD, dpred=None):
    """dim of J(F)_d (parity pD) intersected with the span of monomials satisfying dpred:
    computed as dim J_d - dim(J_d projected...) -- simpler: rank of generators restricted to those columns?
    We compute dim (J_d cap V) = dim J_d + dim V - dim(J_d + V) via ranks."""
    p = R.p
    allidx = Index(d, pD)
    gens = []
    for i in range(NV):
        pi = None if pD is None else (pD - (0 if i < 3 else 1)) % 2
        for e in monomials(d - 5, pi):
            gens.append(Fi[i] * R.mono(e))
    if not gens:
        return 0, 0
    J = build_mat(gens, allidx, p)
    rJ = J.rank()
    if dpred is None:
        return rJ, rJ
    Vm = [R.mono(e) for e in monomials(d, pD, dpred)]
    JV = build_mat(gens + Vm, allidx, p)
    rJV = JV.rank()
    return rJ, rJ + len(Vm) - rJV


# --------------------------------------------------------------------------------------
# normalized E2 system for a fixed (D, n)
# --------------------------------------------------------------------------------------

class E2Builder:
    """Precomputes the F*C columns and Koszul monomial lists for degree d / parity pD."""
    def __init__(self, R, Fp, Fi, d, pD):
        self.R, self.Fp, self.Fi, self.d, self.pD = R, Fp, Fi, d, pD
        prow = None if pD is None else 1
        self.rows = Index(2 * d, prow)
        self.Cmons = monomials(2 * d - 6, prow)
        self.FC = [Fp * R.mono(e) for e in self.Cmons]
        self.kos = []   # (i, j, monomial)
        for i in range(NV):
            for j in range(i + 1, NV):
                pc = None if pD is None else koszul_parity(pD, i, j)
                for e in monomials(d - 4, pc):
                    self.kos.append((i, j, e))

    def columns(self, D):
        R, Fi = self.R, self.Fi
        Di = [D.derivative(i) for i in range(NV)]
        W = {}
        for i in range(NV):
            for j in range(i + 1, NV):
                W[(i, j)] = Di[i] * Fi[j] - Di[j] * Fi[i]
        cols = []
        for (i, j, e) in self.kos:
            mu = R.mono(e)
            col = D * (mu.derivative(i) * Fi[j] - mu.derivative(j) * Fi[i]) - mu * W[(i, j)]
            cols.append(col)
        return cols

    def solve(self, D, n, want_leftnull=False):
        """Decide -B_D(n) in F S + im T_D. Returns dict(consistent, rank, rankaug, rows, cols)."""
        R = self.R
        t0 = time.time()
        cols = self.columns(D) + self.FC
        b = -R.B(D, n)
        M = build_mat(cols, self.rows, R.p, rhs=b)
        t1 = time.time()
        Mr, rank = M.rref()
        piv = pivots_of_rref(Mr, rank)
        C = len(cols)
        consistent = not (piv and piv[-1] == C)
        rankM = rank - (0 if consistent else 1)
        out = dict(consistent=consistent, rank_M=rankM, rank_aug=rank, rows=self.rows.n, cols=C,
                   n_koszul=len(self.kos), n_C=len(self.Cmons), t_build=round(t1 - t0, 1),
                   t_rref=round(time.time() - t1, 1))
        if consistent:
            x = [0] * C
            for i, j in enumerate(piv):
                x[j] = int(Mr[i, C])
            out['x'] = x
        return out

    def unpack(self, x):
        R = self.R
        c = {}
        K = [R.zero() for _ in range(NV)]
        for (i, j, e), v in zip(self.kos, x[:len(self.kos)]):
            if v:
                mu = R.mono(e, v)
                K[i] += mu * self.Fi[j]
                K[j] -= mu * self.Fi[i]
        Cp = R.zero()
        for e, v in zip(self.Cmons, x[len(self.kos):]):
            if v:
                Cp += R.mono(e, v)
        return K, Cp


# --------------------------------------------------------------------------------------
# Jacobian reducer in one degree (normal forms as coordinate vectors in R_deg)
# --------------------------------------------------------------------------------------

class Reducer:
    def __init__(self, R, Fi, deg, par, verbose=True):
        t0 = time.time()
        p = R.p
        self.R, self.deg, self.par = R, deg, par
        self.idx = Index(deg, par)
        gens = []
        for i in range(NV):
            pi = None if par is None else (par - (0 if i < 3 else 1)) % 2
            for e in monomials(deg - 5, pi):
                gens.append(Fi[i] * R.mono(e))
        # matrix with generators as ROWS (so rref gives reduced row basis of J_deg)
        n = self.idx.n
        flat = [0] * (len(gens) * n)
        for r, P in enumerate(gens):
            for e, c in P.to_dict().items():
                flat[r * n + self.idx.idx[e]] = c
        G = nmod_mat(len(gens), n, flat, p)
        Gr, rank = G.rref()
        self.rank = rank
        self.piv = pivots_of_rref(Gr, rank)
        pivset = set(self.piv)
        self.std = [j for j in range(n) if j not in pivset]
        self.ring_dim = n - rank
        # Gr restricted to pivot rows and standard columns: rank x ring_dim
        self.Gstd = nmod_mat(rank, len(self.std), [int(Gr[i, j]) for i in range(rank) for j in self.std], p)
        self.t = round(time.time() - t0, 1)
        if verbose:
            print(f'[Reducer deg={deg} par={par}] mons={n} gens={len(gens)} rank={rank} ring_dim={self.ring_dim} ({self.t}s)',
                  file=sys.stderr, flush=True)

    def nf_matrix(self, polys):
        """Normal forms of a list of polys as an (len(polys) x ring_dim) nmod_mat."""
        p = self.R.p
        k = len(polys)
        vecs = [self.idx.vec(P) for P in polys]
        Vpiv = nmod_mat(k, self.rank, [v[j] for v in vecs for j in self.piv], p)
        Vstd = nmod_mat(k, len(self.std), [v[j] for v in vecs for j in self.std], p)
        return Vstd - Vpiv * self.Gstd


def relaxation_quadrics(R, Fi, fam, d, pD, red=None, verbose=True):
    """Quadratic forms Phi_r(a) = coordinate r of [B_{D(a)}(n(a))] in R_{2d}; returns (list of dict
    {(k,l): coeff} for k<=l, ring_dim, reducer)."""
    basis = fam['basis']
    m = len(basis)
    if red is None:
        red = Reducer(R, Fi, 2 * d, None if pD is None else 1, verbose=verbose)
    t0 = time.time()
    polys = []
    pairs = []
    divn = [R.div(n) for (_, n) in basis]
    for k in range(m):
        Dk = basis[k][0]
        for l in range(m):
            nl = basis[l][1]
            polys.append(Dk * divn[l] - R.apply(nl, Dk))
            pairs.append((k, l))
    NF = red.nf_matrix(polys)
    p = R.p
    quad = [dict() for _ in range(red.ring_dim)]
    for r, (k, l) in enumerate(pairs):
        for c in range(red.ring_dim):
            v = int(NF[r, c])
            if v:
                key = (min(k, l), max(k, l))
                quad[c][key] = (quad[c].get(key, 0) + v) % p
    quad = [{kk: vv for kk, vv in q.items() if vv} for q in quad]
    if verbose:
        print(f'[relaxation d={d} pD={pD}] m={m} pairs={len(pairs)} ring_dim={red.ring_dim} '
              f'nonzero quadrics={sum(1 for q in quad if q)} ({round(time.time()-t0,1)}s)', file=sys.stderr, flush=True)
    return quad, red


def write_msolve(path, quad, m, p, varprefix='a'):
    names = [f'{varprefix}{k}' for k in range(m)]
    polys = []
    for q in quad:
        if not q:
            continue
        terms = []
        for (k, l), v in sorted(q.items()):
            if k == l:
                terms.append(f'{v}*{names[k]}^2')
            else:
                terms.append(f'{v}*{names[k]}*{names[l]}')
        polys.append('+'.join(terms))
    with open(path, 'w') as fh:
        fh.write(','.join(names) + '\n')
        fh.write(f'{p}\n')
        fh.write(',\n'.join(polys) + '\n')
    return len(polys)


def write_singular(path, quad, m, p, varprefix='a', extra=''):
    names = [f'{varprefix}{k}' for k in range(m)]
    polys = []
    for q in quad:
        if not q:
            continue
        terms = []
        for (k, l), v in sorted(q.items()):
            if k == l:
                terms.append(f'{v}*{names[k]}^2')
            else:
                terms.append(f'{v}*{names[k]}*{names[l]}')
        polys.append('+'.join(terms))
    with open(path, 'w') as fh:
        fh.write(f'ring r = {p}, ({",".join(names)}), dp;\n')
        fh.write('ideal I;\n')
        for i, s in enumerate(polys):
            fh.write(f'I[{i+1}] = {s};\n')
        fh.write('option(redSB);\nint t0 = timer;\nideal G = std(I);\n')
        fh.write('print("dimension of S/I (Krull) = " + string(dim(G)));\n')
        fh.write('print("size of reduced GB = " + string(size(G)));\n')
        fh.write('print("vdim (length if 0-dim) = " + string(vdim(G)));\n')
        fh.write('print("time = " + string(timer - t0));\n')
        fh.write(extra)
        fh.write('quit;\n')
    return len(polys)


def e1_lift(R, Fp, Ap, Fi, D):
    """One lift n with n(F) = A D (any parity), or None if A D not in J."""
    p = R.p
    d = D.total_degree()
    target = Ap * D
    rows = Index(d + 6, R.parity(target))
    cols = []
    tags = []
    pt = R.parity(target)
    for i in range(NV):
        pi = None if pt is None else (pt - (0 if i < 3 else 1)) % 2
        for e in monomials(d + 1, pi):
            cols.append(Fi[i] * R.mono(e)); tags.append((i, e))
    M = build_mat(cols, rows, p, rhs=target)
    Mr, rank = M.rref()
    piv = pivots_of_rref(Mr, rank)
    C = len(cols)
    if piv and piv[-1] == C:
        return None
    n = [{} for _ in range(NV)]
    for r, j in enumerate(piv):
        v = int(Mr[r, C])
        if v:
            i, e = tags[j]
            n[i][e] = v
    return [R.ctx.from_dict(x) for x in n]


def load_any(p, F, A):
    R = Ring(p)
    Fp, Ap = R.from_dict(F), R.from_dict(A)
    Fi = [Fp.derivative(i) for i in range(NV)]
    return R, Fp, Ap, Fi
