"""Wave-2 EXCLUSION/RIGOR library for symmetric-sextic-residual (q = 1).

Hessian relaxation (necessary condition, derived in wave2-symmetric-sextic-residual-excl.md):
  Let D be E1-admissible (A D in J(F)) and N0 any polynomial vector with sum_i N0_i F_i = A D.
  If (D, N, Q) is a certificate (q = 1) then, because F is irreducible and the syzygies of grad F are
  Koszul, N = N0 - (H/6) z + K_B and the cleared identity reduces to
        D div N0 - N0 . grad D + (K_B terms) = F Q,
  and every K_B term lies in J(F).  Hence   R0(D) := D div N0 - N0 . grad D  in J(F)_{2d}.
  For a linear family D(a) = sum a_k D_k with N0(a) = sum a_k u_k this is a system of quadrics in a.
  The solution set is a cone, invariant under D -> D + F g (illegal radial directions), so we restrict to a
  complement of F * S_{d-6} and ask whether the only solution is a = 0.

All linear algebra mod p (python-flint).  Rigor: every modular rank used in an exclusion is compared with
its characteristic-0 (complete-intersection / Molien-Koszul) value; equality makes the relevant cokernels
p-torsion-free, which is what the cone argument needs (see the report, section "Rigor conventions").
"""
import sys, time, json, itertools, os, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from fractions import Fraction
import numpy as np
from flint import nmod_mat, fmpq_mat, fmpq
from hodge import core, jacobian

IID = 'symmetric-sextic-residual'
RES = f'/tmp/hodge/results/{IID}'
PK = f'{RES}/wave2-excl'
LOG = f'{RES}/wave2-excl.log.jsonl'
PRIMES = core.PRIMES
P1, P2, P3 = PRIMES[0], PRIMES[1], PRIMES[2]

def log(rec, path=LOG):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('lens', 'excl')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

# ---------------------------------------------------------------- polynomial helpers
def var(i):
    e = [0]*6; e[i] = 1
    return {tuple(e): Fraction(1)}
X = [var(0), var(1), var(2)]; Y = [var(3), var(4), var(5)]
def e_k(vs, k):
    out = {}
    for c in itertools.combinations(vs, k):
        t = {(0,)*6: Fraction(1)}
        for v in c: t = core.mul(t, v)
        out = core.add(out, t)
    return out
def p_k(vs, k):
    out = {}
    for v in vs: out = core.add(out, core.power(v, k))
    return out
def dlt(vs):
    return core.mul(core.mul(core.sub(vs[0], vs[1]), core.sub(vs[0], vs[2])), core.sub(vs[1], vs[2]))
def act(Pol, perm):
    out = {}
    for e, c in Pol.items():
        e2 = [0]*6
        for i in range(6): e2[perm[i]] += e[i]
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}
def divergence(N, p=None):
    out = {}
    for i in range(6): out = core.add(out, core.pderiv(N[i], i, p), p)
    return out
def dot_grad(N, D, p=None):
    out = {}
    for i in range(6): out = core.add(out, core.mul(N[i], core.pderiv(D, i, p), p), p)
    return out
def R0(D, N, p=None):
    """D div N - N . grad D"""
    return core.sub(core.mul(D, divergence(N, p), p), dot_grad(N, D, p), p)

# ---------------------------------------------------------------- Hilbert functions (characteristic 0)
def nS(n, par=None):
    if n < 0: return 0
    if par is None: return math.comb(n + 5, 5)
    return sum(math.comb(n - j + 2, 2) * math.comb(j + 2, 2) for j in range(par, n + 1, 2))

def dimR(n, par=None):
    """dim of the Jacobian ring R = S/J(F) in degree n (and Y-parity) for smooth F (regular sequence)."""
    tot = 0
    for a in range(4):
        for b in range(4):
            m = n - 5 * (a + b)
            if m < 0: continue
            pr = None if par is None else (par - b) % 2
            tot += (-1) ** (a + b) * math.comb(3, a) * math.comb(3, b) * nS(m, pr)
    return tot

def dimJ(n, par=None):
    return nS(n, par) - dimR(n, par)

# Molien-Koszul equivariant Hilbert function for G x <tau>, G = S3(X) x S3(Y)
S3 = list(itertools.permutations(range(3)))
def group_elements():
    """list of (perm tuple (i -> perm[i]), signs tuple, sgnX, sgnY, tau flag)"""
    out = []
    for pX in S3:
        for pY in S3:
            for t in (0, 1):
                perm = tuple(pX) + tuple(3 + j for j in pY)
                signs = tuple(1 if i < 3 or t == 0 else -1 for i in range(6))
                out.append((perm, signs, perm_sign(pX), perm_sign(pY), t))
    return out
def perm_sign(pp):
    inv = sum(1 for i in range(len(pp)) for j in range(i + 1, len(pp)) if pp[i] > pp[j])
    return -1 if inv % 2 else 1
def cycles(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i in seen: continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j); cyc.append(j); j = perm[j]
        out.append(cyc)
    return out
def _series_ratio(perm, signs, nmax):
    """coefficients of det(1 - t^5 g)/det(1 - t g) up to t^nmax (g = signed permutation)."""
    num = [0] * (nmax + 1); num[0] = 1
    den = [0] * (nmax + 1); den[0] = 1
    def mulpoly(a, k, eps):   # multiply truncated series a by (1 - eps t^k)
        b = a[:]
        for i in range(k, nmax + 1): b[i] -= eps * a[i - k]
        return b
    for cyc in cycles(perm):
        l = len(cyc); eps = 1
        for i in cyc: eps *= signs[i]
        num = mulpoly(num, 5 * l, eps)
        den = mulpoly(den, l, eps)
    # series num/den
    out = [0] * (nmax + 1)
    for n in range(nmax + 1):
        s = num[n]
        for k in range(1, n + 1): s -= den[k] * out[n - k]
        out[n] = s   # den[0] = 1
    return out
def dimR_char(n, chi, par, nmax=None):
    """dim R_n^{(chi, parity)}, chi in {'triv','sgnX','sgnY','sgnsgn'} of G = S3 x S3, Y-parity par."""
    tot = Fraction(0)
    for perm, signs, sX, sY, t in group_elements():
        c = {'triv': 1, 'sgnX': sX, 'sgnY': sY, 'sgnsgn': sX * sY}[chi]
        eps = (-1) ** (t * par)
        tot += c * eps * _series_ratio(perm, signs, n)[n]
    tot /= 72
    assert tot.denominator == 1
    return int(tot)
def dimS_char(n, chi, par):
    """dim S_n^{(chi, parity)} (polynomial ring) by Molien."""
    tot = Fraction(0)
    for perm, signs, sX, sY, t in group_elements():
        c = {'triv': 1, 'sgnX': sX, 'sgnY': sY, 'sgnsgn': sX * sY}[chi]
        eps = (-1) ** (t * par)
        den = [0] * (n + 1); den[0] = 1
        for cyc in cycles(perm):
            l = len(cyc); e = 1
            for i in cyc: e *= signs[i]
            b = den[:]
            for i in range(l, n + 1): b[i] -= e * den[i - l]
            den = b
        out = [0] * (n + 1)
        for m in range(n + 1):
            s = 1 if m == 0 else 0
            for k in range(1, m + 1): s -= den[k] * out[m - k]
            out[m] = s
        tot += c * eps * out[n]
    tot /= 72
    assert tot.denominator == 1
    return int(tot)

# ---------------------------------------------------------------- dense matrix helpers
def to_arr(Pol, p):
    if not Pol:
        return np.zeros((0, 6), dtype=np.int64), np.zeros(0, dtype=np.int64)
    items = list(core.reduce_mod(Pol, p).items())
    return np.array([e for e, _ in items], dtype=np.int64), np.array([c for _, c in items], dtype=np.int64)
BASE = 64
def codes_of(E):
    code = np.zeros(E.shape[0], dtype=np.int64)
    for i in range(6): code = code * BASE + E[:, i]
    return code
class MonIndex:
    def __init__(self, mons):
        self.mons = mons
        E = np.array(mons, dtype=np.int64)
        self.codes = codes_of(E)
        self.order = np.argsort(self.codes)
        self.sorted = self.codes[self.order]
        self.n = len(mons)
    def vec(self, Pol, p):
        E, C = to_arr(Pol, p)
        v = np.zeros(self.n, dtype=np.int64)
        if E.shape[0]:
            pos = self.order[np.searchsorted(self.sorted, codes_of(E))]
            np.add.at(v, pos, C)
        return v % p
def mat_from_polys(polys, index, p):
    """nmod_mat with columns = coefficient vectors of polys."""
    R, C = index.n, len(polys)
    M = np.zeros((R, C), dtype=np.int64)
    for j, P in enumerate(polys):
        M[:, j] = index.vec(P, p)
    return nmod_mat(R, C, M.ravel().tolist(), p)
def np_from_nmod(M):
    return np.array(M.tolist(), dtype=object)

# ---------------------------------------------------------------- reducer modulo J in one degree (dense, vectorised)
class Reducer:
    """Normal forms modulo J(F)_deg (parity par) mod p, vectorised: NF(P) = P - P[piv] * Rref."""
    def __init__(self, F, deg, par, p, verbose=False):
        t0 = time.time()
        self.p, self.deg, self.par = p, deg, par
        Fp = core.reduce_mod(F, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        mons = core.monomials(deg, par)
        self.index = MonIndex(mons)
        gens = []
        for i in range(6):
            pm = par if i < 3 else (par + 1) % 2
            for m in core.monomials(deg - 5, pm):
                gens.append(core.shift(Fi[i], m, 1, p))
        self.ngens = len(gens)
        if gens:
            M = mat_from_polys(gens, self.index, p).transpose()    # rows = generators
            Mr, rank = M.rref()
            self.rank = rank
            rows = np.array([[int(Mr[i, j]) for j in range(self.index.n)] for i in range(rank)], dtype=np.int64) if rank else np.zeros((0, self.index.n), dtype=np.int64)
            piv = []
            for i in range(rank):
                nz = np.nonzero(rows[i])[0]
                piv.append(int(nz[0]))
            self.pivots = np.array(piv, dtype=np.int64)
            self.rref = rows
        else:
            self.rank = 0; self.pivots = np.zeros(0, dtype=np.int64); self.rref = np.zeros((0, self.index.n), dtype=np.int64)
        self.standard = np.array([j for j in range(self.index.n) if j not in set(self.pivots.tolist())], dtype=np.int64)
        self.ci_rank = dimJ(deg, par)
        self.ring_dim = self.index.n - self.rank
        self.t_build = time.time() - t0
        if verbose:
            print(f'[Reducer] deg={deg} par={par} p={p} mons={self.index.n} gens={self.ngens} rank={self.rank} (CI {self.ci_rank}) ring={self.ring_dim} {self.t_build:.1f}s', file=sys.stderr, flush=True)
    def coords(self, polys):
        """matrix (len(polys) x ring_dim): coordinates of the normal forms in the standard monomials."""
        p = self.p
        Pm = np.zeros((len(polys), self.index.n), dtype=np.int64)
        for i, P in enumerate(polys): Pm[i] = self.index.vec(P, p)
        if self.rank:
            A = nmod_mat(Pm.shape[0], self.rank, Pm[:, self.pivots].ravel().tolist(), p)
            B = nmod_mat(self.rank, self.index.n, self.rref.ravel().tolist(), p)
            prod = np.array((A * B).tolist(), dtype=np.int64)
            NF = (Pm - prod) % p
        else:
            NF = Pm
        return NF[:, self.standard]

# ---------------------------------------------------------------- E1-admissible family (full space, one parity)
class E1Family:
    """All D of degree d and Y-parity pD with A D in J(F) (mod p), with particular solutions u:
       D_k = sum over basis, u_k polynomial vectors with sum_i u_k[i] F_i = A D_k.
       Also identifies the J-part and the radial (F-multiple) subspace."""
    def __init__(self, F, A, d, pD, p, verbose=False):
        t0 = time.time()
        self.F, self.A, self.d, self.pD, self.p = F, A, d, pD, p
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        tp = (1 + pD) % 2
        rowmons = core.monomials(d + 6, tp)
        index = MonIndex(rowmons)
        Dmons = core.monomials(d, pD)
        Acols = [core.shift(Ap, m, 1, p) for m in Dmons]
        Jinfo, Jcols = [], []
        for i in range(6):
            pm = tp if i < 3 else (tp + 1) % 2
            for m in core.monomials(d + 1, pm):
                Jinfo.append((i, m)); Jcols.append(core.shift(Fi[i], m, 1, p))
        nA, nJ = len(Acols), len(Jcols)
        M = mat_from_polys(Acols + Jcols, index, p)
        Xn, nullity = M.nullspace()
        self.rows, self.nA, self.nJ, self.nullity = index.n, nA, nJ, nullity
        self.rankM = nA + nJ - nullity
        self.rankJ = nmod_mat(index.n, nJ, [int(x) for x in np.array(M.tolist(), dtype=object)[:, nA:].ravel()], p).rank() if False else None
        # rref of the null vectors (rows), D-coords first
        if nullity:
            Nm = nmod_mat(nullity, nA + nJ, [int(Xn[j, k]) for k in range(nullity) for j in range(nA + nJ)], p)
            Nr, rk = Nm.rref()
            rows = np.array([[int(Nr[i, j]) for j in range(nA + nJ)] for i in range(rk)], dtype=np.int64)
            keep = [i for i in range(rk) if np.nonzero(rows[i])[0][0] < nA]
            K = rows[keep]
        else:
            K = np.zeros((0, nA + nJ), dtype=np.int64)
        self.dimK = K.shape[0]
        self.Dmons, self.Jinfo = Dmons, Jinfo
        self.K = K        # rows: [D-coords | J-coords]; A D = -sum_J coords * (m F_i)
        # J_d part
        if d >= 5:
            Jd = []
            for i in range(6):
                pm = pD if i < 3 else (pD + 1) % 2
                for m in core.monomials(d - 5, pm): Jd.append(core.shift(Fi[i], m, 1, p))
            self.dimJd = mat_from_polys(Jd, MonIndex(Dmons), p).rank() if Jd else 0
        else:
            self.dimJd = 0
        self.ker_R = self.dimK - self.dimJd
        # rank bookkeeping for the J_{d+6} generators (CI value)
        self.ci_rankJ_top = dimJ(d + 6, tp)
        self.ci_dimJd = dimJ(d, pD) if d >= 5 else 0
        # radial subspace: F * S_{d-6}^{pD}
        self.radial = []
        if d >= 6:
            for m in core.monomials(d - 6, pD):
                self.radial.append(core.shift(Fp, m, 1, p))
        self.t_build = time.time() - t0
        if verbose:
            print(f'[E1Family] d={d} pD={pD} p={p} rows={index.n} nA={nA} nJ={nJ} nullity={nullity} dimK={self.dimK} dimJd={self.dimJd} ker_R={self.ker_R} (CI dimJd {self.ci_dimJd}) {self.t_build:.1f}s', file=sys.stderr, flush=True)
    def rank_check(self):
        """rank of the J_{d+6} generator block mod p vs CI value (needed for rigor)."""
        p = self.p
        Fp = core.reduce_mod(self.F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        tp = (1 + self.pD) % 2
        index = MonIndex(core.monomials(self.d + 6, tp))
        Jcols = [core.shift(Fi[i], m, 1, p) for (i, m) in self.Jinfo]
        r = mat_from_polys(Jcols, index, p).rank()
        return r, self.ci_rankJ_top
    def pair(self, coords):
        """(D, u) for a coordinate row (length nA + nJ)."""
        p = self.p
        D = {}
        for j, m in enumerate(self.Dmons):
            v = int(coords[j]) % p
            if v: D[m] = v
        u = [{} for _ in range(6)]
        for j, (i, m) in enumerate(self.Jinfo):
            v = int(coords[self.nA + j]) % p
            if v: u[i][m] = (p - v) % p
        return D, u
    def complement_basis(self):
        """basis of K modulo the radial subspace: list of (D, u) rows; also returns radial coords."""
        p = self.p
        Dindex = MonIndex(self.Dmons)
        nr = len(self.radial)
        if nr == 0:
            return [self.pair(self.K[i]) for i in range(self.dimK)], 0
        Rm = np.zeros((nr, self.nA), dtype=np.int64)
        for i, P in enumerate(self.radial): Rm[i] = Dindex.vec(P, p)
        Rr = nmod_mat(nr, self.nA, Rm.ravel().tolist(), p)
        Rrr, rk = Rr.rref()
        assert rk == nr
        rows = np.array([[int(Rrr[i, j]) for j in range(self.nA)] for i in range(rk)], dtype=np.int64)
        piv = [int(np.nonzero(rows[i])[0][0]) for i in range(rk)]
        # u for radial basis rows (F m): A F m = sum_i (A m z_i / 6) F_i  -> u_i = A m z_i / 6 ; we need u for the rref'd rows: solve via K instead:
        # simpler: reduce each K row modulo the radial rows using the D-coords only, then recompute u by combining the K rows' u (linear).
        # Express radial rref rows in terms of K rows: radial rows are in K (A F m in J), so solve R = C * K_D  (C = coefficients).
        KD = nmod_mat(self.dimK, self.nA, self.K[:, :self.nA].ravel().tolist(), p)
        # solve C * KD = rows  <=>  KD^T C^T = rows^T ; use rref of [KD^T | rows^T]
        aug = nmod_mat(self.nA, self.dimK + rk, np.concatenate([self.K[:, :self.nA].T, rows.T], axis=1).ravel().tolist(), p)
        augr, rka = aug.rref()
        assert rka == self.dimK, 'K rows dependent?'
        Cm = np.zeros((rk, self.dimK), dtype=np.int64)
        for i in range(self.dimK):
            # pivot of row i of augr must be column i (KD rows independent)
            for j in range(rk):
                Cm[j, i] = int(augr[i, self.dimK + j])
        # full coordinates (D|J) of the radial rref rows:
        Rfull = (Cm @ self.K) % p
        # reduce K rows modulo Rfull on the pivot columns
        Kred = self.K.copy()
        for i in range(self.dimK):
            for j, pc in enumerate(piv):
                c = int(Kred[i, pc]) % p
                if c:
                    Kred[i] = (Kred[i] - c * Rfull[j]) % p
        # independent subset
        Km = nmod_mat(self.dimK, self.nA + self.nJ, Kred.ravel().tolist(), p)
        Kr, rk2 = Km.rref()
        assert rk2 == self.dimK - nr, (rk2, self.dimK, nr)
        out = []
        for i in range(rk2):
            row = np.array([int(Kr[i, j]) for j in range(self.nA + self.nJ)], dtype=np.int64)
            out.append(self.pair(row))
        return out, nr

# ---------------------------------------------------------------- sector (character) families via charsolve
import charsolve as CS
GROUPS = {'triv': ([1, 1]), 'sgnX': ([-1, 1]), 'sgnY': ([1, -1]), 'sgnsgn': ([-1, -1])}
def group(chi):
    return CS.Group([[0, 1, 2], [3, 4, 5]], GROUPS[chi])

class SectorFamily:
    """chi-semi-invariant D of degree d, parity pD, with A D in J (mod p), via the twisted E1 system."""
    def __init__(self, F, A, d, pD, chi, p, verbose=False):
        t0 = time.time()
        self.F, self.A, self.d, self.pD, self.chi, self.p = F, A, d, pD, chi, p
        G = group(chi); self.G = G
        q = 1
        pA = core.parity_of(A); degN = d + 6 * q - 5
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        tfold, nrows = CS._tfold_factory(G, d + 6 * q, (pA + pD) % 2, p)
        cols = []; colinfo = []
        for s in G.slots:
            FE, FC = CS.to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
            for m in CS.orbit_reps(degN, CS.pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                cols.append(tfold(FE + np.array(m, dtype=np.int64), FC)); colinfo.append(('N', s, m))
        nN = len(cols)
        Dl = CS.live_reps(d, pD, G)
        self.Dreps = [r for r, _ in Dl]; self.Dbasis = [P for _, P in Dl]
        for P in self.Dbasis:
            E_, C_ = CS.to_arr(core.mul(Ap, core.reduce_mod(P, p), p), p)
            cols.append(tfold(E_, C_))
        nD = len(self.Dbasis)
        self.nN, self.nD, self.rows = nN, nD, nrows
        self.colinfo = colinfo
        if nD == 0:
            self.dimK = 0; self.K = np.zeros((0, nN + nD), dtype=np.int64); self.dimJd = 0; self.ker_R = 0; self.rankM = 0
            self.t_build = time.time() - t0
            return
        M = np.stack(cols, axis=1)
        Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
        Xn, nullity = Mn.nullspace()
        self.rankM = nN + nD - nullity
        if nullity:
            # rows = null vectors with D-coords first
            Nm = nmod_mat(nullity, nN + nD, [int(Xn[(nN + j) if j < nD else (j - nD), k]) for k in range(nullity) for j in range(nN + nD)], p)
            Nr, rk = Nm.rref()
            rows = np.array([[int(Nr[i, j]) for j in range(nN + nD)] for i in range(rk)], dtype=np.int64)
            keep = [i for i in range(rk) if np.nonzero(rows[i])[0][0] < nD]
            K = rows[keep]
        else:
            K = np.zeros((0, nN + nD), dtype=np.int64)
        self.K = K; self.dimK = K.shape[0]
        # J part in the sector
        Jpolys = CS.J_char_basis(Fp, d, pD, G, p)
        Jvecs = [[Pj.get(r, 0) % p for r in self.Dreps] for Pj in Jpolys]
        self.dimJd = nmod_mat(len(Jvecs), nD, [v for row in Jvecs for v in row], p).rank() if Jvecs else 0
        self.ker_R = self.dimK - self.dimJd
        # radial: F * (chi-semi-invariant g of degree d-6, parity pD)  == F * signed orbit sums
        self.radial = []
        if d >= 6:
            for r, P in CS.live_reps(d - 6, pD, G):
                self.radial.append(core.reduce_mod(core.mul(Fp, P, p), p))
        self.t_build = time.time() - t0
        if verbose:
            print(f'[SectorFamily] chi={chi} d={d} pD={pD} p={p} rows={nrows} nN={nN} nD={nD} dimK={self.dimK} dimJd={self.dimJd} ker={self.ker_R} radial={len(self.radial)} {self.t_build:.1f}s', file=sys.stderr, flush=True)
    def pair(self, coords):
        """(D, u) for a K row: D = sum coords[j] * Dbasis[j]; u = unfold of the N part with sign: sum N_i F_i = A D."""
        p = self.p; G = self.G
        D = {}
        for j in range(self.nD):
            v = int(coords[j]) % p
            if v: D = core.add(D, core.scale(self.Dbasis[j], v, p), p)
        # columns: sum_N x_N (q F_s m twisted-folded) + sum_D x_D (A P_D) = 0  =>  N(F) = -A D with N from x_N
        N = [{} for _ in range(6)]
        for j, (kind, s, m) in enumerate(self.colinfo):
            v = int(coords[self.nD + j]) % p
            if not v: continue
            mono = {m: (p - v) % p}      # minus sign: u.gradF = A D
            for perm, c in G.perms:
                N[perm[s]] = core.add(N[perm[s]], core.scale(act(mono, perm), c % p, p), p)
        # the fold sums over the whole group: each equivariant unknown contributes |orbit| copies; normalise:
        # The twisted fold of column (s,m) equals tfold(F_s m) and the unfolded N above is sum_g chi(g) g.(m e_s);
        # consistency check below decides the normalisation.
        inv_order = pow(G.order, -1, p)
        N = [core.scale(Ni, inv_order, p) for Ni in N]
        lhs = {}
        Fp = core.reduce_mod(self.F, p)
        for i in range(6): lhs = core.add(lhs, core.mul(N[i], core.pderiv(Fp, i, p), p), p)
        rhs = core.mul(core.reduce_mod(self.A, p), D, p)
        if lhs != rhs: raise RuntimeError('unfold normalisation failed')
        return D, N
    def complement_basis(self):
        p = self.p
        pairs = [self.pair(self.K[i]) for i in range(self.dimK)]
        nr = len(self.radial)
        if nr == 0: return pairs, 0
        # reduce modulo radial using D-coordinates on Dreps
        def dvec(D): return np.array([D.get(r, 0) % p for r in self.Dreps], dtype=np.int64)
        Rm = np.array([dvec(P) for P in self.radial], dtype=np.int64)
        Rr, rk = nmod_mat(nr, self.nD, Rm.ravel().tolist(), p).rref()
        assert rk == nr
        rows = np.array([[int(Rr[i, j]) for j in range(self.nD)] for i in range(rk)], dtype=np.int64)
        piv = [int(np.nonzero(rows[i])[0][0]) for i in range(rk)]
        # radial rows as combinations of K rows (D-coords): solve
        KD = np.array([dvec(D) for D, _ in pairs], dtype=np.int64)
        aug = nmod_mat(self.nD, self.dimK + rk, np.concatenate([KD.T, rows.T], axis=1).ravel().tolist(), p)
        augr, rka = aug.rref()
        assert rka == self.dimK
        Cm = np.zeros((rk, self.dimK), dtype=np.int64)
        for i in range(self.dimK):
            for j in range(rk): Cm[j, i] = int(augr[i, self.dimK + j])
        # pairs for radial rows: combine
        rad_pairs = []
        for j in range(rk):
            D = {}; u = [{} for _ in range(6)]
            for i in range(self.dimK):
                c = int(Cm[j, i]) % p
                if c:
                    D = core.add(D, core.scale(pairs[i][0], c, p), p)
                    for s in range(6): u[s] = core.add(u[s], core.scale(pairs[i][1][s], c, p), p)
            rad_pairs.append((D, u))
        # reduce each pair modulo radial pairs on pivot coordinates
        red = []
        for D, u in pairs:
            D2, u2 = dict(D), [dict(us) for us in u]
            for j, pc in enumerate(piv):
                c = D2.get(self.Dreps[pc], 0) % p
                if c:
                    D2 = core.sub(D2, core.scale(rad_pairs[j][0], c, p), p)
                    for s in range(6): u2[s] = core.sub(u2[s], core.scale(rad_pairs[j][1][s], c, p), p)
            red.append((D2, u2))
        # independent subset by D-coords
        M = np.array([dvec(D) for D, _ in red], dtype=np.int64)
        Mr, rk2 = nmod_mat(len(red), self.nD, M.ravel().tolist(), p).rref()
        assert rk2 == self.dimK - nr, (rk2, self.dimK, nr)
        # pick rows greedily
        out = []; cur = np.zeros((0, self.nD), dtype=np.int64)
        for (D, u), v in zip(red, M):
            test = np.concatenate([cur, v[None, :]], axis=0)
            r = nmod_mat(test.shape[0], self.nD, test.ravel().tolist(), p).rank()
            if r == cur.shape[0] + 1:
                cur = test; out.append((D, u))
        assert len(out) == rk2
        return out, nr

# ---------------------------------------------------------------- relaxation quadrics
def relaxation_quadrics(pairs, F, p, reducer=None, deg=None, verbose=False, folded=None):
    """pairs: list of (D_k, u_k) mod p (same degree d).  Returns (standard-monomial count, quadrics) where
    quadrics is a list of dicts {(k,l): coeff} (k<=l) -- one per coordinate of R_{2d} (or of the folded quotient)."""
    n = len(pairs)
    t0 = time.time()
    polys = []; keys = []
    # R(a) = sum_{k,l} a_k a_l (D_k div u_l - u_l . grad D_k)
    Rkl = {}
    for k in range(n):
        Dk = pairs[k][0]
        gradDk = [core.pderiv(Dk, i, p) for i in range(6)]
        for l in range(n):
            ul = pairs[l][1]
            divu = divergence(ul, p)
            t = core.mul(Dk, divu, p)
            for i in range(6): t = core.sub(t, core.mul(ul[i], gradDk[i], p), p)
            Rkl[(k, l)] = t
    for k in range(n):
        for l in range(k, n):
            if k == l: polys.append(Rkl[(k, k)])
            else: polys.append(core.add(Rkl[(k, l)], Rkl[(l, k)], p))
            keys.append((k, l))
    t1 = time.time()
    if folded is not None:
        C = folded.coords(polys)
    else:
        C = reducer.coords(polys)
    t2 = time.time()
    nq = C.shape[1]
    quadrics = []
    for s in range(nq):
        col = C[:, s]
        nz = np.nonzero(col)[0]
        if len(nz) == 0: continue
        quadrics.append({keys[i]: int(col[i]) for i in nz})
    if verbose:
        print(f'[relax] n={n} pairs={len(polys)} products {t1-t0:.1f}s, reduction {t2-t1:.1f}s, coords={nq} nonzero quadrics={len(quadrics)}', file=sys.stderr, flush=True)
    return nq, quadrics

# ---------------------------------------------------------------- folded (G-invariant) reducer for invariant R
class FoldedReducer:
    """Membership/coordinates modulo J(F)_deg^{G} for G-invariant polynomials of degree deg, parity par.
    Rows = G-orbits of monomials, columns = Sym_G(m F_s)."""
    def __init__(self, F, deg, par, p, verbose=False):
        t0 = time.time()
        self.p, self.deg, self.par = p, deg, par
        G = group('triv'); self.G = G
        Fp = core.reduce_mod(F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        rowmons = np.array(core.monomials(deg, par), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); self.nrows = len(self.orbit_codes)
        def fold(E, C):
            idx = np.searchsorted(self.orbit_codes, G.canon(E))
            v = np.bincount(idx, weights=C.astype(np.float64), minlength=self.nrows)
            return (np.rint(v).astype(np.int64)) % p
        self.fold = fold
        cols = []
        for s in G.slots:
            FE, FC = CS.to_arr(Fi[s], p); sb = G.stab_blocks(s)
            for m in CS.orbit_reps(deg - 5, (par + (1 if s >= 3 else 0)) % 2, lambda E: G.canon(E, sb)):
                cols.append(fold(FE + np.array(m, dtype=np.int64), FC))
        self.ngens = len(cols)
        M = np.stack(cols, axis=1).T   # rows = generators (folded)
        Mr, rank = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rref()
        self.rank = rank
        rows = np.array([[int(Mr[i, j]) for j in range(self.nrows)] for i in range(rank)], dtype=np.int64) if rank else np.zeros((0, self.nrows), dtype=np.int64)
        self.rref = rows
        self.pivots = np.array([int(np.nonzero(rows[i])[0][0]) for i in range(rank)], dtype=np.int64)
        self.standard = np.array([j for j in range(self.nrows) if j not in set(self.pivots.tolist())], dtype=np.int64)
        self.ring_dim = self.nrows - self.rank
        self.ci_ring_dim = dimR_char(deg, 'triv', par)
        self.ci_rank = self.nrows - self.ci_ring_dim
        self.t_build = time.time() - t0
        if verbose:
            print(f'[FoldedReducer] deg={deg} par={par} p={p} orbits={self.nrows} gens={self.ngens} rank={self.rank} (CI {self.ci_rank}) ring={self.ring_dim} (CI {self.ci_ring_dim}) {self.t_build:.1f}s', file=sys.stderr, flush=True)
    def coords(self, polys):
        p = self.p
        Pm = np.zeros((len(polys), self.nrows), dtype=np.int64)
        for i, P in enumerate(polys):
            E, C = CS.to_arr(P, p)
            if E.shape[0]: Pm[i] = self.fold(E, C)
        if self.rank:
            A = nmod_mat(Pm.shape[0], self.rank, Pm[:, self.pivots].ravel().tolist(), p)
            B = nmod_mat(self.rank, self.nrows, self.rref.ravel().tolist(), p)
            NF = (Pm - np.array((A * B).tolist(), dtype=np.int64)) % p
        else:
            NF = Pm
        return NF[:, self.standard]

# ---------------------------------------------------------------- msolve / Singular I/O
def write_msolve(path, n, quadrics, p, varprefix='a'):
    names = [f'{varprefix}{k+1}' for k in range(n)]
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n')
        f.write(f'{p}\n')
        lines = []
        for qd in quadrics:
            terms = []
            for (k, l), c in sorted(qd.items()):
                c = c % p
                if k == l: terms.append(f'{c}*{names[k]}^2')
                else: terms.append(f'{c}*{names[k]}*{names[l]}')
            lines.append('+'.join(terms))
        f.write(',\n'.join(lines) + '\n')
    return names
def write_singular(path, n, quadrics, p, varprefix='a', extra=''):
    names = [f'{varprefix}{k+1}' for k in range(n)]
    with open(path, 'w') as f:
        f.write(f'ring r = {p}, ({",".join(names)}), dp;\n')
        f.write('ideal I;\n')
        for qd in quadrics:
            terms = []
            for (k, l), c in sorted(qd.items()):
                c = c % p
                if k == l: terms.append(f'{c}*{names[k]}^2')
                else: terms.append(f'{c}*{names[k]}*{names[l]}')
            f.write('I[size(I)+1] = ' + '+'.join(terms) + ';\n')
        f.write('option(redSB);\nint t0 = timer;\nideal G = std(I);\n')
        f.write('"dim(std(I)) ="; dim(G); "vdim ="; vdim(G); "time(s) ="; (timer - t0) div 1000;\n')
        f.write(extra)
        f.write('quit;\n')
    return names
def msolve_dimension(outfile):
    """parse msolve output: returns ('empty'|'positive'|'zero-dim', raw head)."""
    txt = open(outfile).read().strip()
    head = txt[:200]
    if txt.startswith('[-1]'): return 'empty', head
    if txt.startswith('[1,'): return 'positive', head
    return 'zero-dim', head
def gb_pure_powers(gbfile, names):
    """from an msolve -g 1 leading-ideal dump, check which variables have a pure power as a leading term."""
    txt = open(gbfile).read()
    found = {}
    import re
    for tok in re.split(r'[,\n\[\]]', txt):
        tok = tok.strip()
        if not tok: continue
        m = re.fullmatch(r'(\d+\*)?([A-Za-z_]\w*)(\^(\d+))?', tok)
        if m:
            v = m.group(2); k = int(m.group(4)) if m.group(4) else 1
            found[v] = min(found.get(v, 99), k)
    return {v: found.get(v) for v in names}

# ---------------------------------------------------------------- exact (rational) twisted system for Farkas certificates
class ExactCharSystem:
    """Exact integer version of charsolve.CharSystem (F scaled to integers).  Rows = H-orbits, columns = equivariant unknowns."""
    def __init__(self, F, A, q, D, G, verbose=False):
        t0 = time.time()
        assert q == 1
        assert G.is_semi_invariant(D)
        # clear denominators of F and D
        lf = 1
        for c in F.values(): lf = math.lcm(lf, Fraction(c).denominator)
        ld = 1
        for c in D.values(): ld = math.lcm(ld, Fraction(c).denominator)
        Fz = {e: int(Fraction(c) * lf) for e, c in F.items()}
        Dz = {e: int(Fraction(c) * ld) for e, c in D.items()}
        Az = {e: int(Fraction(c)) for e, c in A.items()}
        self.lf, self.ld = lf, ld
        self.G = G; self.F, self.A, self.D = Fz, Az, Dz
        d = core.degree(Dz); self.d = d
        pD = core.parity_of(Dz); pA = core.parity_of(Az)
        degN = d + 1; degQ = 2 * d - 6; degE = 2 * d + 6
        Fi = [core.pderiv(Fz, i) for i in range(6)]
        Di = [core.pderiv(Dz, i) for i in range(6)]
        FD = core.mul(Fz, Dz)
        rowmons = np.array(core.monomials(degE, pA), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); self.nrows = len(self.orbit_codes)
        def fold(P):
            v = [0] * self.nrows
            if not P: return v
            E = np.array(list(P.keys()), dtype=np.int64)
            idx = np.searchsorted(self.orbit_codes, G.canon(E))
            for i, c in zip(idx.tolist(), P.values()):
                v[i] += int(c)
            return v
        cols = []; colinfo = []
        for s in G.slots:
            G_s = core.sub(core.scale(core.mul(Dz, Fi[s]), -1), core.mul(Fz, Di[s]))
            sb = G.stab_blocks(s)
            for m in CS.orbit_reps(degN, CS.pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                P = core.shift(G_s, m, 1)
                if m[s]:
                    m2 = list(m); m2[s] -= 1
                    P = core.add(P, core.shift(FD, tuple(m2), m[s]))
                cols.append(fold(P)); colinfo.append(('N', s, m))
        self.nN = len(cols)
        if degQ >= 0:
            F2 = core.scale(core.mul(Fz, Fz), -1)
            for m in CS.orbit_reps(degQ, pA, G.canon):
                cols.append(fold(core.shift(F2, m, 1))); colinfo.append(('Q', -1, m))
        self.nQ = len(cols) - self.nN
        self.rhs = fold(core.scale(core.mul(Az, core.mul(Dz, Dz)), -1))
        self.cols = cols; self.colinfo = colinfo; self.ncols = len(cols)
        self.t_build = time.time() - t0
        if verbose: print(f'[exact] d={d} rows={self.nrows} cols={self.ncols} built {self.t_build:.1f}s', file=sys.stderr, flush=True)
    def matrix(self):
        M = fmpq_mat(self.nrows, self.ncols + 1)
        for j, col in enumerate(self.cols):
            for i, v in enumerate(col):
                if v: M[i, j] = v
        for i, v in enumerate(self.rhs):
            if v: M[i, self.ncols] = v
        return M
    def ranks(self):
        M = self.matrix()
        r_aug = M.rank()
        M0 = fmpq_mat(self.nrows, self.ncols)
        for j, col in enumerate(self.cols):
            for i, v in enumerate(col):
                if v: M0[i, j] = v
        r = M0.rank()
        return r, r_aug
    def farkas(self):
        """integer vector y (length nrows) with y M = 0 and y b != 0, via exact rref of [M | b | I]."""
        n, m = self.nrows, self.ncols
        Aug = fmpq_mat(n, m + 1 + n)
        for j, col in enumerate(self.cols):
            for i, v in enumerate(col):
                if v: Aug[i, j] = v
        for i, v in enumerate(self.rhs):
            if v: Aug[i, m] = v
        for i in range(n): Aug[i, m + 1 + i] = 1
        R, rk = Aug.rref()
        # find a row whose first m entries vanish and whose entry m is nonzero
        for i in range(rk):
            if all(R[i, j] == 0 for j in range(m)) and R[i, m] != 0:
                y = [Fraction(int(R[i, m + 1 + k].p), int(R[i, m + 1 + k].q)) for k in range(n)]
                l = 1
                for c in y: l = math.lcm(l, c.denominator)
                y = [int(c * l) for c in y]
                return y
        return None
    def check_farkas(self, y):
        yM = [sum(y[i] * col[i] for i in range(self.nrows)) for col in self.cols]
        yb = sum(y[i] * self.rhs[i] for i in range(self.nrows))
        return all(v == 0 for v in yM), yb


# ---------------------------------------------------------------- vectorised polynomial products into a fixed monomial index
def arr(P, p):
    return to_arr(P, p)
def mulvec(E1, C1, E2, C2, index, p):
    """coefficient vector (over index) of the product of two polys given as exponent/coefficient arrays (mod p)."""
    if E1.shape[0] == 0 or E2.shape[0] == 0:
        return np.zeros(index.n, dtype=np.int64)
    E = (E1[:, None, :] + E2[None, :, :]).reshape(-1, 6)
    C = ((C1[:, None] * C2[None, :]) % p).reshape(-1)
    pos = index.order[np.searchsorted(index.sorted, codes_of(E))]
    # exact: values < 2^31, count per bucket < 2^20  -> float64 sums exact (< 2^53)
    assert E.shape[0] < 4000000
    v = np.bincount(pos, weights=C.astype(np.float64), minlength=index.n)
    return (np.rint(v).astype(np.int64)) % p

def relaxation_matrix(pairs, index, p, verbose=False):
    """rows: for k<=l the coefficient vector (over index = monomials of degree 2d, right parity) of
    R_kl + R_lk (k<l) or R_kk, where R_kl = D_k div(u_l) - u_l . grad D_k.  Returns (keys, matrix)."""
    n = len(pairs); t0 = time.time()
    Darr = [arr(D, p) for D, _ in pairs]
    gradD = [[arr(core.pderiv(D, i, p), p) for i in range(6)] for D, _ in pairs]
    divu = [arr(divergence(u, p), p) for _, u in pairs]
    uarr = [[arr(u[i], p) for i in range(6)] for _, u in pairs]
    def Rkl(k, l):
        v = mulvec(Darr[k][0], Darr[k][1], divu[l][0], divu[l][1], index, p)
        for i in range(6):
            v = v - mulvec(uarr[l][i][0], uarr[l][i][1], gradD[k][i][0], gradD[k][i][1], index, p)
        return v % p
    keys = []; rows = []
    for k in range(n):
        for l in range(k, n):
            if k == l: rows.append(Rkl(k, k))
            else: rows.append((Rkl(k, l) + Rkl(l, k)) % p)
            keys.append((k, l))
    M = np.array(rows, dtype=np.int64) if rows else np.zeros((0, index.n), dtype=np.int64)
    if verbose: print(f'[relax-matrix] n={n} rows={M.shape[0]} mons={index.n} {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
    return keys, M

def coords_from_matrix(reducer, Pm):
    """normal-form coordinates for a dense coefficient matrix Pm (rows x index.n) -- Reducer version."""
    p = reducer.p
    if reducer.rank:
        A = nmod_mat(Pm.shape[0], reducer.rank, (Pm[:, reducer.pivots] % p).ravel().tolist(), p)
        B = nmod_mat(reducer.rank, reducer.index.n, reducer.rref.ravel().tolist(), p)
        NF = (Pm - np.array((A * B).tolist(), dtype=np.int64)) % p
    else:
        NF = Pm % p
    return NF[:, reducer.standard]

def folded_coords_from_matrix(fred, Pm, index):
    """FoldedReducer coordinates for a dense full coefficient matrix over `index` (same degree/parity)."""
    p = fred.p
    E = np.array(index.mons, dtype=np.int64)
    fidx = np.searchsorted(fred.orbit_codes, fred.G.canon(E))
    Fm = np.zeros((Pm.shape[0], fred.nrows), dtype=np.int64)
    for i in range(Pm.shape[0]):
        Fm[i] = (np.rint(np.bincount(fidx, weights=(Pm[i] % p).astype(np.float64), minlength=fred.nrows)).astype(np.int64)) % p
    if fred.rank:
        A = nmod_mat(Fm.shape[0], fred.rank, Fm[:, fred.pivots].ravel().tolist(), p)
        B = nmod_mat(fred.rank, fred.nrows, fred.rref.ravel().tolist(), p)
        NF = (Fm - np.array((A * B).tolist(), dtype=np.int64)) % p
    else:
        NF = Fm
    return NF[:, fred.standard]

def quadrics_from_coords(keys, C):
    quads = []
    for s in range(C.shape[1]):
        col = C[:, s]
        nz = np.nonzero(col)[0]
        if len(nz) == 0: continue
        quads.append({keys[i]: int(col[i]) for i in nz})
    return quads

# ---------------------------------------------------------------- analysis of a leading (monomial) ideal from msolve -g 1
def parse_leading_ideal(leadfile, names):
    import re
    txt = open(leadfile).read()
    body = txt[txt.index('['):] if '[' in txt else ''
    mons = []
    for tok in re.split(r'[,\n\[\]]', body):
        tok = tok.strip().rstrip(':')
        if not tok: continue
        supp = {}
        ok = True
        for fac in tok.split('*'):
            m = re.fullmatch(r'(\d+)?\s*([A-Za-z_]\w*)(\^(\d+))?', fac.strip())
            if not m:
                if re.fullmatch(r'\d+', fac.strip()): continue
                ok = False; break
            v = m.group(2); k = int(m.group(4)) if m.group(4) else 1
            supp[v] = supp.get(v, 0) + k
        if ok and supp: mons.append(supp)
    return mons
def krull_dim_monomial(mons, names):
    """Krull dimension of k[names]/(monomial ideal) = max size of a variable set S such that no generator has support inside S."""
    supports = [frozenset(m.keys()) for m in mons]
    if any(len(s) == 0 for s in supports): return -1
    best = [0]
    names = list(names)
    # branch and bound: choose S greedily with pruning
    def rec(i, S):
        if len(S) + (len(names) - i) <= best[0]: return
        if i == len(names):
            best[0] = max(best[0], len(S)); return
        v = names[i]
        S2 = S | {v}
        if all(not (s <= S2) for s in supports if v in s):
            rec(i + 1, S2)
        rec(i + 1, S)
    rec(0, frozenset())
    return best[0]
def dehomogenize_msolve(msfile, outfile, var_index, p):
    """write the chart a_{var_index} = 1 version of a homogeneous msolve input (substitute the variable by 1)."""
    lines = open(msfile).read().splitlines()
    names = lines[0].split(',')
    v = names[var_index]
    rest = [n for n in names if n != v]
    import re
    polys = lines[2:]
    def sub(poly):
        poly = re.sub(rf'\*{v}\^\d+', '', poly); poly = re.sub(rf'\*{v}(?![\w])', '', poly)
        poly = re.sub(rf'(?<![\w]){v}\^\d+\*', '', poly); poly = re.sub(rf'(?<![\w]){v}\*', '', poly)
        poly = re.sub(rf'(?<![\w]){v}\^\d+(?![\w])', '1', poly); poly = re.sub(rf'(?<![\w]){v}(?![\w])', '1', poly)
        return poly
    with open(outfile, 'w') as f:
        f.write(','.join(rest) + '\n' + f'{p}\n')
        f.write('\n'.join(sub(pl) for pl in polys) + '\n')
    return rest
