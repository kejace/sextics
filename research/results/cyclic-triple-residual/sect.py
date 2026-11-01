"""Character-sector-reduced fixed-D solver (any q) and E1-kernel sweeps, exact over Q or mod p.
Uses the full projective diagonal symmetry group of F (automatically computed) restricted to the subgroup
under which D (and A) are semi-invariant.  Split system (equivalent to the single identity, see notes):
   E1:  q * sum_i N_i F_i = A D          (rows: degree d+6q, character chi_A + chi_D)
   E2:  D div N - N(D) = F^q Q           (rows: degree 2d,   character chi_A + 2 chi_D)
   chi(N_i) = chi_A + chi_D + e_i,  chi(Q) = chi_A + 2 chi_D.
"""
import sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from math import lcm
from flint import fmpz_mat, nmod_mat, fmpq_mat, fmpq, fmpz
from itertools import product

PRIMES = core.PRIMES

def integral(P):
    """Scale a Fraction-coefficient poly to integer coefficients (content not removed)."""
    L = 1
    for c in P.values():
        L = lcm(L, Fraction(c).denominator)
    return {m: int(Fraction(c) * L) for m, c in P.items() if c}

def elems(F):
    """All g in (Z/6)^6 (normalized g5 = 0) with g.m = 0 mod 6 for every monomial m of F."""
    mons = list(F)
    out = []
    for g in product(range(6), repeat=5):
        g = g + (0,)
        if all(sum(a*b for a, b in zip(m, g)) % 6 == 0 for m in mons):
            out.append(g)
    return out

def semi_invariant_subgroup(G, polys):
    """Elements of G under which every poly in polys is semi-invariant (all its monomials share a character)."""
    out = []
    for g in G:
        ok = True
        for P in polys:
            ms = list(P)
            if not ms: continue
            c0 = sum(a*b for a, b in zip(ms[0], g)) % 6
            if any(sum(a*b for a, b in zip(m, g)) % 6 != c0 for m in ms[1:]):
                ok = False; break
        if ok: out.append(g)
    return out

def generators(G):
    """Greedy generating set of the finite abelian group G (list of tuples mod 6, g5 = 0)."""
    gens = []; span = {tuple([0]*6)}
    for g in G:
        if g in span: continue
        gens.append(g)
        new = set()
        for s in span:
            for k in range(6):
                new.add(tuple((a + k*b) % 6 for a, b in zip(s, g)))
        span = new
    assert len(span) == len(G), (len(span), len(G))
    return gens

def make_chi(gens):
    def chi(m):
        return tuple(sum(a*b for a, b in zip(m, g)) % 6 for g in gens)
    return chi

def cadd(*cs):
    return tuple(sum(x) % 6 for x in zip(*cs))

def chi_of(P, chi):
    cs = set(chi(m) for m in P)
    assert len(cs) == 1, 'not semi-invariant'
    return cs.pop()

def imul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a+b for a, b in zip(m1, m2)); r[m] = r.get(m, 0) + c1*c2
    return {m: c for m, c in r.items() if c}

def ideriv(p, i):
    r = {}
    for m, c in p.items():
        if m[i]:
            mm = list(m); mm[i] -= 1; r[tuple(mm)] = r.get(tuple(mm), 0) + c*m[i]
    return r

def ipow(p, k):
    r = {(0,)*6: 1}
    for _ in range(k): r = imul(r, p)
    return r

def mons(deg):
    if deg < 0: return
    for c in itertools.combinations(range(deg+5), 5):
        prev = -1; e = []
        for x in c: e.append(x-prev-1); prev = x
        e.append(deg+4-prev); yield tuple(e)

E6 = [tuple(int(j == i) for j in range(6)) for i in range(6)]

class Sys:
    def __init__(self, F, A, q, D, verbose=False):
        """F, A, D with Fraction/int coefficients; scaled internally to integers."""
        t0 = time.time()
        self.q = q
        Fz, Az, Dz = integral(F), integral(A), integral(D)
        self.Fz, self.Az, self.Dz = Fz, Az, Dz
        self.D_orig = D
        d = max(sum(m) for m in Dz); assert all(sum(m) == d for m in Dz), 'D not homogeneous'
        self.d = d
        G = elems(Fz)
        G = semi_invariant_subgroup(G, [Az, Dz])
        self.G_order = len(G)
        gens = generators(G); self.gens = gens
        chi = make_chi(gens); self.chi = chi
        chiA, chiD = chi_of(Az, chi), chi_of(Dz, chi)
        cN = [cadd(chiA, chiD, chi(E6[i])) for i in range(6)]; cQ = cadd(chiA, chiD, chiD)
        self.cN, self.cQ = cN, cQ
        dF = [ideriv(Fz, i) for i in range(6)]; dD = [ideriv(Dz, i) for i in range(6)]
        degN = d + 6*q - 5; degQ = 2*d - 6
        r1 = [m for m in mons(d + 6*q) if chi(m) == cadd(chiA, chiD)]
        r2 = [m for m in mons(2*d) if chi(m) == cQ]
        ri = {('1', m): k for k, m in enumerate(r1)}; n1 = len(r1)
        for k, m in enumerate(r2): ri[('2', m)] = n1 + k
        self.nrows = len(r1) + len(r2); self.n1 = n1
        cols = []; meta = []
        for i in range(6):
            for m in mons(degN):
                if chi(m) != cN[i]: continue
                col = {}
                for e, c in imul({m: 1}, dF[i]).items():
                    k = ri[('1', e)]; col[k] = col.get(k, 0) + q*c
                t2 = imul(Dz, ideriv({m: 1}, i))
                for e, c in imul({m: 1}, dD[i]).items(): t2[e] = t2.get(e, 0) - c
                for e, c in t2.items():
                    if c:
                        k = ri[('2', e)]; col[k] = col.get(k, 0) + c
                cols.append(col); meta.append(('N', i, m))
        self.nN = len(cols)
        Fq = ipow(Fz, q)
        for m in mons(degQ):
            if chi(m) != cQ: continue
            col = {}
            for e, c in imul({m: 1}, Fq).items():
                k = ri[('2', e)]; col[k] = col.get(k, 0) - c
            cols.append(col); meta.append(('Q', -1, m))
        self.nQ = len(cols) - self.nN
        rhs = {}
        for e, c in imul(Az, Dz).items():
            k = ri[('1', e)]; rhs[k] = rhs.get(k, 0) + c
        self.cols, self.meta, self.rhs = cols, meta, rhs
        self.ncols = len(cols)
        self.nnz = sum(len(c) for c in cols)
        self.t_build = round(time.time() - t0, 2)
        if verbose:
            print(f'[sect] d={d} q={q} |G|={self.G_order} rows={self.nrows} (E1 {n1}) cols={self.ncols} (N {self.nN}, Q {self.nQ}) nnz={self.nnz} built {self.t_build}s', file=sys.stderr, flush=True)

    def dense(self, p=None, aug=True):
        R, C = self.nrows, self.ncols
        W = C + (1 if aug else 0)
        flat = [0] * (R * W)
        for j, col in enumerate(self.cols):
            for r, v in col.items(): flat[r*W + j] = v if p is None else v % p
        if aug:
            for r, v in self.rhs.items(): flat[r*W + C] = v if p is None else v % p
        return fmpz_mat(R, W, flat) if p is None else nmod_mat(R, W, flat, p)

    def ranks(self, p=None):
        """(rank, rank_aug) exactly (p None) or mod p."""
        t0 = time.time()
        M = self.dense(p, aug=True)
        ra = M.rank()
        M0 = self.dense(p, aug=False)
        r0 = M0.rank()
        return r0, ra, round(time.time() - t0, 2)

    def solve_exact(self):
        """Exact rational particular solution via fmpq_mat (full rows). Returns (N, Q) as Fraction polys, or None."""
        R, C = self.nrows, self.ncols
        M = self.dense(None, aug=False)
        b = fmpz_mat(R, 1, [self.rhs.get(r, 0) for r in range(R)])
        # least-squares-free approach: use rref of augmented matrix over Q
        Maug = fmpq_mat(self.dense(None, aug=True))
        Mr, rank = Maug.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and Mr[i, j] == 0: j += 1
            pivots.append(j)
        if pivots and pivots[-1] == C: return None
        x = [Fraction(0)] * C
        for i, j in enumerate(pivots):
            v = Mr[i, C]; x[j] = Fraction(int(v.p), int(v.q))
        N = [{} for _ in range(6)]; Q = {}
        for (kind, i, m), v in zip(self.meta, x):
            if v == 0: continue
            if kind == 'N': N[i][m] = v
            else: Q[m] = v
        # we solved with integral Fz, Az, Dz: identity q sum N_i Fz_i = Az Dz ; B(N) = Fz^q Q.
        # Convert to original F, A, D: F = Fz/LF, A = Az/LA, D = Dz/LD.
        LF = Fraction(next(iter(self.Fz.values()))) / Fraction(next(iter(core.reduce_mod(self.Fz, 10**9) and self.Fz.values()))) if False else None
        return N, Q

def rescale_solution(F, A, D, Fz, Az, Dz, N, Q, q):
    """Given a solution of the integral system, return (N', Q') for the original (F, A, D)."""
    def ratio(P, Pz):
        m = next(iter(Pz)); return Fraction(P[m]) / Fraction(Pz[m])   # P = ratio * Pz
    rF, rA, rD = ratio(F, Fz), ratio(A, Az), ratio(D, Dz)
    # q sum N_i Fz_i = Az Dz  =>  q sum (N_i * rA*rD/rF) F_i = A D
    s = rA * rD / rF
    N2 = [core.scale(n, s) for n in N]
    # B(N) = Fz^q Q ; B is linear in N and in D: B_{D}(N) with D = rD Dz: B_D(N2) = rD * s * B_Dz(N) = rD*s*Fz^q*Q = rD*s*(F/rF)^q Q
    Q2 = core.scale(Q, rD * s / rF**q)
    return N2, Q2

def test(F, A, q, D, exact=True, primes=(PRIMES[0], PRIMES[1]), verbose=True, lift=True):
    S = Sys(F, A, q, D, verbose=verbose)
    out = {'d': S.d, 'G_order': S.G_order, 'rows': S.nrows, 'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ, 'nnz': S.nnz, 't_build': S.t_build}
    if exact and S.nrows * S.ncols <= 4_000_000:
        r0, ra, t = S.ranks(None)
        out.update(rank=r0, rank_aug=ra, nullity=S.ncols - r0, consistent=(r0 == ra), exact_rank=True, t_rank=t)
    else:
        r0, ra, t = S.ranks(primes[0])
        out.update(rank=r0, rank_aug=ra, nullity=S.ncols - r0, consistent=(r0 == ra), exact_rank=False, prime=primes[0], t_rank=t)
        if not out['consistent'] and len(primes) > 1:
            r0b, rab, tb = S.ranks(primes[1])
            out.update(rank2=r0b, rank_aug2=rab, consistent2=(r0b == rab), prime2=primes[1])
            out['consistent'] = out['consistent'] or out['consistent2']
    if out['consistent'] and lift:
        sol = S.solve_exact()
        if sol is not None:
            N, Q = rescale_solution(F, A, D, S.Fz, S.Az, S.Dz, sol[0], sol[1], q)
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            out['verified'] = bool(ok)
            out['N'] = [core.to_str(n) for n in N]; out['Q'] = core.to_str(Q)
        else:
            out['verified'] = False; out['note'] = 'exact lift failed although ranks agree'
    if verbose:
        print({k: v for k, v in out.items() if k not in ('N', 'Q')}, file=sys.stderr, flush=True)
    return out

# ---------------- E1-kernel sweep by sectors ----------------
def e1_kernel(F, A, d, p=PRIMES[0], want_vectors=True, G=None):
    """dim {D in S_d : A D in J(F)} / J(F)_d, split by character of the symmetry group of F (restricted
    so that A is semi-invariant).  Returns dict chi -> (kerA_in_S, dimJ, ker_mod_J, kernel basis polys mod p)."""
    Fz, Az = integral(F), integral(A)
    if G is None:
        G = semi_invariant_subgroup(elems(Fz), [Az])
    gens = generators(G); chi = make_chi(gens)
    chiA = chi_of(Az, chi)
    dF = [ideriv(Fz, i) for i in range(6)]
    # group degree-d monomials by character
    bych = {}
    for m in mons(d): bych.setdefault(chi(m), []).append(m)
    res = {}
    for c, Dm in sorted(bych.items()):
        t0 = time.time()
        target = cadd(chiA, c)
        cols = []; kinds = []
        for m in Dm:
            cols.append({tuple(a+b for a, b in zip(m, e)): (v % p) for e, v in Az.items()}); kinds.append('A')
        for i in range(6):
            ci = cadd(target, chi(E6[i]))   # chi(m') + chi(F_i) = target, chi(F_i) = -e_i
            for m in mons(d + 1):
                if chi(m) != ci: continue
                col = {}
                for e, v in dF[i].items():
                    k = tuple(a+b for a, b in zip(m, e)); col[k] = (col.get(k, 0) + v) % p
                cols.append(col); kinds.append('J')
        rows = {}
        for col in cols:
            for e in col: rows.setdefault(e, len(rows))
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, col in enumerate(cols):
            for e, v in col.items(): flat[rows[e]*C + j] = v
        M = nmod_mat(R, C, flat, p)
        X, nullity = M.nullspace()
        nA = len(Dm)
        vecs = []
        if nullity:
            proj = nmod_mat(nA, nullity, [int(X[j, k]) for j in range(nA) for k in range(nullity)], p)
            # column-reduce to get a basis of the projection
            Pr, kA = proj.transpose().rref()
            for i in range(kA):
                vecs.append({Dm[j]: int(Pr[i, j]) for j in range(nA) if int(Pr[i, j])})
        else:
            kA = 0
        # J_d in this character: columns F_i * S_{d-5}
        dimJ = 0
        if d >= 5:
            cols2 = []
            for i in range(6):
                ci = cadd(c, chi(E6[i]))
                for m in mons(d - 5):
                    if chi(m) != ci: continue
                    col = {}
                    for e, v in dF[i].items():
                        k = tuple(a+b for a, b in zip(m, e)); col[k] = (col.get(k, 0) + v) % p
                    cols2.append(col)
            if cols2:
                rows2 = {}
                for col in cols2:
                    for e in col: rows2.setdefault(e, len(rows2))
                flat2 = [0] * (len(rows2) * len(cols2))
                for j, col in enumerate(cols2):
                    for e, v in col.items(): flat2[rows2[e]*len(cols2) + j] = v
                dimJ = nmod_mat(len(rows2), len(cols2), flat2, p).rank()
        res[c] = {'dim_S': nA, 'kerA_in_S': kA, 'dimJ': dimJ, 'ker_mod_J': kA - dimJ, 'rows': R, 'cols': C,
                  'seconds': round(time.time() - t0, 2), 'vectors': vecs if want_vectors else None}
    return res, gens

if __name__ == '__main__':
    # self-test against known exact results: fqlx cone X2^6+Y0^6+Y1^6+Y2^6 -> cols 440, rank 436, aug 437
    F, A, q, rec = core.load_instance('cyclic-triple-residual')
    D = core.parse('X2**6 + Y0**6 + Y1**6 + Y2**6')
    print(test(F, A, q, D))
    # control: D = F (illegal) must be consistent and lift (checks the lifting/rescaling code path)
    print(test(F, A, q, F)['verified'])
