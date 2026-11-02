"""Helper library for the off-paired attack: character bookkeeping, reduced solvers, logging."""
import sys, time, json, itertools, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
import flint
from flint import nmod_mat

IID = 'off-paired'
LOG = '/tmp/hodge/results/off-paired/log.jsonl'
WG = (0, 2, 0, 1, 3, 3)   # generator g of the diagonal symmetry group (sixth-root weights)
WH = (0, 0, 1, 0, 0, 5)   # generator h
P0 = core.PRIMES[0]

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def char(e):
    return (sum(a*b for a, b in zip(WG, e)) % 6, sum(a*b for a, b in zip(WH, e)) % 6)

def char_of(P):
    cs = {char(e) for e in P}
    return cs.pop() if len(cs) == 1 else None

def cadd(c1, c2):
    return ((c1[0]+c2[0]) % 6, (c1[1]+c2[1]) % 6)

def cneg(c):
    return ((-c[0]) % 6, (-c[1]) % 6)

VARCHAR = [char(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def monos_char(deg, c, par=None):
    """Monomials of degree deg with (g,h)-character c (par ignored; tau is in G)."""
    return [m for m in monomials(deg) if char(m) == c]

class ReducedSystem:
    """Fixed-D certificate system restricted to the character sector of an eigenvector D."""
    def __init__(self, F, A, q, D, p=P0, verbose=False):
        t0 = time.time()
        self.p = p; self.q = q
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
        cD = char_of(Dp); cA = char_of(Ap)
        assert cD is not None, 'D not an eigenvector'
        self.cD = cD
        d = core.degree(D); self.d = d
        degN = d + 6*q - 5; degQ = 2*d - 6
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        G = [scale(add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p), p-1, p) for i in range(NV)]  # -(q D F_i + F D_i)  [engine v2 sign fix]
        Fq1 = scale(power(Fp, q+1, p), p-1, p)
        cols = []; colpolys = []
        for i in range(NV):
            ci = cadd(cadd(cD, VARCHAR[i]), cA)
            for m in monos_char(degN, ci):
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        self.nN = len(cols)
        if degQ >= 0:
            cQ = cadd(cA, cadd(cD, cD))
            for m in monos_char(degQ, cQ):
                cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
        self.nQ = len(cols) - self.nN
        rhs = scale(mul(Ap, mul(Dp, Dp, p), p), p-1, p)
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows: rows[e] = len(rows)
        for e in rhs:
            if e not in rows: rows[e] = len(rows)
        self.cols = cols; self.colpolys = colpolys; self.rows = rows; self.rhs = rhs
        self.nrows = len(rows); self.ncols = len(cols)
        self.t_build = time.time() - t0
        if verbose:
            print(f'[reduced] d={d} char={cD} rows={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}) built {self.t_build:.1f}s', file=sys.stderr)

    def solve(self, verbose=False):
        t0 = time.time()
        R, C, p = self.nrows, self.ncols, self.p
        W = C + 1
        flat = [0] * (R * W)
        rows = self.rows
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                flat[rows[e]*W + j] = v
        for e, v in self.rhs.items():
            flat[rows[e]*W + C] = v
        M = nmod_mat(R, W, flat, p)
        t1 = time.time()
        Mr, rank = M.rref()
        t2 = time.time()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0]*C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        # rank of the coefficient matrix alone
        rank0 = rank - (0 if consistent else 1)
        res = {'consistent': consistent, 'rank': rank0, 'nullity': C - rank0, 'x': x, 'pivots': pivots,
               'rows': R, 'cols': C, 't_build': round(t1-t0, 2), 't_rref': round(t2-t1, 2)}
        if verbose:
            print(f'[reduced solve] consistent={consistent} rank={rank0} nullity={C-rank0} rref {t2-t1:.1f}s', file=sys.stderr)
        return res

    def unpack(self, x):
        N = [{} for _ in range(NV)]; Q = {}
        for (kind, i, m), v in zip(self.cols, x):
            if not v: continue
            if kind == 'N': N[i][m] = v
            else: Q[m] = v
        return N, Q


def reduced_solve(F, A, q, D, p=P0, verbose=False):
    S = ReducedSystem(F, A, q, D, p, verbose)
    r = S.solve(verbose)
    r['S'] = S
    return r


def lift_solution(F, A, q, D, primes=None, verbose=False, max_primes=10):
    """Multi-prime reduced solve + CRT + rational reconstruction + exact verification."""
    primes = primes or core.PRIMES
    acc = None; piv0 = None; S = None
    for k, p in enumerate(primes[:max_primes]):
        S = ReducedSystem(F, A, q, D, p, verbose)
        r = S.solve(verbose)
        if not r['consistent']:
            return 'inconsistent', r
        if piv0 is None: piv0 = r['pivots']
        elif piv0 != r['pivots']:
            acc = None; piv0 = r['pivots']
        x = r['x']
        if acc is None: acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m*p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, i, mon), v in zip(S.cols, rat):
                if v == 0: continue
                if kind == 'N': N[i][mon] = v
                else: Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                return 'solved', {'N': N, 'Q': Q, 'primes_used': k+1, 'rank': r['rank'], 'nullity': r['nullity']}
    return 'lift-failed', {'primes_used': max_primes}


def e1_reduced(F, A, q, D, p=P0):
    """Necessary test A*D in J(F) using the character sector only."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    cD = char_of(Dp); cA = char_of(Ap)
    d = core.degree(D); degN = d + 6*q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    for i in range(NV):
        ci = cadd(cadd(cD, VARCHAR[i]), cA)
        for m in monos_char(degN, ci):
            cols.append(shift(Fi[i], m, q, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*(C+1))
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*(C+1)+j] = v
    for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = v
    M = nmod_mat(R, C+1, flat, p)
    ra = M.rank()
    M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    r0 = M0.rank()
    return ra == r0, {'rank': r0, 'rank_aug': ra, 'rows': R, 'cols': C}


def colon_space(F, A, q, d, c, p=P0, verbose=False):
    """Space of D of degree d, character c, with q*N(F) = A*D solvable.  Returns (dim, basis mod p, info).
    Basis vectors are polys (dict) mod p in the degree-d character-c monomials."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    cA = char_of(Ap)
    degN = d + 6*q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []; kinds = []
    for i in range(NV):
        ci = cadd(cadd(c, VARCHAR[i]), cA)
        for m in monos_char(degN, ci):
            cols.append(shift(Fi[i], m, q, p)); kinds.append(('N', i, m))
    Dmons = monos_char(d, c)
    for m in Dmons:
        cols.append(scale(shift(Ap, m, 1, p), p-1, p)); kinds.append(('D', -1, m))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    if C == 0 or len(Dmons) == 0:
        return 0, [], {'rows': R, 'cols': C, 'nD': len(Dmons)}
    flat = [0]*(R*C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*C+j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    didx = [j for j, k in enumerate(kinds) if k[0] == 'D']
    if nullity == 0:
        return 0, [], {'rows': R, 'cols': C, 'nD': len(Dmons)}
    proj = nmod_mat(len(didx), nullity, [int(X[j, k]) for j in didx for k in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(didx):
            v = int(Tr[i, jj])
            if v: poly[kinds[j][2]] = v
        basis.append(poly)
    return rk, basis, {'rows': R, 'cols': C, 'nD': len(Dmons), 'nullity': nullity}


def equivariant_hilbert(deg, c):
    """dim of (S/J(F))_{deg, char c} assuming J(F) is a complete intersection of six quintics
    (Koszul exactness): coefficient extraction from prod(1 - t^5 chi_i^{-1}) / prod(1 - t chi_i)."""
    # represent series as dict {(deg, char): coeff} truncated at deg
    from collections import defaultdict
    ser = defaultdict(int); ser[(0, (0, 0))] = 1
    # multiply by 1/(1 - t chi_i) for each variable
    for i in range(NV):
        wi = VARCHAR[i]
        new = defaultdict(int)
        for (dg, ch), v in ser.items():
            k = 0
            while dg + k <= deg:
                chk = ((ch[0] + k*wi[0]) % 6, (ch[1] + k*wi[1]) % 6)
                new[(dg + k, chk)] += v
                k += 1
        ser = new
    # multiply by (1 - t^5 chi_i^{-1})
    for i in range(NV):
        wi = cneg(VARCHAR[i])
        new = defaultdict(int)
        for (dg, ch), v in ser.items():
            new[(dg, ch)] += v
            if dg + 5 <= deg:
                new[(dg + 5, cadd(ch, wi))] -= v
        ser = new
    return ser[(deg, c)]


def rand_combo(basis, p=P0, seed=None):
    rng = random.Random(seed)
    out = {}
    for b in basis:
        cco = rng.randrange(1, p)
        out = add(out, scale(b, cco, p), p)
    return out


def lift_poly_small(P, p):
    """Rational reconstruction of a mod-p poly's coefficients (small rationals)."""
    out = {}
    for e, v in P.items():
        r = core.ratrecon(v, p)
        if r is None: return None
        if r: out[e] = r
    return out


def served_targets(F, q, D, p=P0, cA=(3, 0), verbose=False):
    """Space of degree-6q targets of character cA served by eigenvector D (mod p), and its dimension modulo J(F).
    Returns dict with dim_served, dim_J_in_sector, dim_nontrivial, basis (list of polys mod p)."""
    Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p)
    cD = char_of(Dp)
    d = core.degree(D)
    degN = d + 6*q - 5; degQ = 2*d - 6
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    FD = mul(Fp, Dp, p)
    G = [scale(add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p), p-1, p) for i in range(NV)]
    Fq1 = scale(power(Fp, q+1, p), p-1, p)
    D2 = mul(Dp, Dp, p)
    cols = []; colpolys = []
    for i in range(NV):
        ci = cadd(cadd(cD, VARCHAR[i]), cA)
        for m in monos_char(degN, ci):
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); colpolys.append(col)
    if degQ >= 0:
        for m in monos_char(degQ, cadd(cA, cadd(cD, cD))):
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
    Amons = monos_char(6*q, cA)
    for m in Amons:
        cols.append(('A', -1, m)); colpolys.append(shift(D2, m, 1, p))
    rows = {}
    for col in colpolys:
        for e in col:
            if e not in rows: rows[e] = len(rows)
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): flat[rows[e]*C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    aidx = [j for j, c in enumerate(cols) if c[0] == 'A']
    if nullity == 0:
        served = []
    else:
        proj = nmod_mat(len(aidx), nullity, [int(X[j, k]) for j in aidx for k in range(nullity)], p)
        Tr, rk = proj.transpose().rref()
        served = []
        for i in range(rk):
            poly = {}
            for jj, j in enumerate(aidx):
                v = int(Tr[i, jj])
                if v: poly[cols[j][2]] = v
            served.append(poly)
    # Jacobian part of the sector
    Jvecs = []
    for i in range(NV):
        for j in range(NV):
            if cadd(VARCHAR[i], cneg(VARCHAR[j])) == cA:
                Jvecs.append(shift(Fi[j], tuple(1 if k == i else 0 for k in range(NV)), 1, p))
    idx = {m: k for k, m in enumerate(Amons)}
    def vec(P):
        v = [0]*len(Amons)
        for e, val in P.items(): v[idx[e]] = val % p
        return v
    rowsJ = [vec(P) for P in Jvecs]
    rkJ = nmod_mat(len(rowsJ), len(Amons), [x for r in rowsJ for x in r], p).rank() if rowsJ else 0
    rowsS = rowsJ + [vec(P) for P in served]
    rkS = nmod_mat(len(rowsS), len(Amons), [x for r in rowsS for x in r], p).rank() if rowsS else 0
    return {'dim_served': len(served), 'dim_J': rkJ, 'dim_nontrivial': rkS - rkJ, 'basis': served, 'rows': R, 'cols': C}


# ---------------------------------------------------------------------------
# group switching: 'G' = full diagonal group (chars in Z/6 x Z/6); 'Gp' = index-2 subgroup
# ker(chi_(3,0)) (chars (2*(g mod 3), h mod 6) embedded in Z/6 x Z/6)
# ---------------------------------------------------------------------------
def _char_G(e):
    return (sum(a*b for a, b in zip(WG, e)) % 6, sum(a*b for a, b in zip(WH, e)) % 6)

def _char_Gp(e):
    return (2 * (sum(a*b for a, b in zip(WG, e)) % 3), sum(a*b for a, b in zip(WH, e)) % 6)

def set_group(name):
    global char, VARCHAR, GROUP
    GROUP = name
    char = _char_G if name == 'G' else _char_Gp
    VARCHAR = [char(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
    monos_char.cache_clear() if hasattr(monos_char, 'cache_clear') else None
    return char

GROUP = 'G'
