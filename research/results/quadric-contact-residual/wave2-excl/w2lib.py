"""Wave-2 EXCLUSION library for quadric-contact-residual (independent re-implementation; new prime).

Conventions: variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5; polys are dicts {exp6: coeff}.
G_aff character of a monomial (fqlx): chi(e) = (e0+e3 mod 2, e2+e5 mod 2, e1-e4 mod 6, e3+e4+e5 mod 2).
Default prime P2 = 2147483629 (core.PRIMES[1]); wave 1 used 2147483647, the dossier 1000000009, Singular 32003.
"""
import sys, os, json, time, itertools, subprocess, re
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
import flint
from flint import nmod_mat

IID = 'quadric-contact-residual'
LANE = '/tmp/hodge/results/quadric-contact-residual'
PK = LANE + '/wave2-excl'
LOG = os.environ.get('HODGE_W2LOG', LANE + '/wave2-excl.log.jsonl')
P2 = int(os.environ.get('HODGE_P', core.PRIMES[1]))

F, A, q, rec = core.load_instance(IID)
assert q == 1
p = P2
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
Hess = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
CHI0 = (0, 0, 0, 0)

def chi(e):
    return ((e[0] + e[3]) % 2, (e[2] + e[5]) % 2, (e[1] - e[4]) % 6, (e[3] + e[4] + e[5]) % 2)
def chi_add(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2, (a[2] + b[2]) % 6, (a[3] + b[3]) % 2)
def chi_neg(a):
    return ((-a[0]) % 2, (-a[1]) % 2, (-a[2]) % 6, (-a[3]) % 2)
def chi_of(P):
    cs = {chi(e) for e in P}
    return cs.pop() if len(cs) == 1 else None
CHI_A = chi_of(A)
CHI_Z = [chi(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
assert chi_of(F) == CHI0 and CHI_A == (0, 0, 0, 1)

def log(exp, **kw):
    rec = {'exp': exp, 'prime': p, 'time': time.strftime('%Y-%m-%dT%H:%M:%S')}
    rec.update(kw)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print('[LOG]', json.dumps(rec)[:400], file=sys.stderr, flush=True)
    return rec

def vfield(N, P):
    """N(P) = sum N_i d_i P"""
    out = {}
    for i in range(NV):
        if N[i]:
            out = add(out, mul(N[i], pderiv(P, i, p), p), p)
    return out
def divg(N):
    out = {}
    for i in range(NV):
        out = add(out, pderiv(N[i], i, p), p)
    return out
def BD(D, N):
    """B_D(N) = D div N - N(D)"""
    return sub(mul(D, divg(N), p), vfield(N, D), p)

# --------------------------------------------------------------------------------------------
class Sector:
    """Degree-m polynomials of G_aff character c; Jacobian subspace J_m[c] (rows m'*F_i) in rref; normal forms."""
    def __init__(self, m, c, verbose=False):
        t0 = time.time()
        self.m, self.c = m, c
        self.mons = [e for e in monomials(m) if chi(e) == c]
        self.idx = {e: i for i, e in enumerate(self.mons)}
        self.n = len(self.mons)
        self.gens = []
        rows = []
        if m >= 5:
            for mm in monomials(m - 5):
                for i in range(NV):
                    if chi_add(chi(mm), chi_neg(CHI_Z[i])) == c:
                        rows.append(shift(Fi[i], mm, 1, p)); self.gens.append((i, mm))
        self.nJ = len(rows)
        if rows:
            flat = [0] * (self.nJ * self.n)
            for r, P in enumerate(rows):
                for e, v in P.items():
                    flat[r * self.n + self.idx[e]] = v
            M = nmod_mat(self.nJ, self.n, flat, p)
            self.Mgen = M
            Mr, rk = M.rref()
            self.rank = rk
            piv = []
            for i in range(rk):
                j = 0
                while int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            self.piv = piv
            self.Rb = nmod_mat(rk, self.n, [int(Mr[i, j]) for i in range(rk) for j in range(self.n)], p) if rk else None
        else:
            self.Mgen, self.rank, self.piv, self.Rb = None, 0, [], None
        ps = set(self.piv)
        self.nonpiv = [j for j in range(self.n) if j not in ps]
        self.dimR = self.n - self.rank
        if verbose:
            print(f'[Sector m={m} c={c}] n={self.n} Jrows={self.nJ} rank={self.rank} dimR={self.dimR} ({time.time()-t0:.1f}s)', file=sys.stderr)
    def vec(self, P):
        v = [0] * self.n
        for e, cv in P.items():
            v[self.idx[e]] = cv % p
        return v
    def nf_matrix(self, vecs):
        """normal forms restricted to non-pivot columns (k x dimR nmod_mat)"""
        k = len(vecs)
        if k == 0:
            return nmod_mat(0, len(self.nonpiv), [], p)
        V = nmod_mat(k, self.n, [x for v in vecs for x in v], p)
        if self.rank:
            Pm = nmod_mat(k, self.rank, [v[j] for v in vecs for j in self.piv], p)
            V = V - Pm * self.Rb
        return nmod_mat(k, len(self.nonpiv), [int(V[i, j]) for i in range(k) for j in self.nonpiv], p)
    def nf(self, P):
        N = self.nf_matrix([self.vec(P)])
        return [int(N[0, j]) for j in range(N.ncols())]
    def in_J(self, P):
        return all(x == 0 for x in self.nf(P))
    def lift(self, P):
        """coefficients x over self.gens with sum x_g * gen = P, or None"""
        v = self.vec(P)
        n, nJ = self.n, self.nJ
        flat = []
        for j in range(n):
            flat.extend(int(self.Mgen[r, j]) for r in range(nJ)); flat.append(v[j])
        Aug = nmod_mat(n, nJ + 1, flat, p)
        R, rk = Aug.rref()
        x = [0] * nJ
        for i in range(rk):
            j = 0
            while int(R[i, j]) == 0: j += 1
            if j == nJ: return None
            x[j] = int(R[i, nJ])
        return x
    def lift_field(self, P):
        x = self.lift(P)
        if x is None: return None
        N = [{} for _ in range(NV)]
        for (i, mm), cv in zip(self.gens, x):
            if cv: N[i] = add(N[i], {mm: cv}, p)
        assert sub(vfield(N, Fp), P, p) == {}, 'lift check failed'
        return N

_SEC = {}
def sector(m, c, verbose=False):
    if (m, c) not in _SEC:
        _SEC[(m, c)] = Sector(m, c, verbose=verbose)
    return _SEC[(m, c)]

# --------------------------------------------------------------------------------------------
def e1_profile(d, verbose=True):
    """per character: dim S_d[c], dim J_d[c], dim{D in S_d[c]: A D in J}, extras (basis polys mod J_d[c])."""
    out, extras = {}, {}
    for c in sorted({chi(e) for e in monomials(d)}):
        mons = [e for e in monomials(d) if chi(e) == c]
        S6 = sector(d + 6, chi_add(c, CHI_A))
        NF = S6.nf_matrix([S6.vec(shift(Ap, e, 1, p)) for e in mons])
        X, nul = NF.transpose().nullspace()
        Sd = sector(d, c)
        extra = nul - Sd.rank
        out[str(c)] = {'dimS': len(mons), 'dimJ': Sd.rank, 'kernel': nul, 'extra': extra}
        if extra > 0:
            kers = [{e: int(X[i, j]) for i, e in enumerate(mons) if int(X[i, j])} for j in range(nul)]
            NFk = Sd.nf_matrix([Sd.vec(P) for P in kers])
            chosen, cur = [], []
            for i, P in enumerate(kers):
                test = nmod_mat(len(cur) + 1, NFk.ncols(), [int(NFk[t, j]) for t in cur + [i] for j in range(NFk.ncols())], p)
                if test.rank() == len(cur) + 1:
                    cur.append(i); chosen.append(P)
                if len(cur) == extra: break
            extras[str(c)] = chosen
        if verbose: print(d, c, out[str(c)], file=sys.stderr, flush=True)
    return out, extras

# --------------------------------------------------------------------------------------------
def quadrics_from_family(Ds, Ns, d, chars=None, verbose=True):
    """Hessian relaxation: for D = sum c_k D_k, N = sum c_k N_k (E1 lifts, N_k(F) = A D_k), the mod-J_{2d}
    necessary condition B_D(N) in J_{2d} gives quadrics sum_{k<=l} q_{kl} c_k c_l = 0, one per coordinate of
    R_{2d} (all characters that occur).  Returns (quads, info) with quads = list of {(k,l): coeff}."""
    t0 = time.time()
    n = len(Ds)
    # group B_{D_k}(N_l) by character of the result; D_k must be eigenvectors (possibly different chars)
    cD = [chi_of(D) for D in Ds]
    assert all(c is not None for c in cD), 'family members must be G_aff eigenvectors'
    divN = [divg(N) for N in Ns]
    quads = []
    by_char = {}
    for k in range(n):
        for l in range(n):
            cc = chi_add(CHI_A, chi_add(cD[k], cD[l]))
            by_char.setdefault(cc, []).append((k, l))
    nq_tot = 0
    for cc, pairs in by_char.items():
        S2d = sector(2 * d, cc, verbose=verbose)
        vecs = []
        for (k, l) in pairs:
            B = sub(mul(Ds[k], divN[l], p), vfield(Ns[l], Ds[k]), p)
            vecs.append(S2d.vec(B))
        NF = S2d.nf_matrix(vecs)
        m = NF.ncols()
        for j in range(m):
            coeffs = {}
            for r, (k, l) in enumerate(pairs):
                v = int(NF[r, j])
                if v:
                    key = (min(k, l), max(k, l))
                    coeffs[key] = (coeffs.get(key, 0) + v) % p
            coeffs = {kk: v for kk, v in coeffs.items() if v}
            if coeffs:
                quads.append(coeffs)
        nq_tot += m
        if verbose: print(f'  char {cc}: {len(pairs)} pairs, dimR={m} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    info = {'n_params': n, 'n_quadrics': len(quads), 'n_target_coords': nq_tot, 'seconds': round(time.time() - t0, 1)}
    return quads, info

def quad_rank_and_squares(quads, n):
    pairs = [(k, l) for k in range(n) for l in range(k, n)]
    pidx = {pr: i for i, pr in enumerate(pairs)}
    if not quads:
        return 0, [], pairs
    flat = [0] * (len(quads) * len(pairs))
    for r, cf in enumerate(quads):
        for pr, v in cf.items():
            flat[r * len(pairs) + pidx[pr]] = v
    QM = nmod_mat(len(quads), len(pairs), flat, p)
    full = QM.rank()
    forced = []
    for k in range(n):
        j = pidx[(k, k)]
        cols = [i for i in range(len(pairs)) if i != j]
        sub_ = nmod_mat(len(quads), len(cols), [int(QM[r, i]) for r in range(len(quads)) for i in cols], p)
        if sub_.rank() < full:
            forced.append(k)
    return full, forced, pairs

def write_systems(quads, n, prefix, drop=()):
    """write msolve (.ms) and Singular (.sing) inputs for the quadric ideal in unknowns c_k, k not in drop."""
    keep = [k for k in range(n) if k not in set(drop)]
    kidx = {k: i for i, k in enumerate(keep)}
    names = [f'c{k}' for k in keep]
    eqs = []
    for cf in quads:
        terms = [f'{v}*c{k}*c{l}' for (k, l), v in cf.items() if k in kidx and l in kidx]
        if terms:
            eqs.append('+'.join(terms))
    with open(prefix + '.ms', 'w') as f:
        f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(eqs) + '\n')
    with open(prefix + '.sing', 'w') as f:
        f.write(f'ring r = {p}, ({",".join(names)}), dp;\nideal I = {",".join(eqs) if eqs else "0"};\n')
        f.write('ideal G = std(I);\nprint("SDIM " + string(dim(G)));\nif (dim(G) == 0) { print("SVDIM " + string(vdim(G))); }\n')
        f.write('if (dim(G) > 0) { LIB "primdec.lib"; list L = minAssGTZ(G); int i; for (i=1;i<=size(L);i++){ print("COMPONENT"); print(L[i]); } }\nquit;\n')
    return len(eqs), len(names)

def singular_verdict(sing_file, timeout=3600):
    """Run Singular; returns dict(sdim, svdim, rc, out_tail).  sdim==0 <=> only the origin (over the algebraic closure of F_p)."""
    t0 = time.time()
    r = subprocess.run(['timeout', str(timeout), 'Singular', '-q', sing_file], capture_output=True, text=True)
    sd = re.search(r'SDIM (-?\d+)', r.stdout); sv = re.search(r'SVDIM (\d+)', r.stdout)
    return {'sdim': int(sd.group(1)) if sd else None, 'svdim': int(sv.group(1)) if sv else None, 'rc': r.returncode,
            'seconds': round(time.time() - t0, 1), 'out_tail': (r.stdout + r.stderr)[-1500:]}

def msolve_verdict(ms_file, threads=4, timeout=3600):
    """msolve on a homogeneous system: output '[-1]:' = inconsistent (impossible for homogeneous), '[0, n, ...' = zero-dim,
    '[1, n, -1, []]:' = positive dimension.  We add nothing; interpret ourselves."""
    t0 = time.time()
    out = ms_file + '.out'
    r = subprocess.run(['timeout', str(timeout), 'msolve', '-v', '1', '-t', str(threads), '-f', ms_file, '-o', out], capture_output=True, text=True)
    txt = open(out).read() if os.path.exists(out) else ''
    return {'rc': r.returncode, 'seconds': round(time.time() - t0, 1), 'head': txt[:200].replace('\n', ' '), 'log_tail': (r.stdout + r.stderr)[-800:]}

# --------------------------------------------------------------------------------------------
class SectorSystem:
    """Sector-restricted fixed-D certificate system (q=1), D a G_aff eigenvector. Independent re-implementation."""
    def __init__(self, D, verbose=False):
        t0 = time.time()
        cD = chi_of(D)
        assert cD is not None
        d = core.degree(D)
        Dp = core.reduce_mod(D, p)
        self.D, self.d, self.cD = Dp, d, cD
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        # unknown N_i = sum x m: contributes  -D F_i m - F D_i m + F D d_i(m)
        G = [scale(add(mul(Dp, Fi[i], p), mul(Fp, Di[i], p), p), p - 1, p) for i in range(NV)]
        F2neg = scale(power(Fp, 2, p), p - 1, p)
        cols, colpolys = [], []
        cA2D = chi_add(CHI_A, chi_add(cD, cD))
        for i in range(NV):
            ci = chi_add(CHI_A, chi_add(cD, CHI_Z[i]))
            for m in monomials(d + 1):
                if chi(m) != ci: continue
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        if 2 * d - 6 >= 0:
            for m in monomials(2 * d - 6):
                if chi(m) != cA2D: continue
                cols.append(('Q', -1, m)); colpolys.append(shift(F2neg, m, 1, p))
        rhs = scale(mul(Ap, mul(Dp, Dp, p), p), p - 1, p)   # move A D^2 to the right: sum cols = -A D^2
        rows = {}
        for col in colpolys:
            for e in col: rows.setdefault(e, len(rows))
        for e in rhs: rows.setdefault(e, len(rows))
        assert all(chi(e) == cA2D for e in rows)
        self.cols, self.colpolys, self.rows, self.rhs = cols, colpolys, rows, rhs
        self.nrows, self.ncols = len(rows), len(cols)
        self.t_build = round(time.time() - t0, 1)
        if verbose: print(f'[SectorSystem d={d} chi={cD}] rows={self.nrows} cols={self.ncols} ({self.t_build}s)', file=sys.stderr, flush=True)
    def solve(self):
        t0 = time.time()
        R, C = self.nrows, self.ncols
        flat = [0] * (R * (C + 1))
        for j, col in enumerate(self.colpolys):
            for e, v in col.items(): flat[self.rows[e] * (C + 1) + j] = v
        for e, v in self.rhs.items(): flat[self.rows[e] * (C + 1) + C] = v
        M = nmod_mat(R, C + 1, flat, p)
        Mr, rank = M.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        rank0 = rank - (0 if consistent else 1)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
        return {'consistent': consistent, 'rank': rank0, 'nullity': C - rank0, 'rows': R, 'cols': C, 'x': x, 't_rref': round(time.time() - t0, 1)}
    def unpack(self, x):
        N = [{} for _ in range(NV)]; Q = {}
        for (kind, i, m), v in zip(self.cols, x):
            if not v: continue
            if kind == 'N': N[i][m] = v
            else: Q[m] = v
        return N, Q

def sector_e1(D):
    cD = chi_of(D); d = core.degree(D)
    S = sector(d + 6, chi_add(cD, CHI_A))
    return S.in_J(mul(Ap, core.reduce_mod(D, p), p))

def test_D(D, name, family, note=''):
    """E1 (sector) then full sector system; verify a consistent solution exactly mod p; log."""
    t0 = time.time()
    d = core.degree(D)
    cD = chi_of(D)
    if cD is None:
        r = log('fixedD', D=name, family=family, d=d, result='skipped-mixed-character', note=note)
        return r
    if not sector_e1(D):
        r = log('fixedD', D=name, family=family, d=d, chi=str(cD), e1=False, result='E1-fail', seconds=round(time.time() - t0, 1), note=note)
        return r
    S = SectorSystem(D, verbose=True)
    res = S.solve()
    extra = {}
    if res['consistent']:
        N, Q = S.unpack(res['x'])
        Dp = core.reduce_mod(D, p)
        E = core.defect(Fp, Ap, Dp, 1, N, p)
        extra['verified_mod_p'] = (sub(E, mul(power(Fp, 2, p), Q, p), p) == {})
    r = log('fixedD', D=name, family=family, d=d, chi=str(cD), e1=True, result='CONSISTENT' if res['consistent'] else 'inconsistent',
            rank=res['rank'], nullity=res['nullity'], rows=res['rows'], cols=res['cols'], seconds=round(time.time() - t0, 1), note=note, **extra)
    if res['consistent']:
        r['N'], r['Q'] = N, Q
        with open(PK + f'/CONSISTENT-{name.replace(" ", "_").replace("*", "s").replace("/", "_")}.json', 'w') as f:
            json.dump({'D': core.to_str(D), 'N': [core.to_str(x) for x in N], 'Q': core.to_str(Q), 'p': p}, f)
    return r

# building blocks (exact)
def X(i):
    return {tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)}
X0, X1, X2, Y0, Y1, Y2 = [X(i) for i in range(NV)]
ONE = {(0,) * NV: Fraction(1)}
S0 = add(mul(X0, X0), mul(Y0, Y0)); P0 = mul(X0, Y0); DELTA = sub(mul(X0, X0), mul(Y0, Y0))
S1 = add(mul(X1, X1), mul(Y1, Y1)); P1 = mul(X1, Y1); S2 = add(mul(X2, X2), mul(Y2, Y2)); P2 = mul(X2, Y2)
Q1 = add(add(mul(P1, P1), mul(P1, P2)), mul(P2, P2))
Q2 = mul(P1, add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
GS = add(add(scale(mul(S0, S0), 3), scale(mul(P0, P0), 2)), scale(Q1, 5))
GP = add(scale(mul(P0, S0), 4), scale(Q2, Fraction(5, 7)))
D8A = sub(scale(GS, 2), GP); D8B = add(scale(GS, 2), GP)
D8 = mul(D8A, D8B)
FX = [pderiv(F, i) for i in range(NV)]
