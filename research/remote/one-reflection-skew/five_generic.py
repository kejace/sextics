"""Five-variable reduction for F = f(X0,X1,X2,Y0,Y1) + Y2^6, q = 1.

A certificate (D,N,Q) for (F,A,1) with A free of Y2 exists iff the strict 5-variable equation
    A D^2 - D N(f) + (6/5) f (D div N - N(D)) = 0
has a solution with D = D(x) homogeneous of degree d, N_i of degree d+1 (5 components).
(Direction 5->6: (D, (N, Y2 div N/5), Q=0) is a strict 6-variable certificate.)
Polynomials are stored as 6-tuple dicts with Y2-exponent 0 so core helpers apply.
rho-parity = (e3+e4) mod 2 plays the role of Y-parity; D rho-even, A rho-odd,
N_0..N_2 rho-odd, N_3,N_4 rho-even.
"""
import sys, time, json, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from fractions import Fraction
from hodge import core
from hodge.core import mul, add, sub, scale, pderiv, shift, power
from flint import nmod_mat
import os, json
INST = os.environ.get('HODGE_INST', 'one-reflection-skew')
F, A, q, _rec = core.load_instance(INST)
def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    with open(f'/tmp/hodge/results/one-reflection-skew/control_{INST}.jsonl', 'a') as fh:
        fh.write(json.dumps({"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent, "rank": rank,
                             "nullity": nullity, "seconds": round(seconds, 2), "note": note, "instance": INST}) + '\n')
class _O: pass
orslib = _O(); orslib.save_witness = lambda D, N, Q: core.verify_witness(F, A, q, D, N, Q)[0]

P1, P2 = core.PRIMES[0], core.PRIMES[1]
f = {e: c for e, c in F.items() if e[5] == 0}          # f = F - Y2^6
assert len(f) == len(F) - 1 and all(sum(e) == 6 and e[5] == 6 for e in F if e[5])
A_plus = A
A_minus = A
X = [core.parse(v) for v in core.VARS]

def rpar(e):
    return (e[3] + e[4]) & 1

def spar(e):  # sigma: (X1,Y1) -> -(X1,Y1)
    return (e[1] + e[4]) & 1

def sparity_of(P):
    ps = {spar(e) for e in P}
    return ps.pop() if len(ps) == 1 else None

def rparity_of(P):
    ps = {rpar(e) for e in P}
    return ps.pop() if len(ps) == 1 else None

def mons5(deg, par=None):
    if deg < 0: return []
    out = []
    for c in itertools.combinations(range(deg + 4), 4):
        e = []; prev = -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + 3 - prev)
        e = tuple(e) + (0,)
        if par is None or rpar(e) == par:
            out.append(e)
    return out

def fi(i, p=None):
    return pderiv(core.reduce_mod(f, p) if p else f, i, p)

INV5 = {}
def inv(a, p):
    return pow(a, -1, p)

class Sys5:
    """Columns for N (5 components) and optionally for targets (rho-odd degree-6 monomials)."""
    def __init__(self, D, p, with_targets=True, verbose=False, use_sigma=True):
        t0 = time.time()
        self.p = p
        fp = core.reduce_mod(f, p); Dp = core.reduce_mod(D, p)
        d = core.degree(D); self.d = d
        pD = rparity_of(D)
        assert pD is not None, 'D has mixed rho-parity'
        pA = rparity_of(A); assert pA is not None
        pars = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(5)]
        self.pars = pars
        sD = sparity_of(D)
        self.use_sigma = bool(use_sigma and sD is not None)
        if self.use_sigma:
            # sigma-character of N_i: chi(A) + chi(D) + chi(f) + chi(Z_i); chi(A)=0 (X0^3Y0^3), chi(f)=0
            self.spars = [(0 + sD + 0 + (1 if i in (1, 4) else 0)) % 2 for i in range(5)]
            self.sT = 0  # targets: sigma-even only (A's sector)
        else:
            self.spars = None
        c65 = (6 * inv(5, p)) % p
        fD = scale(mul(fp, Dp, p), c65, p)                          # (6/5) f D
        G = []
        for i in range(5):
            Di = pderiv(Dp, i, p)
            G.append(add(mul(Dp, fi(i, p), p), scale(mul(fp, Di, p), c65, p), p))  # D f_i + (6/5) f D_i
        G = [scale(g, p - 1, p) for g in G]
        cols = []; colpolys = []
        for i in range(5):
            for m in mons5(d + 1, pars[i]):
                if self.use_sigma and spar(m) != self.spars[i]: continue
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(fD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        self.nN = len(cols)
        self.D2 = mul(Dp, Dp, p)
        self.tmons = mons5(6, 1) if with_targets else []
        if self.use_sigma: self.tmons = [m for m in self.tmons if spar(m) == 0]
        for m in self.tmons:
            cols.append(('A', -1, m)); colpolys.append(shift(self.D2, m, 1, p))
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows: rows[e] = len(rows)
        self.cols, self.colpolys, self.rows = cols, colpolys, rows
        self.nrows, self.ncols = len(rows), len(cols)
        self.nnz = sum(len(c) for c in colpolys)
        self.t_build = time.time() - t0
        if verbose:
            print(f'[sys5] d={d} rows={self.nrows} cols={self.ncols} (N {self.nN}, T {len(self.tmons)}) nnz={self.nnz} build {self.t_build:.1f}s', file=sys.stderr)

    def matrix(self, rhs=None):
        R, C, p = self.nrows, self.ncols, self.p
        W = C + (1 if rhs is not None else 0)
        flat = [0] * (R * W)
        rows = self.rows
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                flat[rows[e] * W + j] = v
        if rhs is not None:
            for e, v in rhs.items():
                if e in rows:
                    flat[rows[e] * W + C] = v
                else:
                    raise ValueError('rhs monomial outside row space (cannot happen for A D^2)')
        return nmod_mat(R, W, flat, p)

def target_vec(Aq, tmons, p):
    Ap = core.reduce_mod(Aq, p)
    idx = {m: k for k, m in enumerate(tmons)}
    v = [0] * len(tmons)
    for e, c in Ap.items():
        assert e in idx, f'target monomial {e} not rho-odd degree 6'
        v[idx[e]] = c
    return v

def jacobian6_odd(p, sigma_even=False):
    """Basis (mod p) of J_f in degree 6, rho-odd part (optionally sigma-even part): Z_j f_i."""
    out = []
    for i in range(5):
        for j in range(5):
            P = mul(core.reduce_mod(X[j], p), fi(i, p), p)
            if P and rparity_of(P) == 1 and (not sigma_even or sparity_of(P) == 0):
                out.append(P)
    return out

def analyze(D, p=P1, verbose=False, targets=(('A', A), ('A+', A_plus), ('A-', A_minus)), use_sigma=False):
    """Full served-target analysis of D mod p.  Returns dict.  With use_sigma and sigma-pure D only the
    sigma-even target sector (A's sector) is computed; A+/A- (sigma-mixed) are then tested via their
    sigma-even projection, which is A/2 -- so they are reported as served iff A is."""
    t0 = time.time()
    S = Sys5(D, p, with_targets=True, verbose=verbose, use_sigma=use_sigma)
    if S.use_sigma:
        targets = tuple((n, (A if n in ('A+', 'A-') else T)) for n, T in targets)
    M = S.matrix()
    t1 = time.time()
    Xn, nullity = M.nullspace()
    t2 = time.time()
    nT = len(S.tmons); C = S.ncols; nN = S.nN
    res = {'d': S.d, 'rows': S.nrows, 'cols': C, 'nN': nN, 'nT': nT, 'nullity_full': nullity,
           'p': p, 't_build': round(S.t_build, 1), 't_null': round(t2 - t1, 1)}
    res['sigma_sector'] = S.use_sigma
    if nullity == 0:
        res.update({'served_dim': 0, 'served_mod_J': 0, 'rank_N': nN, 'nullity_N': 0, 'served_basis': []})
        for name, _ in targets: res['serves_' + name] = False
        res['seconds'] = time.time() - t0
        return res
    proj = nmod_mat(nT, nullity, [int(Xn[nN + j, k]) for j in range(nT) for k in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for j in range(nT):
            v = int(Tr[i, j])
            if v: poly[S.tmons[j]] = v
        basis.append(poly)
    res['served_dim'] = rk
    res['rank_N'] = nN - (nullity - rk)
    res['nullity_N'] = nullity - rk
    # modulo J_f (degree 6, rho-odd part)
    Jb = jacobian6_odd(p, sigma_even=S.use_sigma)
    res['sigma_sector'] = S.use_sigma
    idx = {m: k for k, m in enumerate(S.tmons)}
    def vec(P):
        v = [0] * nT
        for e, c in P.items(): v[idx[e]] = c
        return v
    Jrows = [vec(P) for P in Jb]
    MJ = nmod_mat(len(Jrows), nT, [x for r in Jrows for x in r], p)
    rJ = MJ.rank()
    MJS = nmod_mat(len(Jrows) + rk, nT, [x for r in Jrows for x in r] + [x for b in basis for x in vec(b)], p)
    rJS = MJS.rank()
    res['dim_J6_odd'] = rJ
    res['served_mod_J'] = rJS - rJ
    for name, Aq in targets:
        va = target_vec(Aq, S.tmons, p)
        MA = nmod_mat(len(Jrows) + rk + 1, nT, [x for r in Jrows for x in r] + [x for b in basis for x in vec(b)] + va, p)
        res['serves_' + name] = (MA.rank() == rJS)
        # also strictly (not mod J): A in served space itself
        MA2 = nmod_mat(rk + 1, nT, [x for b in basis for x in vec(b)] + va, p)
        res['serves_strict_' + name] = (MA2.rank() == rk)
    res['served_basis'] = basis
    res['seconds'] = time.time() - t0
    return res

def solve_exact(D, Aq=A, primes=(P1, P2, core.PRIMES[2], core.PRIMES[3]), verbose=False):
    """Fixed-D strict 5-var solve for target Aq with rational lift; verifies as a 6-var certificate."""
    acc = None; piv0 = None
    for k, p in enumerate(primes):
        S = Sys5(D, p, with_targets=False, verbose=verbose, use_sigma=False)
        rhs = scale(mul(core.reduce_mod(Aq, p), S.D2, p), p - 1, p)
        M = S.matrix(rhs)
        Mr, rank = M.rref()
        C = S.ncols
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            pivots.append(j)
        if pivots and pivots[-1] == C:
            return 'inconsistent', {'rank': rank, 'rows': S.nrows, 'cols': C, 'p': p}
        x = [0] * C
        for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
        if piv0 is None or piv0 != pivots:
            acc = (x, p); piv0 = pivots
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(6)]
            for (kind, i, mon), v in zip(S.cols, rat):
                if v: N[i][mon] = v
            # 6-variable lift: u = (N/D, Y2*div(N/D)/5) has denominator D^2:
            # D6 = D^2, N6_i = N_i*D (i<5), N6_5 = Y2*(D div N - N(D))/5.  Strict (Q = 0).
            divN = {}; ND = {}
            for i in range(5):
                divN = add(divN, pderiv(N[i], i)); ND = add(ND, mul(N[i], pderiv(D, i)))
            BD = sub(mul(D, divN), ND)
            N6 = [mul(N[i], D) for i in range(5)] + [scale(mul(X[5], BD), Fraction(1, 5))]
            D6 = mul(D, D)
            ok, info = core.verify_witness(F, Aq, 1, D6, N6, {})
            if ok:
                return 'solved', {'N': N, 'N6': N6, 'D6': D6, 'Q': {}, 'rank': rank, 'nullity': C - rank, 'primes_used': k + 1}
            elif verbose:
                print('[lift] verified false, continuing', file=sys.stderr)
    return 'lift-failed', {}

def e1_kernel5(d, p=P1, par=0):
    """dim{D of rho-parity par, degree d : A D in J_f}, dim J_f,d (same parity), and a basis of kernel (mod p)."""
    t0 = time.time()
    Ap = core.reduce_mod(A, p)
    Dm = mons5(d, par)
    cols = [shift(Ap, m, 1, p) for m in Dm]
    # N_i monomials of degree d+1 with parity so that N_i f_i has parity of A D = 1 + par
    for i in range(5):
        pi = (1 + par - (0 if i < 3 else 1)) % 2
        for m in mons5(d + 1, pi):
            cols.append(shift(fi(i, p), m, 1, p))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    Xn, nullity = M.nullspace()
    nD = len(Dm)
    if nullity:
        proj = nmod_mat(nD, nullity, [int(Xn[j, k]) for j in range(nD) for k in range(nullity)], p)
        Tr, kerD = proj.transpose().rref()
        basis = []
        for i in range(kerD):
            poly = {}
            for j in range(nD):
                v = int(Tr[i, j])
                if v: poly[Dm[j]] = v
            basis.append(poly)
    else:
        kerD, basis = 0, []
    # dim J_d with parity par
    Jc = []
    if d >= 5:
        for i in range(5):
            pi = (par - (0 if i < 3 else 1)) % 2
            for m in mons5(d - 5, pi):
                Jc.append(shift(fi(i, p), m, 1, p))
    if Jc:
        ridx = {m: k for k, m in enumerate(Dm)}
        flatJ = [0] * (nD * len(Jc))
        for j, col in enumerate(Jc):
            for e, v in col.items(): flatJ[ridx[e] * len(Jc) + j] = v
        dimJ = nmod_mat(nD, len(Jc), flatJ, p).rank()
    else:
        dimJ = 0
    return {'d': d, 'par': par, 'nD': nD, 'kerD': kerD, 'dimJ': dimJ, 'ker_mod_J': kerD - dimJ,
            'rows': R, 'cols': C, 'seconds': round(time.time() - t0, 1)}, basis

def test5(D_expr, family, note='', Dpoly=None, verbose=False, do_exact=True):
    """Analyze D (served targets, A/A+/A- membership), log, and lift exactly when served."""
    D = Dpoly if Dpoly is not None else core.parse(D_expr)
    assert all(e[5] == 0 for e in D), 'D must be free of Y2'
    r = analyze(D, verbose=verbose)
    d = r['d']
    served_any = r['serves_A'] or r['serves_A+'] or r['serves_A-']
    status = 'inconsistent'
    if r['serves_A'] and do_exact:
        st, data = solve_exact(D, A, verbose=verbose)
        status = st
        if st == 'solved':
            print('!!!!!!!! WITNESS FOUND (5-var) !!!!!!!!', D_expr, file=sys.stderr)
            orslib.save_witness(data['D6'], data['N6'], data['Q'])
            r['N'] = data['N']
    elif r['serves_A']:
        status = 'consistent-modp'
    if (r['serves_A+'] or r['serves_A-']) and not r['serves_A']:
        status = 'half-served'
        print('!!!! HALF TARGET SERVED (mod p) !!!!', D_expr, r['serves_A+'], r['serves_A-'], file=sys.stderr)
    note2 = (note + f' | sigma_sector={r["sigma_sector"]} 5var served_dim={r["served_dim"]} served_mod_J={r["served_mod_J"]} '
             f'A={r["serves_A"]} A+={r["serves_A+"]} A-={r["serves_A-"]} rows={r["rows"]} cols={r["cols"]} status={status}')
    log(D_expr, d, '5var:' + family, True, bool(r['serves_A']), r['rank_N'], r['nullity_N'], r['seconds'], note2)
    r['status'] = status; r['D_expr'] = D_expr; r['family'] = family
    print(f'{D_expr[:70]:70s} d={d} served={r["served_dim"]}/{r["served_mod_J"]} A={r["serves_A"]} A+={r["serves_A+"]} A-={r["serves_A-"]} '
          f'rank={r["rank_N"]} null={r["nullity_N"]} {r["seconds"]:.1f}s', flush=True)
    return r
