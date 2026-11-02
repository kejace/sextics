"""Independent exact/modular machinery for the off-paired exclusion lens (wave 2).

Independent of wave-1 code (bilinear.py/oplib.py): own character bookkeeping, own fixed-D matrix builder
(sparse, any prime, or exact over Q with Fractions), rational inconsistency certificates (left-kernel vectors y
with y^T M = 0, y^T b = 1, reconstructed by CRT over several primes and verified EXACTLY over Q), and own
bilinear (D,N) sector systems exported to msolve / Singular (mod p and over QQ).
Engine functions used only for polynomial arithmetic (hodge.core) -- which is itself re-checked by core.defect.
"""
import sys, os, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, power
from flint import nmod_mat, fmpq_mat, fmpq, fmpz_mat

IID = 'off-paired'
LANE = '/tmp/hodge/results/off-paired'
LOG = LANE + '/wave2-excl.log.jsonl'
WG = (0, 2, 0, 1, 3, 3); WH = (0, 0, 1, 0, 0, 5)
PRIMES = [2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483497, 2147483489, 2147483477,
          2147483423, 2147483399, 2147483353, 2147483323, 2147483269, 2147483249, 2147483237, 2147483179]
# NOTE: wave-1 used 1000003 (fixed-D) and 1073741827 (msolve); dossier used 1000003 / 1009 / 101.  All of ours differ.

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def char(e):
    return (sum(a*b for a, b in zip(WG, e)) % 6, sum(a*b for a, b in zip(WH, e)) % 6)
def char_of(P):
    cs = {char(e) for e in P}
    return cs.pop() if len(cs) == 1 else None
def cadd(*cs):
    a = b = 0
    for c in cs: a += c[0]; b += c[1]
    return (a % 6, b % 6)
VARCHAR = [char(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
def monos_char(deg, c):
    return [m for m in monomials(deg) if char(m) == c] if deg >= 0 else []

def log(rec):
    rec = dict(rec); rec.setdefault('engine', 2); rec.setdefault('lens', 'excl'); rec['time'] = time.strftime('%Y-%m-%d %H:%M:%S')
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
    return rec

# ---------------------------------------------------------------- fixed-D system (own builder)
def fixedD_columns(F, A, q, D, p=None, sector=True):
    """Columns of the linear map (N,Q) -> -q D N(F) + F (D div N - N(D)) - F^{q+1} Q  and rhs = -A D^2,
    i.e. the identity  A D^2 + sum cols = 0  <=>  M x = b with b = -A D^2.
    p=None: exact over Q (Fraction coefficients).  sector=True restricts N_i, Q to the characters forced by D
    (lossless for eigenvector D).  Returns (cols:list of (kind,i,m), colpolys, rhs)."""
    red = (lambda P: core.reduce_mod(P, p)) if p else (lambda P: {e: Fraction(v) for e, v in P.items()})
    Fp, Ap, Dp = red(F), red(A), red(D)
    d = core.degree(D); degN = d + 6*q - 5; degQ = 2*d - 6
    cD, cA = char_of(Dp), char_of(Ap)
    if sector: assert cD is not None, 'D must be an eigenvector for sector restriction'
    Fi = [pderiv(Fp, i, p) for i in range(NV)]; Di = [pderiv(Dp, i, p) for i in range(NV)]
    FD = mul(Fp, Dp, p)
    neg = (lambda P: scale(P, p-1, p)) if p else (lambda P: scale(P, -1))
    G = [neg(add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p)) for i in range(NV)]   # -(q D F_i + F D_i)
    Fq1 = neg(power(Fp, q+1, p))
    parA = core.parity_of(A)
    cols = []; colpolys = []
    for i in range(NV):
        if sector: ms = monos_char(degN, cadd(cD, cA, VARCHAR[i]))
        else: ms = monomials(degN, (parA + core.parity_of(D) + (1 if i >= 3 else 0)) % 2)
        for m in ms:
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); colpolys.append(col)
    if degQ >= 0:
        ms = monos_char(degQ, cadd(cA, cD, cD)) if sector else monomials(degQ, parA % 2)
        for m in ms:
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
    rhs = neg(mul(Ap, mul(Dp, Dp, p), p))
    return cols, colpolys, rhs

def to_sparse(colpolys, rhs):
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    return rows

def dense_mod(colpolys, rhs, rows, p):
    R, C = len(rows), len(colpolys)
    flat = [0]*(R*(C+1))
    for j, col in enumerate(colpolys):
        for e, v in col.items(): flat[rows[e]*(C+1)+j] = int(v) % p
    for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = int(v) % p
    return nmod_mat(R, C+1, flat, p)

def ranks_mod(colpolys, rhs, p):
    rows = to_sparse(colpolys, rhs)
    Mb = dense_mod(colpolys, rhs, rows, p)
    R, C = Mb.nrows(), Mb.ncols()-1
    M = nmod_mat(R, C, [int(Mb[i, j]) for i in range(R) for j in range(C)], p)
    return M.rank(), Mb.rank(), R, C

# ---------------------------------------------------------------- rational inconsistency certificate
def _lin_indep_rows(M):
    """indices of a maximal set of linearly independent rows of nmod_mat M (pivot columns of rref(M^T))."""
    T, rk = M.transpose().rref()
    piv = []; r = 0
    for j in range(T.ncols()):
        if r < rk and int(T[r, j]) != 0:
            piv.append(j); r += 1
    return piv, rk

def inconsistency_certificate(colpolys, rhs, primes=PRIMES, max_primes=12, verbose=False):
    """Exact proof that  M x = b  has no solution over Q (hence over C): a rational vector y with
    y^T M = 0 and y^T b = 1.  Construction (mod p1): r = rank M; S = r independent rows of M; extra row i with
    rank [M|b]_{S+i} = r+1; y = the unique left-kernel vector of M_{S+i} normalized by y.b = 1.  The same canonical
    y is computed mod further primes, CRT-combined, rationally reconstructed and verified EXACTLY over Q.
    Returns dict(ok, y (sparse {row_monomial: Fraction}), r, rows_used, primes_used) or ok=False."""
    rows = to_sparse(colpolys, rhs); R = len(rows); C = len(colpolys)
    inv = {i: e for e, i in rows.items()}
    p1 = primes[0]
    Mb = dense_mod(colpolys, rhs, rows, p1)
    M = nmod_mat(R, C, [int(Mb[i, j]) for i in range(R) for j in range(C)], p1)
    S, r = _lin_indep_rows(M)
    rkb = Mb.rank()
    if rkb != r + 1:
        return {'ok': False, 'reason': f'rank M={r}, rank [M|b]={rkb} mod {p1}: consistent mod p1', 'r': r, 'rank_aug': rkb}
    # find extra row
    sub = nmod_mat(r, C+1, [int(Mb[i, j]) for i in S for j in range(C+1)], p1)
    extra = None
    for i in range(R):
        if i in set(S): continue
        test = nmod_mat(r+1, C+1, [int(Mb[k, j]) for k in S + [i] for j in range(C+1)], p1)
        if test.rank() == r + 1:
            extra = i; break
    assert extra is not None
    used = S + [extra]
    def y_mod(p):
        Mb_p = dense_mod(colpolys, rhs, rows, p)
        Msub = nmod_mat(r+1, C, [int(Mb_p[k, j]) for k in used for j in range(C)], p)
        K, k = Msub.transpose().nullspace()      # left kernel of Msub: vectors of length r+1
        if k != 1: return None
        y = [int(K[i, 0]) for i in range(r+1)]
        yb = sum(y[t] * int(Mb_p[used[t], C]) for t in range(r+1)) % p
        if yb == 0: return None
        s = pow(yb, -1, p)
        return [(v * s) % p for v in y]
    # CRT over primes
    acc = None; mod = 1; used_primes = []
    for p in primes[:max_primes]:
        yp = y_mod(p)
        if yp is None:
            if verbose: print('bad prime', p); 
            continue
        if acc is None: acc, mod = yp, p
        else:
            acc = [core.crt_pair(a, mod, b, p) for a, b in zip(acc, yp)]; mod *= p
        used_primes.append(p)
        rec = [core.ratrecon(a, mod) for a in acc]
        if all(v is not None for v in rec):
            # exact verification over Q
            y = {used[t]: rec[t] for t in range(r+1) if rec[t] != 0}
            ok = True
            for j, col in enumerate(colpolys):
                s = Fraction(0)
                for e, v in col.items():
                    i = rows[e]
                    if i in y: s += y[i] * Fraction(v)
                if s != 0: ok = False; break
            if ok:
                yb = sum(y[rows[e]] * Fraction(v) for e, v in rhs.items() if rows[e] in y)
                if yb == 1:
                    return {'ok': True, 'y': {core.to_str({inv[i]: 1}): str(v) for i, v in y.items()}, 'r': r, 'rank_aug': rkb,
                            'rows': R, 'cols': C, 'nrows_used': r+1, 'primes_used': used_primes, 'max_height': max(max(abs(v.numerator), v.denominator) for v in y.values())}
    return {'ok': False, 'reason': 'rational reconstruction failed within prime budget', 'r': r, 'rank_aug': rkb, 'primes_used': used_primes}

def verify_certificate_exact(colpolys, rhs, ycert):
    """Re-verify a certificate dict (monomial-string -> Fraction-string) exactly over Q from scratch."""
    y = {tuple(next(iter(core.parse(k)))): Fraction(v) for k, v in ycert.items()}
    for col in colpolys:
        if sum(y[e] * Fraction(v) for e, v in col.items() if e in y) != 0: return False
    return sum(y[e] * Fraction(v) for e, v in rhs.items() if e in y) == 1

# ---------------------------------------------------------------- exact rank over Q (small systems)
def exact_ranks_Q(colpolys, rhs):
    rows = to_sparse(colpolys, rhs); R = len(rows); C = len(colpolys)
    M = fmpq_mat(R, C+1)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): M[rows[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    for e, v in rhs.items(): M[rows[e], C] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    rkb = M.rank()
    M0 = fmpq_mat(R, C)
    for j, col in enumerate(colpolys):
        for e, v in col.items(): M0[rows[e], j] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
    return M0.rank(), rkb, R, C

# ---------------------------------------------------------------- bilinear (D,N) sector system, exact over Q or mod p
def colon_sector_and_lifts(F, A, q, d, c, p=None):
    """Basis of (J:A)_{d,c} modulo the Euler gauge F*S_{d-6,c}, with first-equation lifts N (q N(F) = A D),
    and the Koszul-type kernel {N : N(F) = 0} in the matching characters.  Exact over Q when p is None.
    Returns (chosen [(D,N)], syz [N], info)."""
    red = (lambda P: core.reduce_mod(P, p)) if p else (lambda P: {e: Fraction(v) for e, v in P.items()})
    Fp, Ap = red(F), red(A)
    cA = char_of(Ap); degN = d + 6*q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    neg = (lambda P: scale(P, p-1, p)) if p else (lambda P: scale(P, -1))
    kinds = []; cols = []
    Dmons = monos_char(d, c)
    for m in Dmons:
        cols.append(neg(shift(Ap, m, 1, p))); kinds.append(('D', -1, m))
    for i in range(NV):
        for m in monos_char(degN, cadd(c, VARCHAR[i], cA)):
            cols.append(shift(Fi[i], m, q, p)); kinds.append(('N', i, m))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    if p:
        flat = [0]*(R*C)
        for j, col in enumerate(cols):
            for e, v in col.items(): flat[rows[e]*C+j] = int(v)
        M = nmod_mat(R, C, flat, p)
        X, n = M.nullspace()
        T = nmod_mat(n, C, [int(X[j, k]) for k in range(n) for j in range(C)], p)
        Tr, rk = T.rref()
        getv = lambda i, j: int(Tr[i, j])
    else:
        M = fmpq_mat(R, C)
        for j, col in enumerate(cols):
            for e, v in col.items(): M[rows[e], j] = fmpq(v.numerator, v.denominator)
        # nullspace over Q: rref of M, free columns
        Rr, rk0 = M.rref()
        pivs = []; r = 0
        for j in range(C):
            if r < rk0 and Rr[r, j] != 0: pivs.append(j); r += 1
        free = [j for j in range(C) if j not in set(pivs)]
        n = len(free)
        # kernel vectors: for each free f: x_f = 1, x_piv[i] = -Rr[i,f]
        T = fmpq_mat(n, C)
        for k, f in enumerate(free):
            T[k, f] = 1
            for i, pj in enumerate(pivs): T[k, pj] = -Rr[i, f]
        Tr, rk = T.rref()
        getv = lambda i, j: Fraction(int(Tr[i, j].p), int(Tr[i, j].q))
    nD = len(Dmons)
    part = []; syz = []
    for i in range(rk):
        D = {}; N = [{} for _ in range(NV)]
        for j in range(C):
            v = getv(i, j)
            if not v: continue
            kind, comp, m = kinds[j]
            if kind == 'D': D[m] = v
            else: N[comp][m] = v
        if D: part.append((D, N))
        else: syz.append(N)
    # gauge complement
    gauge = [mul(Fp, {m: (1 if p else Fraction(1))}, p) for m in monos_char(d-6, c)]
    idx = {m: i for i, m in enumerate(Dmons)}
    def vec(D):
        v = [0]*nD
        for e, val in D.items(): v[idx[e]] = val
        return v
    def rank_of(rowsl):
        if not rowsl: return 0
        if p: return nmod_mat(len(rowsl), nD, [int(x) for r in rowsl for x in r], p).rank()
        Mq = fmpq_mat(len(rowsl), nD)
        for i, r in enumerate(rowsl):
            for j, x in enumerate(r):
                if x: Mq[i, j] = fmpq(Fraction(x).numerator, Fraction(x).denominator)
        return Mq.rank()
    base = [vec(g) for g in gauge]
    cur = rank_of(base)
    chosen = []
    for (D, N) in part:
        if rank_of(base + [vec(D_) for D_, _ in chosen] + [vec(D)]) > cur + len(chosen):
            chosen.append((D, N))
    return chosen, syz, {'colon_dim': len(part), 'gauge_dim': cur, 'syz_dim': len(syz), 'nD': nD, 'e1_rows': R, 'e1_cols': C}

def B_poly(D, N, p=None):
    divN = {}; ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

def cokernel_F(F, deg, sector, p=None):
    """Functionals on S_{deg,sector} vanishing on F*S_{deg-6,sector} (left nullspace of multiplication by F)."""
    red = (lambda P: core.reduce_mod(P, p)) if p else (lambda P: {e: Fraction(v) for e, v in P.items()})
    Fp = red(F)
    tgt = monos_char(deg, sector); src = monos_char(deg-6, sector)
    idx = {m: i for i, m in enumerate(tgt)}
    R, C = len(tgt), len(src)
    if C == 0:
        return idx, None, R
    if p:
        flat = [0]*(R*C)
        for j, m in enumerate(src):
            for e, v in shift(Fp, m, 1, p).items(): flat[idx[e]*C+j] = int(v)
        M = nmod_mat(R, C, flat, p)
        L, k = M.transpose().nullspace()
        Lrows = [[int(L[i, j]) for i in range(R)] for j in range(k)]
    else:
        M = fmpq_mat(C, R)   # transpose directly
        for j, m in enumerate(src):
            for e, v in shift(Fp, m, 1).items(): M[j, idx[e]] = fmpq(v.numerator, v.denominator)
        Rr, rk0 = M.rref()
        pivs = []; r = 0
        for j in range(R):
            if r < rk0 and Rr[r, j] != 0: pivs.append(j); r += 1
        free = [j for j in range(R) if j not in set(pivs)]
        Lrows = []
        for f in free:
            v = [Fraction(0)]*R; v[f] = Fraction(1)
            for i, pj in enumerate(pivs): v[pj] = -Fraction(int(Rr[i, f].p), int(Rr[i, f].q))
            Lrows.append(v)
        k = len(Lrows)
    return idx, Lrows, k

def build_bilinear(F, A, q, d, c, p=None, verbose=True):
    """Equations (quadratic forms in (t,s)) of the sector bilinear system; exact over Q if p None."""
    t0 = time.time()
    chosen, syz, info = colon_sector_and_lifts(F, A, q, d, c, p)
    cA = char_of(core.reduce_mod(A, p) if p else A)
    sector = cadd(cA, c, c)
    nt, ns = len(chosen), len(syz)
    if nt == 0: return None, info
    idx, Lrows, k = cokernel_F(F, 2*d + 6*q - 6 + 0, sector, p) if False else cokernel_F(F, 2*d + 6*q - 6, sector, p)
    # degree of B_D(N) = d + (d+6q-5) - 1 = 2d + 6q - 6  (q=1: 2d)
    Ds = [D for D, N in chosen]; Ns = [N for D, N in chosen]
    R = len(idx)
    keys = []; vecs = []
    for a in range(nt):
        for b in range(nt):
            P = B_poly(Ds[a], Ns[b], p)
            v = [0]*R
            for e, val in P.items(): v[idx[e]] = val
            keys.append((('t', a), ('t', b)) if a <= b else (('t', b), ('t', a))); vecs.append(v)
        for j in range(ns):
            P = B_poly(Ds[a], syz[j], p)
            v = [0]*R
            for e, val in P.items(): v[idx[e]] = val
            keys.append((('t', a), ('s', j))); vecs.append(v)
    eqs = [{} for _ in range(k)]
    if p:
        Mv = nmod_mat(len(vecs), R, [int(x) for v in vecs for x in v], p)
        if Lrows is None: Pm = Mv
        else:
            Lt = nmod_mat(R, k, [Lrows[j][i] for i in range(R) for j in range(k)], p)
            Pm = Mv * Lt
        for rowi, key in enumerate(keys):
            for r in range(k):
                val = int(Pm[rowi, r])
                if val: eqs[r][key] = (eqs[r].get(key, 0) + val) % p
    else:
        Mv = fmpq_mat(len(vecs), R)
        for i, v in enumerate(vecs):
            for j, x in enumerate(v):
                if x: Mv[i, j] = fmpq(x.numerator, x.denominator)
        if Lrows is None: Pm = Mv
        else:
            Lt = fmpq_mat(R, k)
            for j in range(k):
                for i in range(R):
                    x = Lrows[j][i]
                    if x: Lt[i, j] = fmpq(x.numerator, x.denominator)
            Pm = Mv * Lt
        for rowi, key in enumerate(keys):
            for r in range(k):
                val = Pm[rowi, r]
                if val != 0:
                    eqs[r][key] = eqs[r].get(key, Fraction(0)) + Fraction(int(val.p), int(val.q))
    eqs = [e for e in eqs if e]
    info.update({'nt': nt, 'ns': ns, 'neq': len(eqs), 'neq_raw': k, 'build_s': round(time.time()-t0, 1), 'p': p})
    if verbose: print(f'[bilinear d={d} c={c} p={p}] {info}', file=sys.stderr)
    return {'eqs': eqs, 'nt': nt, 'ns': ns, 'Ds': Ds, 'Ns': Ns, 'syz': syz, 'p': p, 'd': d, 'c': c}, info

def vname(v): return f't{v[1]}' if v[0] == 't' else f's{v[1]}'

def write_msolve(sysd, path, p, seed=1):
    """msolve input mod p (system must be mod p or exact; exact coefficients are reduced mod p)."""
    nt, ns = sysd['nt'], sysd['ns']
    names = [f't{i}' for i in range(nt)] + [f's{j}' for j in range(ns)]
    def cf(x):
        if isinstance(x, Fraction): return (x.numerator * pow(x.denominator, -1, p)) % p
        return int(x) % p
    polys = []
    for eq in sysd['eqs']:
        terms = [f'{cf(v)}*{vname(a)}*{vname(b)}' for (a, b), v in eq.items() if cf(v)]
        if terms: polys.append('+'.join(terms))
    rng = random.Random(seed)
    polys.append('+'.join(f'{rng.randrange(1, p)}*t{i}' for i in range(nt)) + f'+{p-1}')
    with open(path, 'w') as f: f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
    return len(polys)

def write_singular(sysd, path, p=0, seed=1, extra_cmds=''):
    """Singular script: ring over QQ (p=0) or ZZ/p; prints 'UNIT' if std == 1 else 'NONUNIT dim=...'."""
    nt, ns = sysd['nt'], sysd['ns']
    names = [f't{i}' for i in range(nt)] + [f's{j}' for j in range(ns)]
    def cf(x):
        if p:
            if isinstance(x, Fraction): return str((x.numerator * pow(x.denominator, -1, p)) % p)
            return str(int(x) % p)
        x = Fraction(x); return f'({x.numerator}/{x.denominator})' if x.denominator != 1 else str(x.numerator)
    polys = []
    for eq in sysd['eqs']:
        terms = [f'{cf(v)}*{vname(a)}*{vname(b)}' for (a, b), v in eq.items() if Fraction(v) != 0]
        if terms: polys.append('+'.join(terms))
    rng = random.Random(seed)
    polys.append('+'.join(f'{rng.randrange(1, 10**6)}*t{i}' for i in range(nt)) + '-1')
    with open(path, 'w') as f:
        f.write(f'ring r = {p},({",".join(names)}),dp;\noption(redSB);\nideal I = {",\n".join(polys)};\n')
        f.write(extra_cmds)
        f.write('int t0 = timer;\nideal G = std(I);\nif (size(G)==1 and deg(G[1])==0) { print("UNIT"); } else { print("NONUNIT dim=" + string(dim(G)) + " size=" + string(size(G))); }\nprint("seconds=" + string(timer-t0));\nquit;\n')
    return len(polys)

def run_msolve(path, out, threads=1, timeout=3600):
    import subprocess
    t0 = time.time()
    try:
        subprocess.run(['msolve', '-f', path, '-o', out, '-t', str(threads), '-v', '0'], capture_output=True, text=True, timeout=timeout)
        res = open(out).read().strip() if os.path.exists(out) else ''
    except subprocess.TimeoutExpired:
        return 'timeout', time.time()-t0, ''
    dt = time.time()-t0
    if res.startswith('[-1]'): return 'no-solution', dt, res
    if res.startswith('[1,'): return 'positive-dim', dt, res[:2000]
    return 'finite-or-other', dt, res[:2000]
