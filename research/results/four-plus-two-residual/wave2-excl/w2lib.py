"""wave-2 EXCLUSION lens library for four-plus-two-residual (q=1, A = X0^3 Y0^3).

Everything here is exact linear algebra (python-flint) either mod p (nmod_mat) or over Q (fmpz_mat, fraction-free).
Conventions (engine): variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5; polys are dicts {exponent tuple: coeff}.
Core = indices (0,1,3,4): H(X0,X1,Y0,Y1); binary block = (2,5): B(X2,Y2); F = H + B; a = A = X0^3 Y0^3 (core).

Rigor conventions used by every job (see README in the packet):
 * [MODP]  a statement computed modulo one prime; evidence only.
 * [RANK]  modular-rank argument: for an integer matrix, rank mod p <= rank over Q; combined with an EXACT
           upper bound on the rational rank (an explicit rational solution space of the right dimension) it
           proves rational (hence complex) inconsistency.  Used for fixed-D exclusions (koszul_nullity_exact).
 * [EXACT] computed over Q with fmpz_mat (fraction-free); field-independent (ranks do not change under Q -> C).
"""
import sys, json, time, os, random, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, symmetry
from hodge.core import monomials, shift, mul, add, sub, scale, pderiv, reduce_mod, power, NV
from flint import nmod_mat, fmpz_mat, fmpq_mat
from fractions import Fraction

IID = 'four-plus-two-residual'
LANE = f'/tmp/hodge/results/{IID}'
PACKET = f'{LANE}/wave2-excl'
REMOTE_OUT = f'{PACKET}/remote'
F, A, q, rec = core.load_instance(IID)
assert q == 1
CORE = (0, 1, 3, 4)
BIN = (2, 5)
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
B = {e: c for e, c in F.items() if e[2] + e[5] == 6}
assert add(H, B) == F
P1 = 2147483647          # wave-1 prime
P2 = 1000000007
P3 = 998244353
P4 = 2147483629
PRIMES_NEW = [P2, P3, P4]


def log(job, ent):
    """Append one JSON line to the job's own output file (never collides with other jobs / rsync)."""
    os.makedirs(REMOTE_OUT, exist_ok=True)
    ent = dict(ent); ent.setdefault('job', job); ent.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(f'{REMOTE_OUT}/{job}.jsonl', 'a') as f:
        f.write(json.dumps(ent) + '\n')
    print(json.dumps(ent), flush=True)


def save_json(job, name, obj):
    os.makedirs(REMOTE_OUT, exist_ok=True)
    with open(f'{REMOTE_OUT}/{job}-{name}.json', 'w') as f:
        json.dump(obj, f)


# ----------------------------------------------------------------------------------------------
# monomials restricted to a variable subset
# ----------------------------------------------------------------------------------------------
def mons(deg, idx, par=None):
    """Exponent tuples of total degree deg supported on variable indices idx, optional Y-parity."""
    idx = tuple(idx)
    out = []
    for m in monomials(deg, par):
        if all(m[j] == 0 for j in range(NV) if j not in idx):
            out.append(m)
    return out


def ci_hilbert(n, e=5):
    """Hilbert function of a complete intersection of n forms of degree e in n variables: coefficients of (1+t+..+t^{e-1})^n."""
    poly = [1]
    for _ in range(n):
        new = [0] * (len(poly) + e - 1)
        for i, c in enumerate(poly):
            for j in range(e):
                new[i + j] += c
        poly = new
    return poly


# ----------------------------------------------------------------------------------------------
# matrices from column polynomials
# ----------------------------------------------------------------------------------------------
def build_modp(cols, p, rows=None, extra=None):
    """nmod_mat whose columns are the polys in cols (mod p); rows indexed by exponent tuples."""
    rows = {} if rows is None else dict(rows)
    allc = list(cols) + (list(extra) if extra else [])
    for col in allc:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(allc)
    flat = [0] * (R * C)
    for j, col in enumerate(allc):
        for e, v in col.items():
            flat[rows[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p), rows


def build_exact(cols, rows=None, extra=None):
    """fmpz_mat whose columns are the polys in cols scaled to integers (per-column lcm of denominators);
    column scaling preserves rank and the nullspace up to per-coordinate scaling.  Returns (M, rows, scales)."""
    rows = {} if rows is None else dict(rows)
    allc = list(cols) + (list(extra) if extra else [])
    for col in allc:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(allc)
    flat = [0] * (R * C)
    scales = []
    for j, col in enumerate(allc):
        den = 1
        for v in col.values():
            v = Fraction(v); den = den * v.denominator // math.gcd(den, v.denominator)
        scales.append(den)
        for e, v in col.items():
            v = Fraction(v) * den
            assert v.denominator == 1
            flat[rows[e] * C + j] = int(v)
    return fmpz_mat(R, C, flat), rows, scales


def rank_modp(cols, p, extra=None):
    M, _ = build_modp(cols, p, extra=extra)
    return int(M.rank())


def rank_exact(cols, extra=None):
    M, _, _ = build_exact(cols, extra=extra)
    return int(M.rank())


# ----------------------------------------------------------------------------------------------
# Jacobian rings: Hilbert functions (regular-sequence check) and multiplication kernels
# ----------------------------------------------------------------------------------------------
def jacobian_dims(P, idx, kmax, p=None, par=None):
    """dim J(P)_k and dim R(P)_k for k = 0..kmax on variables idx (mod p, or exact over Q if p is None).
    If par is given (0/1) only the Y-parity-par part is counted."""
    Pi = [pderiv(P, i) for i in idx]
    out = []
    for k in range(kmax + 1):
        ms = mons(k, idx, par)
        if k < 5:
            out.append((len(ms), 0, len(ms))); continue
        cols = []
        for j, i in enumerate(idx):
            pj = 0 if i < 3 else 1          # parity of P_i for Y-even P
            for m in mons(k - 5, idx, None if par is None else (par - pj) % 2):
                cols.append(shift(Pi[j], m, 1))
        r = rank_exact(cols) if p is None else rank_modp([reduce_mod(c, p) for c in cols], p)
        out.append((len(ms), r, len(ms) - r))
    return out


def regular_sequence_check(P, idx, p=None):
    """Return (ok, detail): ok iff the partials of P on idx form a regular sequence (mod p or over Q), i.e.
    the Jacobian ring is Artinian with the complete-intersection Hilbert function (1+t+..+t^4)^n."""
    n = len(idx)
    ci = ci_hilbert(n)
    socle = 4 * n
    dims = jacobian_dims(P, idx, socle + 1, p)
    hf = [d[2] for d in dims]
    ok = hf[:socle + 1] == ci and hf[socle + 1] == 0
    return ok, {'hilbert': hf, 'ci': ci}


def mult_kernel(P, idx, a, k, par, p=None):
    """Kernel of multiplication by a on R(P)_k (variables idx, Y-parity par):
       ker_full = {D in R_k^par : a*D in J(P)_{k+deg a}};  returns dict with
       dim_full, dim_J (= dim J(P)_k^par), dim_modJ, basis_full (list of polys, exact ints if p None else mod p),
       reps (basis elements independent modulo J_k).
    Exact (p None): integer kernel vectors from fmpz_mat.nullspace (fraction-free)."""
    da = core.degree(a)
    Pi = [pderiv(P, i) for i in idx]
    msD = mons(k, idx, par)
    tpar = None if par is None else (par + core.parity_of(a)) % 2
    Acols = [mul(a, {m: 1}) for m in msD]
    Jcols = []
    for j, i in enumerate(idx):
        pj = 0 if i < 3 else 1
        for m in mons(k + da - 5, idx, None if tpar is None else (tpar - pj) % 2):
            Jcols.append(shift(Pi[j], m, 1))
    nA = len(msD)
    if p is None:
        M, rows, scales = build_exact(Acols + Jcols)
        X, nul = M.nullspace()
        vecs = []
        for kk in range(nul):
            v = [int(X[j, kk]) for j in range(nA)]
            if any(v): vecs.append(v)
    else:
        M, rows = build_modp([reduce_mod(c, p) for c in Acols + Jcols], p)
        X, nul = M.nullspace()
        vecs = []
        for kk in range(nul):
            v = [int(X[j, kk]) for j in range(nA)]
            if any(v): vecs.append(v)
    # the projection of the nullspace onto the D-coordinates is ker_full; row-reduce it
    if not vecs:
        basis = []
    elif p is None:
        Mv = fmpz_mat(len(vecs), nA, [x for v in vecs for x in v])
        R, den, rk = Mv.rref()
        basis = []
        for i in range(rk):
            row = [int(R[i, j]) for j in range(nA)]
            g = 0
            for x in row: g = math.gcd(g, abs(x))
            basis.append({msD[j]: Fraction(row[j] // g) for j in range(nA) if row[j]})
    else:
        Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
        R, rk = Mv.rref()
        basis = [{msD[j]: int(R[i, j]) for j in range(nA) if int(R[i, j])} for i in range(rk)]
    # J_k^par
    Jk = []
    if k >= 5:
        for j, i in enumerate(idx):
            pj = 0 if i < 3 else 1
            for m in mons(k - 5, idx, None if par is None else (par - pj) % 2):
                Jk.append(shift(Pi[j], m, 1))
    if p is None:
        rJ = rank_exact(Jk) if Jk else 0
    else:
        Jk = [reduce_mod(c, p) for c in Jk]
        rJ = rank_modp(Jk, p) if Jk else 0
    reps = []
    cur = list(Jk); currank = rJ
    for b in basis:
        r = rank_exact(cur + [b]) if p is None else rank_modp(cur + [b], p)
        if r > currank:
            reps.append(b); cur.append(b); currank += 1
    return {'k': k, 'par': par, 'p': p, 'n_mons': nA, 'dim_full': len(basis), 'dim_J': rJ,
            'dim_modJ': len(basis) - rJ, 'basis_full': basis, 'reps': reps}


def poly_to_json(P):
    return {' '.join(map(str, e)): str(c) for e, c in P.items()}


def poly_from_json(d):
    out = {}
    for k, v in d.items():
        e = tuple(int(x) for x in k.split())
        out[e] = Fraction(v) if '/' in v or not v.lstrip('-').isdigit() or True else int(v)
    return out


def small_height_reps(reps, msD=None):
    """LLL-reduce an integer basis (list of polys with integer coeffs) to get smaller coefficients; same span over Q."""
    if not reps: return reps
    keys = sorted({e for r in reps for e in r})
    M = fmpz_mat(len(reps), len(keys), [int(Fraction(r.get(e, 0))) for r in reps for e in keys])
    L = M.lll()
    out = []
    for i in range(L.nrows()):
        row = {keys[j]: Fraction(int(L[i, j])) for j in range(len(keys)) if int(L[i, j])}
        if row: out.append(row)
    return out


# ----------------------------------------------------------------------------------------------
# fixed-D certificate systems (general builder: optional parity-free unknowns, optional sector filter)
# ----------------------------------------------------------------------------------------------
def numerator_parities(D):
    pD = core.parity_of(D); pA = core.parity_of(A)
    if pD is None: return None
    return [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]


def build_fixed_D(D, p, allowed=None, parity_free=False, target_sector=False):
    """Columns of the cleared identity for fixed D (q=1):
         sum_i N_i-columns + Q-columns  = -A D^2     (mod p if p, exact otherwise)
       N_i column for monomial m: -(D F_i + F D_i) m + F D d_i(m);  Q column: -F^2 m.
       parity_free: no Y-parity restriction on N, Q (needed for mixed-parity D).
       allowed(kind,i,m): sector filter (symmetry.sector_filter) applied to N and Q monomials.
       target_sector: also add columns + m D^2 for Y-odd m of degree 6 in the character sector of A
                      (served-target computation).  Returns dict(cols, colpolys, rhs, nN, nQ, nT)."""
    Fp = reduce_mod(F, p) if p else F
    Ap = reduce_mod(A, p) if p else A
    Dp = reduce_mod(D, p) if p else D
    d = core.degree(D)
    pars = None if parity_free else numerator_parities(D)
    if pars is None and not parity_free:
        raise ValueError('mixed parity D: use parity_free=True')
    degN, degQ = d + 1, 2 * d - 6
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    FD = mul(Fp, Dp, p)
    G = [scale(add(mul(Dp, Fi[i], p), mul(Fp, Di[i], p), p), (p - 1) if p else -1, p) for i in range(NV)]
    F2 = scale(mul(Fp, Fp, p), (p - 1) if p else -1, p)
    cols, colpolys = [], []
    for i in range(NV):
        ms = monomials(degN, None if parity_free else pars[i])
        if allowed is not None:
            ms = symmetry.filter_monomials(ms, 'N', i, allowed)
        for m in ms:
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); colpolys.append(col)
    nN = len(cols)
    if degQ >= 0:
        ms = monomials(degQ, None if parity_free else core.parity_of(A))
        if allowed is not None:
            ms = symmetry.filter_monomials(ms, 'Q', -1, allowed)
        for m in ms:
            cols.append(('Q', -1, m)); colpolys.append(shift(F2, m, 1, p))
    nQ = len(cols) - nN
    D2 = mul(Dp, Dp, p)
    rhs = scale(mul(Ap, D2, p), (p - 1) if p else -1, p)
    nT = 0
    if target_sector:
        grp = GROUP
        ms = [m for m in monomials(6, 1) if symmetry.character({m: 1}, grp.generators[0], grp.order) is not None
              and all(symmetry.character({m: 1}, w, grp.order) == symmetry.character(A, w, grp.order) for w in grp.elements)]
        for m in ms:
            cols.append(('A', -1, m)); colpolys.append(shift(D2, m, 1, p))
        nT = len(ms)
    return {'cols': cols, 'colpolys': colpolys, 'rhs': rhs, 'nN': nN, 'nQ': nQ, 'nT': nT, 'd': d, 'p': p}


GROUP = symmetry.diagonal_group(F)
TARGET_SECTOR = [m for m in monomials(6, 1)
                 if all(symmetry.character({m: 1}, w, GROUP.order) == symmetry.character(A, w, GROUP.order) for w in GROUP.elements)]


def modp_consistency(S, compressed=False, k_extra=40, seed=1):
    """rank(M), rank([M|b]) mod p for a system dict S.  compressed: random row compression (lower bounds on ranks)."""
    p = S['p']
    t0 = time.time()
    if not compressed:
        M, rows = build_modp(S['colpolys'], p, extra=[S['rhs']])
        C = len(S['colpolys'])
        rMb = int(M.rank())
        # rank of M alone: drop last column
        R = M.nrows()
        flat = [0] * (R * C)
        rows2 = {}
        for j, col in enumerate(S['colpolys']):
            for e, v in col.items():
                flat[rows[e] * C + j] = v % p
        M0 = nmod_mat(R, C, flat, p)
        rM = int(M0.rank())
        return {'rank_M': rM, 'rank_Mb': rMb, 'consistent': rM == rMb, 'rows': R, 'cols': C,
                'nullity': C - rM, 'seconds': round(time.time() - t0, 2), 'compressed': False}
    else:
        import numpy as np, scipy.sparse as sp
        rows = {}
        for col in S['colpolys']:
            for e in col: rows.setdefault(e, len(rows))
        for e in S['rhs']: rows.setdefault(e, len(rows))
        C = len(S['colpolys']); R = len(rows)
        ri, ci, vals = [], [], []
        for j, col in enumerate(S['colpolys']):
            for e, v in col.items():
                ri.append(rows[e]); ci.append(j); vals.append(v % p)
        for e, v in S['rhs'].items():
            ri.append(rows[e]); ci.append(C); vals.append(v % p)
        Msp = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(R, C + 1))
        k = min(R, C + 1 + k_extra)
        Cm = fixed_d._compress(Msp, p, k, seed)
        Mn = nmod_mat(k, C + 1, Cm.ravel().tolist(), p)
        rMb = int(Mn.rank())
        Mn0 = nmod_mat(k, C, Cm[:, :C].copy().ravel().tolist(), p)
        rM = int(Mn0.rank())
        return {'rank_M': rM, 'rank_Mb': rMb, 'consistent': rM == rMb, 'rows': R, 'cols': C,
                'nullity': C - rM, 'seconds': round(time.time() - t0, 2), 'compressed': True, 'k': k}


def koszul_nullity_exact(D, allowed=None, parity_free=False, p=None):
    """EXACT (over Q, or mod p if p given) dimension of the homogeneous solution space of the fixed-D system:
         {(N,Q): -D N(F) + F (D div N - N(D)) = F^2 Q}  =  {(m z/6 + K, Q)}  with
         m in R_d (parity/sector of the N-sector), K = sum_{i<j} k_ij (F_j d_i - F_i d_j) Koszul,
         D div K - K(D) = F Q.   Valid (bijective parametrisation) for d <= 8 (k_ij of degree d-4 <= 4 < 5).
       Returns dict(n_m, n_k, koszul_rank, koszul_nullity, nullity_total).  nullity_total = n_m + koszul_nullity
       is the exact nullity of the (restricted) fixed-D coefficient matrix M."""
    d = core.degree(D)
    assert d <= 8, 'parametrisation injective only for d <= 8'
    pars = None if parity_free else numerator_parities(D)
    pA = core.parity_of(A)
    # m: N_i = m z_i has parity pars[i] -> m parity = pars[i] - parity(z_i) = pars[0]
    par_m = None if parity_free else pars[0]
    ms_m = monomials(d, par_m)
    if allowed is not None:
        z0 = tuple(1 if j == 0 else 0 for j in range(NV))
        ms_m = [m for m in ms_m if allowed('N', 0, tuple(a + b for a, b in zip(m, z0)))]
    n_m = len(ms_m)
    Fi = [pderiv(F, i) for i in range(NV)]
    Di = [pderiv(D, i) for i in range(NV)]
    cols, info = [], []
    for i in range(NV):
        for j in range(i + 1, NV):
            # k_ij parity: K_i = k_ij F_j has parity pars[i] -> par(k) = pars[i] - par(F_j)
            park = None if parity_free else (pars[i] - (0 if j < 3 else 1)) % 2
            ms = monomials(d - 4, park)
            if allowed is not None:
                mFj = next(iter(Fi[j]))
                ms = [m for m in ms if allowed('N', i, tuple(a + b for a, b in zip(m, mFj)))]
            for m in ms:
                # K = m (F_j d_i - F_i d_j):  D div K - K(D) = D (d_i m F_j - d_j m F_i) - m (F_j D_i - F_i D_j)
                dm_i = pderiv({m: Fraction(1)}, i); dm_j = pderiv({m: Fraction(1)}, j)
                t1 = mul(D, sub(mul(dm_i, Fi[j]), mul(dm_j, Fi[i])))
                t2 = mul({m: Fraction(1)}, sub(mul(Fi[j], Di[i]), mul(Fi[i], Di[j])))
                cols.append(sub(t1, t2)); info.append(('k', i, j, m))
    n_k = len(cols)
    degQ = 2 * d - 6
    if degQ >= 0:
        msQ = monomials(degQ, None if parity_free else pA)
        if allowed is not None:
            msQ = symmetry.filter_monomials(msQ, 'Q', -1, allowed)
        for m in msQ:
            cols.append(scale(shift(F, m, 1), -1)); info.append(('Q', m))
    t0 = time.time()
    if p is None:
        r = rank_exact(cols) if cols else 0
    else:
        r = rank_modp([reduce_mod(c, p) for c in cols], p) if cols else 0
    nul = len(cols) - r
    return {'n_m': n_m, 'n_k': n_k, 'n_Q': len(cols) - n_k, 'koszul_rank': r, 'koszul_nullity': nul,
            'nullity_total': n_m + nul, 'seconds': round(time.time() - t0, 2), 'exact': p is None}


def fixed_D_test(D, p, allowed=None, parity_free=False, compressed=False, exact_nullity=True):
    """Full test for one D at prime p: E1 (A D in J?), E2 consistency mod p (ranks), and (if d <= 8 and
    exact_nullity) the exact rational nullity of M for the [RANK] argument.  Returns the record."""
    t0 = time.time()
    d = core.degree(D)
    rec_ = {'d': d, 'p': p, 'sector': allowed is not None, 'parity_free': parity_free}
    S = build_fixed_D(D, p, allowed=allowed, parity_free=parity_free)
    r = modp_consistency(S, compressed=compressed)
    rec_.update({'consistent_modp': r['consistent'], 'rank_M_modp': r['rank_M'], 'rank_Mb_modp': r['rank_Mb'],
                 'rows': r['rows'], 'cols': r['cols'], 'nullity_modp': r['nullity'], 'compressed': r['compressed']})
    if exact_nullity and d <= 8 and not r['consistent']:
        kz = koszul_nullity_exact(D, allowed=allowed, parity_free=parity_free)
        rec_['nullity_exact_Q'] = kz['nullity_total']
        rec_['koszul'] = {k: v for k, v in kz.items() if k != 'seconds'}
        rec_['rigorous'] = (kz['nullity_total'] == r['nullity']) and not r['compressed']
        rec_['result'] = 'inconsistent over C [RANK]' if rec_['rigorous'] else 'inconsistent mod p (nullity mismatch: inconclusive)'
    elif not r['consistent']:
        rec_['result'] = 'inconsistent mod p [MODP]' + (' (compressed)' if r['compressed'] else '')
    else:
        rec_['result'] = 'CONSISTENT mod p -- LIFT NEEDED'
    rec_['seconds'] = round(time.time() - t0, 2)
    return rec_


def e1_test(D, p):
    """A*D in J(F) mod p?"""
    ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
    return bool(ok), info


def sector_allowed(D):
    allowed, sub = symmetry.sector_filter(F, A, q, D, GROUP)
    return allowed, sub


# ----------------------------------------------------------------------------------------------
# served-target sector (sum principle)
# ----------------------------------------------------------------------------------------------
def served_sector(D, p, allowed=None):
    """Basis (mod p) of the space of targets A' in the character sector of A (TARGET_SECTOR monomials)
    served by denominator D: nullspace of [N-cols | Q-cols | A'-cols] projected to A'-coordinates.
    If allowed is None the sector filter of D is used (lossless when D is a G-eigenvector)."""
    if allowed is None and core.parity_of(D) is not None:
        allowed, sub = sector_allowed(D)
    S = build_fixed_D(D, p, allowed=allowed, parity_free=(core.parity_of(D) is None), target_sector=True)
    M, rows = build_modp(S['colpolys'], p)
    X, nul = M.nullspace()
    a_idx = [j for j, c in enumerate(S['cols']) if c[0] == 'A']
    if nul == 0:
        return [], S
    proj = nmod_mat(len(a_idx), nul, [int(X[j, k]) for j in a_idx for k in range(nul)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[S['cols'][j][2]] = v
        basis.append(poly)
    return basis, S


def in_span_modp(target, polys, p):
    """Is target (poly mod p) in the span of polys (mod p)?"""
    polys = [reduce_mod(P, p) if any(isinstance(v, Fraction) for v in P.values()) else P for P in polys]
    r0 = rank_modp(polys, p) if polys else 0
    r1 = rank_modp(polys + [reduce_mod(target, p)], p)
    return r1 == r0


def jacobian_target_cols(p):
    """J(F)_6 Y-odd columns mod p (targets in J are trivially served)."""
    Fi = [pderiv(reduce_mod(F, p), i, p) for i in range(NV)]
    cols = []
    for i in range(NV):
        for m in monomials(1, (1 - (0 if i < 3 else 1)) % 2):
            cols.append(shift(Fi[i], m, 1, p))
    return cols


# ----------------------------------------------------------------------------------------------
# random members of the E1-admissible family (Thom-Sebastiani description)
# ----------------------------------------------------------------------------------------------
def random_poly(ms, p, rng):
    return {m: rng.randrange(1, p) for m in ms}


def random_combo(polys, p, rng):
    out = {}
    for P in polys:
        out = add(out, scale(reduce_mod(P, p) if any(isinstance(v, Fraction) for v in P.values()) else P, rng.randrange(1, p), p), p)
    return out


# ----------------------------------------------------------------------------------------------
# small-height kernel representatives via modular lattice reduction (exactly verified)
# ----------------------------------------------------------------------------------------------
BIGP = [2305843009213693951, 2305843009213693921, 2305843009213693907]   # 61-bit primes (first is 2^61-1)


def small_kernel_reps(P, idx, a, k, par, nprimes=3, verify=True):
    """Short integer representatives of ker_full = {D in R_k^par : a D in J(P)} (variables idx) and of its
    quotient modulo J(P)_k.  Method: compute the kernel modulo m = product of `nprimes` 61-bit primes (CRT of
    the nullspaces; pivot structure must agree), LLL-reduce the lattice {x in Z^n : x mod m in ker_m}, then
    keep the LLL vectors that are EXACTLY in the kernel over Q (verified by exact rank) and independent.
    Returns dict(reps_full, reps_modJ, max_abs_coeff, verified, dim_full, dim_J)."""
    da = core.degree(a)
    Pi = [pderiv(P, i) for i in idx]
    msD = mons(k, idx, par)
    nA = len(msD)
    tpar = (par + core.parity_of(a)) % 2
    Acols = [mul(a, {m: 1}) for m in msD]
    Jcols = []
    for j, i in enumerate(idx):
        pj = 0 if i < 3 else 1
        for m in mons(k + da - 5, idx, (tpar - pj) % 2):
            Jcols.append(shift(Pi[j], m, 1))
    # modular kernels (projected to D-coordinates, rref form)
    bases = []
    for p in BIGP[:nprimes]:
        M, rows = build_modp([reduce_mod(c, p) for c in Acols + Jcols], p)
        X, nul = M.nullspace()
        vecs = [[int(X[j, kk]) for j in range(nA)] for kk in range(nul)]
        vecs = [v for v in vecs if any(v)]
        if not vecs:
            bases.append(([], [])); continue
        Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
        R, rk = Mv.rref()
        rowsR = [[int(R[i, j]) for j in range(nA)] for i in range(rk)]
        piv = [next(j for j in range(nA) if row[j]) for row in rowsR]
        bases.append((rowsR, piv))
    dim_full = len(bases[0][0])
    if dim_full == 0:
        return {'reps_full': [], 'reps_modJ': [], 'max_abs_coeff': 0, 'verified': True, 'dim_full': 0, 'dim_J': 0}
    if any(len(b[0]) != dim_full or b[1] != bases[0][1] for b in bases):
        raise RuntimeError('pivot structure differs between primes')
    # CRT lift of the rref rows
    mod = 1; lifted = [[0] * nA for _ in range(dim_full)]
    for (rowsR, piv), p in zip(bases, BIGP[:nprimes]):
        for i in range(dim_full):
            for j in range(nA):
                lifted[i][j] = core.crt_pair(lifted[i][j], mod, rowsR[i][j], p)[0] if mod > 1 else rowsR[i][j]
        mod *= p
    # lattice {x : x mod `mod` in span(lifted)} generated by lifted rows and mod*e_j
    gens = [row[:] for row in lifted] + [[mod if jj == j else 0 for jj in range(nA)] for j in range(nA)]
    L = fmpz_mat(len(gens), nA, [x for g in gens for x in g]).lll()
    cand = []
    for i in range(L.nrows()):
        row = [int(L[i, j]) for j in range(nA)]
        if any(row): cand.append(row)
    cand.sort(key=lambda r: max(abs(x) for x in r))
    # exact verification: a*D in J over Q
    rJ6 = rank_exact(Jcols) if Jcols else 0
    reps_full = []
    cur = list(Jcols); currank = rJ6
    checked = 0
    for row in cand:
        if len(reps_full) >= dim_full: break
        g = 0
        for x in row: g = math.gcd(g, abs(x))
        D = {msD[j]: Fraction(row[j] // g) for j in range(nA) if row[j]}
        checked += 1
        if verify:
            if rank_exact(Jcols + [mul(a, D)]) != rJ6:
                continue   # not in the exact kernel
        # independence among the chosen reps (projected): rank of reps so far
        if rank_exact(reps_full + [D]) == len(reps_full) + 1:
            reps_full.append(D)
    # modulo J_k
    Jk = []
    if k >= 5:
        for j, i in enumerate(idx):
            pj = 0 if i < 3 else 1
            for m in mons(k - 5, idx, (par - pj) % 2):
                Jk.append(shift(Pi[j], m, 1))
    rJ = rank_exact(Jk) if Jk else 0
    reps_modJ = []; cur = list(Jk); currank = rJ
    for D in reps_full:
        if rank_exact(cur + [D]) > currank:
            reps_modJ.append(D); cur.append(D); currank += 1
    mx = max((abs(int(c)) for D in reps_full for c in D.values()), default=0)
    return {'reps_full': reps_full, 'reps_modJ': reps_modJ, 'max_abs_coeff': mx, 'verified': len(reps_full) == dim_full,
            'dim_full': dim_full, 'dim_J': rJ, 'n_checked': checked}
