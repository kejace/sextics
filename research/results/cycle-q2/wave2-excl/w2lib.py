"""Wave-2 EXCLUSION/RIGOR library for lane cycle-q2.

Conventions: engine v2 (hodge.core / hodge.symmetry), certificate identity
  E = A D^2 - q D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) = F^(q+1) Q.
All systems here are built with F2 = 2F (integer coefficients); consistency for (2F,A,q,D)
is equivalent to consistency for (F,A,q,D) (N -> N/2, Q -> Q/2^(q+1)).

Farkas certificate of inconsistency: an exact rational vector y supported on r+1 rows of the
sector-restricted system with y^T M = 0 and y^T b = 1 (verified with exact arithmetic); it proves
that no (N,Q) over any field of characteristic 0 solves the sector system, hence (averaging) that D
admits no certificate at all.
"""
import sys, os, json, time, hashlib, math
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
from fractions import Fraction
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, symmetry, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat, fmpz_mat, fmpq_mat

LANE = '/tmp/hodge/results/cycle-q2'
import socket
_HOST = socket.gethostname()
# remote jobs write to a remote-only file (merged into the canonical local log by merge.py); this avoids the rsync race
# in which `rjob up` pushed the local copy of wave2-excl.log.jsonl over the remote one and lost lines appended there.
LOG = f'{LANE}/wave2-excl.log.jsonl' if _HOST == 'brunos' else f'{LANE}/wave2-excl.remote.log.jsonl'
PACKET = f'{LANE}/wave2-excl'
P1, P2, P3 = 1000000007, 998244353, 2147483629   # wave-2 primes (wave-1 used 2^31-1, 4194301, 4194287; dossier 10007, 65521)

def log(**kw):
    kw.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    kw.setdefault('engine', 2)
    kw.setdefault('lens', 'excl')
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw) + '\n')

def load():
    F, A, q, rec = core.load_instance('cycle-q2')
    G = symmetry.diagonal_group(F)
    return F, A, q, G

def integral(D):
    """scale D by a positive integer so that all coefficients are integers; returns (D', factor)."""
    den = 1
    for v in D.values():
        den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    return {e: int(Fraction(v) * den) for e, v in D.items()}, den

def sector_system(F, A, q, D, p, G=None):
    """engine-v2 character-restricted system mod p (lossless by averaging). Returns (S, sub)."""
    G = G or symmetry.diagonal_group(F)
    return symmetry.restricted_system(F, A, q, D, p=p, group=G)

def solve_modp(F, A, q, D, p, G=None, compressed=False):
    S, sub = sector_system(F, A, q, D, p, G)
    r = fixed_d.solve_compressed(S) if compressed else S.solve()
    out = {'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows,
           'rank_note': 'engine rank = rank of [M|b]; for inconsistent systems rank(M) = rank-1',
           'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ, 'full_cols': S.full_cols, 'stab': len(sub), 'p': p}
    return out, S, r

# ---------------------------------------------------------------------------------------------
# exact integer system with the same columns as a mod-p restricted system

def exact_columns(F, A, q, D, cols):
    """exact integer column polys (for F2 = 2F, D integral) and rhs, same column order as `cols`."""
    F2 = scale(F, 2)
    Dz, _ = integral(D)
    d = core.degree(D)
    Fi = [pderiv(F2, i) for i in range(NV)]
    Di = [pderiv(Dz, i) for i in range(NV)]
    FD = mul(F2, Dz)
    G = [scale(add(scale(mul(Dz, Fi[i]), q), mul(F2, Di[i])), -1) for i in range(NV)]
    Fq1 = scale(power(F2, q + 1), -1)
    colpolys = []
    for (kind, i, m) in cols:
        if kind == 'N':
            col = shift(G[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
        else:
            col = shift(Fq1, m, 1)
        colpolys.append({e: int(v) for e, v in col.items()})
    rhs = scale(mul(A, mul(Dz, Dz)), -1)
    rhs = {e: int(v) for e, v in rhs.items()}
    return colpolys, rhs, Dz

def farkas(F, A, q, D, p, G=None, save=None, Dstr=None, verbose=False, max_store=60_000_000):
    """Exact Farkas certificate of inconsistency of the sector system for D (if inconsistent mod p).

    Returns dict(status, ...).  status in {'consistent-mod-p', 'certified', 'failed'}."""
    t0 = time.time()
    out, S, r = solve_modp(F, A, q, D, p, G)
    res = dict(out); res['D'] = Dstr or core.to_str(D); res['d'] = S.d
    if out['consistent']:
        res['status'] = 'consistent-mod-p'; res['seconds'] = round(time.time() - t0, 2)
        return res
    n = S.ncols; m = S.nrows
    piv = r['pivots']
    assert piv[-1] == n, 'inconsistent system must have the rhs column as last pivot'
    C = piv[:-1]; rk = len(C)
    assert rk == r['rank'] - 1, (rk, r['rank'])   # engine 'rank' is the augmented rank when inconsistent
    res['rank_M'] = rk; res['rank_aug'] = rk + 1; res['nullity_M'] = n - rk
    # exact system
    colpolys, rhs, Dz = exact_columns(F, A, q, D, S.cols)
    rows = {}
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    mrows = len(rows)
    # transpose mod p to get pivot rows of [M|b] (as pivot columns of the transpose)
    flatT = [0] * ((n + 1) * mrows)
    for j, col in enumerate(colpolys):
        base = j * mrows
        for e, v in col.items():
            flatT[base + rows[e]] = v % p
    base = n * mrows
    for e, v in rhs.items():
        flatT[base + rows[e]] = v % p
    MT = nmod_mat(n + 1, mrows, flatT, p)
    MTr, rkT = MT.rref()
    assert rkT == rk + 1, (rkT, rk)
    Srows = []
    for i in range(rkT):
        j = 0
        while j < mrows and int(MTr[i, j]) == 0:
            j += 1
        Srows.append(j)
    assert len(Srows) == rk + 1
    # B = M[S, C]  ((rk+1) x rk), exact integers; T = B^T (rk x (rk+1)); want y with T y = 0
    rowpos = {ri: k for k, ri in enumerate(Srows)}
    Bt = [[0] * (rk + 1) for _ in range(rk)]
    for jj, j in enumerate(C):
        for e, v in colpolys[j].items():
            ri = rows[e]
            if ri in rowpos:
                Bt[jj][rowpos[ri]] = v
    # mod-p kernel of T to pick j0
    Tp = nmod_mat(rk, rk + 1, [v % p for row in Bt for v in row], p)
    X, nullity = Tp.nullspace()
    assert nullity == 1, nullity
    j0 = next(k for k in range(rk + 1) if int(X[k, 0]) != 0)
    others = [k for k in range(rk + 1) if k != j0]
    Tsq = fmpz_mat(rk, rk, [Bt[i][k] for i in range(rk) for k in others])
    rhsv = fmpz_mat(rk, 1, [-Bt[i][j0] for i in range(rk)])
    t1 = time.time()
    ysol = Tsq.solve(rhsv)          # exact rational
    t2 = time.time()
    y = [None] * (rk + 1)
    y[j0] = Fraction(1)
    for idx, k in enumerate(others):
        y[k] = Fraction(int(ysol[idx, 0].p), int(ysol[idx, 0].q))
    # exact verification over all columns and rhs
    yrow = {Srows[k]: y[k] for k in range(rk + 1)}
    bad = 0
    for j, col in enumerate(colpolys):
        s = Fraction(0)
        for e, v in col.items():
            ri = rows[e]
            if ri in yrow:
                s += yrow[ri] * v
        if s != 0:
            bad += 1
    yb = Fraction(0)
    for e, v in rhs.items():
        ri = rows[e]
        if ri in yrow:
            yb += yrow[ri] * v
    ok = (bad == 0 and yb != 0)
    inv_rows = {v: k for k, v in rows.items()}
    cert_rows = [list(inv_rows[ri]) for ri in Srows]
    ytxt = [str(v) for v in y]
    h = hashlib.sha256(('\n'.join(ytxt) + '\n' + json.dumps(cert_rows)).encode()).hexdigest()
    res.update({'status': 'certified' if ok else 'failed', 'farkas_rows': rk + 1, 'bad_cols': bad,
                'yb': str(yb), 'sha256_y': h, 'bits_max': max(max(v.numerator.bit_length(), v.denominator.bit_length()) for v in y),
                't_solve': round(t2 - t1, 2), 'seconds': round(time.time() - t0, 2), 'F_scaled': '2F', 'D_scale': core.to_str(Dz) if len(core.to_str(Dz)) < 400 else 'see D (integralized)'})
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        size = sum(len(s) for s in ytxt)
        cert = {'instance': 'cycle-q2', 'F': '2F', 'q': q, 'D': core.to_str(Dz), 'p_used_for_pivots': p,
                'columns_note': 'character-sector restricted columns of engine v2 symmetry.restricted_system',
                'row_monomials': cert_rows, 'y': ytxt if size <= max_store else None, 'y_omitted_size': None if size <= max_store else size,
                'sha256_y': h, 'y_dot_b': str(yb), 'verified_exact': ok, 'rank_mod_p': rk,
                'pivot_cols': [list(S.cols[j][2]) + [S.cols[j][0], S.cols[j][1]] for j in C] if size <= max_store else None,
                'j0': j0}
        with open(save, 'w') as f:
            json.dump(cert, f)
    if verbose:
        print(res, file=sys.stderr)
    return res

def cert_from_file(path, F, A, q):
    """Re-verify a stored certificate file exactly (independent of pivots)."""
    c = json.load(open(path))
    D = core.parse(c['D'])
    G = symmetry.diagonal_group(F)
    S, sub = symmetry.restricted_system(F, A, q, D, p=c['p_used_for_pivots'], group=G)
    colpolys, rhs, Dz = exact_columns(F, A, q, D, S.cols)
    y = {tuple(r): Fraction(v) for r, v in zip(c['row_monomials'], c['y'])}
    bad = sum(1 for col in colpolys if sum((y.get(e, 0) * v for e, v in col.items()), Fraction(0)) != 0)
    yb = sum((y.get(e, 0) * v for e, v in rhs.items()), Fraction(0))
    return bad == 0 and yb != 0, {'bad': bad, 'yb': str(yb)}
