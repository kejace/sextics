"""Exact fixed-denominator system (UV or XY frame) + rigorous Farkas certificate of inconsistency.

System (q=1):  for unknown N_i (deg d+1) and Q (deg 2d-6):
    -D N(F) + F (D div N - N(D)) - F^2 Q  =  -A D^2         (columns: N-monomials, Q-monomials; rhs = -A D^2)
Optional sector filter `allowed(kind,i,m)` (symmetry.sector_filter) restricts unknown monomials; optional parity
restriction (numerator_parities) as in the engine.  Columns are scaled to integers; rhs scaled to integers.

farkas(): mod-p rank test; if inconsistent, find a row basis R and column basis C of [M|b] mod p with b in C,
solve the square exact system [M|b]_{R,C}^T y = e_b over Q (flint Dixon), and verify EXACTLY that
y^T M = 0 and y^T b != 0.  Such a y proves inconsistency over Q, hence over C (rational linear system)."""
import sys, time, json, random
sys.path.insert(0, '.')
from k3lib import *
from hodge import fixed_d
from math import gcd

def build_exact(F, A, q, D, allowed=None, parity=True, verbose=False):
    t0 = time.time()
    d = core.degree(D)
    degN = d + 6 * q - 5; degQ = 2 * d - 6
    Fi = [pderiv(F, i) for i in range(NV)]; Di = [pderiv(D, i) for i in range(NV)]
    FD = mul(F, D)
    G = [scale(add(scale(mul(D, Fi[i]), q), mul(F, Di[i])), -1) for i in range(NV)]
    Fq1 = scale(power(F, q + 1), -1)
    pars = fixed_d.numerator_parities(D, A) if parity else [None] * NV
    pQ = core.parity_of(A) if parity else None
    cols, colpolys = [], []
    for i in range(NV):
        for m in monomials(degN, pars[i]):
            if allowed is not None and not allowed('N', i, m): continue
            col = shift(G[i], m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
            cols.append(('N', i, m)); colpolys.append(col)
    if degQ >= 0:
        for m in monomials(degQ, pQ):
            if allowed is not None and not allowed('Q', -1, m): continue
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1))
    rhs = scale(mul(A, mul(D, D)), -1)
    rows = {}
    for col in colpolys:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    icols, _ = polys_to_int_rows(colpolys, rows)
    irhs, _ = polys_to_int_rows([rhs], rows); irhs = irhs[0]
    S = dict(d=d, cols=cols, icols=icols, irhs=irhs, rows=rows, nrows=len(rows), ncols=len(cols),
             nN=sum(1 for c in cols if c[0] == 'N'), nQ=sum(1 for c in cols if c[0] == 'Q'),
             nnz=sum(len(c) for c in icols), build_time=round(time.time() - t0, 1))
    if verbose: print(f"[exact system] d={d} rows={S['nrows']} cols={S['ncols']} (N {S['nN']}, Q {S['nQ']}) nnz={S['nnz']} {S['build_time']}s", flush=True)
    return S

def dense_modp(S, p, rows_subset=None, with_rhs=True):
    rs = rows_subset if rows_subset is not None else list(range(S['nrows']))
    rpos = {r: i for i, r in enumerate(rs)}
    R, C = len(rs), S['ncols'] + (1 if with_rhs else 0)
    flat = [0] * (R * C)
    for j, col in enumerate(S['icols']):
        for r, v in col.items():
            i = rpos.get(r)
            if i is not None: flat[i * C + j] = v % p
    if with_rhs:
        for r, v in S['irhs'].items():
            i = rpos.get(r)
            if i is not None: flat[i * C + S['ncols']] = v % p
    return nmod_mat(R, C, flat, p)

def pivot_cols(Mr, rank, C):
    piv = []
    for i in range(rank):
        j = 0
        while j < C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    return piv

def farkas(S, p=1000000007, verbose=True, seed=0):
    """Returns dict(status, rank, rank_aug, ..., certificate y (dict row->int) if inconsistent and verified)."""
    t0 = time.time()
    nr, nc = S['nrows'], S['ncols']
    Mb = dense_modp(S, p)
    rank_aug = Mb.rank()
    M0 = dense_modp(S, p, with_rhs=False)
    rank = M0.rank()
    out = dict(prime=p, rows=nr, cols=nc, rank_modp=rank, rank_aug_modp=rank_aug, nullity_modp=nc - rank)
    if rank_aug == rank:
        out['status'] = 'consistent-modp'
        if verbose: print('[farkas] consistent mod p (rank', rank, ')', flush=True)
        return out
    out['status'] = 'inconsistent-modp'
    t1 = time.time()
    # row basis of [M|b]: pivots of rref of transpose
    Tr, rk = Mb.transpose().rref()
    assert rk == rank_aug
    R = pivot_cols(Tr, rk, nr)
    MbR = dense_modp(S, p, rows_subset=R)
    MbRr, rk2 = MbR.rref(); assert rk2 == rank_aug
    Cset = pivot_cols(MbRr, rk2, nc + 1)
    assert Cset[-1] == nc, 'b column must be a pivot for an inconsistent system'
    t2 = time.time()
    if verbose: print(f'[farkas] inconsistent mod p: rank {rank} aug {rank_aug}; row basis {len(R)} col basis {len(Cset)} ({t2-t1:.1f}s)', flush=True)
    # exact square system: Sq[i][j] = [M|b][R[i], C[j]]
    n = len(R)
    Sq = fmpz_mat(n, n)
    rpos = {r: i for i, r in enumerate(R)}
    for jj, j in enumerate(Cset):
        col = S['icols'][j] if j < nc else S['irhs']
        for r, v in col.items():
            i = rpos.get(r)
            if i is not None: Sq[i, jj] = v
    e = fmpz_mat(n, 1); e[n - 1, 0] = 1    # the b column is last in Cset
    y = Sq.transpose().solve(e)   # fmpq_mat, Sq^T y = e
    t3 = time.time()
    # scale to integers
    den = 1
    ys = [y[i, 0] for i in range(n)]
    for v in ys:
        den = den * int(v.q) // gcd(den, int(v.q))
    yint = {R[i]: int(ys[i].p) * (den // int(ys[i].q)) for i in range(n) if ys[i] != 0}
    # exact verification: y^T M == 0 for all columns, y^T b != 0
    bad = 0
    for j, col in enumerate(S['icols']):
        s = 0
        for r, v in col.items():
            yv = yint.get(r)
            if yv is not None: s += yv * v
        if s != 0:
            bad += 1
    yb = sum(yint.get(r, 0) * v for r, v in S['irhs'].items())
    t4 = time.time()
    out.update(dict(certificate_support=len(yint), certificate_bits=max(abs(v).bit_length() for v in yint.values()),
                    verified=(bad == 0 and yb != 0), bad_columns=bad, yTb_nonzero=(yb != 0),
                    t_modp=round(t2 - t0, 1), t_solve=round(t3 - t2, 1), t_verify=round(t4 - t3, 1)))
    if out['verified']:
        out['status'] = 'inconsistent-Q-certified'
        out['certificate'] = yint
    if verbose: print(f"[farkas] exact certificate: support {len(yint)} bits {out['certificate_bits']} verified={out['verified']} (solve {t3-t2:.1f}s verify {t4-t3:.1f}s)", flush=True)
    return out

def save_certificate(path, S, out, meta):
    rec = dict(meta=meta, prime=out['prime'], rows=out['rows'], cols=out['cols'], rank_modp=out['rank_modp'],
               rank_aug_modp=out['rank_aug_modp'], status=out['status'], verified=out.get('verified'),
               certificate_support=out.get('certificate_support'), certificate_bits=out.get('certificate_bits'))
    if 'certificate' in out:
        inv = {i: e for e, i in S['rows'].items()}
        rec['certificate'] = {str(inv[r]): str(v) for r, v in out['certificate'].items()}
        rec['column_scaling_note'] = 'columns of the system are each scaled to integers (polys_to_int_rows); the certificate y satisfies y^T M_int = 0, y^T b_int != 0, which is scaling-invariant'
    json.dump(rec, open(path, 'w'))
    return rec
