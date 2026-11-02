"""Wave-2 EXCLUSION/RIGOR library for lane off-paired (q = 1, A = X0^3 Y0^3).

Rigor conventions used everywhere in this packet
------------------------------------------------
* A modular rank is a LOWER bound for the rational rank (rank_p(M) <= rank_Q(M)); a random
  row-compressed rank is a lower bound for the modular rank.
* Fixed-D system  M x = b  (x = coefficients of N_0..N_5, Q; b = -A D^2), restricted to the
  character sector forced by an eigenvector D (lossless by averaging).  The kernel of M over Q
  contains the EXPLICIT rational subspace (q = 1)
      {(L Z, 0) : L in S_d[chi_A + chi_D]}            (Euler gauge; identity B_D(LZ) = 6 L D)
    + {(F W, B_D(W)) : W_i in S_{d-5}[chi_A + chi_D + w_i]}   (F-multiples)
  whose dimension is exactly  e = s_d(L) + sum_i s_{d-5}(W_i) - s_{d-6}(L)  (intersection L = F L').
  Hence rank_Q(M) <= C - e.  If rank_p(M) = C - e and rank_p([M|b]) = C - e + 1 then
      rank_Q([M|b]) >= C - e + 1 > C - e >= rank_Q(M),
  so the rational system is inconsistent: PROVEN exclusion of the denominator D (and by the
  prescribed-support theorem of every D' with the same reduced support, at q = 1).
* E1 test  A D in J(F)_{d+6}: the rational rank of the Jacobian generator matrix in a character
  sector equals dim J_{d+6}[chi] computed EXACTLY from the equivariant Koszul Hilbert series of
  the complete intersection J(F) (F is smooth: coefficient-load bound).  If rank_p = that number
  and rank_p(aug) = +1, then A D not in J(F) over Q: PROVEN.
"""
import sys, os, json, time, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from collections import defaultdict
from hodge import core, fixed_d, symmetry, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, power, reduce_mod
import numpy as np
from flint import nmod_mat

LANE = '/tmp/hodge/results/off-paired'
PK = LANE + '/wave2-excl'
LOG = LANE + '/wave2-excl.log.jsonl'
IID = 'off-paired'
P1, P2 = 2147483629, 2147483587      # wave-2 primes; wave 1 used 1000003 / 2147483647 / 1073741827 (msolve)
WG = (0, 2, 0, 1, 3, 3)              # generator g (sixth-root weights)
WH = (0, 0, 1, 0, 0, 5)              # generator h
N6 = 6

_cache = {}

def load():
    if 'inst' not in _cache:
        F, A, q, rec = core.load_instance(IID)
        _cache['inst'] = (F, A, q)
    return _cache['inst']

def group(F):
    if 'G' not in _cache:
        _cache['G'] = symmetry.diagonal_group(F)
    return _cache['G']

def log(rec, path=LOG):
    rec = dict(rec)
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', fixed_d.ENGINE_VERSION)
    with open(path, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec

# ---------------------------------------------------------------- characters
def char(m):
    return (sum(a * b for a, b in zip(WG, m)) % N6, sum(a * b for a, b in zip(WH, m)) % N6)

def char_of(P):
    cs = {char(e) for e in P}
    return cs.pop() if len(cs) == 1 else None

def cadd(a, b):
    return ((a[0] + b[0]) % N6, (a[1] + b[1]) % N6)

def cneg(a):
    return ((-a[0]) % N6, (-a[1]) % N6)

VARCHAR = [char(tuple(1 if j == i else 0 for j in range(NV))) for i in range(NV)]
_mc = {}

def monos_char(deg, c):
    key = (deg, c)
    if key not in _mc:
        _mc[key] = [m for m in monomials(deg) if char(m) == c] if deg >= 0 else []
    return _mc[key]

def equivariant_hilbert(deg, c):
    """dim (S/J(F))_{deg, c} for the complete intersection J(F) of six quintics (Koszul exact):
    coefficient of t^deg [c] in prod_i (1 - t^5 [-w_i]) / (1 - t [w_i])."""
    ser = defaultdict(int); ser[(0, (0, 0))] = 1
    for i in range(NV):
        wi = VARCHAR[i]
        new = defaultdict(int)
        for (dg, ch), v in ser.items():
            k = 0
            while dg + k <= deg:
                new[(dg + k, ((ch[0] + k * wi[0]) % N6, (ch[1] + k * wi[1]) % N6))] += v
                k += 1
        ser = new
    for i in range(NV):
        wi = cneg(VARCHAR[i])
        new = defaultdict(int)
        for (dg, ch), v in ser.items():
            new[(dg, ch)] += v
            if dg + 5 <= deg:
                new[(dg + 5, cadd(ch, wi))] -= v
        ser = new
    return ser[(deg, c)]

def dim_J_sector(deg, c):
    """Exact dim J(F)_{deg}[c] over Q (Koszul)."""
    return len(monos_char(deg, c)) - equivariant_hilbert(deg, c)

# ---------------------------------------------------------------- polynomial blocks
_blocks = None
def blocks():
    """Named building blocks as expression strings (parenthesized) for core.parse."""
    global _blocks
    if _blocks is None:
        F, A, q = load()
        def s(P): return '(' + core.to_str(P) + ')'
        f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6')
        Qx = sub(sub(F, f6), g6)
        H = core.parse('X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1/2')
        C = core.parse('X1**3*Y1/2+X0**2*X1*Y0/3')
        D0 = core.parse('X0**4+X0*X1*Y0*Y1/4+X1*X2*Y0*Y2/9')
        R = core.parse('X0**3*Y0*Y1/2+X0**2*X2*Y0*Y2/3+X1**2*X2*Y1*Y2/2')
        B = {'f6': f6, 'g6': g6, 'Qx': Qx, 'H': H, 'C': C, 'D0': D0, 'R': R, 'F': F, 'A': A}
        for i in range(NV):
            B['F' + core.VARS[i]] = pderiv(F, i)
        B['P2'] = core.parse('X2*Y2'); B['M'] = core.parse('X2**6-Y2**6'); B['S'] = core.parse('X2**6+Y2**6')
        B['M2'] = core.parse('X2**2-Y2**2'); B['S2'] = core.parse('X2**2+Y2**2')
        _blocks = {k: s(v) for k, v in B.items()}
        _blocks['_polys'] = B
    return _blocks

def expand(expr):
    """Expand an expression over block names and variables into a core.parse string."""
    B = blocks()
    out = expr
    # longest names first
    for k in sorted((k for k in B if k != '_polys'), key=len, reverse=True):
        out = _replace_token(out, k, B[k])
    return out

def _replace_token(s, name, rep):
    import re
    return re.sub(r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])', rep, s)

def parse_expr(expr):
    return core.parse(expand(expr))

# ---------------------------------------------------------------- sector counting
def _count(deg, req, sub, n):
    if deg < 0:
        return 0
    mons = monomials(deg)
    if not sub:
        return len(mons)
    E = np.array(mons, dtype=np.int64)
    W = np.array(sub, dtype=np.int64).T
    T = (E @ W) % n
    r = np.array(req, dtype=np.int64)
    return int(np.sum(np.all(T == r[None, :], axis=1)))

def explicit_kernel_dim(F, A, q, D, G=None):
    """Exact dimension of the explicit rational kernel of the restricted fixed-D system (q = 1)."""
    G = G or group(F)
    allowed, sub = symmetry.sector_filter(F, A, q, D, G)
    n = allowed.n; reqN = allowed.reqN
    d = core.degree(D)
    if not sub:
        raise ValueError('D not an eigenvector of the diagonal group')
    reqL = [(reqN[0][g] - sub[g][0]) % n for g in range(len(sub))]
    for i in range(NV):
        assert all((reqN[i][g] - sub[g][i]) % n == reqL[g] for g in range(len(sub)))
    degL = d + 6 * q - 6
    e = _count(degL, reqL, sub, n)
    if q == 1:
        e += sum(_count(d - 5, reqN[i], sub, n) for i in range(NV)) - _count(d - 6, reqL, sub, n)
    return e

# ---------------------------------------------------------------- receipts
def fixedD_receipt(F, A, q, D, p, label='', family='', compressed=None, verbose=False):
    """Restricted fixed-D solve mod p with the explicit-kernel rank sandwich."""
    t0 = time.time()
    G = group(F)
    S, sub = symmetry.restricted_system(F, A, q, D, p, G, verbose=verbose)
    C = S.ncols; Rr = S.nrows
    if compressed is None:
        compressed = Rr * (C + 1) > 40_000_000
    r = fixed_d.solve_compressed(S, verbose=verbose) if compressed else S.solve(verbose=verbose)
    consistent = r['consistent']
    rank_aug = r['rank']
    rank_M = rank_aug - (0 if consistent else 1)
    e = explicit_kernel_dim(F, A, q, D, G)
    proven = (not consistent) and (rank_M == C - e)
    rec = {'exp': 'fixedD', 'D': label, 'family': family, 'd': core.degree(D), 'char': char_of(reduce_mod(D, p)),
           'prime': p, 'rows': Rr, 'cols': C, 'nN': S.nN, 'nQ': S.nQ, 'subgroup': len(sub),
           'rank_M': rank_M, 'rank_aug': rank_aug, 'nullity': C - rank_M, 'explicit_kernel': e,
           'consistent': consistent, 'compressed': bool(compressed),
           'result': ('solvable-modp' if consistent else ('PROVEN-inconsistent' if proven else 'inconsistent-modp(nullity>explicit)')),
           'seconds': round(time.time() - t0, 1)}
    if consistent:
        rec['x'] = r['x']; rec['S'] = S
    return rec

def e1_receipt(F, A, q, D, p, label='', family=''):
    """E1 test A*D in J(F) in D's sector, with the exact Koszul dimension of J_{d+6}[chi]."""
    t0 = time.time()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
    cD = char_of(Dp); cA = char_of(Ap)
    assert cD is not None
    d = core.degree(D); degN = d + 6 * q - 5
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
    Rn, Cn = len(rows), len(cols)
    flat = [0] * (Rn * (Cn + 1))
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * (Cn + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (Cn + 1) + Cn] = v
    M = nmod_mat(Rn, Cn + 1, flat, p)
    ra = M.rank()
    M0 = nmod_mat(Rn, Cn, [flat[r * (Cn + 1) + j] for r in range(Rn) for j in range(Cn)], p)
    r0 = M0.rank()
    cT = cadd(cA, cD)
    dJ = dim_J_sector(d + 6 * q, cT)
    passes = (ra == r0)
    proven_fail = (not passes) and (r0 == dJ)
    rec = {'exp': 'e1', 'D': label, 'family': family, 'd': d, 'char': cD, 'prime': p, 'rows': Rn, 'cols': Cn,
           'rank_J': r0, 'rank_aug': ra, 'dimJ_exact': dJ, 'e1': passes,
           'result': ('e1-pass-modp' if passes else ('PROVEN-E1-fail' if proven_fail else 'e1-fail-modp(rank<exact)')),
           'seconds': round(time.time() - t0, 1)}
    return rec

def full_receipt(F, A, q, D, primes, label, family, do_full=True, compressed=None, verbose=False):
    """E1 (first prime) then, if E1 passes, the full fixed-D solve on each prime.  Returns list of records."""
    out = []
    r1 = e1_receipt(F, A, q, D, primes[0], label, family)
    out.append(r1)
    if not r1['e1']:
        if len(primes) > 1:
            out.append(e1_receipt(F, A, q, D, primes[1], label, family))
        return out
    if do_full:
        for p in primes:
            r = fixedD_receipt_certified(F, A, q, D, p, label, family, compressed=compressed, verbose=verbose)
            out.append({k: v for k, v in r.items() if k not in ('x', 'S')})
            if r['consistent']:
                break
    return out

# ---------------------------------------------------------------- kernel-lift certification
def _build_modular(S, with_rhs=True):
    """Dense nmod_mat of the restricted system (rows x cols[+1])."""
    R, C, p = S.nrows, S.ncols, S.p
    W = C + (1 if with_rhs else 0)
    flat = [0] * (R * W)
    for j, col in enumerate(S.colpolys):
        for e, v in col.items():
            flat[S.rows[e] * W + j] = v
    if with_rhs:
        for e, v in S.rhs.items():
            flat[S.rows[e] * W + C] = v
    return nmod_mat(R, W, flat, p)

def _compressed_matrix(S, k_extra=40, seed=1):
    M = fixed_d._sparse_matrix(S)
    C = S.ncols
    k = min(S.nrows, C + 1 + k_extra)
    Cm = fixed_d._compress(M, S.p, k, seed)
    return nmod_mat(k, C + 1, Cm.ravel().tolist(), S.p), k

def kernel_vectors(S, Maug):
    """Nullspace of the coefficient part of Maug (mod p), as rref'd rows (lists of ints)."""
    R, C, p = Maug.nrows(), S.ncols, S.p
    M0 = nmod_mat(R, C, [int(Maug[i, j]) for i in range(R) for j in range(C)], p)
    X, nul = M0.nullspace()
    if nul == 0:
        return []
    T = nmod_mat(nul, C, [int(X[j, k]) for k in range(nul) for j in range(C)], p)
    Tr, rk = T.rref()
    return [[int(Tr[i, j]) for j in range(C)] for i in range(rk)]

def certify_kernel(F, A, q, D, S, vecs, p, verbose=False):
    """Rationally reconstruct modular kernel vectors (single prime) and verify EXACTLY that each is a
    kernel vector of the rational fixed-D system:  -q D N(F) + F B_D(N) - F^{q+1} Q = 0.
    Returns (n_verified, n_total, max_height, n_unreconstructable)."""
    t0 = time.time()
    Fq1 = power(F, q + 1)
    AD2 = mul(A, mul(D, D))
    nver = 0; maxh = 0; bad = 0
    for v in vecs:
        rat = []
        ok = True
        for a in v:
            r = core.ratrecon(a, p) if a else Fraction(0)
            if r is None:
                ok = False; break
            rat.append(r)
        if not ok:
            bad += 1; continue
        N = [{} for _ in range(NV)]; Q = {}
        for (kind, i, m), x in zip(S.cols, rat):
            if x == 0: continue
            if kind == 'N': N[i][m] = x
            else: Q[m] = x
            maxh = max(maxh, abs(x.numerator), x.denominator)
        E = sub(core.defect(F, A, D, q, N), AD2)
        E = sub(E, mul(Fq1, Q))
        if not E:
            nver += 1
    if verbose:
        print(f'[certify] {nver}/{len(vecs)} kernel vectors verified exactly, {bad} unreconstructable, max height {maxh}, {time.time()-t0:.1f}s', file=sys.stderr)
    return nver, len(vecs), maxh, bad

def fixedD_receipt_certified(F, A, q, D, p, label='', family='', compressed=None, certify=True, verbose=False):
    """Restricted fixed-D solve mod p; rank sandwich via explicit kernel, else via exact kernel-lift certification."""
    t0 = time.time()
    G = group(F)
    S, sub = symmetry.restricted_system(F, A, q, D, p, G, verbose=verbose)
    C = S.ncols; Rr = S.nrows
    if compressed is None:
        compressed = Rr * (C + 1) > 40_000_000
    if compressed:
        Maug, k = _compressed_matrix(S)
    else:
        Maug = _build_modular(S); k = Rr
    t1 = time.time()
    Mr, rank_aug = Maug.rref()
    pivots = []
    for i in range(rank_aug):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    rank_M = rank_aug - (0 if consistent else 1)
    e = explicit_kernel_dim(F, A, q, D, G)
    rec = {'exp': 'fixedD', 'D': label, 'family': family, 'd': core.degree(D), 'char': char_of(reduce_mod(D, p)),
           'prime': p, 'rows': Rr, 'cols': C, 'nN': S.nN, 'nQ': S.nQ, 'subgroup': len(sub), 'k_rows': k,
           'rank_M': rank_M, 'rank_aug': rank_aug, 'nullity': C - rank_M, 'explicit_kernel': e,
           'consistent': consistent, 'compressed': bool(compressed), 't_rref': round(time.time() - t1, 1)}
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
        rec['x'] = x; rec['S'] = S; rec['result'] = 'solvable-modp'
    elif rank_M == C - e:
        rec['result'] = 'PROVEN-inconsistent'; rec['rigor'] = 'explicit-kernel'
    elif certify:
        vecs = kernel_vectors(S, Maug)
        nver, ntot, maxh, bad = certify_kernel(F, A, q, D, S, vecs, p, verbose=verbose)
        rec.update({'kernel_verified': nver, 'kernel_total': ntot, 'kernel_max_height': maxh, 'kernel_unrecon': bad})
        if nver == ntot == C - rank_M:
            rec['result'] = 'PROVEN-inconsistent'; rec['rigor'] = 'kernel-certified'
        else:
            rec['result'] = 'inconsistent-modp(kernel-lift-incomplete)'
    else:
        rec['result'] = 'inconsistent-modp(nullity>explicit)'
    rec['seconds'] = round(time.time() - t0, 1)
    return rec
