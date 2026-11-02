"""Degree-7 Hessian-relaxation closure for mixtures of parity sectors (extends fqlx's single-sector result).

Unknowns: C in S_2^6 (126), gauge (Z_k Y2 in C_5) dropped -> 120 unknowns c_a, a = (m, j), D = sum c_a m F_j.
Necessary condition for a certificate with denominator D (q = 1): C^T Hess(F) C in (J:A)_8 = J_8 + <D8*>.
Quadrics: coordinates of NF_{J_8}(C^T Hess F C) per G_aff character block, the D8* coordinate removed in the trivial block.
For a subset T of the 8 parity sectors (b0, b2, y) we restrict to the unknowns supported in T (others = 0) and certify
"only c = 0" by showing I_3 = S_3 on those unknowns (graded degree-3 Macaulay rank, per character block), which transfers
from F_p to Q-bar (the system has Z_(p)-coefficients because the Jacobian ideal has full rank mod p in degrees 8 and 14).
Fallback when dim (S/I)_3 > 0: Singular dim (and msolve -g 1) on the restricted system.
Usage: python3 d7subsets.py <sizes, e.g. 1 2 3 4> [--singular] [--maxsub N]"""
import sys, os, json, time, itertools, subprocess
from w2lib import *
OUT = PK + '/d7subsets'
os.makedirs(OUT, exist_ok=True)
args = [a for a in sys.argv[1:] if not a.startswith('--')]
sizes = [int(a) for a in args] or [1, 2, 3, 4]
use_singular = '--singular' in sys.argv
maxsub = None
for a in sys.argv[1:]:
    if a.startswith('--maxsub='): maxsub = int(a.split('=')[1])

t0 = time.time()
names = [(m, j) for m in monomials(2) for j in range(NV) if not (j == 5 and m[5] >= 1)]
n = len(names); assert n == 120
chi_u = [chi_add(chi(m), chi_neg(CHI_Z[j])) for m, j in names]
par_u = [(c[0], c[1], c[3]) for c in chi_u]
sectors = sorted(set(par_u))
assert len(sectors) == 8
# Jacobian rank certification in degrees 8 and 14 (needed for the transfer to Q): sum of sector ranks vs Hilbert function
HF = [1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, 1686, 1506, 1246, 951, 666, 426, 246, 126, 56, 21, 6, 1]
jac_ok = {}
for mdeg in (8, 14):
    chars = sorted({chi(e) for e in monomials(mdeg)})
    rk = sum(sector(mdeg, c).rank for c in chars)
    jac_ok[mdeg] = (rk, len(monomials(mdeg)) - HF[mdeg], rk == len(monomials(mdeg)) - HF[mdeg])
print('Jacobian rank mod p (deg: rank, expected, ok):', jac_ok, file=sys.stderr, flush=True)
assert all(v[2] for v in jac_ok.values())
# (J:A)_8 = J_8 + <D8*>: certify mod p in the trivial sector (dim kernel of A: R_8[0] -> R_14[chi_A]) and D8* not in J_8
S8_0 = sector(8, CHI0)
D8p = core.reduce_mod(D8, p)
v8 = S8_0.nf(D8p); assert any(v8), 'D8* in J_8?!'
# full quadric system (per character block)
by_char = {}
for a in range(n):
    for b in range(a, n):
        by_char.setdefault(chi_add(chi_u[a], chi_u[b]), []).append((a, b))
quads = []        # list of (char, {(a,b): v})
for cc, pairs in sorted(by_char.items()):
    S8 = sector(8, cc)
    vecs = []
    for (a, b) in pairs:
        ma, ia = names[a]; mb, ib = names[b]
        P = shift(Hess[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 2 if a != b else 1, p)
        vecs.append(S8.vec(P))
    NF = S8.nf_matrix(vecs)
    m = NF.ncols()
    cols = [[int(NF[r, j]) for r in range(len(pairs))] for j in range(m)]
    if cc == CHI0:
        j0 = max(range(m), key=lambda j: 1 if v8[j] else 0)
        inv = pow(v8[j0], p - 2, p)
        cols = [[(cols[j][r] - (v8[j] * inv % p) * cols[j0][r]) % p for r in range(len(pairs))] for j in range(m) if j != j0]
    for col in cols:
        cf = {pairs[r]: v for r, v in enumerate(col) if v}
        if cf: quads.append((cc, cf))
print(f'built {len(quads)} quadrics in {n} unknowns ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
# linearized rank of the full system
pairs_all = [(a, b) for a in range(n) for b in range(a, n)]
pidx = {pr: i for i, pr in enumerate(pairs_all)}
flat = [0] * (len(quads) * len(pairs_all))
for r, (cc, cf) in enumerate(quads):
    for pr, v in cf.items(): flat[r * len(pairs_all) + pidx[pr]] = v
rk2 = nmod_mat(len(quads), len(pairs_all), flat, p).rank()
log('d7-full-system', d=7, n_unknowns=n, n_quadrics=len(quads), quadric_rank=rk2, jacobian_rank_check={str(k): v for k, v in jac_ok.items()},
    note='Hessian relaxation C^T Hess F C in J_8 + <D8*>; 120 gauge-fixed unknowns; quadrics per G_aff block; dossier: rank 950')

def restricted(T):
    keep = [a for a in range(n) if par_u[a] in T]
    kset = set(keep)
    qs = []
    for cc, cf in quads:
        cf2 = {(a, b): v for (a, b), v in cf.items() if a in kset and b in kset}
        if cf2: qs.append((cc, cf2))
    return keep, qs

def degree3_closure(keep, qs):
    """dim (S/I)_3 on the unknowns `keep` via per-character Macaulay blocks. Returns (quotient_dim, blocks)."""
    kset = set(keep)
    blocks = {}
    for r, (cq, cf) in enumerate(qs):
        for v in keep:
            blocks.setdefault(chi_add(cq, chi_u[v]), []).append((v, r))
    quot, binfo = 0, {}
    for cc, lst in blocks.items():
        cols3 = {}
        rows = []
        for (v, r) in lst:
            row = {}
            for (k, l), val in qs[r][1].items():
                key = tuple(sorted((v, k, l)))
                if key not in cols3: cols3[key] = len(cols3)
                row[cols3[key]] = (row.get(cols3[key], 0) + val) % p
            rows.append(row)
        allm = [tuple(sorted((a, b, c))) for a in keep for b in keep for c in keep if a <= b <= c and chi_add(chi_add(chi_u[a], chi_u[b]), chi_u[c]) == cc]
        for key in allm:
            if key not in cols3: cols3[key] = len(cols3)
        nc = len(cols3)
        flat = [0] * (len(rows) * nc)
        for r, row in enumerate(rows):
            for j, val in row.items(): flat[r * nc + j] = val
        rk3 = nmod_mat(len(rows), nc, flat, p).rank() if rows else 0
        quot += nc - rk3
        binfo[str(cc)] = (len(rows), nc, rk3)
    # degree-3 monomials of characters not hit by any row are counted too
    hit = set(blocks.keys())
    for a in keep:
        for b in keep:
            for c in keep:
                if a <= b <= c:
                    cc = chi_add(chi_add(chi_u[a], chi_u[b]), chi_u[c])
                    if cc not in hit:
                        quot += 1
    return quot, binfo

done = set()
try:
    for line in open(LOG):
        rec = json.loads(line)
        if rec.get('exp') == 'd7-parity-subset' and rec.get('prime') == p:
            done.add(rec['subset'])
except FileNotFoundError:
    pass
for size in sizes:
    subs = list(itertools.combinations(sectors, size))
    if maxsub: subs = subs[:maxsub]
    for T in subs:
        key = ';'.join(''.join(map(str, s)) for s in T)
        if key in done: continue
        t1 = time.time()
        keep, qs = restricted(T)
        if not keep:
            log('d7-parity-subset', d=7, subset=key, size=size, n_unknowns=0, result='EXCLUDED (no unknowns)'); continue
        # quick: square-column test on the restricted system
        quads_r = [{(keep.index(a), keep.index(b)): v for (a, b), v in cf.items()} for cc, cf in qs]
        rk, forced, prs = quad_rank_and_squares(quads_r, len(keep))
        quot3, binfo = degree3_closure(keep, qs)
        res = {'d': 7, 'subset': key, 'size': size, 'n_unknowns': len(keep), 'n_quadrics': len(qs), 'quadric_rank': rk, 'n_forced_by_squares': len(forced),
               'dim_S3': sum(v[1] for v in binfo.values()), 'dim_I3': sum(v[2] for v in binfo.values()), 'dim_quotient3': quot3}
        if quot3 == 0:
            res['result'] = 'EXCLUDED (I_3 = S_3)'
        else:
            res['result'] = 'OPEN at degree 3'
            if use_singular:
                tag = OUT + '/sub_' + key.replace(';', '_')
                write_systems(quads_r, len(keep), tag)
                v = singular_verdict(tag + '.sing', timeout=int(os.environ.get('SING_TIMEOUT', '1800')))
                res['singular'] = {k: v[k] for k in ('sdim', 'svdim', 'rc', 'seconds')}
                if v['sdim'] == 0: res['result'] = 'EXCLUDED (Singular dim 0)'
                elif v['sdim'] is not None and v['sdim'] > 0: res['result'] = 'POSITIVE-DIMENSIONAL (Singular)'
                else: res['result'] = 'OPEN (Singular timeout/failed)'
        res['seconds'] = round(time.time() - t1, 1)
        log('d7-parity-subset', **res, note='mixture of the listed (b0,b2,y) parity sectors; EXCLUDED = only C=0 over F_p-bar hence over Q-bar')
print('DONE', file=sys.stderr)
