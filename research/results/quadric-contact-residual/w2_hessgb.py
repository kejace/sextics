#!/usr/bin/env python3
"""EXACT (mod p) Hessian-congruence closure for character-mixed Jacobian-type denominators of pure Y-parity.

Family: D = sum_k a_k D_k over the Jacobian generators D_k = m F_i of degree d with the parity pattern forced by
Y-parity `parity` of D (plus the known non-Jacobian extras when d = 8, 9, with exact mod-p lifts).  Necessary E2
condition (fqlx's lemma, Koszul fields and F*Q lie in J):   T(a,a) = NF_{J_{2d}}( sum_{k,l} a_k a_l B_{D_k}(N_l) ) = 0.
The gauge (D = F s, which lies in the radical of T) is fixed by setting the coefficients of X0^2, X0 X1, ... (the
monomials s Z_0 / 6 of C_0) to zero.  Output: msolve / Singular input files + Singular dim/vdim run (timeout).
dim 0 (affine) means the only solution is a = 0 over the algebraic closure of F_p, hence (upper semicontinuity)
over Q-bar: no denominator of this degree/parity passes the Hessian congruence -> degree d, parity excluded.

  python3 w2_hessgb.py --d 7 --parity 0 [--timeout 3600] [--p 2147483647]
"""
import sys, os, json, time, argparse, subprocess
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import w2lib as W
from hodge import core, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, reduce_mod
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--parity', type=int, default=0)
ap.add_argument('--p', type=int, default=W.P0)
ap.add_argument('--timeout', type=int, default=3600)
ap.add_argument('--tag', default='hgb')
ap.add_argument('--no-run', action='store_true')
a = ap.parse_args()
d, par, p = a.d, a.parity, a.p
MIXED = par < 0          # --parity -1: mixed Y-parity (both parity patterns, rows of both parities)
t0 = time.time()
Fp = reduce_mod(W.F, p); Ap = reduce_mod(W.A, p)
Fip = [pderiv(Fp, i, p) for i in range(NV)]

# members: Jacobian generators with parity pattern + extras (exact mod-p lifts via linear solve)
members = []
for i in range(NV):
    if MIXED:
        ms = monomials(d - 5)
    else:
        parsN_i = (1 - par - (0 if i < 3 else 1)) % 2      # parity of N_i = A C_i
        ms = monomials(d - 5, (parsN_i + 1) % 2)           # parity of C_i
    for m in ms:
        Dk = shift(Fip[i], m, 1, p)
        Nk = [{} for _ in range(NV)]
        Nk[i] = shift(Ap, m, 1, p)
        members.append((Dk, Nk, ('J', i, m)))
nJ = len(members)
# extras: exact E1 kernel elements beyond J (from the exact Q computation), lifted mod p
fam = W.e1_family(d, None, parity=None if MIXED else par)
extras = fam['extras']
if extras:
    gens6 = W.jgens(d + 6, None, None if MIXED else (par + 1) % 2)
    gp = [shift(Fip[i], m, 1, p) for i, m in gens6]
    rows = {}
    for P in gp:
        for e in P:
            rows.setdefault(e, len(rows))
    for E in extras:
        T = mul(Ap, reduce_mod(E, p), p)
        for e in T:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(gp)
    flat = [0] * (R * (C + len(extras)))
    Wd = C + len(extras)
    for j, P in enumerate(gp):
        for e, v in P.items():
            flat[rows[e] * Wd + j] = v
    for k, E in enumerate(extras):
        T = mul(Ap, reduce_mod(E, p), p)
        for e, v in T.items():
            flat[rows[e] * Wd + C + k] = v
    M = nmod_mat(R, Wd, flat, p)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < Wd and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    # solve for each extra: columns C+k must be non-pivot
    for k, E in enumerate(extras):
        assert (C + k) not in piv, 'extra not liftable?!'
        x = [0] * C
        for i, j in enumerate(piv):
            if j < C:
                x[j] = (-int(Mr[i, C + k])) % p
        Nk = [{} for _ in range(NV)]
        for (i, m), v in zip(gens6, x):
            if v:
                Nk[i] = add(Nk[i], {m: v}, p)
        chk = {}
        for i in range(NV):
            chk = add(chk, mul(Nk[i], Fip[i], p), p)
        assert sub(chk, mul(Ap, reduce_mod(E, p), p), p) == {}, 'lift check failed'
        members.append((reduce_mod(E, p), Nk, ('E', k)))
K = len(members)
print(f'[hessgb] d={d} parity={par} members {K} (J {nJ}, extras {len(extras)})', flush=True)

# reducer modulo J_{2d} on the parity of A D^2 (odd)
red = jacobian._Reducer(Fp, 2 * d, p, None if MIXED else 1)
mons2 = red.mons; idx2 = red.idx
pivots = red.pivots
nonpiv = [j for j in range(len(mons2)) if j not in set(pivots)]
dimR = len(nonpiv)
Rb = nmod_mat(len(pivots), len(mons2), [0] * (len(pivots) * len(mons2)), p)
flat = [0] * (len(pivots) * len(mons2))
for i, row in enumerate(red.rows):
    for e, v in row.items():
        flat[i * len(mons2) + idx2[e]] = v
Rb = nmod_mat(len(pivots), len(mons2), flat, p)
print(f'[hessgb] R_{2*d}[odd] dim {dimR} (J rank {len(pivots)} of {len(mons2)}) {time.time()-t0:.0f}s', flush=True)

# pair vectors B_{D_k}(N_l), all K*K, as rows
vecs = []
for l, (Dl, Nl, _) in enumerate(members):
    divN = {}
    for i in range(NV):
        if Nl[i]:
            divN = add(divN, pderiv(Nl[i], i, p), p)
    for k, (Dk, Nk, _) in enumerate(members):
        Bv = mul(Dk, divN, p)
        for i in range(NV):
            if Nl[i]:
                Bv = sub(Bv, mul(Nl[i], pderiv(Dk, i, p), p), p)
        v = [0] * len(mons2)
        for e, c in Bv.items():
            v[idx2[e]] = c
        vecs.append(v)
Vm = nmod_mat(K * K, len(mons2), [x for v in vecs for x in v], p)
Pm = nmod_mat(K * K, len(pivots), [v[j] for v in vecs for j in pivots], p)
NF = Vm - Pm * Rb
print(f'[hessgb] normal forms done {time.time()-t0:.0f}s', flush=True)
# quadrics: for each nonpivot column j, Q_j(a) = sum_{k<=l} q_{kl} a_k a_l
quads = []
for j in nonpiv:
    cf = {}
    for k in range(K):
        for l in range(k, K):
            v = int(NF[k * K + l, j]) + (int(NF[l * K + k, j]) if l != k else 0)
            v %= p
            if v:
                cf[(k, l)] = v
    if cf:
        quads.append(cf)
# gauge fixing: unknowns (J, 0, X0*s) for s in S_1 of the right parity (C_0 has parity parsN[0]-1)
gauge_idx = []
gauge_s = (monomials(d - 6) if MIXED else monomials(d - 6, par)) if d >= 6 else []     # D = F s, parity(s) = parity(D)
for s in gauge_s:
    m = tuple(s[i] + (1 if i == 0 else 0) for i in range(NV))   # the monomial s*X0 of C_0
    for k, (_, _, tag) in enumerate(members):
        if tag == ('J', 0, m):
            gauge_idx.append(k)
assert len(gauge_idx) == len(gauge_s), (gauge_idx, len(gauge_s))
keep = [k for k in range(K) if k not in set(gauge_idx)]
kidx = {k: i for i, k in enumerate(keep)}
lines = []      # msolve syntax  a1*a2
lines_s = []    # Singular syntax a(1)*a(2)
for cf in quads:
    pairs_kept = [((kidx[k], kidx[l]), v) for (k, l), v in cf.items() if k in kidx and l in kidx]
    if pairs_kept:
        lines.append('+'.join(f'{v}*a{i+1}*a{j+1}' for (i, j), v in pairs_kept))
        lines_s.append('+'.join(f'{v}*a({i+1})*a({j+1})' for (i, j), v in pairs_kept))
nvar = len(keep)
base = f'{W.RES}/hessgb-d{d}-' + ('mixed' if MIXED else f'par{par}')
import re
with open(base + '.ms', 'w') as f:
    f.write(','.join(f'a{i+1}' for i in range(nvar)) + '\n' + str(p) + '\n' + ',\n'.join(lines) + '\n')
with open(base + '.sing', 'w') as f:
    f.write(f'ring r = {p}, (a(1..{nvar})), dp;\nideal I;\n')
    for ln in lines_s:
        f.write('I = I, ' + ln + ';\n')
    f.write('I = simplify(I, 2);\nint t = timer;\noption(redSB);\nideal G = std(I);\n')
    f.write('print("HESSGB dim = " + string(dim(G)) + "  vdim = " + string(vdim(G)) + "  secs = " + string(timer - t));\nquit;\n')
# linear-algebra pre-check: rank of the quadric coefficient matrix on the pair monomials of the kept unknowns
pairs = [(k, l) for k in keep for l in keep if k <= l]
pidx = {pr: i for i, pr in enumerate(pairs)}
flat = [0] * (len(quads) * len(pairs))
nq = 0
for r, cf in enumerate(quads):
    for (k, l), v in cf.items():
        if (k, l) in pidx:
            flat[r * len(pairs) + pidx[(k, l)]] = v
QM = nmod_mat(len(quads), len(pairs), flat, p)
qrank = QM.rank()
rec = {'task': 'hess-gb', 'd': d, 'chi': 'mixed' if MIXED else f'par{par}', 'tag': a.tag, 'p': p, 'members': K, 'nJ': nJ, 'n_extra': len(extras),
       'unknowns': nvar, 'dimR_2d': dimR, 'n_quadrics': len(lines), 'n_pair_monomials': len(pairs), 'quadric_rank': qrank,
       'full_rank_degree2': qrank == len(pairs), 'seconds_build': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
if qrank == len(pairs):
    rec['status'] = 'EXCLUDED-linear (all pair monomials in the ideal: only a = 0)'
    W.log_line(rec)
    print(rec['status'], flush=True)
    sys.exit(0)
if a.no_run:
    W.log_line({**rec, 'status': 'built'})
    sys.exit(0)
# run Singular with a timeout
t1 = time.time()
try:
    r = subprocess.run(['Singular', '-q', base + '.sing'], capture_output=True, text=True, timeout=a.timeout)
    out = (r.stdout + r.stderr)[-2000:]
    rec['singular'] = out.strip()[-500:]
    m = re.search(r'HESSGB dim = (-?\d+)\s+vdim = (-?\d+)', out)
    if m:
        rec['dim'] = int(m.group(1)); rec['vdim'] = int(m.group(2))
        rec['status'] = 'EXCLUDED (affine dim 0: only a=0 mod p, hence over Q-bar)' if rec['dim'] == 0 else f'positive-dimensional solution set (dim {rec["dim"]}): relaxation does not exclude'
    else:
        rec['status'] = 'singular-unparsed'
except subprocess.TimeoutExpired:
    rec['status'] = f'singular-timeout {a.timeout}s'
rec['seconds_gb'] = round(time.time() - t1, 1)
rec['seconds'] = round(time.time() - t0, 1)
W.log_line(rec)
print(json.dumps({k: v for k, v in rec.items() if k != 'singular'}), flush=True)
