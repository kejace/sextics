#!/usr/bin/env python3
"""Sub-systems of a homogeneous quadric system (msolve format) by zeroing variable classes, and cheap
Macaulay-degree-k emptiness tests (VM excl agent).

Variables of quadsys_d7_{even,odd}_p*.ms are c{k} with k indexing receipt['basis'][k] = [component i, monomial m].
Classes: 'X0' = X-component with Y-free monomial (w_X^{(0)}), 'X2' = X-component with Y-quadric monomial (w_X^{(2)}),
         'Y1' = Y-component with X_aY_b monomial (w_Y^{(1)}) [even sector];
         'Y0' = Y-component Y-free (w_Y^{(0)}), 'Y2' = Y-component Y-quadric (w_Y^{(2)}), 'X1' = X-component X_aY_b (w_X^{(1)}) [odd sector].
Usage: python3 w2x_subsys.py file.ms --zero CLASS[,CLASS] [--macaulay k] [--write out.ms]
Emptiness certificate: if the degree-k Macaulay matrix of the restricted system has full column rank mod p
(dim (S/I)_k = 0) then the restricted ideal contains all degree-k forms over Q as well (its degree-k part reduces
into the mod-p one: (I_p)_k is contained in the reduction of the saturated lattice (I_Q)_k ∩ Z[c]_k), so the restricted
system has no nonzero complex solution."""
import sys, re, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2x_log import *
msf = sys.argv[1]
args = sys.argv[2:]
zero = []; k = None; outms = None
for i, a in enumerate(args):
    if a == '--zero': zero = args[i + 1].split(',')
    if a == '--macaulay': k = int(args[i + 1])
    if a == '--write': outms = args[i + 1]
rec = json.load(open(msf.replace('.ms', '_receipt.json')))
basis = rec['basis']
def cls(kvar):
    i, m = basis[kvar]; ydeg = m[3] + m[4] + m[5]
    return ('X' if i < 3 else 'Y') + str(ydeg)
lines = open(msf).read().strip().split('\n')
names = lines[0].split(','); p = int(lines[1])
eqs_txt = [l.rstrip(',').strip() for l in lines[2:] if l.strip()]
vidx = {n: j for j, n in enumerate(names)}
kill = set(n for n in names if cls(int(n[1:])) in zero)
keep = [n for n in names if n not in kill]
kpos = {n: j for j, n in enumerate(keep)}
term_re = re.compile(r'(\d+)\*(c\d+)(?:\*(c\d+)|\^2)?')
quad = []   # list of dict {(a,b): coeff} with a<=b indices into keep
for eq in eqs_txt:
    q = {}
    for c, v1, v2 in term_re.findall(eq):
        if v1 in kill or (v2 and v2 in kill): continue
        a = kpos[v1]; b = kpos[v2] if v2 else a
        if a > b: a, b = b, a
        q[(a, b)] = (q.get((a, b), 0) + int(c)) % p
    q = {key: v for key, v in q.items() if v}
    if q: quad.append(q)
n = len(keep)
# independent quadrics
mon2 = [(a, b) for a in range(n) for b in range(a, n)]; m2idx = {m: j for j, m in enumerate(mon2)}
M2 = nmod_mat(len(quad), len(mon2), p)
for r, q in enumerate(quad):
    for key, v in q.items(): M2[r, m2idx[key]] = v
rk2 = M2.rank()
info = {'exp': 'subsystem', 'ms_file': msf, 'zeroed_classes': zero, 'n_vars': n, 'n_quadrics_nonzero': len(quad), 'rank_quadrics': rk2,
        'dim_S2': len(mon2), 'dim_quotient_2': len(mon2) - rk2, 'prime': p}
print(json.dumps(info), flush=True)
if outms:
    with open(outms, 'w') as f:
        f.write(','.join(keep) + '\n'); f.write(f'{p}\n')
        f.write(',\n'.join('+'.join((f'{v}*{keep[a]}*{keep[b]}' if a != b else f'{v}*{keep[a]}^2') for (a, b), v in q.items()) for q in quad) + '\n')
    info['written'] = outms
if k:
    t0 = time.time()
    # degree-k Macaulay matrix: rows = m * q for monomials m of degree k-2 (as sorted index tuples), cols = degree-k monomials
    monk = list(itertools.combinations_with_replacement(range(n), k)); kidx = {m: j for j, m in enumerate(monk)}
    monkm2 = list(itertools.combinations_with_replacement(range(n), k - 2))
    rows = len(monkm2) * len(quad); cols = len(monk)
    print(f'Macaulay degree {k}: {rows} x {cols} ({time.time()-t0:.1f}s)', flush=True)
    M = nmod_mat(rows, cols, p)
    r = 0
    for mm in monkm2:
        for q in quad:
            for (a, b), v in q.items():
                col = kidx[tuple(sorted(mm + (a, b)))]
                M[r, col] = (int(M[r, col]) + v) % p
            r += 1
    rk = M.rank()
    info.update({'macaulay_degree': k, 'rows': rows, 'cols': cols, 'rank': rk, 'dim_quotient_k': cols - rk,
                 'projectively_empty_certified': rk == cols, 'seconds': round(time.time() - t0, 1)})
    print(json.dumps({kk: vv for kk, vv in info.items() if kk != 'ms_file'}), flush=True)
log2(**info)
