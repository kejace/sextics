"""Certify the Q-nullity of the fixed-D operator [N|Q] for D = F60 (degree 6), needed to make the served-target
claim rigorous.  Strategy: (1) mod-p nullspace (expect 262); (2) check the 244 explicit exact kernels (Euler, F*V);
(3) choose a row subset of full p-rank, compute its EXACT integer nullspace with fmpz_mat, verify every vector
against ALL rows exactly.  If the exact nullity >= p-nullity, the Q-rank equals the p-rank and T_Q <= T_p follows.
Usage: python3 nullNQ_v.py p name"""
import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
from fractions import Fraction
p = int(sys.argv[1]); name = sys.argv[2]
F, A, q = parse_rec()
G = {e: int(c * 9973) for e, c in F.items()}          # integer form G = 9973 F
blocks = {0: {}, 2: {}, 4: {}, 6: {}}
for e, c in G.items(): blocks[e[3] + e[4] + e[5]][e] = c
cands = {'F60': blocks[0], 'F06': blocks[6], 'F60+F06': padd(blocks[0], blocks[6])}
D = cands[name]            # integer, = 9973 * block
d = 6
Gi = [pderiv(G, i) for i in range(6)]
Di = [pderiv(D, i) for i in range(6)]
G2 = pmul(G, G)
t0 = time.time()
# integer columns (exact):  -D N_i G_i + G (D d_i N_i - N_i D_i)   and  -G^2 Q   (scaling by 9973 powers does not change the kernel)
cols = []; names = []
pars = [1, 1, 1, 0, 0, 0]
for i in range(6):
    for m in monomials(d + 1, pars[i]):
        mm = {m: 1}
        col = pscale(pmul(D, pmul(mm, Gi[i])), -1)
        col = padd(col, pmul(G, padd(pmul(D, pderiv(mm, i)), pmul(mm, Di[i]), None, -1)))
        cols.append(col); names.append(('N', i, m))
nN = len(cols)
for m in monomials(2 * d - 6, 1):
    cols.append(pscale(pshift(G2, m), -1)); names.append(('Q', m))
C = len(cols)
rows = monomials(2 * d + 6, 1); ridx = {m: k for k, m in enumerate(rows)}
print(f'[{name}] integer system {len(rows)} x {C} built in {time.time()-t0:.1f}s', flush=True)
# mod-p matrix and nullspace
Mp = nmod_mat(len(rows), C, p)
for j, c in enumerate(cols):
    for e, v in c.items(): Mp[ridx[e], j] = v % p
X, nul = Mp.nullspace()
print(f'[{name}] mod {p}: rank {C - nul}, nullity {nul} ({time.time()-t0:.1f}s)', flush=True)
# pivot rows: rref of transpose
Mt = Mp.transpose()
Rt, rk = Mt.rref()
pivrows = []
for r in range(rk):
    j = 0
    while int(Rt[r, j]) == 0: j += 1
    pivrows.append(j)
assert len(pivrows) == C - nul
print(f'[{name}] {len(pivrows)} independent rows selected ({time.time()-t0:.1f}s)', flush=True)
# exact integer submatrix
sub = fmpz_mat(len(pivrows), C)
rset = {r: k for k, r in enumerate(pivrows)}
for j, c in enumerate(cols):
    for e, v in c.items():
        r = ridx[e]
        if r in rset: sub[rset[r], j] = v
t1 = time.time()
Xq, nulq = sub.nullspace()      # exact over Z/Q
print(f'[{name}] exact nullspace of {len(pivrows)} x {C} integer submatrix: nullity {nulq} ({time.time()-t1:.1f}s)', flush=True)
# verify each nullspace vector against ALL rows exactly
ok = True; bad = 0
colvecs = [[0] * C for _ in range(len(rows))]
for j, c in enumerate(cols):
    for e, v in c.items(): colvecs[ridx[e]][j] = v
for k in range(nulq):
    vec = [int(Xq[j, k]) for j in range(C)]
    for r in range(len(rows)):
        row = colvecs[r]
        s = 0
        for j in range(C):
            if row[j] and vec[j]: s += row[j] * vec[j]
        if s != 0:
            ok = False; bad += 1; break
print(f'[{name}] exact verification against all {len(rows)} rows: all {nulq} vectors in kernel: {ok} (bad {bad}) ({time.time()-t0:.1f}s)', flush=True)
rec = {'p': p, 'D': name, 'rows': len(rows), 'cols': C, 'rank_p': C - nul, 'nullity_p': nul, 'exact_nullity_lower_bound': nulq if ok else None,
       'Q_rank_equals_p_rank': ok and nulq == nul, 'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
with open('/tmp/hodge/results/dense-rational/verify/nullNQ_v.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
if ok and nulq == nul:
    log(f'served-target rigor for D = {name}: Q-rank of [N|Q] equals its mod-{p} rank (exact integer nullspace of dimension {nulq} verified against all rows)', 'CONFIRMED',
        f'rank_p = {C-nul}, exact nullity {nulq}', 'with this, dim T_Q(D) <= dim T_p(D); combined with dim T_p = 18 at two primes the "no primitive target served" statement becomes rigorous over C for this D')
else:
    log(f'served-target rigor for D = {name}: could not certify Q-rank = p-rank', 'UNVERIFIED', json.dumps(rec), '')
