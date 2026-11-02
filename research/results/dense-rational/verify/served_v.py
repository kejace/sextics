"""Independent served-target computation for a fixed degree-6 Y-even denominator D (q=1).
System: unknowns N (6 x S_7 with parities), Q (S_6^-), A' (S_6^-);  equation
  A' D^2 - D N.gradF + F (D div N - N.gradD) - F^2 Q = 0   in S_18^-.
dim T_p(D) = 226 + rank[N|Q] - rank[N|Q|A'].  Also records nullity of [N|Q] and tests explicit exact kernels
(Euler N = R Z, R in S_6^-: 226;  N = F V, V linear with parities: 18;  Q accordingly).
Usage: python3 served_v.py p name"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
p = int(sys.argv[1]); name = sys.argv[2]
F, A, q = parse_rec()
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
blocks = {0: {}, 2: {}, 4: {}, 6: {}}
for e, c in Fp.items(): blocks[e[3] + e[4] + e[5]][e] = c
cands = {'F60': blocks[0], 'F06': blocks[6], 'F42': blocks[2], 'F24': blocks[4], 'F60+F06': padd(blocks[0], blocks[6], p),
         'F60-F06': padd(blocks[0], blocks[6], p, -1), 'F42+F24': padd(blocks[2], blocks[4], p),
         'fermatX': {(6,0,0,0,0,0): 1, (0,6,0,0,0,0): 1, (0,0,6,0,0,0): 1}}
D = cands[name]
d = 6
Fi = [pderiv(Fp, i, p) for i in range(6)]
Di = [pderiv(D, i, p) for i in range(6)]
F2 = pmul(Fp, Fp, p); D2 = pmul(D, D, p)
t0 = time.time()
cols = []; names = []
pars = [1, 1, 1, 0, 0, 0]  # D Y-even: N_X Y-odd, N_Y Y-even
for i in range(6):
    for m in monomials(d + 1, pars[i]):
        mm = {m: 1}
        col = pscale(pmul(D, pmul(mm, Fi[i], p), p), p - 1, p)
        col = padd(col, pmul(Fp, padd(pmul(D, pderiv(mm, i, p), p), pmul(mm, Di[i], p), p, -1), p), p)
        cols.append(col); names.append(('N', i, m))
nN = len(cols)
for m in monomials(2 * d - 6, 1):
    cols.append(pscale(pshift(F2, m, 1, p), p - 1, p)); names.append(('Q', m))
nNQ = len(cols)
for m in monomials(6, 1):
    cols.append(pshift(D2, m, 1, p)); names.append(('A', m))
rows = monomials(2 * d + 6, 1)
print(f'[{name}] rows {len(rows)} cols {len(cols)} (N {nN}, Q {nNQ-nN}, A 226) built in {time.time()-t0:.1f}s', flush=True)
M = build_mat(cols, rows, p)
r_all = M.rank()
MNQ = build_mat(cols[:nNQ], rows, p)
r_nq = MNQ.rank()
dimT = 226 + r_nq - r_all
print(f'[{name}] rank[N|Q] = {r_nq}  rank[N|Q|A] = {r_all}  dim T_p(D) = {dimT}  nullity[N|Q] = {nNQ - r_nq}  ({time.time()-t0:.1f}s)', flush=True)
# explicit exact kernels of [N|Q]: Euler N = R Z (R in S_6^-), Q = 0 ; and N = F V (V linear, parities), Q = D div V - V.grad D
# count them and check independence mod p (they are exact over Q by the symbolic identities)
kcols = []
for Rm in monomials(6, 1):
    v = {}
    for i in range(6):
        m = [0]*6; m[i] = 1
        v[('N', i, tuple(a + b for a, b in zip(Rm, m)))] = 1
    kcols.append(v)
for i in range(6):
    for lm in monomials(1, (pars[i] - 0) % 2):   # V_i linear with parity pars[i] (F is Y-even)
        v = {}
        # N_i = F * lm
        for e, c in pshift(Fp, lm, 1, p).items(): v[('N', i, e)] = c
        # Q = D div V - V.grad D  = D * d_i(lm) - lm * D_i
        Qp = padd(pscale(D, pderiv({lm: 1}, i, p).get((0,)*6, 0), p), pshift(Di[i], lm, 1, p), p, -1)
        for e, c in Qp.items(): v[('Q', e)] = (v.get(('Q', e), 0) + c) % p
        kcols.append(v)
cidx = {nm: j for j, nm in enumerate(names[:nNQ])}
K = nmod_mat(nNQ, len(kcols), p)
for j, v in enumerate(kcols):
    for nm, c in v.items(): K[cidx[nm], j] = c % p
rk_known = K.rank()
# verify they are in the kernel mod p
Z = MNQ * K
iszero = all(int(Z[i, j]) == 0 for i in range(0, Z.nrows(), 97) for j in range(Z.ncols()))
print(f'[{name}] explicit kernel candidates: {len(kcols)} (226 Euler + 18 F*V), independent rank {rk_known}, in kernel (sampled rows): {iszero}', flush=True)
rec = {'p': p, 'D': name, 'rows': len(rows), 'cols': len(cols), 'rank_NQ': r_nq, 'rank_NQA': r_all, 'dimT_p': dimT,
       'nullity_NQ': nNQ - r_nq, 'explicit_exact_kernel_dim': rk_known, 'explicit_in_kernel': iszero, 'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
with open('/tmp/hodge/results/dense-rational/verify/served_v.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
# save the nullspace of [N|Q] for later exact-lifting attempts
if nNQ - r_nq > rk_known:
    X, nul = MNQ.nullspace()
    import pickle
    pickle.dump({'p': p, 'names': names[:nNQ], 'null': [[int(X[i, j]) for i in range(nNQ)] for j in range(nul)]},
                open(f'/tmp/hodge/results/dense-rational/verify/nullNQ_{name}_p{p}.pkl', 'wb'))
    print('saved nullspace', nul, flush=True)
