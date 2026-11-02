import sys, time, json
from w2lib import *
t0 = time.time()
names = [(m, j) for m in monomials(2) for j in range(NV) if not (j == 5 and m[5] >= 1)]
n = len(names)
chi_u = [chi_add(chi(m), chi_neg(CHI_Z[j])) for m, j in names]
by_char = {}
for a in range(n):
    for b in range(a, n):
        by_char.setdefault(chi_add(chi_u[a], chi_u[b]), []).append((a, b))
quads = []
D8p = core.reduce_mod(D8, p)
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
        v8 = S8.nf(D8p); j0 = max(range(m), key=lambda j: 1 if v8[j] else 0); inv = pow(v8[j0], p - 2, p)
        cols = [[(cols[j][r] - v8[j] * inv % p * cols[j0][r]) % p for r in range(len(pairs))] for j in range(m) if j != j0]
    for col in cols:
        cf = {pairs[r]: v for r, v in enumerate(col) if v}
        if cf: quads.append(cf)
print('built', len(quads), 'quadrics', f'{time.time()-t0:.0f}s', file=sys.stderr)
json.dump({'p': p, 'names': [[list(m), j] for m, j in names], 'chi': [list(c) for c in chi_u], 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads]}, open('d7mixed.quads.local.json', 'w'))
pairs_all = [(a, b) for a in range(n) for b in range(a, n)]
pidx = {pr: i for i, pr in enumerate(pairs_all)}
flat = [0] * (len(quads) * len(pairs_all))
for r, cf in enumerate(quads):
    for pr, v in cf.items(): flat[r * len(pairs_all) + pidx[pr]] = v
QM = nmod_mat(len(quads), len(pairs_all), flat, p)
R2, rk2 = QM.rref()
mono2 = []
for i in range(rk2):
    j = 0
    while int(R2[i, j]) == 0: j += 1
    if all(int(R2[i, jj]) == 0 for jj in range(j + 1, len(pairs_all))): mono2.append(pairs_all[j])
nterms = [len(cf) for cf in quads]
from collections import Counter
sz = Counter(len(v) for v in by_char.values())
print(json.dumps({'n_quadrics': len(quads), 'rank': rk2, 'pairs': len(pairs_all), 'mean_terms': round(sum(nterms)/len(nterms),1), 'max_terms': max(nterms),
                  'n_I2_monomials': len(mono2), 'squares': [a for a,b in mono2 if a==b], 'products': [m for m in mono2 if m[0]!=m[1]][:40],
                  'pair_block_sizes': sorted(len(v) for v in by_char.values()), 'seconds': round(time.time()-t0,1)}))
