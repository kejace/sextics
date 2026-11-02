"""(a)/(b): degree 6 exclusion replicated exactly over Q with the L-relaxation (thread's DEGREE_SIX):
D = sum_i a_i P_i, a_i linear (36 params C), N' = T a + Koszul S (quadratic q_jk); B_D(S) in L := (G)+J(G)^2+(R_ijk),
so necessarily B_D(Ta) in L_12.  We compute exactly K = {P in Q^666 : T.P in L_12} (P = the 666 products C_ij C_kl)
blockwise in V-degree mod 6, and analyse it."""
import sys, time, json
sys.path.insert(0, '.')
from relax import *
from logutil import log
t0 = time.time()
Rx = Relaxation(6)
keys = Rx.quad_keys(); nq = len(keys); kidx = {k: i for i, k in enumerate(keys)}
print('quadratic monomials:', nq)
gens = Rx.L_generators(12)
print('L_12 generators:', len(gens))
# blocks by V-degree mod 6
blocks = {}
def vd(e): return (e[3] + e[4] + e[5]) % 6
rows_all = {}
for e in Rx.E2: rows_all.setdefault(e, None)
for _, _, g in gens:
    for e in g: rows_all.setdefault(e, None)
for e in rows_all: blocks.setdefault(vd(e), {})[e] = len(blocks.get(vd(e), {}))
Z_rows = []   # exact constraint rows (dict over quad index)
info = {}
for r in sorted(blocks):
    rows = blocks[r]; nr = len(rows)
    Lcols = []
    for _, _, g in gens:
        col = {rows[e]: c for e, c in g.items() if vd(e) == r}
        if col: Lcols.append(col)
    Tcols = []
    for k in keys:
        Tcols.append({rows[e]: v[k] for e, v in Rx.E2.items() if vd(e) == r and k in v})
    L = sparse_cols_to_fmpz(Lcols, nr); Tm = sparse_cols_to_fmpz(Tcols, nr)
    t1 = time.time()
    rL = L.rank(); rLT = sparse_cols_to_fmpz(Lcols + Tcols, nr).rank()
    # left null space of L: nullspace of L^T
    Y, ny = L.transpose().nullspace()
    assert ny == nr - rL, (ny, nr, rL)
    Ysub = fmpz_mat(nr, ny)
    for i in range(nr):
        for j in range(ny):
            Ysub[i, j] = Y[i, j]
    Z = Ysub.transpose() * Tm   # ny x nq  exact constraints
    rZ = Z.rank()
    info[r] = dict(rows=nr, Lcols=len(Lcols), rankL=rL, rank_LT=rLT, constraints=ny, rank_constraints=rZ, seconds=round(time.time() - t1, 1))
    print('block', r, info[r], flush=True)
    for i in range(ny):
        row = {j: int(Z[i, j]) for j in range(nq) if Z[i, j] != 0}
        if row: Z_rows.append(row)
# stack and compute kernel K exactly
Zs = fmpz_mat(len(Z_rows), nq)
for i, row in enumerate(Z_rows):
    for j, v in row.items(): Zs[i, j] = v
rankZ = Zs.rank()
K, nk = Zs.nullspace()
print('total constraint rows', len(Z_rows), 'rank', rankZ, 'dim K =', nk, flush=True)
# analyse K: which quadratic monomials are forced to zero (coordinate functional e_k in row space <=> all K vectors have 0 there)
forced = [keys[j] for j in range(nq) if all(K[j, i] == 0 for i in range(nk))]
pname = lambda pid: f'C{Rx.params[pid][1]}{Rx.params[pid][2].index(1)}'
print('forced-zero products:', len(forced), 'of', nq)
free = [keys[j] for j in range(nq) if not all(K[j, i] == 0 for i in range(nk))]
print('free products:', [(pname(a), pname(b)) for a, b in free][:60])
# expected: K spanned by (sum_i C_ii)^2-type? print kernel vectors
kv = []
for i in range(nk):
    v = {f'{pname(keys[j][0])}*{pname(keys[j][1])}': int(K[j, i]) for j in range(nq) if K[j, i] != 0}
    kv.append(v); print('kernel vector', i, v if len(v) < 80 else f'{len(v)} entries')
log(exp='a5_deg6_relaxation', rigor='exact over Q (fmpz_mat), L-relaxation, parity-free, UV frame', n_quadratic=nq,
    blocks=info, rankZ=rankZ, dimK=nk, forced_zero=len(forced), free_products=[(pname(a), pname(b)) for a, b in free],
    kernel_vectors=kv if nk <= 5 else f'{nk} vectors', seconds=round(time.time() - t0, 1))
json.dump(dict(info=info, rankZ=rankZ, dimK=nk, free=[(pname(a), pname(b)) for a, b in free], kernel=kv), open('a5_deg6.json', 'w'), indent=1, default=str)
