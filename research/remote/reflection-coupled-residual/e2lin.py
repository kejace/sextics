"""Second-order (E2 mod J) screen of the whole E1 family D = D0 + sum g_i F_i, g_i in S_{d-5}.
ob(D) := [D div N - N(D)] in R_{2d}^odd is well-defined (independent of the E1 solution N) and affine-linear in g.
Solve ob(D0) + L(g) = 0 mod p.  Usage: python3 e2lin.py <kernel-json> [index]"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv, mul, add, sub, scale
from flint import nmod_mat
P = core.PRIMES[0]
F, A, q, rec = core.load_instance('reflection-coupled-residual')
D0 = core.parse(json.load(open(sys.argv[1]))[int(sys.argv[2]) if len(sys.argv) > 2 else 0])
d = core.degree(D0)
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P); Dp = core.reduce_mod(D0, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]
par = lambda i: 0 if i < 3 else 1
t0 = time.time()
# ---- step 1: particular solution N0 of sum N_i F_i = A D0, deg N_i = d+1, parity: N_i Y-odd for i<3, Y-even for i>=3
monsT = monomials(d + 6, 1); idxT = {m: i for i, m in enumerate(monsT)}
cols = []  # (i, m)
for i in range(NV):
    for m in monomials(d + 1, (1 - par(i)) % 2):
        cols.append((i, m))
R, C = len(monsT), len(cols)
flat = [0] * (R * (C + 1))
for j, (i, m) in enumerate(cols):
    for e, v in shift(Fi[i], m, 1, P).items():
        flat[idxT[e] * (C + 1) + j] = v
AD = mul(Ap, Dp, P)
for e, v in AD.items(): flat[idxT[e] * (C + 1) + C] = v
M = nmod_mat(R, C + 1, flat, P)
Mr, rk = M.rref()
N0 = [dict() for _ in range(NV)]
consistent = True
for r in range(rk):
    j = 0
    while int(Mr[r, j]) == 0: j += 1
    if j == C: consistent = False; break
    val = int(Mr[r, C])
    if val:
        i, m = cols[j]; N0[i][m] = val
print(f'E1 solve: rank {rk}, consistent={consistent} ({time.time()-t0:.0f}s)', flush=True)
assert consistent
# check
chk = sub(sum((mul(N0[i], Fi[i], P) for i in range(NV)), {}), AD, P) if False else None
acc = {}
for i in range(NV): acc = add(acc, mul(N0[i], Fi[i], P), P)
assert sub(acc, AD, P) == {}, 'E1 check failed'
# ---- step 2: rho(D, N) = D div N - N(D)
def rho(D, N):
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N[i], i, P), P)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N[i], pderiv(D, i, P), P), P)
    return sub(mul(D, divN, P), ND, P)
rho0 = rho(Dp, N0)
# ---- step 3: quotient map S_{2d}^odd -> R_{2d}^odd via W = nullspace(Jrows)
mons2 = monomials(2 * d, 1); idx2 = {m: i for i, m in enumerate(mons2)}
Jrows = []
for i in range(NV):
    for m in monomials(2 * d - 5, (1 - par(i)) % 2):
        Jrows.append(shift(Fi[i], m, 1, P))
RJ, C2 = len(Jrows), len(mons2)
flat = [0] * (RJ * C2)
for r, poly in enumerate(Jrows):
    base = r * C2
    for e, v in poly.items(): flat[base + idx2[e]] = v
t1 = time.time()
MJ = nmod_mat(RJ, C2, flat, P)
W, nul = MJ.nullspace()   # C2 x nul (padded); columns span annihilator of rowspace(J)
print(f'J_{2*d}^odd: {RJ} rows x {C2} cols, dim R_{2*d}^odd = {nul} ({time.time()-t1:.0f}s)', flush=True)
Wm = nmod_mat(C2, nul, [int(W[i, k]) for i in range(C2) for k in range(nul)], P)
def cls(poly):
    v = nmod_mat(1, C2, [0] * C2, P)
    row = [0] * C2
    for e, c in poly.items(): row[idx2[e]] = c % P
    return nmod_mat(1, C2, row, P) * Wm   # 1 x nul
ob0 = cls(rho0)
print('ob(D0) zero?', all(int(ob0[0, k]) == 0 for k in range(nul)), flush=True)
# ---- step 4: linear map L(g): g = m e_i  (D' = D + m F_i, N' = N0 + A m e_i)
gbasis = []
for i in range(NV):
    for m in monomials(d - 5, (0 - par(i)) % 2):   # g_i F_i must be Y-even: g_i parity = par(i)
        gbasis.append((i, m))
Lcols = []
t2 = time.time()
for (i, m) in gbasis:
    g = {m: 1}
    T = shift(Fi[i], m, 1, P)                    # T = m F_i
    Ag = mul(Ap, g, P)
    # L = D * d_i(A g) - N0(T) - A g d_i D   (terms linear in g, modulo J)
    term = mul(Dp, pderiv(Ag, i, P), P)
    N0T = {}
    for j in range(NV): N0T = add(N0T, mul(N0[j], pderiv(T, j, P), P), P)
    term = sub(term, N0T, P)
    term = sub(term, mul(Ag, pderiv(Dp, i, P), P), P)
    Lcols.append(cls(term))
print(f'L columns: {len(gbasis)} ({time.time()-t2:.0f}s)', flush=True)
# ---- step 5: solve ob0 + sum g_k L_k = 0  ->  [L | -ob0] rref
ng = len(gbasis)
flat = [0] * (nul * (ng + 1))
for k, Lc in enumerate(Lcols):
    for r in range(nul): flat[r * (ng + 1) + k] = int(Lc[0, r])
for r in range(nul): flat[r * (ng + 1) + ng] = (-int(ob0[0, r])) % P
S = nmod_mat(nul, ng + 1, flat, P)
Sr, rs = S.rref()
rankL = nmod_mat(nul, ng, [flat[r * (ng + 1) + k] for r in range(nul) for k in range(ng)], P).rank()
print(f'E2-mod-J system: {nul} equations x {ng} unknowns; rank L = {rankL}; rank [L|ob] = {rs}; SOLVABLE = {rankL == rs}; solution family dim = {ng - rankL}', flush=True)
out = {'d': d, 'dimR2d_odd': nul, 'ng': ng, 'rankL': rankL, 'rank_aug': rs, 'solvable': rankL == rs, 'ob0_zero': all(int(ob0[0, k]) == 0 for k in range(nul))}
if rankL == rs:
    g = [0] * ng
    for r in range(rs):
        j = 0
        while int(Sr[r, j]) == 0: j += 1
        if j < ng: g[j] = int(Sr[r, ng])
    Dnew = dict(Dp)
    for k, (i, m) in enumerate(gbasis):
        if g[k]: Dnew = add(Dnew, scale(shift(Fi[i], m, 1, P), g[k], P), P)
    out['g'] = {f'{i}:{m}': g[k] for k, (i, m) in enumerate(gbasis) if g[k]}
    out['Dnew_modp'] = core.to_str(Dnew)
    print('particular g:', out['g'])
json.dump(out, open(f'e2lin_d{d}.json', 'w'), indent=1)
print(f'total {time.time()-t0:.0f}s')
