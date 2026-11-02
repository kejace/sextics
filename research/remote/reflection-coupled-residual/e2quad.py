"""Exact E2-mod-J obstruction on the E1 family D(g) = D0 + sum g_i F_i (mod p):
   ob(g) = [ D(g) div N(g) - N(g)(D(g)) ] in R_{2d}^odd,  N(g) = N0 + A g   (exact, quadratic in g).
Reports ob(0), the linearised system, and runs a Newton heuristic (random square projections, random starts)
restricted to the (Z/2)^2-invariant sector.  Usage: python3 e2quad.py <kernel-json> <idx> <n_starts> <n_iter> [seed]"""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv, mul, add, sub, scale
from flint import nmod_mat
P = core.PRIMES[0]
F, A, q, rec = core.load_instance('reflection-coupled-residual')
D0 = core.parse(json.load(open(sys.argv[1]))[int(sys.argv[2])])
nstarts = int(sys.argv[3]); niter = int(sys.argv[4]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
random.seed(seed)
d = core.degree(D0)
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P); Dp = core.reduce_mod(D0, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]
par = lambda i: 0 if i < 3 else 1
t0 = time.time()
# E1 particular solution N0
monsT = monomials(d + 6, 1); idxT = {m: i for i, m in enumerate(monsT)}
cols = [(i, m) for i in range(NV) for m in monomials(d + 1, (1 - par(i)) % 2)]
R, C = len(monsT), len(cols)
flat = [0] * (R * (C + 1))
for j, (i, m) in enumerate(cols):
    for e, v in shift(Fi[i], m, 1, P).items(): flat[idxT[e] * (C + 1) + j] = v
AD = mul(Ap, Dp, P)
for e, v in AD.items(): flat[idxT[e] * (C + 1) + C] = v
Mr, rk = nmod_mat(R, C + 1, flat, P).rref()
N0 = [dict() for _ in range(NV)]
for r in range(rk):
    j = 0
    while int(Mr[r, j]) == 0: j += 1
    assert j < C, 'E1 inconsistent'
    val = int(Mr[r, C])
    if val: N0[cols[j][0]][cols[j][1]] = val
del Mr, flat
acc = {}
for i in range(NV): acc = add(acc, mul(N0[i], Fi[i], P), P)
assert sub(acc, AD, P) == {}
print(f'E1 solved ({time.time()-t0:.0f}s)', flush=True)
# quotient map S_{2d}^odd -> R_{2d}^odd
mons2 = monomials(2 * d, 1); idx2 = {m: i for i, m in enumerate(mons2)}; C2 = len(mons2)
Jrows = [shift(Fi[i], m, 1, P) for i in range(NV) for m in monomials(2 * d - 5, (1 - par(i)) % 2)]
flat = [0] * (len(Jrows) * C2)
for r, poly in enumerate(Jrows):
    base = r * C2
    for e, v in poly.items(): flat[base + idx2[e]] = v
W, nul = nmod_mat(len(Jrows), C2, flat, P).nullspace()
del flat
Wm = nmod_mat(C2, nul, [int(W[i, k]) for i in range(C2) for k in range(nul)], P); del W
print(f'dim R_{2*d}^odd = {nul} ({time.time()-t0:.0f}s)', flush=True)
def cls(poly):
    row = [0] * C2
    for e, c in poly.items(): row[idx2[e]] = c % P
    v = nmod_mat(1, C2, row, P) * Wm
    return [int(v[0, k]) for k in range(nul)]
# sector-restricted g basis: g_i F_i in sector (Y2-even, Y0+Y1-even)
def sector(m): return (m[5] % 2, (m[3] + m[4]) % 2)
Fsec = [sector(next(iter(Fi[i]))) for i in range(NV)]
gbasis = [(i, m) for i in range(NV) for m in monomials(d - 5, par(i))
          if ((sector(m)[0] + Fsec[i][0]) % 2, (sector(m)[1] + Fsec[i][1]) % 2) == (0, 0)]
ng = len(gbasis)
print(f'g basis (invariant sector): {ng} of {sum(1 for i in range(NV) for m in monomials(d-5, par(i)))}', flush=True)
def polys_of(g):
    gi = [dict() for _ in range(NV)]
    for k, (i, m) in enumerate(gbasis):
        if g[k]: gi[i][m] = g[k] % P
    return gi
def ob(g):
    gi = polys_of(g)
    T = {}
    for i in range(NV): T = add(T, mul(gi[i], Fi[i], P), P)
    D = add(Dp, T, P)
    N = [add(N0[i], mul(Ap, gi[i], P), P) for i in range(NV)]
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N[i], i, P), P)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N[i], pderiv(D, i, P), P), P)
    return cls(sub(mul(D, divN, P), ND, P))
ob0 = ob([0] * ng)
print('ob(D0) == 0 ?', not any(ob0), ' nonzero coords:', sum(1 for x in ob0 if x), flush=True)
# linearisation at 0: columns L_k = ob(e_k) - ob(0) - Q(e_k,e_k); get via two evaluations: ob(e_k), ob(-e_k): L_k = (ob(e_k)-ob(-e_k))/2
def jac(g):
    base = ob(g)
    colsJ = []
    for k in range(ng):
        gp = list(g); gp[k] = (gp[k] + 1) % P
        gm = list(g); gm[k] = (gm[k] - 1) % P
        a = ob(gp); b = ob(gm)
        colsJ.append([((a[r] - b[r]) * pow(2, -1, P)) % P for r in range(nul)])
    return base, colsJ
t1 = time.time()
base, L = jac([0] * ng)
print(f'Jacobian at 0 computed ({time.time()-t1:.0f}s)', flush=True)
Lm = nmod_mat(nul, ng, [L[k][r] for r in range(nul) for k in range(ng)], P)
rL = Lm.rank()
Laug = nmod_mat(nul, ng + 1, [L[k][r] if k < ng else (-ob0[r]) % P for r in range(nul) for k in range(ng + 1)], P)
print(f'linearised: rank L = {rL}, rank [L|-ob0] = {Laug.rank()}, linear-solvable = {rL == Laug.rank()}', flush=True)
json.dump({'d': d, 'nul': nul, 'ng': ng, 'ob0_zero': not any(ob0), 'rankL': rL, 'rank_aug': Laug.rank()},
          open(f'e2quad_d{d}_info.json', 'w'))
# Newton heuristic on random square projections
found = []
for s in range(nstarts):
    g = [random.randrange(P) for _ in range(ng)]
    Rp = [[random.randrange(P) for _ in range(nul)] for _ in range(ng)]   # ng x nul projection
    for it in range(niter):
        base, J = jac(g)
        nz = sum(1 for x in base if x)
        print(f'  start {s} iter {it}: ||ob|| nonzero coords = {nz}', flush=True)
        if nz == 0:
            found.append(g); print('!!!!!!!! FULL SOLUTION mod p FOUND', g, flush=True); break
        # square system: Rp * J * delta = -Rp * base
        RJ = [[sum(Rp[a][r] * J[k][r] for r in range(nul)) % P for k in range(ng)] for a in range(ng)]
        rb = [(-sum(Rp[a][r] * base[r] for r in range(nul))) % P for a in range(ng)]
        M = nmod_mat(ng, ng + 1, [RJ[a][k] if k < ng else rb[a] for a in range(ng) for k in range(ng + 1)], P)
        Mr, rr = M.rref()
        if rr < ng:
            print('   singular projected Jacobian, rank', rr); break
        delta = [int(Mr[a, ng]) for a in range(ng)]
        g = [(g[k] + delta[k]) % P for k in range(ng)]
json.dump({'found': found}, open(f'e2quad_d{d}_newton_seed{seed}.json', 'w'))
print(f'done, {len(found)} solutions, total {time.time()-t0:.0f}s')
