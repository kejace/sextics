"""Numerical (float64) screening of the E2-mod-J obstruction on the full E1-admissible family at degree d:
   D = sum_k lambda_k K_k + w(F),  w in S_{d-5}^6 modulo F-multiples;  ob(D) = [B_D(N_D)] in (S/J)_{2d}.
Builds the quadratic map r: R^nU -> R^{dim (S/J)_{2d}} (projection onto the orthogonal complement of J_{2d}, per Y-parity),
then runs Gauss-Newton on the unit sphere from random starts.  SCREENING ONLY -- a small residual is a candidate to be
confirmed exactly; a large minimum residual proves nothing.
Usage: python3 obnum_direct.py <d> [--kernel-rational f.json ...] [--starts N] [--seed s] [--sector a,b]"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
from common import *
import numpy as np
from fractions import Fraction
args = sys.argv[1:]
d = int(args[0]); kfiles = []; nstart = 100; seed = 0; sector = None
i = 1
while i < len(args):
    if args[i] == '--kernel-rational':
        i += 1
        while i < len(args) and not args[i].startswith('--'): kfiles.append(args[i]); i += 1
        continue
    if args[i] == '--starts': nstart = int(args[i + 1]); i += 2; continue
    if args[i] == '--seed': seed = int(args[i + 1]); i += 2; continue
    if args[i] == '--sector': sector = tuple(int(c) for c in args[i + 1].split(',')); i += 2; continue
    i += 1
F, A, q, rec = load()
def fl(P): return {e: float(Fraction(v)) for e, v in P.items()}
Ff = fl(F); Af = fl(A)
def fmul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(x + y for x, y in zip(e1, e2)); out[e] = out.get(e, 0.0) + c1 * c2
    return out
def fadd(P, Q, s=1.0):
    out = dict(P)
    for e, c in Q.items(): out[e] = out.get(e, 0.0) + s * c
    return out
def fpd(P, i):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; out[tuple(e2)] = out.get(tuple(e2), 0.0) + c * e[i]
    return out
def fshift(P, m, c=1.0):
    return {tuple(x + y for x, y in zip(e, m)): v * c for e, v in P.items()}
Fi = [fpd(Ff, i) for i in range(NV)]
def chi_w(i, m):
    return ((m[5] - (1 if i == 5 else 0)) % 2, (m[3] + m[4] - (1 if i in (3, 4) else 0)) % 2)
U = []
for i in range(NV):
    for m in monomials(d - 5):
        if sector is not None and chi_w(i, m) != sector: continue
        if i == 5 and m[5] >= 1: continue
        U.append(('w', i, m))
nW = len(U)
Ks = []; Ns = []
for kf in kfiles:
    for s in json.load(open(kf)):
        K = fl(core.parse(s))
        if sector is not None and klein(next(iter(K))) != sector: continue
        # E1 solve in floats: sum N_i F_i = A K
        parA = 1; parK = core.parity_of(core.parse(s))
        cols = []; polys = []
        for i in range(NV):
            pi = 0 if i < 3 else 1
            for m in monomials(d + 1, (parA + parK - pi) % 2):
                cols.append((i, m)); polys.append(fshift(Fi[i], m))
        rhs = fmul(Af, K)
        rows = {}
        for P in polys:
            for e in P: rows.setdefault(e, len(rows))
        for e in rhs: rows.setdefault(e, len(rows))
        M = np.zeros((len(rows), len(cols)))
        for j, P in enumerate(polys):
            for e, v in P.items(): M[rows[e], j] = v
        b = np.zeros(len(rows))
        for e, v in rhs.items(): b[rows[e]] = v
        x, res, rk, sv = np.linalg.lstsq(M, b, rcond=None)
        rel = np.linalg.norm(M @ x - b) / np.linalg.norm(b)
        print(f'[obnum] E1 float solve for kernel class: rel residual {rel:.2e}, rank {rk}', flush=True)
        N = [{} for _ in range(NV)]
        for (i, m), v in zip(cols, x):
            if abs(v) > 1e-14: N[i][m] = v
        Ks.append(K); Ns.append(N); U.append(('K', len(Ks) - 1))
nU = len(U)
print(f'[obnum d={d}] unknowns {nU} = {nW} J-part + {len(Ks)} kernel', flush=True)
Ds = []; Nsol = []
for u in U:
    if u[0] == 'w':
        _, i, m = u; D = fshift(Fi[i], m); N = [{} for _ in range(NV)]; N[i] = fshift(Af, m)
    else:
        D = Ks[u[1]]; N = Ns[u[1]]
    Ds.append(D); Nsol.append(N)
divN = []
for N in Nsol:
    s = {}
    for i in range(NV): s = fadd(s, fpd(N[i], i))
    divN.append(s)
dD = [[fpd(D, i) for i in range(NV)] for D in Ds]
def Bab(a, b):
    out = fmul(Ds[a], divN[b])
    for i in range(NV):
        if Nsol[b][i]: out = fadd(out, fmul(Nsol[b][i], dD[a][i]), -1.0)
    return out
# orthonormal bases of the left null space of J_{2d} per parity
t0 = time.time()
Z = {}; monsP = {}; idxP = {}
for par in (0, 1):
    mons = monomials(2 * d, par); idx = {e: k for k, e in enumerate(mons)}
    gens = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(2 * d - 5, (par - pi) % 2): gens.append(fshift(Fi[i], m))
    Jm = np.zeros((len(mons), len(gens)))
    for j, P in enumerate(gens):
        for e, v in P.items(): Jm[idx[e], j] = v
    Qf, Rf = np.linalg.qr(Jm, mode='complete')
    diag = np.abs(np.diag(Rf)); rank = int((diag > 1e-9 * diag.max()).sum())
    Z[par] = Qf[:, rank:].T.copy()     # (ring_dim x |mons|)
    monsP[par] = mons; idxP[par] = idx
    print(f'[obnum] J_{2*d} parity {par}: {len(mons)} mons, {len(gens)} gens, rank {rank}, ring dim {len(mons)-rank} (expected {expected_ring_dim(2*d, par)}) ({time.time()-t0:.0f}s)', flush=True)
nV = Z[0].shape[0] + Z[1].shape[0]; off = {0: 0, 1: Z[0].shape[0]}
T = np.zeros((nU, nU, nV))
t0 = time.time()
for a in range(nU):
    for b in range(a, nU):
        Q = Bab(a, b)
        if a != b: Q = fadd(Q, Bab(b, a))
        if not Q: continue
        par = (sum(next(iter(Q))[3:]) % 2)
        v = np.zeros(len(monsP[par]))
        for e, c in Q.items(): v[idxP[par][e]] = c
        r = Z[par] @ v
        T[a, b, off[par]:off[par] + len(r)] = r
        if a != b: T[b, a, off[par]:off[par] + len(r)] = r
    if a % 20 == 0: print(f'[obnum] row {a}/{nU} ({time.time()-t0:.0f}s)', flush=True)
# scale: symmetric T with T[a,b] = T[b,a] = q_ab; ob(w) = sum_{a<=b} w_a w_b q_ab = 1/2 sum_{a,b} w_a w_b T[a,b] + 1/2 sum_a w_a^2 T[a,a]
for a in range(nU): T[a, a] *= 2.0
T *= 0.5
# normalise columns of the residual
sc = np.sqrt((T ** 2).sum(axis=(0, 1))); sc[sc == 0] = 1.0; T /= sc
def resid(w): return np.einsum('abk,a,b->k', T, w, w)
def jac(w): return 2 * np.einsum('abk,b->ka', T, w)
rng = np.random.default_rng(seed)
results = []
t0 = time.time()
for s in range(nstart):
    w = rng.standard_normal(nU); w /= np.linalg.norm(w)
    for it in range(300):
        r = resid(w); Jm = jac(w)
        Aa = np.vstack([Jm, 2 * w[None, :]]); bb = np.concatenate([-r, [0.0]])
        step, *_ = np.linalg.lstsq(Aa, bb, rcond=None)
        w2 = w + step; w2 /= np.linalg.norm(w2)
        if np.linalg.norm(resid(w2)) < np.linalg.norm(r): w = w2
        else:
            w = w + 0.2 * step; w /= np.linalg.norm(w)
        if np.linalg.norm(resid(w)) < 1e-12: break
    f = float(np.linalg.norm(resid(w)))
    results.append((f, s, w.tolist() if f < 1e-6 else None))
    if s % 10 == 0 or f < 1e-6:
        print(f'[obnum] start {s}: resid {f:.3e} best {min(x[0] for x in results):.3e} ({time.time()-t0:.0f}s)', flush=True)
results.sort(key=lambda x: x[0])
tag = f'obnum_d{d}' + (f'_sec{sector[0]}{sector[1]}' if sector else '') + f'_seed{seed}'
json.dump({'d': d, 'nU': nU, 'nV': nV, 'sector': sector, 'nstart': nstart, 'seed': seed, 'kernel_files': kfiles,
           'best': [(f, s) for f, s, _ in results[:20]], 'candidates': [(f, s, w) for f, s, w in results if w is not None][:5],
           'unknowns': [list(u) if u[0] == 'K' else [u[0], u[1], list(u[2])] for u in U]}, open(f'{PK}/{tag}.json', 'w'))
log({'exp': 'obnum-screen', 'd': d, 'sector': sector, 'nU': nU, 'nV': nV, 'nstart': nstart, 'seed': seed,
     'best_resid': results[0][0], 'n_below_1e-6': sum(1 for x in results if x[0] < 1e-6), 'rigor': 'none (numerical screening)'})
