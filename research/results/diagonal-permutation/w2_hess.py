"""Hessian relaxation for Jacobian-type denominators D = w(F) (Y-even), degree d:  A * (w^T Hess(F) w) in J(F).
Exact structure mod p: quadratic map  w -> [w^T H w] in R_{2d-2}^+  composed with mult by [A] into R_{2d+4}^-.
Numerics: random-restart Levenberg-Marquardt on the sphere for the (real) quadratic system; thousands of starts.
usage: python3 w2_hess.py d nstarts seed [--exact-dims]    (appends to w2-hess.jsonl)"""
import os, sys, time, json
os.environ['OPENBLAS_NUM_THREADS'] = '2'
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import monomials, NV, pderiv, shift, mul, add, scale
from flint import nmod_mat
LANE = os.path.dirname(os.path.abspath(__file__))
d = int(sys.argv[1]); nstarts = int(sys.argv[2]); seed = int(sys.argv[3])
PAR = [0, 0, 0, 1, 1, 1]
F, A, q, rec = core.load_instance('diagonal-permutation')
p = core.PRIMES[0]
t0 = time.time()
e = d - 5
wbasis = [(i, m) for i in range(NV) for m in monomials(e, PAR[i])]
nw = len(wbasis)
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
deg2 = 2 * e + 4           # degree of w^T H w  (= 2d - 6 ... careful: H has degree 4, w degree e: 2e+4 = 2d-6)
# reducer at degree deg2 (Y-even) and the multiplication-by-A map to degree deg2+6 (Y-odd)
red2 = jacobian._Reducer(Fp, deg2, p, 0)
red3 = jacobian._Reducer(Fp, deg2 + 6, p, 1)
std2 = red2.standard; sidx2 = {m: i for i, m in enumerate(std2)}
std3 = red3.standard; sidx3 = {m: i for i, m in enumerate(std3)}
print(f'd={d}: nw={nw}, R_{deg2}^+ dim {len(std2)}, R_{deg2+6}^- dim {len(std3)}  ({time.time()-t0:.1f}s)', flush=True)
# matrix of [A]: R_deg2^+ -> R_{deg2+6}^-  in standard monomial bases
MA = np.zeros((len(std3), len(std2)), dtype=np.int64)
for j, m in enumerate(std2):
    nf = red3.normal_form(shift(Ap, m, 1, p))
    for mm, v in nf.items():
        MA[sidx3[mm], j] = v
rankA = nmod_mat(len(std3), len(std2), MA.ravel().tolist(), p).rank()
print(f'rank [A]: R_{deg2}^+ -> R_{deg2+6}^- = {rankA} (kernel {len(std2)-rankA})  ({time.time()-t0:.1f}s)', flush=True)
# quadratic coefficients: for pairs (k,l): NF(m_k m_l H_{i_k i_l}) in std2 coordinates, then apply MA (mod p -> symmetric lift -> float)
# we work over floats with symmetric lifts of the mod-p normal forms: valid because the normal form of an integer polynomial
# mod p equals the reduction of the rational normal form; floats of symmetric lifts are only meaningful if the rational
# NF has small height -- instead, build the quadratic map exactly mod p and SOLVE NUMERICALLY only the float version of
# the *original* (unreduced) equations: Q_m(w) = coefficient extraction via a float complement of J.  So: float route.
import scipy.linalg as sla
Ff = F; Fif = [pderiv(F, i) for i in range(NV)]; Hf = [[pderiv(Fif[i], j) for j in range(NV)] for i in range(NV)]
mons3 = monomials(deg2 + 6, 1); idx3 = {m: j for j, m in enumerate(mons3)}
# J_{deg2+6}^- generators (float) and an orthonormal basis of the complement (rows), via SVD
gens = []
for i in range(NV):
    for m in monomials(deg2 + 1, (1 - PAR[i]) % 2):
        gens.append(shift(Fif[i], m))
G = np.zeros((len(gens), len(mons3)))
for r, g in enumerate(gens):
    for mm, v in g.items(): G[r, idx3[mm]] = float(v)
U, s, Vt = sla.svd(G, full_matrices=True, check_finite=False)
rk = int((s > 1e-9 * s[0]).sum())
Cc = Vt[rk:]                      # complement rows: (dim R_{deg2+6}^-) x n
print(f'float: dim J = {rk}, complement {Cc.shape[0]} (expect {len(std3)})  ({time.time()-t0:.1f}s)', flush=True)
# tensor T[m, k, l] = Cc . vec(A * m_k m_l H_{ik il})
T = np.zeros((Cc.shape[0], nw, nw))
for k in range(nw):
    ik, mk = wbasis[k]
    for l in range(k, nw):
        il, ml = wbasis[l]
        P = shift(mul(A, Hf[ik][il]), tuple(a + b for a, b in zip(mk, ml)))
        if not P: continue
        v = np.zeros(len(mons3))
        for mm, c in P.items(): v[idx3[mm]] = float(c)
        cv = Cc @ v
        T[:, k, l] += cv
        if l != k: T[:, l, k] += cv
print(f'tensor built ({time.time()-t0:.1f}s); equations {T.shape[0]}, unknowns {nw}', flush=True)
# gauge (spurious) directions: w = l*Euler, l in S_{d-6}^+  -> remove by projecting onto the orthogonal complement
idx = {wb: k for k, wb in enumerate(wbasis)}
sp = []
for l in (monomials(d - 6, 0) if d >= 6 else []):
    v = np.zeros(nw)
    for i in range(NV):
        ee = list(l); ee[i] += 1; v[idx[(i, tuple(ee))]] = 1.0
    sp.append(v)
if sp:
    S = np.array(sp).T; Us, ss, _ = np.linalg.svd(S, full_matrices=True); Pm = Us[:, S.shape[1]:]
else:
    Pm = np.eye(nw)
TP = np.einsum('mkl,ka,lb->mab', T, Pm, Pm)
na = Pm.shape[1]
rng = np.random.default_rng(seed)
def resid(u):
    un = u / np.linalg.norm(u)
    return np.einsum('mab,a,b->m', TP, un, un)
def jac(u):
    n = np.linalg.norm(u); un = u / n
    Jc = 2 * np.einsum('mab,b->ma', TP, un)
    # chain rule of normalization: d un / du = (I - un un^T)/n
    return (Jc - np.outer(Jc @ un, un)) / n
best = []
hits = 0
for s_ in range(nstarts):
    u0 = rng.standard_normal(na)
    res = least_squares(resid, u0, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
    val = float(np.linalg.norm(res.fun))
    best.append(val)
    if val < 1e-10:
        hits += 1
        np.save(os.path.join(LANE, f'w2-hess-hit-d{d}-{seed}-{s_}.npy'), Pm @ (res.x / np.linalg.norm(res.x)))
        print(f'!!! hit at start {s_}: residual {val:.3e}', flush=True)
    if (s_ + 1) % 200 == 0:
        print(f'{s_+1} starts: best {min(best):.3e} median {np.median(best):.3e} ({time.time()-t0:.0f}s)', flush=True)
best.sort()
# typical scale of the quadratic map on the sphere for reference
ref = float(np.median([np.linalg.norm(resid(rng.standard_normal(na))) for _ in range(50)]))
out = {'exp': 'hessian-relaxation', 'd': d, 'nw': nw, 'gauge': nw - na, 'neq': int(T.shape[0]), 'rankA_modp': int(rankA),
       'nstarts': nstarts, 'seed': seed, 'hits': hits, 'best10': best[:10], 'median': float(np.median(best)), 'random_scale': ref,
       'seconds': round(time.time() - t0, 1), 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
with open(os.path.join(LANE, 'w2-hess.jsonl'), 'a') as fh: fh.write(json.dumps(out) + '\n')
print(json.dumps(out), flush=True)
