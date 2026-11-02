#!/usr/bin/env python3
"""w2_realkernel.py -- real (float64) orthonormal basis of the non-Jacobian part of K_d^even = {D : A D in J(F)} with E1 lifts,
for the joint stage-2 search (w2_jvp.py --extra).  Bases are those of w2_jvp.Basis (monomials sorted by code).

  Z   = orthonormal complement of col(G) in S_{d+6}^odd,  G = [F_i * S_{d+1}^{par_i}]  (left singular vectors of G)
  Kb  = null(Z^T M_A)  (M_A = multiplication by A : S_d^even -> S_{d+6}^odd)        -> dim K_d
  B   = orthonormal basis of Kb projected off J_d^even = span(m F_i)                 -> dim K_d / J_d
  n_k = least-squares E1 lift:  G n_k = M_A B_k   (exact up to float error, residual reported)
Output: w2_realkernel_d<d>.npz (B, n0..n5, singular-value diagnostics) + a JSON line in w2_realkernel.jsonl.
"""
import sys, os, json, time, argparse
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.linalg as sla, scipy.sparse as sp
from fractions import Fraction
from hodge import core
from hodge.core import NV, pderiv, mul
import w2_jvp as W

ap = argparse.ArgumentParser(); ap.add_argument('--d', type=int, default=10); ap.add_argument('--rJ', type=int, default=None)
a = ap.parse_args()
d = a.d
F, A, q, _ = core.load_instance(W.LANE)
Fi = [pderiv(F, i) for i in range(NV)]
t0 = time.time()
BD = W.Basis(d, 0); Bout = W.Basis(d + 6, 1)
parN = [1 if i < 3 else 0 for i in range(NV)]
BN = [W.Basis(d + 1, parN[i]) for i in range(NV)]
blocks = [W.mult_matrix(Fi[i], BN[i], Bout) for i in range(NV)]
G = sp.hstack(blocks).tocsr()
MA = W.mult_matrix(A, BD, Bout).toarray()
print(f'[realkernel d={d}] G {G.shape} nnz {G.nnz}, MA {MA.shape}; built {time.time()-t0:.1f}s', flush=True)
t1 = time.time()
U, s, Vt = sla.svd(G.toarray(), full_matrices=False, lapack_driver='gesdd', check_finite=False)
print(f'[realkernel] svd {time.time()-t1:.1f}s; s[0]={s[0]:.3e}', flush=True)
rJ = a.rJ or W.JVP._hilbert_odd(d + 6)
n = G.shape[0]
rank = n - rJ
print(f'[realkernel] rJ={rJ} rank={rank}: s[rank-1]={s[rank-1]:.3e} s[rank]={s[rank]:.3e} s[-1]={s[-1]:.3e}', flush=True)
Z = U[:, rank:]
ZM = Z.T @ MA                      # rJ x nD
Kb = sla.null_space(ZM, rcond=1e-9)
print(f'[realkernel] dim K_d = {Kb.shape[1]} (ZM sv: {np.linalg.svd(ZM, compute_uv=False)[-5:]})', flush=True)
# J_d^even inside S_d^even
Bw = [W.Basis(d - 5, 0 if l < 3 else 1) for l in range(NV)]
Jcols = []
for l in range(NV):
    for mu in Bw[l].mons:
        Jcols.append(BD.vec(mul({mu: Fraction(1)}, Fi[l])))
Jm = np.array(Jcols).T
UJ, sJ, _ = sla.svd(Jm, full_matrices=False)
rkJ = int(np.sum(sJ > 1e-9 * sJ[0]))
YJ = UJ[:, :rkJ]
Kc = Kb - YJ @ (YJ.T @ Kb)
Uc, sc, _ = sla.svd(Kc, full_matrices=False)
rc = int(np.sum(sc > 1e-7 * max(1.0, sc[0])))
B = Uc[:, :rc]
print(f'[realkernel] dim J_d = {rkJ}; nonJ part dim = {rc}; sc around cut: {sc[max(0,rc-3):rc+3]}', flush=True)
# E1 lifts via the SVD pseudo-inverse on the rank part
rhs = MA @ B                       # n x rc
Sinv = 1.0 / s[:rank]
Nlift = Vt[:rank, :].T @ (Sinv[:, None] * (U[:, :rank].T @ rhs))   # ncols x rc
res = G @ Nlift - rhs
relres = np.linalg.norm(res, axis=0) / np.linalg.norm(rhs, axis=0)
print(f'[realkernel] E1 lift residuals: max {relres.max():.2e} median {np.median(relres):.2e}', flush=True)
out = {'B': B}
off = 0
for i in range(NV):
    out[f'n{i}'] = Nlift[off:off + BN[i].n, :]; off += BN[i].n
out['s_gap'] = np.array([s[rank - 1], s[rank], s[-1]]); out['sc'] = sc; out['sJ'] = sJ; out['relres'] = relres
np.savez(f'w2_realkernel_d{d}.npz', **out)
rec = {'lane': W.LANE, 'item': 'realkernel', 'd': d, 'dimK': int(Kb.shape[1]), 'dimJ': rkJ, 'dim_nonJ': rc,
       's_rank_minus1': float(s[rank - 1]), 's_rank': float(s[rank]), 's_min': float(s[-1]),
       'lift_relres_max': float(relres.max()), 'seconds': round(time.time() - t0, 1), 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
open('w2_realkernel.jsonl', 'a').write(json.dumps(rec) + '\n')
print(json.dumps(rec), flush=True)
