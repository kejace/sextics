"""Job: full reduced-E2 LM search over the whole E1-admissible family (kernel + J(F)_d) at a non-Jacobian degree
(d >= 10 here), Y-parity pD, F-multiples projected out.  Float basis of {D : A D in J(F)} via SVD of the E1 system.

usage: python3 w2_nonjac_e2.py D PD NRESTARTS SEED [--maxiter 40] [--tag T]
Writes w2-nonjac-d{D}-p{PD}-{tag}.jsonl and caches the admissible basis in w2-admissible-d{D}-p{PD}.npz."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, time, json, argparse, os
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np, scipy.linalg as sla, scipy.sparse as sp
import w2lib as L
from fractions import Fraction as Fr
from hodge import core, fixed_d
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, monomials

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('pD', type=int); ap.add_argument('n', type=int); ap.add_argument('seed', type=int)
ap.add_argument('--maxiter', type=int, default=40); ap.add_argument('--tag', default='a'); ap.add_argument('--kernel-only', action='store_true')
args = ap.parse_args()
d, pD = args.d, args.pD
out = f'w2-nonjac-d{d}-p{pD}-{args.tag}.jsonl'
T0 = time.time()
F, A, q = L.load()
Fi = [pderiv(F, i) for i in range(NV)]
pars = [(1 - pD - (1 if i >= 3 else 0)) % 2 for i in range(NV)]
S_D = L.space(d, pD); S_B6 = L.space(d + 6, 1)
cache = f'w2-admissible-d{d}-p{pD}.npz'
if os.path.exists(cache):
    z = np.load(cache); Dbasis = z['Dbasis']; Nlift = z['Nlift']; ncols_info = z['info'].tolist()
    print('loaded cached admissible basis', Dbasis.shape, flush=True)
else:
    t0 = time.time()
    Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, pars[i])]
    nN = len(Ncols); nD = S_D.n
    M = np.zeros((S_B6.n, nN + nD))
    for j, (i, m) in enumerate(Ncols):
        for e, c in shift(Fi[i], m, 1).items():
            M[S_B6.idx[e], j] = float(c)
    for k, mD in enumerate(S_D.mons):
        for e, c in shift(A, mD, 1).items():
            M[S_B6.idx[e], nN + k] = -float(c)
    M *= 7.0
    print(f'E1 matrix {M.shape} built in {time.time()-t0:.1f}s; SVD...', flush=True)
    t1 = time.time()
    U_, s_, Vt = sla.svd(M, full_matrices=False, check_finite=False, lapack_driver='gesdd')
    tol = s_[0] * 1e-10
    rank = int(np.sum(s_ > tol))
    gap = (s_[rank - 1] / s_[rank]) if rank < len(s_) else float('inf')
    # null space: rows of Vt beyond rank, plus the (nN+nD - min(m,n)) trailing part if n > m
    NS = Vt[rank:, :].T                            # (nN+nD) x (min(m,n)-rank)
    if M.shape[1] > M.shape[0]:
        # complete: null_space handles this properly
        NS = sla.null_space(M, rcond=1e-10)
    print(f'SVD {time.time()-t1:.1f}s rank={rank} nullity={NS.shape[1]} gap={gap:.2e} s_min_kept={s_[rank-1]:.3e} s_max_dropped={(s_[rank] if rank < len(s_) else 0):.3e}', flush=True)
    Dpart = NS[nN:, :]; Npart = NS[:nN, :]
    Ud, sd, Vdt = sla.svd(Dpart, full_matrices=False, check_finite=False)
    rD = int(np.sum(sd > sd[0] * 1e-10))
    Dbasis = Ud[:, :rD]                                                   # nD x rD orthonormal
    # lifts: Dbasis = Dpart @ C  ->  C = Vdt[:rD].T / sd[:rD]  (minimum norm)
    C = (Vdt[:rD, :].T) / sd[:rD][None, :]
    Nlift = Npart @ C                                                      # nN x rD
    chk = np.linalg.norm(M[:, :nN] @ Nlift - (-M[:, nN:]) @ Dbasis) / np.linalg.norm(M[:, nN:] @ Dbasis)
    ncols_info = {'nN': nN, 'nD': nD, 'rank': rank, 'nullity': int(NS.shape[1]), 'admissible_dim': rD, 'lift_rel_err': float(chk), 'gap': float(gap)}
    np.savez(cache, Dbasis=Dbasis, Nlift=Nlift, info=np.array(json.dumps(ncols_info)))
    print('admissible dim', rD, 'lift check', chk, flush=True)
    L.log_line(out, {'experiment': 'nonjac-admissible-basis', 'd': d, 'pD': pD, **ncols_info, 'seconds': round(time.time() - t0, 1)})
Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, pars[i])]
nN = len(Ncols)
rD = Dbasis.shape[1]
# exact J(F)_d part (to measure the "kernel" complement) and F-multiples
Jb = L.jac_basis(d, pD)
Jm = np.array([S_D.vec(Dk) for _, _, Dk in Jb]).T
QJ = sla.orth(Jm)
L.log_line(out, {'experiment': 'nonjac-family-structure', 'd': d, 'pD': pD, 'admissible_dim': rD, 'J_d_dim': int(QJ.shape[1]),
                 'kernel_mod_J_dim': int(rD - QJ.shape[1]), 'J_subset_check': float(np.linalg.norm(QJ - Dbasis @ (Dbasis.T @ QJ)))})
Fm = np.array([S_D.vec(shift(F, h, 1)) for h in monomials(d - 6, pD)]).T
aF = Dbasis.T @ Fm                           # coordinates of F-multiples in the admissible basis
P = L.complement(aF, rD)
if args.kernel_only:
    # restrict to the orthogonal complement of J(F)_d inside the admissible space (pure kernel representatives)
    aJ = Dbasis.T @ QJ
    P = L.complement(np.concatenate([aF, aJ], axis=1), rD)
dimP = P.shape[1]
# family for ReducedE2: polys with float coefficients
family = []
for k in range(rD):
    Dk = S_D.poly(Dbasis[:, k])
    N0 = [{} for _ in range(NV)]
    for j, (i, m) in enumerate(Ncols):
        v = Nlift[j, k]
        if abs(v) > 1e-14:
            N0[i][m] = v
    family.append((Dk, N0))
S = L.ReducedE2(d, pD, family, use_gauge=False)
S.P = P; S.dimP = dimP
u0 = np.random.default_rng(0).standard_normal(dimP); u0 /= np.linalg.norm(u0)
t1 = time.time(); asm = S.assemble(S.P @ u0); r0 = S.residual(S.P @ u0, asm, with_jac=True); t_eval = time.time() - t1
L.log_line(out, {'experiment': 'nonjac-e2-setup', 'd': d, 'pD': pD, 'K': S.K, 'dimP': dimP, 'kernel_only': args.kernel_only, 'rows': S.S_B.n,
                 'red_rows': S.nred, 'nc': S.nc, 'nQ': S.S_Q.n, 'build_seconds': round(S.build_time, 1), 'eval_seconds': round(t_eval, 1),
                 'probe_f': r0['f'], 'probe_rel': r0['rel'], 'probe_rank': r0['rank'], 'maxiter': args.maxiter, 'seed': args.seed})
print(json.dumps({'K': S.K, 'dimP': dimP, 'red_rows': S.nred, 'nc': S.nc, 'eval_seconds': round(t_eval, 1), 'probe_rel': r0['rel']}), flush=True)
def resfun(u):
    a = S.P @ u
    asm = S.assemble(a); res = S.residual(a, asm, with_jac=True)
    return res['f'], res['r'], res['J'] @ S.P
rng = np.random.default_rng(args.seed)
best = None
for it in range(args.n):
    t0 = time.time()
    u0 = rng.standard_normal(dimP); u0 /= np.linalg.norm(u0)
    res = L.lm_sphere(resfun, u0, maxiter=args.maxiter, ftol=1e-30)
    u = res['u']; a = S.P @ u
    asm = S.assemble(a); fin = S.residual(a, asm)
    rel = fin['rel']
    Dv = Dbasis @ a
    kernel_frac = float(np.linalg.norm(Dv - QJ @ (QJ.T @ Dv)) / np.linalg.norm(Dv))
    rec = {'experiment': 'nonjac-e2-restart', 'd': d, 'pD': pD, 'seed': args.seed, 'restart': it, 'f': fin['f'], 'rel': rel, 'bnorm2': fin['bnorm2'],
           'rank_Mc': fin['rank'], 'iters': res['iters'], 'nfev': res['nfev'], 'kernel_component_frac': kernel_frac,
           'seconds': round(time.time() - t0, 1), 'hist': [float(x) for x in res['hist']]}
    if best is None or rel < best: best = rel
    rec['best_rel_so_far'] = best
    if rel < 1e-6:
        rec['D_coeffs'] = [float(x) for x in Dv]
    L.log_line(out, rec)
    print(json.dumps({k: rec[k] for k in ('restart', 'f', 'rel', 'iters', 'nfev', 'kernel_component_frac', 'seconds', 'best_rel_so_far')}), flush=True)
    if rel < 1e-10:
        np.save(f'w2-nonjac-d{d}-p{pD}-{args.tag}-candidate{it}.npy', Dv)
        print('LOW RESIDUAL CANDIDATE saved', flush=True)
print('done', round(time.time() - T0, 1), 's; best rel', best, flush=True)
