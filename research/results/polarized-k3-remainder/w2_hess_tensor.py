"""Hessian relaxation at a Jacobian-type degree in the (U,V)-tensor frame (exact quotient, no big SVD):
LM random restarts on the sphere of cofactor coordinates modulo the 21 Euler directions (common to both parity sectors).
usage: python3 w2_hess_tensor.py D NRESTARTS SEED [--deflate-P8] [--maxiter 60] [--tag T] [--time-budget S]"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2lib as L, w2tensor
ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('n', type=int); ap.add_argument('seed', type=int)
ap.add_argument('--deflate-P8', action='store_true'); ap.add_argument('--maxiter', type=int, default=60); ap.add_argument('--tag', default='a')
ap.add_argument('--time-budget', type=float, default=1e9); ap.add_argument('--mixed', action='store_true', help='parity-free: independent cofactors for f_i(U) and f_i(V)')
args = ap.parse_args()
d = args.d
out = f'w2-hess-tensor-d{d}-{"mixed-" if args.mixed else ""}{args.tag}.jsonl'
T0 = time.time()
TR = w2tensor.TensorRelaxation(d, verbose=False, mixed=args.mixed)
# cofactor pieces: pure sectors -> by V-degree of the cofactor monomial; mixed -> by (block, "purity" degree)
if args.mixed:
    # block U: natural grading by V-degree; block V: by U-degree (pure means degree 0 in the other block's variables)
    vdeg = np.array([(m[3] + m[4] + m[5]) if b == 'U' else (m[0] + m[1] + m[2]) for (_, m, b) in TR.basis])
else:
    vdeg = np.array([m[3] + m[4] + m[5] for (_, m) in TR.basis])
pieces = {int(v): np.where(vdeg == v)[0] for v in sorted(set(vdeg.tolist()))}
# P_8-type subspace (cofactors in C[U]) in quotient coordinates
import scipy.linalg as sla
Pu = np.zeros((TR.K, len(pieces[0]))); Pu[pieces[0], np.arange(len(pieces[0]))] = 1.0
Qp8 = sla.orth(TR.P.T @ Pu, rcond=1e-9); r8 = Qp8.shape[1]
# natural scale of f on the unit sphere (mean over random points) for a relative residual
_rng = np.random.default_rng(12345); _fs = []
for _ in range(32):
    _u = _rng.standard_normal(TR.dimP); _u /= np.linalg.norm(_u); _fs.append(TR.value(_u)[0])
FSCALE = float(np.mean(_fs))
def dist_p8(u):
    w = u - Qp8 @ (Qp8.T @ u); return float(np.linalg.norm(w)), w
ap_power = 2.0  # deflation exponent: rho^power (f vanishes to 4th order in the normal directions of P8, so use dist^2 -> f/dist^4)
def deflate(u):
    rho, w = dist_p8(u); rho = max(rho, 1e-12)
    return rho ** ap_power, ap_power * rho ** (ap_power - 1) * (w / rho)
L.log_line(out, {'experiment': 'hessian-tensor-setup', 'd': d, 'mixed': args.mixed, 'K': TR.K, 'dimP': TR.dimP, 'n_coords': TR.nquot, 'euler_dirs': int(TR.Euler.shape[1]), 'f_scale': FSCALE,
                 'P8_quotient_dim': r8, 'pieces': {str(k): int(len(v)) for k, v in pieces.items()}, 'build_seconds': round(TR.build_time, 1),
                 'deflate_P8': args.deflate_P8, 'seed': args.seed, 'maxiter': args.maxiter})
rng = np.random.default_rng(args.seed)
fs = []; sols = 0; t_start = time.time()
for it in range(args.n):
    if time.time() - T0 > args.time_budget: break
    u0 = rng.standard_normal(TR.dimP); u0 /= np.linalg.norm(u0)
    res = L.lm_sphere(TR.value, u0, maxiter=args.maxiter, deflate=deflate if args.deflate_P8 else None, ftol=1e-28)
    u = res['u']; f = res['f']; fs.append(f)
    c = TR.P @ u
    prof = {str(k): float('%.3g' % float(np.sum(c[idx] ** 2))) for k, idx in pieces.items()}
    f_, r_, J_ = TR.value(u); gnorm = float(np.linalg.norm(J_.T @ r_ - (u @ (J_.T @ r_)) * u))
    rho, _ = dist_p8(u)
    rec = {'experiment': 'hessian-tensor-restart', 'd': d, 'mixed': args.mixed, 'seed': args.seed, 'restart': it, 'f': f, 'rel': f / FSCALE, 'f_eff': res['f_eff'], 'iters': res['iters'],
           'nfev': res['nfev'], 'grad_norm': gnorm, 'dist_to_P8': rho, 'cofactor_vdegree_profile': prof, 'deflate_P8': args.deflate_P8, 'hist_tail': [float('%.3g' % x) for x in res['hist'][-4:]]}
    if f < 1e-18:
        sols += 1
        sv = np.linalg.svd(J_, compute_uv=False)
        rec['tangent_nullity'] = int(np.sum(sv < 1e-8 * sv[0]))
        rec['c'] = [round(float(x), 10) for x in c]
        rec['class'] = 'inside P8 (H in J(f)_10 component)' if rho < 1e-6 else 'OFF P8: NEW'
    L.log_line(out, rec)
fs = np.array(fs)
summ = {'experiment': 'hessian-tensor-batch', 'd': d, 'mixed': args.mixed, 'seed': args.seed, 'deflate_P8': args.deflate_P8, 'restarts': int(len(fs)),
        'zero_residual': int(np.sum(fs < 1e-18)), 'f_quantiles': [float(x) for x in np.quantile(fs, [0, 0.1, 0.5, 0.9, 1])] if len(fs) else [],
        'seconds': round(time.time() - T0, 1)}
L.log_line(out, summ); print(json.dumps(summ), flush=True)
