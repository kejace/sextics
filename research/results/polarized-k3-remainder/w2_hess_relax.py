"""Job: Hessian relaxation search  A*Hess_F(a,a) in J(F)_{2d}  over Jacobian-type D = sum a_i F_i.

usage: python3 w2_hess_relax.py D PD NRESTARTS SEED [--deflate] [--tag T]
  PD = 0 (Y-even), 1 (Y-odd), -1 (parity-free: both sectors together)
Writes w2-hess-d{D}-p{PD}-{tag}.jsonl (one line per batch + one line per zero-residual solution)."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2lib as L
from fractions import Fraction as Fr
from hodge import core
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, monomials

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('pD', type=int); ap.add_argument('n', type=int); ap.add_argument('seed', type=int)
ap.add_argument('--deflate', action='store_true'); ap.add_argument('--tag', default='a'); ap.add_argument('--maxiter', type=int, default=150)
ap.add_argument('--batch', type=int, default=100); ap.add_argument('--uvbideg', type=int, default=None)
args = ap.parse_args()
d, pD, n_restarts, seed = args.d, args.pD, args.n, args.seed
par_free = pD < 0
out = f'w2-hess-d{d}-p{pD}-{args.tag}.jsonl'
t0 = time.time()
F, A, q = L.load()
if args.uvbideg is None:
    H = L.HessianRelaxation(d, 0 if par_free else pD, par_free=par_free)
else:
    fam_uv = L.uv_bidegree_family(d, pD, args.uvbideg)
    H = L.HessianRelaxation(d, pD, cofactors=[c for (_, _, _, c) in fam_uv])
    H.uvpolys = [Dk for (_, Dk, _, _) in fam_uv]
    out = f'w2-hess-d{d}-p{pD}-uv{args.uvbideg}-{args.tag}.jsonl'
basis = H.basis; K = H.K
L.log_line(out, {'experiment': 'hessian-relaxation-setup', 'd': d, 'pD': pD, 'uvbideg': args.uvbideg, 'par_free': par_free, 'K': K,
                 'dimP': H.dimP, 'quotient_dim': H.nquot, 'jac_rank_2d': H.jac_rank, 'build_seconds': round(H.build_time, 1)})

# ---- known spurious family:  D = g * h  (g = f(U) - f(V) in J(F)_6), h quadric of the right parity
def g_poly():
    f = L.f_ternary()
    fU = L.uv_to_xy(f); fV = L.uv_to_xy(L.swapUV(f))
    return sub(fU, fV)
g = g_poly()
assert core.parity_of(g) == 1
S_D = L.space(d, None if par_free else pD)
if basis is not None:
    Mb = np.array([S_D.vec(Dk) for _, _, Dk in basis]).T   # monomials x K
else:
    Mb = np.array([S_D.vec(Dk) for Dk in H.uvpolys]).T
spur = []
hpars = [None] if par_free else [(pD + 1) % 2]        # g odd => h parity = pD + 1
for hp in (hpars if basis is not None else []):
    for h in monomials(d - 6, hp):
        v = S_D.vec(mul(g, {h: Fr(1)}))
        a, *_ = np.linalg.lstsq(Mb, v, rcond=None)
        assert np.linalg.norm(Mb @ a - v) < 1e-8 * max(1, np.linalg.norm(v)), 'g*h not in the family?'
        spur.append(H.P.T @ a)
Sspur = np.array(spur).T if spur else np.zeros((H.dimP, 0))
# orthonormalise
if Sspur.shape[1]:
    Qs, _ = np.linalg.qr(Sspur)
else:
    Qs = np.zeros((H.dimP, 0))
def dist_spur(u):
    w = u - Qs @ (Qs.T @ u)
    return float(np.linalg.norm(w)), w
# sanity: spurious points are exact solutions of the relaxation
for k in range(min(3, Qs.shape[1])):
    f_, r_, J_ = H.value(Qs[:, k])
    L.log_line(out, {'experiment': 'hessian-relaxation-spurious-check', 'd': d, 'pD': pD, 'k': k, 'f': f_,
                     'tangent_nullity': int(np.sum(np.linalg.svd(J_, compute_uv=False) < 1e-8 * max(1e-300, np.linalg.norm(J_))))})
def deflate(u):
    rho, w = dist_spur(u)
    rho = max(rho, 1e-12)
    return rho, w / rho

rng = np.random.default_rng(seed)
def resfun(u):
    return H.value(u)
sols = []
batch = {'experiment': 'hessian-relaxation-batch', 'd': d, 'pD': pD, 'uvbideg': args.uvbideg, 'par_free': par_free, 'deflate': args.deflate,
         'seed': seed, 'restarts': 0, 'zero_residual': 0, 'zero_residual_spurious': 0, 'zero_residual_new': 0,
         'best_f': None, 'best_f_nonspurious': None, 'fs': []}
def flush_batch():
    b = dict(batch); fs = np.array(b.pop('fs'))
    if fs.size:
        b['f_quantiles'] = [float(x) for x in np.quantile(fs, [0, 0.1, 0.5, 0.9, 1])]
        b['n_f_lt_1e-20'] = int(np.sum(fs < 1e-20)); b['n_f_lt_1e-10'] = int(np.sum(fs < 1e-10))
    b['seconds'] = round(time.time() - t0, 1)
    L.log_line(out, b)
for it in range(n_restarts):
    u0 = rng.standard_normal(H.dimP)
    res = L.lm_sphere(resfun, u0, maxiter=args.maxiter, deflate=deflate if args.deflate else None)
    u = res['u']; f = res['f']
    batch['restarts'] += 1; batch['fs'].append(f)
    batch['best_f'] = f if batch['best_f'] is None else min(batch['best_f'], f)
    rho, _ = dist_spur(u)
    if rho > 1e-6:
        batch['best_f_nonspurious'] = f if batch['best_f_nonspurious'] is None else min(batch['best_f_nonspurious'], f)
    if f < 1e-20:
        batch['zero_residual'] += 1
        f_, r_, J_ = H.value(u)
        sv = np.linalg.svd(J_, compute_uv=False)
        tn = int(np.sum(sv < 1e-8 * sv[0]))
        rec = {'experiment': 'hessian-relaxation-solution', 'd': d, 'pD': pD, 'uvbideg': args.uvbideg, 'par_free': par_free, 'deflate': args.deflate,
               'seed': seed, 'restart': it, 'f': f, 'dist_to_spurious': rho, 'tangent_nullity': tn, 'iters': res['iters']}
        if rho < 1e-6:
            batch['zero_residual_spurious'] += 1
            rec['class'] = 'spurious g*h'
        else:
            batch['zero_residual_new'] += 1
            rec['class'] = 'NEW'
            a = H.P @ u
            rec['a'] = [float(x) for x in a]
            # does D(a) have a known factor? test divisibility by g (deg 6) and by the partials F_i (deg 5), and F
            Dpoly = L.a_to_D(basis, a) if basis is not None else S_D.poly(Mb @ a)
            tests = {}
            for name, G in [('g', g)] + [(f'F_{i}', pderiv(F, i)) for i in range(NV)]:
                dg = L.core.degree(G)
                Gm = monomials(d - dg, None)
                cols = [S_D.vec(shift(G, m, 1)) for m in Gm]
                Mg = np.array(cols).T
                v = S_D.vec(Dpoly)
                x, *_ = np.linalg.lstsq(Mg, v, rcond=None)
                tests[name] = float(np.linalg.norm(Mg @ x - v) / max(np.linalg.norm(v), 1e-300))
            rec['divisibility_rel_residuals'] = tests
            try:
                rec['uv_bidegree_profile'] = L.uv_bidegree_profile({e: float(c) for e, c in Dpoly.items()}, d)
            except Exception as ex:
                rec['uv_bidegree_profile'] = str(ex)
            rec['D_sample'] = core.to_str({e: Fr(round(float(c), 6)) for e, c in sorted(Dpoly.items())[:12]})
        sols.append(rec)
        L.log_line(out, rec)
    if (it + 1) % args.batch == 0:
        flush_batch(); batch['fs'] = []; batch['restarts'] = 0; batch['zero_residual'] = batch['zero_residual_spurious'] = batch['zero_residual_new'] = 0
        batch['best_f'] = None; batch['best_f_nonspurious'] = None
if batch['restarts']:
    flush_batch()
print('done', len(sols), 'solutions', round(time.time() - t0, 1), 's', flush=True)
