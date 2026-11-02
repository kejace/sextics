"""Reproduce a given restart of w2pt_search (same point seed and RNG) and refine it with many more iterations, logging the
objective trajectory; then run the candidate pipeline (independent points, rational reconstruction, exact confirmation).
usage: python3 w2_refine.py --d 9 --pD 1 --extra kerfam --sector 1 --seed 43 --wid 4 --restart 0 --tag j9op [--vanish X3]
       [--joint-iters 4000] [--maxiter 1500] [--margin 100]"""
import os, sys, time, json, argparse
os.environ['OPENBLAS_NUM_THREADS'] = '2'; os.environ['OMP_NUM_THREADS'] = '2'
import numpy as np
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2pt import *
from w2sym import w_sector_matrix, c_sector_matrix, symmetrize_pair
import w2pt_search as drv
from scipy.optimize import minimize
from scipy.linalg import null_space as sla_null_space
LANE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--pD', type=int, default=0); ap.add_argument('--extra', default='')
ap.add_argument('--sector', type=int, default=0); ap.add_argument('--vanish', default=''); ap.add_argument('--seed', type=int, required=True)
ap.add_argument('--wid', type=int, required=True); ap.add_argument('--restart', type=int, default=0); ap.add_argument('--tag', required=True)
ap.add_argument('--margin', type=int, default=100); ap.add_argument('--joint-iters', type=int, default=4000); ap.add_argument('--maxiter', type=int, default=1500)
args = ap.parse_args()
F, A, q, rec = core.load_instance('diagonal-permutation')
d = args.d
drv.args_d[0] = d; drv.args_pD[0] = args.pD
extra = drv.load_extra(F, A, args.extra.split(',')) if args.extra else None
if extra and args.sector:
    extra = [symmetrize_pair(D, N, args.sector) for (D, N) in extra]; extra = [(D, N) for (D, N) in extra if D]
wb0 = [(i, m) for i in range(NV) for m in monomials(d - 5, (args.pD + PAR[i]) % 2)]
K0 = len(wb0); Kx = len(extra) if extra else 0
pA = core.parity_of(A)
cols0 = [(i, j, m) for i in range(NV) for j in range(i + 1, NV) for m in monomials(d - 4, ((pA - args.pD - PAR[i]) % 2 - PAR[j]) % 2)]
nK = len(cols0)
if args.sector:
    dim_a = w_sector_matrix(wb0, args.sector).shape[1] + Kx; dim_c = c_sector_matrix(cols0, -args.sector).shape[1]
else:
    dim_a, dim_c = K0 + Kx, nK
r = (dim_a + dim_c + 1) // 2 + args.margin
pts = Points(F, r, seed=args.seed * 7919 + args.wid)
fam = jacobian_family(F, A, d, pts, extra=extra, pD=args.pD)
kb = KoszulBlock(F, A, d, pts, pD=args.pD)
subspace = None
allm = {}
for m in fam.members:
    for e in m['Dpoly']: allm.setdefault(e, len(allm))
Dmat = np.zeros((fam.K, len(allm)))
for k, m in enumerate(fam.members):
    for e, v in m['Dpoly'].items(): Dmat[k, allm[e]] = float(v)
ns = sla_null_space(Dmat.T)
if ns.shape[1]:
    fam.spurious = np.hstack([fam.spurious, ns]) if fam.spurious.size else ns
if args.vanish:
    spec = args.vanish
    if spec.startswith('X'): idxs, k = [0, 1, 2], int(spec[1:])
    elif spec.startswith('Y'): idxs, k = [3, 4, 5], int(spec[1:])
    else: raise ValueError(spec)
    bad = [j for e, j in allm.items() if sum(e[i] for i in idxs) < k]
    subspace = sla_null_space(Dmat[:, bad].T)
csector = None
if args.sector:
    wb = [m['w'] for m in fam.members if m['w'] is not None]
    Ws = w_sector_matrix(wb, args.sector)
    Wfull = np.zeros((fam.K, Ws.shape[1] + (fam.K - len(wb)))); Wfull[:len(wb), :Ws.shape[1]] = Ws
    for t in range(fam.K - len(wb)): Wfull[len(wb) + t, Ws.shape[1] + t] = 1.0
    if subspace is None: subspace = Wfull
    else:
        Nn = sla_null_space(np.hstack([Wfull, -subspace])); subspace = Wfull @ Nn[:Wfull.shape[1], :]
    csector = c_sector_matrix(kb.cols, -args.sector)
vp = VarPro(fam, kb, pts, csector=csector)
if subspace is not None:
    S = fam.spurious; Wn = subspace
    if S.size and S.shape[1]:
        Qs, _ = np.linalg.qr(S); Wn = Wn - Qs @ (Qs.T @ Wn)
    U, s_, _ = np.linalg.svd(Wn, full_matrices=False); rk = int(np.sum(s_ > 1e-10 * max(1.0, s_[0])))
    vp.P = U[:, :rk]; vp.dimP = rk
print(f'[refine] r={r} K={fam.K} dimP={vp.dimP} nK={kb.nK} nKc={vp.nKc}', flush=True)
# reproduce the start of restart `args.restart`
rng = np.random.default_rng(args.seed * 1000003 + args.wid * 101)
for it in range(args.restart + 1):
    u0 = rng.standard_normal(vp.dimP)
u0 /= np.linalg.norm(u0)
f0, rel0, _, _ = vp.value(vp.P @ u0)
print(f'[refine] start: f={f0:.4e} rel={rel0:.4e}', flush=True)
OUT = os.path.join(LANE, f'w2-refine-{args.tag}-w{args.wid}-r{args.restart}.jsonl')
def log(e):
    e['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    with open(OUT, 'a') as fh: fh.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:300], flush=True)
t0 = time.time()
# stage 1: the original 300 varpro iterations (to reproduce), with trajectory logging
traj = []
def cb(xk):
    if len(traj) % 25 == 0:
        f, g = vp.value_and_grad(xk); traj.append((len(traj), f)); log({'stage': 'varpro', 'iter': len(traj) * 25, 'f': f, 'gnorm': float(np.linalg.norm(g)), 'elapsed': round(time.time() - t0, 1)})
    else:
        traj.append((len(traj), None))
res = minimize(vp.value_and_grad, u0, jac=True, method='L-BFGS-B', callback=cb,
               options={'maxiter': args.maxiter, 'gtol': 1e-14, 'ftol': 1e-20, 'maxcor': 50})
u1 = res.x
a = vp.P @ u1; a /= np.linalg.norm(a)
f, rel, c, rho = vp.value(a)
log({'stage': 'varpro-done', 'iters': int(res.nit), 'f': f, 'rel': rel, 'rel_abs': f / vp.bnorm2(vp.P @ u0), 'message': str(res.message), 'elapsed': round(time.time() - t0, 1)})
np.save(os.path.join(LANE, f'w2-refine-{args.tag}-w{args.wid}-r{args.restart}-a.npy'), a)
# stage 2: joint continuation from (a, c) for many cheap iterations, then polish
cp0 = c if csector is None else np.linalg.lstsq(csector, c, rcond=None)[0]
x0 = np.concatenate([u1 / np.linalg.norm(u1), cp0])
res2 = minimize(vp.value_and_grad_joint, x0, jac=True, method='L-BFGS-B', options={'maxiter': args.joint_iters, 'gtol': 1e-16, 'ftol': 1e-22, 'maxcor': 50})
u2 = res2.x[:vp.dimP]
res3 = minimize(vp.value_and_grad, u2 / np.linalg.norm(u2), jac=True, method='L-BFGS-B', options={'maxiter': 300, 'gtol': 1e-16, 'ftol': 1e-22, 'maxcor': 50})
a2 = vp.P @ res3.x; a2 /= np.linalg.norm(a2)
f2, rel2, c2, rho2 = vp.value(a2)
log({'stage': 'joint+polish-done', 'joint_iters': int(res2.nit), 'f_joint': float(res2.fun), 'polish_iters': int(res3.nit), 'f': f2, 'rel': rel2, 'elapsed': round(time.time() - t0, 1)})
np.save(os.path.join(LANE, f'w2-refine-{args.tag}-w{args.wid}-r{args.restart}-a2.npy'), a2)
# independent point set check for the better of the two
abest = a2 if f2 < f else a
pts2 = Points(F, r, seed=args.seed * 31 + 999 + args.wid)
fam2 = jacobian_family(F, A, d, pts2, extra=extra, pD=args.pD); kb2 = KoszulBlock(F, A, d, pts2, pD=args.pD)
vp2 = VarPro(fam2, kb2, pts2, csector=csector)
f3, rel3, _, _ = vp2.value(abest)
log({'stage': 'independent-points', 'f': f3, 'rel': rel3})
# D polynomial (float) of the best point: coefficient profile and rational reconstruction attempts
Df = {}
for k, m in enumerate(fam.members):
    for e, v in m['Dpoly'].items(): Df[e] = Df.get(e, 0.0) + abest[k] * float(v)
vals = sorted((abs(v) for v in Df.values()), reverse=True)
log({'stage': 'D-profile', 'terms': len(Df), 'top5': vals[:5], 'median': vals[len(vals) // 2], 'min': vals[-1]})
from fractions import Fraction
emax = max(Df, key=lambda e: abs(Df[e])); sc = Df[emax]
for md in (10, 100, 1000, 10 ** 4, 10 ** 6):
    Dr = {e: Fraction(v / sc).limit_denominator(md) for e, v in Df.items()}
    err = max(abs(float(Dr.get(e, 0)) - v / sc) for e, v in Df.items())
    log({'stage': 'ratrecon', 'max_den': md, 'max_err': err})
    if err < 1e-6:
        Dr = {e: v for e, v in Dr.items() if v}
        status, data = confirm_exact(F, A, q, Dr)
        log({'stage': 'exact-confirm', 'max_den': md, 'status': status, 'info': {k: v for k, v in data.items() if k in ('rank', 'nullity', 'note', 'e1', 'error')}})
        if status == 'solved':
            w = {'d': d, 'D': core.to_str(Dr), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
            json.dump(w, open(os.path.join(LANE, 'witness-candidate.json'), 'w'), indent=1); log({'stage': 'WITNESS', 'file': 'witness-candidate.json'})
            break
print('REFINE DONE', flush=True)
