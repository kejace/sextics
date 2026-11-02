#!/usr/bin/env python3
"""Hessian relaxation (necessary E2 condition modulo J_{2d}) for CHARACTER-MIXED denominators, numerically.

For D(a) = sum a_k D_k (E1-admissible family with lifts N_k) every E2 solution satisfies
    T(a,a) := NF_{J_{2d}}( sum_{k,l} a_k a_l B_{D_k}(N_l) ) = 0        (Koszul fields and F*Q lie in J_{2d}).
This is fqlx's Hessian congruence (for Jacobian-type D it is A C^T Hess(F) C = 0 mod J).  We minimize
|T(a,a)|^2 on the unit sphere of the gauge-free parameter space (gauge F*S_{d-6} and D(a)=0 directions
projected out) with thousands of random L-BFGS restarts; every near-zero is then polished with the FULL
Koszul variable-projection system (w2lib.KoszulVarPro, same mode) and, if the full residual is tiny,
rationally reconstructed and re-solved exactly.

  python3 w2_hess.py --d 7 --parity 0 --restarts 2000 --procs 8 --tag hessd7e
  python3 w2_hess.py --d 7 --mixed  --restarts 1000 --procs 8 --tag hessd7m
  python3 w2_hess.py --d 8 --parity 0 --restarts 1000 --procs 8 --tag hessd8e
"""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
import multiprocessing as mp
import w2lib as W
from hodge.core import monomials, mul, add, sub, pderiv, shift

HS = None


class HessRelax:
    def __init__(self, members, d, parity=None, verbose=False, nthreads=8):
        t0 = time.time()
        self.d, self.parity = d, parity
        self.members = members
        K = len(members); self.K = K
        # pair vectors B_{D_k}(N_l)
        pair = {}
        rows = {}
        for l, (Dl, Nl, _) in enumerate(members):
            divN = {}
            for i in range(W.NV):
                if Nl[i]:
                    divN = add(divN, pderiv(Nl[i], i))
            for k, (Dk, Nk, _) in enumerate(members):
                Bv = mul(Dk, divN)
                for i in range(W.NV):
                    if Nl[i]:
                        Bv = sub(Bv, mul(Nl[i], pderiv(Dk, i)))
                if Bv:
                    pair[(k, l)] = Bv
                    for e in Bv:
                        rows.setdefault(e, len(rows))
        # J_{2d} generators on the same parity (rows are degree 2d)
        deg = 2 * d
        if parity is None:
            par_rows = None
        else:
            par_rows = 1  # A D^2 is Y-odd when D has pure parity
        gens = W.jgens(deg, None, par_rows)
        gp = [shift(W.Fif[i], m, 1.0) for i, m in gens]
        for P in gp:
            for e in P:
                rows.setdefault(e, len(rows))
        R = len(rows); self.rows = R
        G = np.zeros((R, len(gp)))
        for j, P in enumerate(gp):
            for e, v in P.items():
                G[rows[e], j] = v
        t1 = time.time()
        # left null space of G (orthonormal basis of (S_{2d}/J_{2d})^* in these coordinates)
        os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
        U, s, Vt = sla.svd(G, full_matrices=True, lapack_driver='gesdd', check_finite=False, overwrite_a=True)
        tol = s[0] * 1e-11
        rk = int(np.sum(s > tol))
        Uperp = U[:, rk:]
        self.dimR = Uperp.shape[1]
        self.rankJ = rk
        t2 = time.time()
        del G, U
        # T tensor (K, K, dimR) symmetrized
        T = np.zeros((K, K, self.dimR))
        for (k, l), Bv in pair.items():
            v = np.zeros(R)
            for e, c in Bv.items():
                v[rows[e]] = c
            w = Uperp.T @ v
            T[k, l] += 0.5 * w
            T[l, k] += 0.5 * w
        self.T = T
        self.Tflat = T.reshape(K * K, self.dimR)
        self.Fam, self.fmons = W.family_matrix(members)
        self.P, self.excluded = W.parameter_complement(self.Fam, self.fmons, d, None, parity)
        self.dimP = self.P.shape[1]
        # scale: typical |T| for normalization of reported values
        self.tscale = float(np.sqrt((T ** 2).sum() / max(K * K, 1)))
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[hess-relax] d={d} parity={parity} K={K} rows={R} J-gens={len(gp)} rank J_2d={rk} dimR={self.dimR} '
                  f'dimP={self.dimP} excluded={self.excluded} svd {t2-t1:.0f}s built {self.build_time:.0f}s', file=sys.stderr, flush=True)

    def value_and_grad(self, u):
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        M = np.tensordot(self.T, a, axes=([1], [0]))      # (K, dimR): M[k] = sum_l a_l T_kl
        v = a @ M                                          # T(a,a)
        f = float(v @ v)
        g = 4.0 * (M @ v)
        g_full = (g - float(g @ a) * a) / na
        self._nfev += 1
        return f, self.P.T @ g_full

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        M = np.tensordot(self.T, a, axes=([1], [0]))
        v = a @ M
        return float(v @ v)

    def search(self, n_restarts, seed, maxiter=400, gtol=1e-14, verbose=False):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-20, 'maxcor': 30})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            f = self.value(a)
            out.append({'f': f, 'frel': f / (self.tscale ** 2), 'a': a, 'iters': int(res.nit), 'restart': it, 'seed': seed})
            if verbose:
                print(f'[hess seed={seed}] restart {it}: f={f:.3e} frel={f/self.tscale**2:.3e} iters={res.nit} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
        out.sort(key=lambda c: c['f'])
        return out


def _worker(args):
    seed, n, maxiter = args
    out = HS.search(n, seed, maxiter=maxiter)
    return [{'f': c['f'], 'frel': c['frel'], 'a': c['a'].tolist(), 'iters': c['iters'], 'restart': c['restart'], 'seed': seed} for c in out]


def main():
    global HS
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--parity', type=int, default=None)
    ap.add_argument('--mixed', action='store_true')
    ap.add_argument('--restarts', type=int, default=500)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--maxiter', type=int, default=400)
    ap.add_argument('--tag', default='hess')
    ap.add_argument('--polish', type=float, default=1e-16, help='frel threshold for full-E2 polishing')
    ap.add_argument('--npolish', type=int, default=5)
    a = ap.parse_args()
    t0 = time.time()
    parity = None if a.mixed else a.parity
    base = {'task': 'hess-relax', 'd': a.d, 'chi': 'mixed' if parity is None else f'par{parity}', 'tag': a.tag}
    members, info, fam = W.build_family(a.d, None, parity=parity, mode='gens', verbose=True)
    rec = {**base, 'dimS': info['dimS'], 'dim_kernel': info['dim_kernel'], 'dimJ': info['dimJ'], 'n_extra': info['n_extra'],
           'n_members': info['n_members'], 'max_lift_res': info['max_lift_res']}
    HS = HessRelax(members, a.d, parity=parity, verbose=True, nthreads=a.procs)
    rec.update({'rows': HS.rows, 'rankJ2d': HS.rankJ, 'dimR_2d': HS.dimR, 'dimP': HS.dimP, 'excluded': HS.excluded,
                't_build': round(HS.build_time, 1), 'tscale': HS.tscale})
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    procs = max(1, min(a.procs, a.restarts))
    per = [a.restarts // procs + (1 if i < a.restarts % procs else 0) for i in range(procs)]
    jobs = [(a.seed + 1000 * i, n, a.maxiter) for i, n in enumerate(per) if n > 0]
    cands = []
    ctx = mp.get_context('fork')
    with ctx.Pool(processes=len(jobs)) as pool:
        for res in pool.imap_unordered(_worker, jobs):
            cands.extend(res)
    cands.sort(key=lambda c: c['f'])
    fr = np.array([c['frel'] for c in cands])
    rec.update({'restarts_done': len(cands), 'best_frel': float(fr.min()), 'median_frel': float(np.median(fr)),
                'n_frel_lt_1e-8': int((fr < 1e-8).sum()), 'n_frel_lt_1e-12': int((fr < 1e-12).sum()),
                'n_frel_lt_1e-16': int((fr < 1e-16).sum()),
                'hist_log10': {str(k): int(v) for k, v in zip(*np.unique(np.floor(np.log10(np.maximum(fr, 1e-300))).astype(int), return_counts=True))},
                'mean_iters': float(np.mean([c['iters'] for c in cands])), 'seconds_search': round(time.time() - t0, 1)})
    outp = f'{W.RES}/hess-{a.tag}-d{a.d}-{rec["chi"]}.json'
    json.dump({'rec': rec, 'members': [m[2] for m in members], 'best': cands[:20]}, open(outp, 'w'))
    rec['out'] = os.path.basename(outp)
    # polish the best relaxation minima with the full E2 system
    polished = []
    todo = [c for c in cands if c['frel'] < a.polish][:a.npolish]
    if todo:
        vp = W.KoszulVarPro(members, a.d, None, parity=parity, verbose=True)
        for c in todo:
            tp = time.time()
            out = vp.search(n_restarts=1, seed=0, maxiter=300, init=np.array(c['a']))
            best = out[0]
            pr = {'frel': c['frel'], 'full_rel': best['rel'], 'spurious': best['spurious'], 'iters': best['iters'],
                  'seconds': round(time.time() - tp, 1)}
            if best['rel'] < 1e-10 and not best['spurious']:
                pr['recon'] = W.try_reconstruct_and_confirm(vp, best['a'], verbose=True, tag=f'{a.tag}-d{a.d}')
                if any(x.get('witness_path') for x in pr['recon']):
                    rec['WITNESS'] = True
            polished.append(pr)
            print('[polish]', pr, flush=True)
    rec['polished'] = polished
    rec['seconds'] = round(time.time() - t0, 1)
    rec['status'] = 'done'
    W.log_line(rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('hist_log10',)}), flush=True)


if __name__ == '__main__':
    main()
