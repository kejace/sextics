"""Hessian relaxation (E2 mod J) for Jacobian-type denominators D = w.grad F at degree d:
    necessary condition   A * (w^T Hess(F) w)  in  J(F)_{2d}   (both parity sectors, or mixed),
equivalently E @ vec(w^T H w) = 0 with E = C_{2d} @ ShiftA, C_{2d} = orthonormal complement of J_{2d} in S_{2d}.
Quadratic system T[e,k,l] w_k w_l = 0; Euler gauge w -> w + ell z removed (w restricted to the orthogonal complement
of {ell z}); normalization |w| = 1 (real) or L(w) = 1 (complex).  Levenberg-Marquardt with analytic Jacobian,
many random restarts, multiprocessing.  Residual reported scale-free: |T(w,w)| / |w|^2 (real) ; also relative to
the mean over random points.
Usage: python3 w2_relax.py --d 7 --sector even|odd|mixed --restarts N --seed S --workers W [--complex] [--batch 200]
"""
import sys, os, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--sector', default='even')
ap.add_argument('--restarts', type=int, default=1000); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--workers', type=int, default=4); ap.add_argument('--complex', action='store_true')
ap.add_argument('--batch', type=int, default=200); ap.add_argument('--maxnfev', type=int, default=400)
ap.add_argument('--threads', type=int, default=1)
args = ap.parse_args()
os.environ['OMP_NUM_THREADS'] = str(args.threads); os.environ['OPENBLAS_NUM_THREADS'] = str(args.threads)
import numpy as np
from scipy.optimize import least_squares
import multiprocessing as mp
from w2lib import *

d = args.d; sector = args.sector
pD = {'even': 0, 'odd': 1, 'mixed': None}[sector]
F, A, q = load()
G = scale(F, 9973)                      # integer form, same Jacobian ideal
Gi = [pderiv(G, i) for i in range(NV)]
H = [[pderiv(Gi[i], j) for j in range(NV)] for i in range(NV)]
tag = f'relax_d{d}_{sector}' + ('_cplx' if args.complex else '')
TENSOR = f'{LANE}/w2_{tag}_tensor.npz'


def build():
    t0 = time.time()
    S2 = MonSpace(2 * d, None if pD is None else 1)        # A * (w^T H w) lives here
    # J_{2d} generators (same parity) as rows
    gens = []
    for i in range(NV):
        par = None if pD is None else (1 - ZPAR[i]) % 2
        for m in monomials(2 * d - 5, par):
            gens.append(shift(Gi[i], m, 1))
    Jm = np.zeros((len(gens), S2.n))
    for r, g in enumerate(gens):
        ks, vs = poly_arrays(g); Jm[r, S2.index(ks)] = vs
    U, s, Vt = np.linalg.svd(Jm, full_matrices=True)
    rk = int((s > 1e-9 * s[0]).sum())
    C = Vt[rk:]                                            # (nS2 - rk) x nS2
    print(f'[relax] d={d} sector={sector}: dim S_2d={S2.n} dim J_2d={rk} quotient={C.shape[0]} ({time.time()-t0:.1f}s)', flush=True)
    # w basis (gauge removed later by projection)
    basis = []
    for i in range(NV):
        par = None if pD is None else (pD + ZPAR[i]) % 2
        for m in monomials(d - 5, par):
            basis.append((i, m))
    nw = len(basis)
    kA = key(next(iter(A)))
    T = np.zeros((C.shape[0], nw, nw))
    for k in range(nw):
        ik, mk = basis[k]
        for l in range(k, nw):
            il, ml = basis[l]
            P = H[ik][il]
            if not P:
                continue
            ks, vs = poly_arrays(P)
            ks = ks + key(mk) + key(ml) + kA
            v = np.zeros(S2.n); v[S2.index(ks)] = vs
            cv = C @ v
            T[:, k, l] += cv
            if l != k:
                T[:, l, k] += cv
    # gauge directions: w = ell z, ell in S_{d-6}^{pD}
    Gsp = MonSpace(d - 6, pD)
    bidx = {b: j for j, b in enumerate(basis)}
    Gm = np.zeros((nw, Gsp.n))
    for c, ell in enumerate(Gsp.mons):
        for i in range(NV):
            e = list(ell); e[i] += 1
            Gm[bidx[(i, tuple(e))], c] = 1.0
    if Gsp.n:
        Ug, sg, _ = np.linalg.svd(Gm, full_matrices=True)
        ng = int((sg > 1e-9 * sg[0]).sum()); Pp = Ug[:, ng:]
    else:
        ng = 0; Pp = np.eye(nw)
    # tensor in reduced coordinates
    Tr = np.einsum('ekl,ka,lb->eab', T, Pp, Pp)
    print(f'[relax] nw={nw} gauge={ng} reduced={Pp.shape[1]} equations={Tr.shape[0]} tensor built ({time.time()-t0:.1f}s)', flush=True)
    np.savez(TENSOR, T=Tr, Pp=Pp, ng=ng, nw=nw, neq=Tr.shape[0])
    return Tr, Pp


def worker_init():
    global TT, NEQ, NV_
    z = np.load(TENSOR)
    TT = z['T']; NEQ = TT.shape[0]; NV_ = TT.shape[1]


def run_seed(seed):
    rng = np.random.default_rng(seed)
    n = NV_
    t0 = time.time()
    if args.complex:
        L = rng.standard_normal(n) + 1j * rng.standard_normal(n)

        def resid(x):
            w = x[:n] + 1j * x[n:]
            r = np.einsum('ekl,k,l->e', TT, w, w)
            lin = L @ w - 1.0
            return np.concatenate([r.real, r.imag, [lin.real, lin.imag]])

        def jac(x):
            w = x[:n] + 1j * x[n:]
            Jc = 2 * np.einsum('ekl,l->ek', TT, w)
            top = np.hstack([Jc.real, -Jc.imag]); bot = np.hstack([Jc.imag, Jc.real])
            lt = np.hstack([L.real, -L.imag]); lb = np.hstack([L.imag, L.real])
            return np.vstack([top, bot, lt[None, :], lb[None, :]])
        w0 = rng.standard_normal(n) + 1j * rng.standard_normal(n); w0 /= (L @ w0)
        x0 = np.concatenate([w0.real, w0.imag])
        res = least_squares(resid, x0, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=args.maxnfev)
        w = res.x[:n] + 1j * res.x[n:]
        rr = np.einsum('ekl,k,l->e', TT, w, w)
        val = float(np.linalg.norm(rr) / max(np.linalg.norm(w) ** 2, 1e-300))
    else:
        def resid(w):
            r = np.einsum('ekl,k,l->e', TT, w, w)
            return np.concatenate([r, [w @ w - 1.0]])

        def jac(w):
            return np.vstack([2 * np.einsum('ekl,l->ek', TT, w), 2 * w[None, :]])
        w0 = rng.standard_normal(n); w0 /= np.linalg.norm(w0)
        res = least_squares(resid, w0, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=args.maxnfev)
        w = res.x
        rr = np.einsum('ekl,k,l->e', TT, w, w)
        val = float(np.linalg.norm(rr) / max(w @ w, 1e-300))
    out = {'seed': int(seed), 'res': val, 'nfev': int(res.nfev), 'time': round(time.time() - t0, 2)}
    if val < 1e-8:
        # local dimension of the solution set: rank of the Jacobian of the quadratic map at w
        Jc = 2 * np.einsum('ekl,l->ek', TT, w)
        sv = np.linalg.svd(Jc, compute_uv=False)
        rk = int((sv > 1e-8 * sv[0]).sum())
        out['jac_rank'] = rk; out['local_dim'] = n - rk - 1
        np.save(f'{LANE}/w2_{tag}_hit_{seed}.npy', w)
    return out


if __name__ == '__main__':
    if not os.path.exists(TENSOR):
        build()
    z = np.load(TENSOR); neq = int(z['neq']); nred = int(z['T'].shape[1]); ng = int(z['ng'])
    # random-point residual scale for the relative report
    worker_init()
    rng = np.random.default_rng(12345)
    scale0 = float(np.median([np.linalg.norm(np.einsum('ekl,k,l->e', TT, w, w)) / (w @ w)
                              for w in [rng.standard_normal(nred) for _ in range(20)]]))
    seeds = list(range(args.seed, args.seed + args.restarts))
    t0 = time.time()
    allres = []
    with mp.Pool(args.workers, initializer=worker_init) as pool:
        for b0 in range(0, len(seeds), args.batch):
            chunk = seeds[b0:b0 + args.batch]
            outs = pool.map(run_seed, chunk)
            for o in outs:
                append_json(f'{LANE}/w2_{tag}.jsonl', o)
            vals = np.array([o['res'] for o in outs]); allres += list(vals)
            hits = [o for o in outs if o['res'] < 1e-8]
            rec = {'batch': 'hessian-relaxation', 'family': 'jacobian-type w-space (E2 mod J)', 'd': d, 'sector': sector,
                   'complex': bool(args.complex), 'unknowns': nred, 'gauge_removed': ng, 'equations': neq,
                   'restarts': len(chunk), 'seeds': [chunk[0], chunk[-1]], 'best_res': float(vals.min()),
                   'median_res': float(np.median(vals)), 'q10_res': float(np.quantile(vals, 0.1)),
                   'random_point_res': scale0, 'best_rel_to_random': float(vals.min() / scale0),
                   'hits_below_1e-8': len(hits), 'hit_local_dims': sorted({o.get('local_dim') for o in hits}) if hits else [],
                   'cum_restarts': len(allres), 'cum_best': float(min(allres)), 'seconds': round(time.time() - t0, 1)}
            print(json.dumps(rec), flush=True)
            log_batch(rec)
