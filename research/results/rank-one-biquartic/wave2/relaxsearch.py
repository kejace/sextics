"""Relaxation-first numerical search: minimize |Phi(a (x) a)|^2 where Phi(a (x) a) = projection of B(D(a), N0(a)) onto
(S/J)_{2d} (float, orthonormal complement of J_{2d}); ~100x cheaper per evaluation than the full normalized E2.
Zeros of the relaxation (rel < 1e-8) are refined on the full E2 (L-BFGS from that point) and logged.

  python3 wave2/relaxsearch.py --d 8 --pD 0 --restarts 2000 --seed 0 --workers 8 --tag r8e_a [--family all|kernel]
"""
import sys, os, time, json, argparse
if '--workers' in sys.argv and int(sys.argv[sys.argv.index('--workers') + 1]) > 1:
    for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[_v] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
import numpy as np
from scipy.optimize import minimize
import multiprocessing as mp

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--pD', type=int, default=0)
ap.add_argument('--restarts', type=int, default=200)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--workers', type=int, default=1)
ap.add_argument('--maxiter', type=int, default=500)
ap.add_argument('--tag', required=True)
ap.add_argument('--family', default='all')
ap.add_argument('--rcond', type=float, default=1e-9)
ap.add_argument('--refine-rel', type=float, default=1e-7, help='relaxation rel below which the full E2 refinement runs')
ap.add_argument('--refine-maxiter', type=int, default=300)
args = ap.parse_args()
if os.path.exists(f'{OUT}/STOP_{args.tag}'):
    print('[relaxsearch] STOP file present; skipping', file=sys.stderr); sys.exit(0)

F, A, q = load()
d, pD = args.d, args.pD
t0 = time.time()
prob = E2Float(F, A, d, pD)
fam = E1Family(prob, rcond=args.rcond)
Dbasis, Nbasis = fam.Dbasis, fam.Nbasis
if args.family == 'kernel':
    gens = []
    for i in range(NV):
        par = (pD - par_idx(i)) % 2
        Sm = Space(d - 5, par)
        if Sm.n:
            gens.append(mult_op(prob.Fi[i], Sm, prob.SD).toarray())
    G = np.hstack(gens)
    QJ, sJ, _ = sla.svd(G, full_matrices=False)
    rJ = int((sJ > 1e-9 * sJ[0]).sum()); QJ = QJ[:, :rJ]
    Pperp = fam.Dbasis.T - QJ @ (QJ.T @ fam.Dbasis.T)
    U, s, Vt = sla.svd(Pperp, full_matrices=False)
    rk = int((s > 1e-7 * s[0]).sum())
    Dbasis = U[:, :rk].T
    coords = fam.Dbasis @ Dbasis.T
    Nbasis = (fam.Nbasis.T @ coords).T
K = Dbasis.shape[0]
# orthonormal complement of J_{2d} inside Srow (float)
S2 = prob.Srow
Jm = sp.hstack([mult_op(prob.Fi[i], Space(2 * d - 5, (prob.pA - par_idx(i)) % 2), S2) for i in range(NV)]).toarray()
Uj, sj, _ = sla.svd(Jm, full_matrices=True, check_finite=False, lapack_driver='gesdd')
rJ2 = int((sj > 1e-9 * sj[0]).sum())
Qperp = np.ascontiguousarray(Uj[:, rJ2:])            # Srow.n x ring_dim
ring_dim = Qperp.shape[1]
del Jm, Uj
# spurious complement chart (relaxation is invariant under spurious directions -> optimize on the complement only)
SG = Space(d - 6, (pD - core.parity_of(F)) % 2) if d >= 6 else None
if SG is not None and SG.n:
    FG = mult_op(F, SG, prob.SD).toarray()
    coords = Dbasis @ FG
    Sq, _ = np.linalg.qr(coords)
    Pfull, _ = np.linalg.qr(np.hstack([Sq, np.eye(K)]))
    Pc = Pfull[:, Sq.shape[1]:K]
    n_sp = Sq.shape[1]
else:
    Pc = np.eye(K); n_sp = 0
n_u = Pc.shape[1]
print(f'[relaxsearch] d={d} pD={pD} family={args.family} K={K} spurious={n_sp} n_u={n_u} ring_dim(S/J)_{2*d}={ring_dim} build {time.time()-t0:.1f}s', file=sys.stderr, flush=True)

def phi_and_grad(u):
    nu = np.linalg.norm(u); uh = u / nu
    a = Pc @ uh
    Dv = Dbasis.T @ a; Nv = Nbasis.T @ a
    B = prob.Bvec(Dv, Nv)
    y = Qperp.T @ B                     # relaxation residual (ring coordinates)
    f = float(y @ y)
    w = Qperp @ y
    gD, gN = prob.adjoint_B(w, Dv, Nv)
    ga = 2.0 * (Dbasis @ gD + Nbasis @ gN)
    gu_dir = Pc.T @ ga
    gu = (gu_dir - (gu_dir @ uh) * uh) / nu
    return f, gu

def relax_rel(a):
    Dv = Dbasis.T @ a; Nv = Nbasis.T @ a
    B = prob.Bvec(Dv, Nv)
    y = Qperp.T @ B
    return float(np.linalg.norm(y) / max(np.linalg.norm(B), 1e-300)), float(np.linalg.norm(y)), float(np.linalg.norm(B))

def full_e2(a, maxiter=0):
    """full normalized-E2 residual at a (optionally refined by L-BFGS over the complement chart + spurious s)."""
    Sq_local = Sq if n_sp else None
    def val_grad(z):
        u = z[:n_u]; s = z[n_u:]
        nu = np.linalg.norm(u); uh = u / nu
        aa = Pc @ uh + (Sq_local @ s if n_sp else 0.0)
        Dv = Dbasis.T @ aa; Nv = Nbasis.T @ aa
        B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
        Br = prob.reduce(B); Mr = prob.reduce(M)
        c, *_ = sla.lstsq(Mr, -Br, cond=1e-10, lapack_driver='gelsd', check_finite=False)
        r = Br + Mr @ c
        w = prob.reduce_T(r)
        Ntot = Nv + prob.Kmap @ c
        gD, gN = prob.adjoint_B(w, Dv, Ntot)
        ga = 2.0 * (Dbasis @ gD + Nbasis @ gN)
        gu_dir = Pc.T @ ga
        gu = (gu_dir - (gu_dir @ uh) * uh) / nu
        gs = Sq_local.T @ ga if n_sp else np.zeros(0)
        return float(r @ r), np.concatenate([gu, gs])
    u0 = Pc.T @ a
    z = np.concatenate([u0, np.zeros(n_sp)])
    if maxiter > 0:
        res = minimize(val_grad, z, jac=True, method='L-BFGS-B', options={'maxiter': maxiter, 'gtol': 1e-12, 'ftol': 1e-18})
        z = res.x
    u = z[:n_u]; s = z[n_u:]
    aa = Pc @ (u / np.linalg.norm(u)) + (Sq_local @ s if n_sp else 0.0)
    Dv = Dbasis.T @ aa; Nv = Nbasis.T @ aa
    B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
    Br = prob.reduce(B); Mr = prob.reduce(M)
    c, *_ = sla.lstsq(Mr, -Br, cond=1e-10, lapack_driver='gelsd', check_finite=False)
    r = Br + Mr @ c
    return float(np.linalg.norm(r) / max(prob.Bpieces_norm(Dv, Nv), 1e-300)), aa, Dv

def run(worker_id, seeds):
    recs = []
    path = f'{OUT}/{args.tag}.restarts.jsonl'
    for sd in seeds:
        if os.path.exists(f'{OUT}/STOP_{args.tag}'):
            break
        rng = np.random.default_rng(sd)
        u0 = rng.standard_normal(n_u); u0 /= np.linalg.norm(u0)
        t1 = time.time()
        res = minimize(phi_and_grad, u0, jac=True, method='L-BFGS-B', options={'maxiter': args.maxiter, 'gtol': 1e-13, 'ftol': 1e-20})
        a = Pc @ (res.x / np.linalg.norm(res.x))
        rrel, yn, bn = relax_rel(a)
        rec = {'tag': args.tag, 'd': d, 'pD': pD, 'seed': int(sd), 'relax_f': float(res.fun), 'relax_rel': rrel, 'relax_norm': yn, 'B_norm': bn,
               'iters': int(res.nit), 'time': round(time.time() - t1, 1)}
        if rrel < args.refine_rel:
            t2 = time.time()
            e2_0, _, _ = full_e2(a, 0)
            e2_1, aa, Dv = full_e2(a, args.refine_maxiter)
            rec.update({'e2_rel_at_point': e2_0, 'e2_rel_refined': e2_1, 'refine_time': round(time.time() - t2, 1)})
            if e2_1 < 1e-6:
                rec['Dcoef'] = [float(x) for x in Dv]
                log_line(f'{OUT}/{args.tag}.cands.jsonl', rec)
        log_line(path, rec)
        recs.append(rec)
    return recs

if __name__ == '__main__':
    seeds = list(range(args.seed, args.seed + args.restarts))
    ts = time.time()
    if args.workers > 1:
        chunks = [seeds[i::args.workers] for i in range(args.workers)]
        with mp.get_context('fork').Pool(args.workers) as pool:
            results = pool.starmap(run, [(i, ch) for i, ch in enumerate(chunks) if ch])
        recs = [r for rs in results for r in rs]
    else:
        recs = run(0, seeds)
    rr = sorted(r['relax_rel'] for r in recs)
    e2s = sorted(r['e2_rel_refined'] for r in recs if 'e2_rel_refined' in r)
    summary = {'batch': 'relaxsearch', 'tag': args.tag, 'family': args.family, 'd': d, 'pD': pD, 'K': K, 'n_spurious': n_sp, 'ring_dim_2d': ring_dim,
               'restarts': len(recs), 'seed0': args.seed, 'best_relax_rel': rr[0] if rr else None, 'median_relax_rel': rr[len(rr)//2] if rr else None,
               'n_relax_lt_1e-6': sum(1 for x in rr if x < 1e-6), 'n_refined': len(e2s), 'best_e2_rel_refined': e2s[0] if e2s else None,
               'n_e2_lt_1e-6': sum(1 for x in e2s if x < 1e-6), 'time_s': round(time.time() - ts, 1)}
    log_line(f'{OUT}/{args.tag}.batch.jsonl', summary)
    print(json.dumps(summary), flush=True)
