"""Variable-projection multistart search over an E1-admissible family (normalized E2, q = 1).

  python3 wave2/e2search.py --d 7 --pD 0 --restarts 200 --seed 0 --workers 8 --tag d7a [--maxiter 400]
         [--family J|all|file.json] [--rcond 1e-9] [--cond 1e-10] [--batch-note "..."]

Chart: a = Pc u/|u| + S s with S = coordinates of F*S_{d-6} (spurious F-multiples), u free direction,
s free (unbounded) -> lossless parametrization of denominators nonzero modulo F.
Per restart one JSON line in out/<tag>.restarts.jsonl; candidates (rel < 1e-6) in out/<tag>.cands.jsonl;
rel < 1e-8 triggers rational reconstruction + exact confirmation (fixed_d).  Batch summary in out/<tag>.batch.jsonl.
"""
import sys, os, time, json, argparse, math
# single-threaded BLAS per worker when running several workers (set before numpy is imported)
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
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--workers', type=int, default=1)
ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--tag', required=True)
ap.add_argument('--family', default='all', help="all = full E1-admissible space; J = Jacobian part only; kernel = complement of J_d inside the admissible space (pure E1-kernel directions); path.json = list of D expressions")
ap.add_argument('--rcond', type=float, default=1e-9)
ap.add_argument('--cond', type=float, default=1e-10)
ap.add_argument('--cand-rel', type=float, default=1e-6)
ap.add_argument('--confirm-rel', type=float, default=1e-8)
ap.add_argument('--batch-note', default='')
ap.add_argument('--report-every', type=int, default=25)
ap.add_argument('--instance', default=None, help='override instance id (validation runs)')
ap.add_argument('--method', default='lm', choices=['lm', 'lbfgs'], help='local optimizer: lm = Levenberg-Marquardt with the variable-projection (Kaufman) Jacobian; lbfgs = L-BFGS-B on f')
ap.add_argument('--lm-iters', type=int, default=60)
ap.add_argument('--smax', type=float, default=5.0, help='bound on the spurious (F-multiple) component |s| in the chart |u|=1')
args = ap.parse_args()
# cooperative skip: a file out/STOP_<tag> makes this invocation exit immediately (used to drop steps from running chains)
if os.path.exists(f'{OUT}/STOP_{args.tag}'):
    print(f'[search] STOP file present for tag {args.tag}; skipping this step', file=sys.stderr)
    log_line(f'{OUT}/{args.tag}.batch.jsonl', {'batch': 'e2search', 'tag': args.tag, 'd': args.d, 'pD': args.pD, 'skipped': True,
             'note': 'skipped via STOP file (degree excluded exactly / reprioritized)', 'K': None, 'n_spurious': None, 'n_c': None, 'n_C': None,
             'rows': None, 'rows_red': None, 'nN': None, 'restarts': 0, 'best_rel': None, 'median_rel': None, 'worst_rel': None,
             'n_rel_lt_1e-3': 0, 'n_rel_lt_1e-6': 0, 'n_rel_lt_1e-8': 0, 'witness': False, 'time_s': 0})
    sys.exit(0)

F, A, q = load(args.instance)
IID_RUN = args.instance or IID
d, pD = args.d, args.pD
prob = E2Float(F, A, d, pD)
t0 = time.time()
if args.family in ('all', 'J'):
    fam = E1Family(prob, rcond=args.rcond)
    Dbasis, Nbasis = fam.Dbasis, fam.Nbasis
    if args.family == 'J':
        # restrict to J_d: span of F_i m  with exact lifts A m e_i
        cols, lifts = [], []
        for i in range(NV):
            par = (pD - par_idx(i)) % 2
            Sm = Space(d - 5, par)
            for m in Sm.mons:
                Dm = shift(prob.Fi[i], m, 1)
                Nm = [{} for _ in range(NV)]; Nm[i] = shift(A, m, 1)
                cols.append(prob.SD.vec(Dm)); lifts.append(prob.Nvec(Nm))
        G = np.array(cols); Ls = np.array(lifts)
        # orthonormalize D-part, carry lifts along
        U, s, Vt = sla.svd(G.T, full_matrices=False)
        r = int((s > 1e-9 * s[0]).sum())
        Dbasis = U[:, :r].T
        # lifts for the orthonormal basis: Dbasis = (Vt[:r].T / s) ... solve G^T y = Dbasis^T
        Y = np.linalg.lstsq(G.T, Dbasis.T, rcond=None)[0]   # (ngen x r)
        Nbasis = (Ls.T @ Y).T
    fam_info = {'K': fam.K, 'J_rank': fam.J_rank, 'ring_dim': fam.ring_dim, 'J_gap': fam.J_gap, 'C_gap': fam.C_gap, 'e1_rel': fam.e1_rel}
elif args.family == 'kernel':
    fam = E1Family(prob, rcond=args.rcond)
    gens = []
    for i in range(NV):
        par = (pD - par_idx(i)) % 2
        Sm = Space(d - 5, par)
        if Sm.n:
            gens.append(mult_op(prob.Fi[i], Sm, prob.SD).toarray())
    if gens:
        G = np.hstack(gens)
        QJ, sJ, _ = sla.svd(G, full_matrices=False)
        rJ = int((sJ > 1e-9 * sJ[0]).sum()); QJ = QJ[:, :rJ]
        Pperp = fam.Dbasis.T - QJ @ (QJ.T @ fam.Dbasis.T)      # SD.n x K, columns = family basis minus J-part
        U, s, Vt = sla.svd(Pperp, full_matrices=False)
        rk = int((s > 1e-7 * s[0]).sum()) if s.size else 0
        Dbasis = U[:, :rk].T                                   # rk x SD.n orthonormal, orthogonal to J_d
        coords = fam.Dbasis @ Dbasis.T                         # K x rk
        Nbasis = (fam.Nbasis.T @ coords).T
        fam_info = {'K': int(rk), 'K_full': fam.K, 'dim_J_d': rJ, 'sv_gap': [float(s[rk-1]) if rk else None, float(s[rk]) if rk < len(s) else None]}
    else:
        Dbasis, Nbasis = fam.Dbasis, fam.Nbasis
        fam_info = {'K': fam.K, 'K_full': fam.K, 'dim_J_d': 0}
else:
    exprs = json.load(open(args.family))
    cols, lifts = [], []
    famE = E1Family(prob, rcond=args.rcond)
    for ex in exprs:
        Dm = core.parse(ex)
        v = prob.SD.vec(Dm)
        co = famE.Dbasis @ v
        fit = np.linalg.norm(famE.Dbasis.T @ co - v) / np.linalg.norm(v)
        if fit > 1e-8:
            print(f'[family] member {ex[:60]} not E1-admissible (fit {fit:.2e}); skipped', file=sys.stderr)
            continue
        cols.append(v); lifts.append(famE.Nbasis.T @ co)
    G = np.array(cols); Ls = np.array(lifts)
    U, s, Vt = sla.svd(G.T, full_matrices=False)
    r = int((s > 1e-9 * s[0]).sum())
    Dbasis = U[:, :r].T
    Y = np.linalg.lstsq(G.T, Dbasis.T, rcond=None)[0]
    Nbasis = (Ls.T @ Y).T
    fam = famE
    fam_info = {'K': int(r), 'K_full': famE.K}
K = Dbasis.shape[0]
# spurious subspace coordinates (F multiples)
SG = Space(d - 6, (pD - core.parity_of(F)) % 2) if d >= 6 else None
if SG is not None and SG.n:
    FG = mult_op(F, SG, prob.SD).toarray()
    coords = Dbasis @ FG
    fit_sp = np.linalg.norm(Dbasis.T @ coords - FG) / np.linalg.norm(FG)
    Sq, _ = np.linalg.qr(coords)
    if fit_sp > 1e-6:
        print(f'[family] F-multiples not inside family (fit {fit_sp:.2e}) -> no spurious directions removed', file=sys.stderr)
        Sq = np.zeros((K, 0))
else:
    Sq = np.zeros((K, 0)); fit_sp = 0.0
n_sp = Sq.shape[1]
# complement
if n_sp:
    Pfull, _ = np.linalg.qr(np.hstack([Sq, np.eye(K)]))
    Pc = Pfull[:, n_sp:K]
else:
    Pc = np.eye(K)
n_u = Pc.shape[1]
print(f'[search] d={d} pD={pD} K={K} spurious={n_sp} n_u={n_u} n_c={prob.n_c} rows_red={prob.n_red} rows={prob.Srow.n} build {time.time()-t0:.1f}s', file=sys.stderr, flush=True)

# generic rank statistics of the reduced Koszul matrix
rng0 = np.random.default_rng(12345)
a0 = Pc @ rng0.standard_normal(n_u); a0 /= np.linalg.norm(a0)
Mr0 = prob.reduce(prob.koszul_matrix(Dbasis.T @ a0))
sv = sla.svdvals(Mr0)
rk0 = int((sv > 1e-9 * sv[0]).sum())
rank_stats = {'Mred_shape': list(Mr0.shape), 'rank': rk0, 'nullity': int(Mr0.shape[1] - rk0),
              'sv_min_kept': float(sv[rk0 - 1]), 'sv_next': float(sv[rk0]) if rk0 < len(sv) else 0.0, 'sv_max': float(sv[0])}
print('[search] generic reduced Koszul rank', rank_stats, file=sys.stderr, flush=True)
del Mr0

def unpack(z):
    u = z[:n_u]; s = z[n_u:]
    nu = np.linalg.norm(u)
    if n_sp:
        ns = np.linalg.norm(s)
        if ns > args.smax:
            s = s * (args.smax / ns); z[n_u:] = s
    a = Pc @ (u / nu) + (Sq @ s if n_sp else 0.0)
    return a, u, nu

def value_and_grad(z):
    a, u, nu = unpack(z)
    Dv = Dbasis.T @ a; Nv = Nbasis.T @ a
    B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
    Br = prob.reduce(B); Mr = prob.reduce(M)
    c, *_ = sla.lstsq(Mr, -Br, cond=args.cond, lapack_driver='gelsd', check_finite=False)
    r = Br + Mr @ c
    f = float(r @ r)
    w = prob.reduce_T(r)
    Ntot = Nv + prob.Kmap @ c
    gD, gN = prob.adjoint_B(w, Dv, Ntot)
    ga = 2.0 * (Dbasis @ gD + Nbasis @ gN)
    # chain rule
    gu_dir = Pc.T @ ga                        # d/d(u/|u|)
    uh = u / nu
    gu = (gu_dir - (gu_dir @ uh) * uh) / nu
    gs = Sq.T @ ga if n_sp else np.zeros(0)
    return f, np.concatenate([gu, gs])

def evaluate(a):
    Dv = Dbasis.T @ a; Nv = Nbasis.T @ a
    B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
    Br = prob.reduce(B); Mr = prob.reduce(M)
    c, *_ = sla.lstsq(Mr, -Br, cond=args.cond, lapack_driver='gelsd', check_finite=False)
    r = Br + Mr @ c
    scale = prob.Bpieces_norm(Dv, Nv)
    return float(np.linalg.norm(r)), float(np.linalg.norm(Br)), float(np.linalg.norm(B)), c, Dv, scale

def scale0_of(u):
    """cancellation-free scale of the u-part alone (s = 0): immune to inflation by spurious components."""
    a0 = Pc @ (u / np.linalg.norm(u))
    return prob.Bpieces_norm(Dbasis.T @ a0, Nbasis.T @ a0)

def spurious_rel(Dv):
    """relative residual of least-squares division of D by F."""
    if d < 6 or SG is None or SG.n == 0:
        return 1.0
    x, *_ = sla.lstsq(FG, Dv, check_finite=False)
    return float(np.linalg.norm(FG @ x - Dv) / max(np.linalg.norm(Dv), 1e-300))

def reconstruct_and_confirm(Dv, tag):
    """Try rational reconstruction of D (several scalings / denominators) and exact confirmation."""
    out = []
    order = np.argsort(-np.abs(Dv))
    tried = set()
    for pivot in order[:3]:
        v = Dv / Dv[pivot]
        for max_den in (1, 2, 3, 4, 6, 7, 12, 14, 21, 42, 84, 168, 1000, 10**4, 10**5, 10**6):
            fr = [Fraction(float(x)).limit_denominator(max_den) for x in v]
            err = max(abs(float(f_) - x) for f_, x in zip(fr, v))
            if err > 1e-7:
                continue
            Dq = {m: f_ for m, f_ in zip(prob.SD.mons, fr) if f_ != 0}
            key = core.to_str(Dq)
            if key in tried:
                continue
            tried.add(key)
            rec = {'tag': tag, 'max_den': max_den, 'err': err, 'D': key[:4000], 'nterms': len(Dq)}
            # exact tests
            divF, _ = jacobian.divisible_by(F, Dq)
            rec['F_divides_D'] = bool(divF)
            if divF:
                rec['status'] = 'spurious (F | D)'; out.append(rec); continue
            ok, info = fixed_d.e1_solvable(F, A, q, Dq)
            rec['e1'] = bool(ok); rec['e1_info'] = info
            if not ok:
                rec['status'] = 'e1 fails exactly'; out.append(rec); continue
            try:
                res = fixed_d.solve_fixed_D_fast(F, A, q, Dq, verbose=False)
                rec['consistent'] = bool(res['consistent']); rec['rank'] = res['rank']; rec['nullity'] = res['nullity']
                if res['consistent'] and res.get('exact'):
                    okw, winfo = core.verify_witness(F, A, q, Dq, res['N'], res['Q'])
                    rec['verified'] = bool(okw)
                    if okw:
                        rec['status'] = 'WITNESS'
                        w = {'instance': IID_RUN, 'd': d, 'D': core.to_str(Dq), 'N': [core.to_str(n) for n in res['N']],
                             'Q': core.to_str(res['Q']), 'family': tag}
                        wpath = f'{LANE}/witness-candidate.json' if args.instance is None else f'{OUT}/witness-{IID_RUN}-{tag}.json'
                        with open(wpath, 'w') as fh:
                            json.dump(w, fh, indent=1)
                        print('!!!!!!!!!! WITNESS FOUND !!!!!!!!!!', key[:200], file=sys.stderr, flush=True)
                        out.append(rec)
                        return out, True
                elif res['consistent']:
                    rec['status'] = 'consistent mod p but exact lift failed'
                else:
                    rec['status'] = 'inconsistent exactly (mod p)'
            except Exception as ex:
                rec['status'] = f'solver error {ex!r}'
            out.append(rec)
    return out, False

def residual_and_jacobian(z):
    """Reduced residual r(z) and its Kaufman Jacobian wrt z = (u, s) (holding the optimal c fixed)."""
    a, u, nu = unpack(z)
    uh = u / nu
    Dv = Dbasis.T @ a; Nv = Nbasis.T @ a
    B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
    Br = prob.reduce(B); Mr = prob.reduce(M)
    c, *_ = sla.lstsq(Mr, -Br, cond=args.cond, lapack_driver='gelsd', check_finite=False)
    r = Br + Mr @ c
    Ntot = Nv + prob.Kmap @ c
    Ja = np.empty((prob.n_red, K))
    for m in range(K):
        col = prob.Bvec(Dbasis[m], Ntot) + prob.Bvec(Dv, Nbasis[m])
        Ja[:, m] = prob.reduce(col)
    # chain rule: da/du = Pc (I - uh uh^T)/nu ; da/ds = Sq
    Ju = (Ja @ Pc - np.outer(Ja @ (Pc @ uh), uh)) / nu
    Jz = np.hstack([Ju, Ja @ Sq]) if n_sp else Ju
    return r, Jz, c

def lm_minimize(z0, iters, verbose=False):
    z = z0.copy()
    r, J, c = residual_and_jacobian(z)
    f = float(r @ r)
    lam = 1e-3
    nfev = 1
    for it in range(iters):
        JtJ = J.T @ J; g = J.T @ r
        if np.linalg.norm(g) < 1e-14:
            break
        dscale = np.sqrt(np.maximum(np.diag(JtJ), 1e-12))
        improved = False
        for _ in range(12):
            try:
                step = -np.linalg.solve(JtJ + lam * np.diag(dscale ** 2), g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            z_new = z + step
            a_new, u_new, nu_new = unpack(z_new)
            Dv = Dbasis.T @ a_new; Nv = Nbasis.T @ a_new
            B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
            Br = prob.reduce(B); Mr = prob.reduce(M)
            c_new, *_ = sla.lstsq(Mr, -Br, cond=args.cond, lapack_driver='gelsd', check_finite=False)
            r_new = Br + Mr @ c_new
            f_new = float(r_new @ r_new)
            nfev += 1
            if f_new < f:
                # renormalize the u-part to keep the chart well conditioned
                u_new = z_new[:n_u]; z_new[:n_u] = u_new / np.linalg.norm(u_new)
                z, f = z_new, f_new
                lam = max(lam / 3.0, 1e-12)
                improved = True
                break
            lam *= 4.0
        if not improved:
            break
        r, J, c = residual_and_jacobian(z)
        nfev += 1
        if f < 1e-28:
            break
    return z, f, it + 1, nfev

def run_restarts(worker_id, seeds):
    rng = np.random.default_rng(seeds[0] * 7919 + worker_id)
    recs = []
    path = f'{OUT}/{args.tag}.restarts.jsonl'
    cpath = f'{OUT}/{args.tag}.cands.jsonl'
    for sd in seeds:
        if os.path.exists(f'{OUT}/STOP_{args.tag}'):
            print(f'[search] STOP file for {args.tag}: ending batch early after {len(recs)} restarts in this worker', file=sys.stderr, flush=True)
            break
        rng = np.random.default_rng(sd)
        u0 = rng.standard_normal(n_u); u0 /= np.linalg.norm(u0)
        z0 = np.concatenate([u0, np.zeros(n_sp)])
        t1 = time.time()
        if args.method == 'lm':
            zf, ff, nit, nfev = lm_minimize(z0, args.lm_iters)
            class _R: pass
            res = _R(); res.x = zf; res.fun = ff; res.nit = nit; res.nfev = nfev; res.status = 0
        else:
            res = minimize(value_and_grad, z0, jac=True, method='L-BFGS-B',
                           options={'maxiter': args.maxiter, 'gtol': 1e-12, 'ftol': 1e-18, 'maxcor': 20})
        a, u, nu = unpack(res.x)
        rn, brn, bn, c, Dv, scale = evaluate(a)
        scale0 = scale0_of(u)
        rel = rn / max(scale0, 1e-300)
        sp_rel = spurious_rel(Dv)
        rec = {'tag': args.tag, 'd': d, 'pD': pD, 'seed': int(sd), 'f': float(res.fun), 'rel': rel, 'rel_scaled': rn / max(scale, 1e-300), 'rel_B': rn / max(brn, 1e-300), 'scale': scale, 'scale0': scale0,
               'resnorm': rn, 'Bred_norm': brn, 'iters': int(res.nit), 'nfev': int(res.nfev), 'status': int(res.status),
               's_norm': float(np.linalg.norm(res.x[n_u:])) if n_sp else 0.0, 'spurious_rel': sp_rel, 'method': args.method,
               'time': round(time.time() - t1, 1)}
        log_line(path, rec)
        if rel < args.cand_rel:
            crec = dict(rec); crec['a'] = [float(x) for x in a]; crec['Dcoef'] = [float(x) for x in Dv]
            log_line(cpath, crec)
            if rel < args.confirm_rel and sp_rel > 1e-6:
                conf, found = reconstruct_and_confirm(Dv, args.tag)
                for cr in conf:
                    log_line(f'{OUT}/{args.tag}.confirm.jsonl', cr)
                if found:
                    rec['WITNESS'] = True
        recs.append(rec)
    return recs

if __name__ == '__main__':
    seeds = list(range(args.seed, args.seed + args.restarts))
    t_start = time.time()
    if args.workers > 1:
        chunks = [seeds[i::args.workers] for i in range(args.workers)]
        ctx = mp.get_context('fork')
        with ctx.Pool(args.workers) as pool:
            results = pool.starmap(run_restarts, [(i, ch) for i, ch in enumerate(chunks) if ch])
        recs = [r for rs in results for r in rs]
    else:
        recs = run_restarts(0, seeds)
    rels = sorted(r['rel'] for r in recs)
    summary = {'batch': 'e2search', 'tag': args.tag, 'family': args.family, 'd': d, 'pD': pD, 'K': K, 'n_spurious': n_sp,
               'n_c': prob.n_c, 'n_C': prob.SC.n, 'rows': prob.Srow.n, 'rows_red': prob.n_red, 'nN': prob.nN,
               'restarts': len(recs), 'seed0': args.seed, 'maxiter': args.maxiter, 'method': args.method,
               'best_rel': rels[0] if rels else None, 'median_rel': rels[len(rels)//2] if rels else None,
               'worst_rel': rels[-1] if rels else None,
               'n_rel_lt_1e-3': sum(1 for r in rels if r < 1e-3), 'n_rel_lt_1e-6': sum(1 for r in rels if r < 1e-6),
               'n_rel_lt_1e-8': sum(1 for r in rels if r < 1e-8),
               'witness': any(r.get('WITNESS') for r in recs),
               'fam_info': fam_info, 'rank_stats': rank_stats, 'spurious_fit': float(fit_sp),
               'time_s': round(time.time() - t_start, 1), 'note': args.batch_note}
    log_line(f'{OUT}/{args.tag}.batch.jsonl', summary)
    print(json.dumps(summary), flush=True)
