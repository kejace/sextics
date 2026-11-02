"""Step 2/3 driver: multi-restart reduced variable-projection search in one (d, chi, kind) family.

usage: python3 -u w2_varpro.py --d 8 --chi 0,0 --kind full|jac|ker --hess 0|1 --restarts 1000
                               --seed 100 --workers 8 --maxiter 300 --name <jobname> [--chunk 10]
Per-restart records -> wave2-<name>.jsonl ; batch summary -> wave2-search.log.jsonl ; report rewritten.
Candidates with rel < 1e-5 are polished (Levenberg-Marquardt on the VP residual); rel < 1e-10 ->
rational reconstruction + exact sector-restricted re-solve; a solved candidate is written to
witness-candidate.json and announced with WITNESS in the log.
"""
import sys, os, time, json, argparse, math
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import numpy as np
from scipy.optimize import least_squares
import w2lib as W, w2fast as WF
from hodge import core
import multiprocessing as mp

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--chi', type=str, required=True)       # "a,b"
ap.add_argument('--kind', type=str, default='full')
ap.add_argument('--hess', type=int, default=0)
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--workers', type=int, default=4)
ap.add_argument('--maxiter', type=int, default=300)
ap.add_argument('--chunk', type=int, default=10)
ap.add_argument('--name', type=str, required=True)
ap.add_argument('--polish_tol', type=float, default=1e-5)
ap.add_argument('--cand_tol', type=float, default=1e-7)
ap.add_argument('--qproj', type=int, default=1)
args = ap.parse_args()
chi = tuple(int(x) for x in args.chi.split(','))

F, A, q = W.load()
t0 = time.time()
fam, S = WF.build_family(F, A, q, args.d, chi, args.kind, hessian=bool(args.hess), verbose=True, qproj=bool(args.qproj))
if fam is None:
    W.log_batch({'step': 'varpro', 'name': args.name, 'd': args.d, 'chi': list(chi), 'kind': args.kind,
                 'hess': args.hess, 'restarts': 0, 'note': 'empty family', 'spaces': S.summary()})
    print('empty family'); sys.exit(0)
info = fam.info(); info['spaces'] = S.summary()
# rank / nullity statistics of the linear part at a random point
rng0 = np.random.default_rng(12345)
a0 = rng0.standard_normal(fam.K); a0 /= np.linalg.norm(a0)
M0, b0 = fam.assemble(a0)
sv = np.linalg.svd(M0, compute_uv=False)
rank0 = int(np.sum(sv > 1e-9 * sv[0])) if sv.size else 0
info['rank_M_random'] = rank0; info['nullity_M_random'] = int(M0.shape[1] - rank0)
info['independent_equations_random'] = rank0
info['cokernel_random'] = int(M0.shape[0] - rank0)
print('[info]', json.dumps(info), flush=True)


def residual_vec(a):
    return fam.residual_vec(a)


def polish(a, max_nfev=None):
    """LM/TRF on the VP residual r(a) with finite-difference Jacobian; returns (a, rel)."""
    K = fam.K
    res = least_squares(residual_vec, np.asarray(a, float), method='trf', jac='2-point',
                        max_nfev=max_nfev or 40 * K, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    a1 = res.x / np.linalg.norm(res.x)
    f, x, r = fam.value(a1)
    return a1, fam.rel(f)


def work(task):
    seed, n = task
    t = time.time()
    out = fam.search(n_restarts=n, seed=seed, maxiter=args.maxiter)
    recs = []
    for c in out:
        rec = {'seed': seed, 'restart': c['restart'], 'rel': float(c['rel']), 'value': float(c['value']),
               'iters': int(c['iters']), 'spurious': bool(c['spurious'])}
        if c['rel'] < args.polish_tol:
            a1, rel1 = polish(c['a'])
            rec['polished_rel'] = float(rel1)
            rec['a'] = [float(v) for v in a1]
            rec['rel_before_polish'] = float(c['rel'])
            rec['rel'] = float(min(rel1, c['rel']))
            if rel1 > c['rel']:
                rec['a'] = [float(v) for v in c['a']]
        elif c['rel'] < 1e-3:
            rec['a'] = [float(v) for v in c['a']]
        recs.append(rec)
    return recs, time.time() - t


def D_float_str(a):
    Dv = fam.DB.T @ np.asarray(a, float)
    Dm = W.mons_chi(fam.d, fam.c)
    terms = sorted(zip(np.abs(Dv), range(len(Dv))), reverse=True)
    out = []
    for mag, j in terms[:12]:
        if mag < 1e-6 * terms[0][0]:
            break
        e = Dm[j]
        mon = '*'.join(f'{v}**{k}' if k > 1 else v for v, k in zip(['X0','X1','X2','Y0','Y1','Y2'], e) if k)
        out.append(f'{Dv[j]:+.6f}*{mon}')
    return ' '.join(out)


def try_exact(a, rel):
    """Rational reconstruction with increasing denominators + exact re-solve."""
    tried = []
    for md in (10, 100, 1000, 10000, 100000):
        for tol in (1e-6, 1e-8, 1e-10):
            for fr in fam.reconstruct(a, max_den=md, tol=tol):
                Dx = fam.D_exact(fr)
                Ds = W.D_str(Dx)
                if Ds in tried or not Dx:
                    continue
                tried.append(Ds)
                st, data = fam.confirm(fr)
                print(f'[exact] max_den={md} tol={tol} -> {st} rank={data.get("rank")} nullity={data.get("nullity")} D={Ds[:100]}', flush=True)
                if st == 'solved':
                    return st, data, fr
    return 'none', {'tried': len(tried)}, None


per_job = f'/tmp/hodge/results/off-paired/wave2-{args.name}.jsonl'
cand_path = '/tmp/hodge/results/off-paired/wave2-candidates.jsonl'
tasks = []
nleft = args.restarts; k = 0
while nleft > 0:
    n = min(args.chunk, nleft)
    tasks.append((args.seed * 100003 + k * 7919, n)); nleft -= n; k += 1
rels = []; best = None; tsearch = 0.0; nspur = 0; ncand = 0
with mp.get_context('fork').Pool(args.workers) as pool, open(per_job, 'a') as fj:
    for recs, dt in pool.imap_unordered(work, tasks):
        tsearch += dt
        for r in recs:
            r['d'] = args.d; r['chi'] = list(chi); r['kind'] = args.kind; r['hess'] = args.hess
            fj.write(json.dumps(r) + '\n')
            rels.append(r['rel']); nspur += r['spurious']
            if best is None or r['rel'] < best['rel']:
                best = r
            if r['rel'] < args.cand_tol and not r['spurious']:
                ncand += 1
                with open(cand_path, 'a') as fc:
                    fc.write(json.dumps({'name': args.name, **r, 'D_float': D_float_str(r['a'])}) + '\n')
                if r['rel'] < 1e-10:
                    st, data, fr = try_exact(np.array(r['a']), r['rel'])
                    if st == 'solved':
                        wit = {'d': args.d, 'D': W.D_str(data['D']), 'N': [W.D_str(n) for n in data['N']],
                               'Q': W.D_str(data['Q']), 'name': args.name, 'chi': list(chi), 'kind': args.kind}
                        ok, vinfo = core.verify_witness(F, A, q, data['D'], data['N'], data['Q'])
                        wit['verify'] = ok
                        with open('/tmp/hodge/results/off-paired/witness-candidate.json', 'w') as fw:
                            json.dump(wit, fw, indent=1)
                        W.log_batch({'step': 'WITNESS', 'name': args.name, 'd': args.d, 'chi': list(chi),
                                     'kind': args.kind, 'D': wit['D'], 'verify': ok})
                        print('WITNESS FOUND', wit['D'], flush=True)
        fj.flush()
        done = len(rels)
        print(f'[progress] {done}/{args.restarts} best rel {min(rels):.3e} median {float(np.median(rels)):.3e} '
              f'elapsed {time.time()-t0:.0f}s', flush=True)

rels_arr = np.array(rels)
lr = np.log10(np.maximum(rels_arr, 1e-300))
rec = {'step': 'varpro', 'name': args.name, 'd': args.d, 'chi': list(chi), 'kind': args.kind, 'hess': args.hess,
       'restarts': len(rels), 'maxiter': args.maxiter, 'seed': args.seed, 'workers': args.workers,
       'best_rel': float(rels_arr.min()), 'median_rel': float(np.median(rels_arr)),
       'q10_log10': float(np.quantile(lr, 0.1)), 'q90_log10': float(np.quantile(lr, 0.9)),
       'n_below_1e-4': int(np.sum(rels_arr < 1e-4)), 'n_below_1e-7': int(np.sum(rels_arr < 1e-7)),
       'n_below_1e-10': int(np.sum(rels_arr < 1e-10)), 'n_spurious': int(nspur), 'n_candidates': ncand,
       'seconds_wall': round(time.time() - t0, 1), 'seconds_cpu_search': round(tsearch, 1),
       'family': info, 'best_a': best.get('a') if best and 'a' in best else None,
       'per_restart_file': os.path.basename(per_job)}
W.log_batch(rec)
os.system('python3 /tmp/hodge/results/off-paired/w2_report.py > /dev/null 2>&1')
print('[done]', json.dumps({k: v for k, v in rec.items() if k not in ('family', 'best_a')}), flush=True)
