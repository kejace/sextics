"""Wave-2 joint search driver (sector-restricted Koszul variable projection), multiprocess.

  python3 -u w2_varpro.py --tasks "6:0:0-6,6:1:0-6" --family J --restarts 400 --chunk 50 --seed 1000 \
        --method lm --maxiter 300 --procs 8 --tag vp6 [--reg 1e-8] [--pkl w2_sectors.pkl]

tasks: comma-separated d:parD:k-range items.  family J = J(F)_d sector generators (Jacobian type, complete for
d <= 8 by the E1-kernel theorem); KJ = kernel reps (from the sector pickle) + J generators (complete E1-admissible
sector space for sigma-eigen D of that parity); K = kernel reps only.
Each (task, chunk) runs in a worker process with 1 BLAS thread; after every chunk a JSON line is appended to
wave2-search.log.jsonl and the per-task state file w2-<tag>-d<d>-p<parD>-k<k>.json is rewritten.
Candidates with rel < 1e-10 (non-spurious) are rationally reconstructed (rem_F chart, increasing denominators)
and re-solved exactly; a verified witness is written to witness-candidate.json.
"""
import os, sys, time, json, math, argparse, pickle, traceback
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['MKL_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import w2lib as L
from hodge import core
from multiprocessing import Pool

LANE = L.LANE
F, A, q = L.load()
STORE = None


def parse_tasks(s):
    out = []
    for item in s.split(','):
        item = item.strip()
        if not item:
            continue
        d, parD, ks = item.split(':')
        if '-' in ks:
            a, b = ks.split('-'); kl = list(range(int(a), int(b) + 1))
        else:
            kl = [int(x) for x in ks.split('+')]
        for k in kl:
            out.append((int(d), int(parD), k))
    return out


def build_family(d, parD, k, family, pkl):
    fam = []
    if family in ('KJ', 'K'):
        store = pickle.load(open(pkl, 'rb'))
        info = store.get(str((d, parD, k)))
        if info is None:
            raise RuntimeError(f'sector {(d, parD, k)} missing in {pkl}')
        reps, lifts = info['reps'], info['lifts']
        if (len(reps) < info.get('nontrivial', 0) or any(N is None for N in lifts)) and info.get('nontrivial', 0) > 0:
            # incomplete rational representatives (large height): use the orthonormal floating-point kernel basis
            if info.get('reps_float') and len(info['reps_float']) == info['nontrivial']:
                reps, lifts = info['reps_float'], info['lifts_float']
            else:
                reps, lifts, finfo = L.float_kernel_reps(F, A, d, parD, k)
            if len(reps) != info['nontrivial']:
                raise RuntimeError(f'float kernel reps {len(reps)} != nontrivial {info["nontrivial"]}')
        for idx, (D, N) in enumerate(zip(reps, lifts)):
            if N is None:
                Nf, relf = L.float_e1_lift(F, A, D, restrict=True)
                if relf > 1e-9:
                    raise RuntimeError(f'float lift failed for rep {idx}: rel {relf}')
                N = Nf
            fam.append((D, N, ('K', idx)))
    if family in ('KJ', 'J'):
        fam += L.jacobian_family(F, A, d, parD, k, restrict=True)
    return fam


def try_candidate(kv, cand, d, parD, k, tag, logf):
    """Rational reconstruction + exact confirmation of a numerical candidate.  Returns dict (status etc)."""
    a = cand['a']
    results = []
    tried = set()
    for max_den in (10, 100, 1000, 10**4, 10**5, 10**6):
        tol = 1e-7 if max_den <= 1000 else 1e-9
        for Dc in kv.reconstruct_D(a, max_den=max_den, tol=tol):
            s = core.to_str(Dc)
            if s in tried:
                continue
            tried.add(s)
            t = time.time()
            try:
                st, data = L.confirm_exact(F, A, Dc, restrict=True)
            except Exception as ex:
                st, data = 'error', {'error': repr(ex)}
            rec = {'max_den': max_den, 'status': st, 'rank': data.get('rank'), 'nullity': data.get('nullity'),
                   'seconds': round(time.time() - t, 1), 'terms': len(Dc), 'D': s[:400]}
            results.append(rec)
            L.wlog({'family': f'candidate-confirm:{tag}', 'd': d, 'parD': parD, 'k': k, 'rel': cand['rel'], **rec,
                    'note': 'exact re-solve of a rationally reconstructed varpro candidate (rem_F chart)'})
            if st == 'solved':
                N, Q = data['N'], data['Q']
                ok, info = core.verify_witness(F, A, 1, Dc, N, Q)
                wit = {'d': d, 'D': core.to_str(Dc), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q),
                       'verify': ok, 'info': str(info), 'tag': tag, 'sector': (d, parD, k)}
                path = f'{LANE}/witness-candidate.json'
                json.dump(wit, open(path, 'w'), indent=1)
                L.wlog({'family': 'WITNESS', 'd': d, 'parD': parD, 'k': k, 'D': core.to_str(Dc)[:2000], 'verify': ok,
                        'note': f'WITNESS FOUND by {tag}; saved {path}'})
                print('WITNESS FOUND', path, flush=True)
                return {'status': 'solved', 'results': results, 'witness': path}
    return {'status': 'unconfirmed', 'results': results}


def run_chunk(args):
    (d, parD, k, family, restarts, seed, method, maxiter, reg, tag, pkl, sparse) = args
    t0 = time.time()
    state_path = f'{LANE}/w2-{tag}-d{d}-p{parD}-k{k}.json'
    try:
        fam = build_family(d, parD, k, family, pkl)
        if not fam:
            L.wlog({'family': f'{family}-varpro', 'd': d, 'parD': parD, 'k': k, 'K': 0, 'restarts': 0, 'tag': tag,
                    'note': 'empty family in this sector'})
            return {'task': (d, parD, k), 'empty': True}
        kv = L.KoszulVarPro(F, A, fam, restrict=True, k=k, reg=reg)
        tb = time.time() - t0
        rels = []
        best = None
        near = []
        cands_log = []

        def cb(c):
            rels.append(c['rel'])
            nonlocal best
            if (best is None) or ((c['spurious'], c['rel']) < (best['spurious'], best['rel'])):
                best = c
            if c['rel'] < 1e-6 and not c['spurious']:
                near.append({'rel': c['rel'], 'a': [float(x) for x in c['a']], 'restart': c['restart']})
            if c['rel'] < 1e-10 and not c['spurious']:
                res = try_candidate(kv, c, d, parD, k, tag, None)
                cands_log.append({'rel': c['rel'], 'status': res['status'], 'n_tried': len(res['results'])})
        relax_stats = None
        if method == 'lbfgs':
            out = kv.search(n_restarts=restarts, seed=seed, maxiter=maxiter, callback=cb)
        elif method == 'joint':
            out = kv.search_joint(n_restarts=restarts, seed=seed, max_nfev=maxiter, callback=cb)
        elif method == 'two':
            out, relax_stats = kv.search_two_stage(n_restarts=restarts, seed=seed, max_nfev=maxiter, callback=cb)
        elif method == 'clm':
            out = kv.search_clm(n_restarts=restarts, seed=seed, max_nfev=maxiter, callback=cb)
        else:
            out = kv.search_lm(n_restarts=restarts, seed=seed, max_nfev=maxiter, callback=cb, sparse=sparse)
        rels_np = np.array(rels) if rels else np.array([np.nan])
        nsp = sum(1 for c in out if c['spurious'])
        rec = {'family': f'{family}-varpro', 'd': d, 'parD': parD, 'k': k, 'K': kv.K, 'dimP': kv.dimP, 'rows': kv.rows,
               'cols': kv.ncols, 'nP': kv.nP, 'nQ': kv.nQ, 'excluded_F_dim': kv.excluded_dim, 'restarts': len(rels),
               'best_rel': float(np.nanmin(rels_np)), 'q25_rel': float(np.nanpercentile(rels_np, 25)),
               'median_rel': float(np.nanmedian(rels_np)), 'max_rel': float(np.nanmax(rels_np)),
               'n_lt_1e-3': int(np.sum(rels_np < 1e-3)), 'n_lt_1e-6': int(np.sum(rels_np < 1e-6)),
               'n_lt_1e-10': int(np.sum(rels_np < 1e-10)), 'n_spurious': nsp, 'seed': seed, 'method': method,
               'maxiter': maxiter, 'reg': reg, 'tag': tag, 'sparse': sparse, 'build_s': round(tb, 1), 'seconds': round(time.time() - t0, 1),
               'mean_iters': float(np.mean([c['iters'] for c in out])) if out else None,
               'candidates': cands_log, 'relax': relax_stats,
               'note': f'sector-restricted Koszul varpro over the {family} family (sigma-eigen Y-{"even" if parD == 0 else "odd"} D of weight {k})'}
        L.wlog(rec)
        # per-task state: merge with existing
        st = {}
        if os.path.exists(state_path):
            try:
                st = json.load(open(state_path))
            except Exception:
                st = {}
        st.setdefault('task', [d, parD, k]); st.setdefault('chunks', [])
        st['chunks'].append({'seed': seed, 'restarts': len(rels), 'best_rel': rec['best_rel'], 'median_rel': rec['median_rel'],
                             'seconds': rec['seconds'], 'method': method})
        st['sizes'] = {'K': kv.K, 'dimP': kv.dimP, 'rows': kv.rows, 'cols': kv.ncols, 'nP': kv.nP, 'nQ': kv.nQ}
        tops = st.get('top', [])
        for c in out[:5]:
            tops.append({'rel': c['rel'], 'spurious': c['spurious'], 'a': [float(x) for x in c['a']], 'seed': seed, 'restart': c['restart']})
        tops.sort(key=lambda c: (c['spurious'], c['rel']))
        st['top'] = tops[:10]
        st['near_misses'] = (st.get('near_misses', []) + near)[:50]
        st['total_restarts'] = sum(c['restarts'] for c in st['chunks'])
        st['best_rel'] = min(c['best_rel'] for c in st['chunks'])
        st['hist'] = {}
        allrels = st.get('all_rels', []) + [float(x) for x in rels]
        st['all_rels'] = allrels[-5000:]
        json.dump(st, open(state_path + '.tmp', 'w')); os.replace(state_path + '.tmp', state_path)
        return {'task': (d, parD, k), 'seed': seed, 'best_rel': rec['best_rel'], 'restarts': len(rels), 'seconds': rec['seconds']}
    except Exception as ex:
        tb = traceback.format_exc()
        L.wlog({'family': f'{family}-varpro', 'd': d, 'parD': parD, 'k': k, 'tag': tag, 'seed': seed, 'error': repr(ex),
                'note': 'worker error', 'traceback': tb[-1500:]})
        return {'task': (d, parD, k), 'seed': seed, 'error': repr(ex)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tasks', required=True)
    ap.add_argument('--family', default='J')
    ap.add_argument('--restarts', type=int, default=200)
    ap.add_argument('--chunk', type=int, default=50)
    ap.add_argument('--seed', type=int, default=1000)
    ap.add_argument('--method', default='lm')
    ap.add_argument('--maxiter', type=int, default=300)
    ap.add_argument('--reg', type=float, default=0.0)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--tag', default='vp')
    ap.add_argument('--pkl', default=f'{LANE}/w2_sectors.pkl')
    ap.add_argument('--sparse', type=int, default=0, help='sparse random starts with this many nonzero family coordinates')
    a = ap.parse_args()
    tasks = parse_tasks(a.tasks)
    jobs = []
    for ti, (d, parD, k) in enumerate(tasks):
        nch = max(1, math.ceil(a.restarts / a.chunk))
        for c in range(nch):
            n = min(a.chunk, a.restarts - c * a.chunk)
            jobs.append((d, parD, k, a.family, n, a.seed + 1000 * ti + c, a.method, a.maxiter, a.reg, a.tag, a.pkl, a.sparse))
    # interleave chunks across tasks so every sector progresses
    jobs.sort(key=lambda j: (j[5] % 1000, j[0], j[1], j[2]))
    print(f'{len(tasks)} tasks, {len(jobs)} chunks, procs {a.procs}', flush=True)
    t0 = time.time()
    with Pool(a.procs) as pool:
        for r in pool.imap_unordered(run_chunk, jobs):
            print(json.dumps(r, default=str), f'[{time.time()-t0:.0f}s]', flush=True)
    print('VARPRO_DONE', flush=True)


if __name__ == '__main__':
    main()
