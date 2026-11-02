"""Wave-2 joint search driver (remote).  Random-restart variable projection over linear denominator families, one
character sector at a time, with the fixed-D unknowns sector-restricted (w2lib.SectorVarPro); candidates with small
residual are polished, rationally reconstructed with growing denominators and re-solved exactly.

usage: python3 w2_search.py --job NAME --family SPEC --degrees 6,7 [--sectors all|"0,0,3;3,3,0"] [--restarts 25]
                            [--batches 40] [--procs 8] [--threads 1] [--maxiter 300] [--seed0 1000] [--instance IID]
family SPEC:  kernelJ            full E1-admissible space {D in S_d^sector : A D in J(F)}   (= kernel + J(F)_d)
              J                  J(F)_d in the sector (Jacobian-type / w-space: D = sum w_i F_i), full E2 via varpro
              hess               Hessian relaxation on the w-space (N = A w, no Koszul correction), E2 only
              vars:X0,X1,Y0      E1-admissible D depending only on the listed variables (cones with vertex the complementary subspace)
              prod:<expr>        E1-admissible D divisible by the given polynomial factor (D = factor * c)
Each batch appends one JSON line to wave2-<job>.jsonl; candidates to wave2-<job>-cands.jsonl.  Resumable (done batches skipped)."""
import os, sys, argparse, json, time, math, hashlib
ap = argparse.ArgumentParser()
ap.add_argument('--job', required=True); ap.add_argument('--family', required=True); ap.add_argument('--degrees', required=True)
ap.add_argument('--sectors', default='all'); ap.add_argument('--restarts', type=int, default=25); ap.add_argument('--batches', type=int, default=40)
ap.add_argument('--procs', type=int, default=8); ap.add_argument('--threads', type=int, default=1); ap.add_argument('--maxiter', type=int, default=300)
ap.add_argument('--seed0', type=int, default=1000); ap.add_argument('--instance', default='four-plus-two-residual')
ap.add_argument('--cand-rel', type=float, default=1e-8); ap.add_argument('--max-K', type=int, default=100000)
ap.add_argument('--method', default='lbfgs', choices=['lbfgs', 'lm'])
args = ap.parse_args()
os.environ['OPENBLAS_NUM_THREADS'] = str(args.threads); os.environ['OMP_NUM_THREADS'] = str(args.threads); os.environ['MKL_NUM_THREADS'] = str(args.threads)
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
from w2lib import *
import numpy as np
import multiprocessing as mp

OUT = f'{LANE}/wave2-{args.job}.jsonl'
CANDS = f'{LANE}/wave2-{args.job}-cands.jsonl'
inst = Instance(args.instance)
DEGS = [int(x) for x in args.degrees.split(',')]
_cache = {}


def parse_sectors(spec, d):
    labs = sorted(inst.labels_in_degree(d))
    if spec == 'all':
        return labs
    want = []
    for part in spec.split(';'):
        t = tuple(int(x) for x in part.split(','))
        if t in labs:
            want.append(t)
    return want


def build_family(spec, d, lab):
    """Returns (obj, info) where obj is SectorVarPro / HessRelax or None."""
    key = (spec, d, lab)
    if key in _cache:
        return _cache[key]
    t0 = time.time()
    info = {'family': spec, 'd': d, 'sector': list(lab)}
    obj = None
    try:
        if spec == 'hess':
            hr = HessRelax(inst, d, lab)
            info.update({'K': hr.K, 'dimP': hr.dimP, 'rows': getattr(hr, 'rows', 0), 'nQ': getattr(hr, 'nQ', 0), 'excluded_F_dim': getattr(hr, 'excluded_dim', 0)})
            obj = hr if hr.dimP > 0 else None
        else:
            V = None
            if spec == 'J':
                V = inst.jac_generators_sector(d, lab)
                if not V:
                    info.update({'K': 0, 'dimP': 0}); _cache[key] = (None, info); return _cache[key]
            elif spec.startswith('vars:'):
                vs = spec[5:].split(',')
                idx = [core.VARS.index(v) for v in vs]
                V = [{m: 1} for m in inst.sector_monomials(d, lab) if all(m[i] == 0 for i in range(NV) if i not in idx)]
                if not V:
                    info.update({'K': 0, 'dimP': 0}); _cache[key] = (None, info); return _cache[key]
            elif spec.startswith('prod:'):
                fac = core.parse(spec[5:])
                lf = label_of_poly(inst, fac); df = core.degree(fac)
                if lf is None or df > d:
                    info.update({'K': 0, 'dimP': 0, 'note': 'factor not an eigenvector or too big'}); _cache[key] = (None, info); return _cache[key]
                V = [mul(fac, {m: Fraction(1)}) for m in inst.sector_monomials(d - df, inst.sublab(lab, lf))]
                if not V:
                    info.update({'K': 0, 'dimP': 0}); _cache[key] = (None, info); return _cache[key]
            fam = e1_family(inst, d, lab, V=V)
            if fam is None or fam['K'] == 0:
                info.update({'K': 0 if fam is None else fam['K'], 'dimP': 0, 'K0': 0 if fam is None else fam['K0']})
                _cache[key] = (None, info); return _cache[key]
            if fam['K'] > args.max_K:
                info.update({'K': fam['K'], 'dimP': -1, 'note': 'K above --max-K; skipped'}); _cache[key] = (None, info); return _cache[key]
            svp = SectorVarPro(inst, d, lab, fam['basis'])
            info.update({'K': svp.K, 'K0': fam['K0'], 'dimP': svp.dimP, 'rows': svp.rows, 'cols': svp.ncols, 'nN': svp.nN, 'nQ': svp.nQ,
                         'full_nN': svp.full_nN, 'full_nQ': svp.full_nQ, 'excluded_F_dim': svp.excluded_dim,
                         'sv_kept_max': fam['sv_kept_max'], 'sv_dropped_min': fam['sv_dropped_min'], 'nnz_Mstack': int(svp.MstackT.nnz)})
            obj = svp if svp.dimP > 0 else None
    except Exception as ex:
        import traceback
        info.update({'error': f'{ex.__class__.__name__}: {ex}', 'traceback': traceback.format_exc()[-1500:], 'dimP': 0})
        obj = None
    info['build_seconds'] = round(time.time() - t0, 1)
    _cache[key] = (obj, info)
    return _cache[key]


def hist(rels):
    edges = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1.0]
    h = {}
    for r in rels:
        lab = None
        for e in edges:
            if r < e:
                lab = f'<{e:.0e}'; break
        lab = lab or '>=1'
        h[lab] = h.get(lab, 0) + 1
    return h


def mon_str(m):
    return '*'.join(f'{core.VARS[i]}^{k}' if k > 1 else core.VARS[i] for i, k in enumerate(m) if k) or '1'


def run_task(task):
    spec, d, lab, batch, seed = task
    obj, info = build_family(spec, d, lab)
    rec = {'kind': 'batch', 'job': args.job, 'instance': args.instance, 'batch': batch, 'seed': seed, 'engine': ENGINE}
    rec.update(info)
    if obj is None:
        rec['note'] = info.get('note', 'empty family (or all F-multiples)'); rec['restarts'] = 0
        return rec if batch == 0 else None
    t0 = time.time()
    cands = []
    if obj.dimP <= 2:
        if batch > 0:
            return None
        if spec == 'hess':
            rels = []
            if obj.dimP == 1:
                a = obj.P[:, 0]; f, vv = obj.value(a); cands = [{'rel': f / max(vv, 1e-300), 'a': a, 'spurious': False}]
            else:
                for t in np.linspace(0, np.pi, 720, endpoint=False):
                    a = obj.P @ np.array([math.cos(t), math.sin(t)]); f, vv = obj.value(a)
                    cands.append({'rel': f / max(vv, 1e-300), 'a': a, 'spurious': False})
                cands.sort(key=lambda c: c['rel'])
                best = cands[:8]
                res = obj.search(n_restarts=8, seed=seed)
                cands = best + res
        else:
            scan = obj.grid_scan()
            cands = [{'rel': r, 'a': a, 'spurious': obj.is_spurious(a)[0]} for r, a in scan[:16]]
            if obj.dimP == 2:
                for r, a in scan[:4]:
                    a1, r1 = obj.polish(a, maxiter=500)
                    cands.append({'rel': r1, 'a': a1, 'spurious': obj.is_spurious(a1)[0]})
        cands.sort(key=lambda c: (c.get('spurious', False), c['rel']))
        rec['mode'] = 'grid'
        rec['restarts'] = len(cands)
    else:
        if args.method == 'lm' and spec != 'hess':
            cands = obj.search_lm(n_restarts=args.restarts, seed=seed, max_nfev=args.maxiter)
        else:
            cands = obj.search(n_restarts=args.restarts, seed=seed, maxiter=args.maxiter)
        rec['mode'] = args.method
        rec['restarts'] = len(cands)
    rels = [c['rel'] for c in cands if not c.get('spurious', False)]
    rec['n_spurious'] = sum(1 for c in cands if c.get('spurious', False))
    if rels:
        rels_sorted = sorted(rels)
        rec['best_rel'] = rels_sorted[0]; rec['median_rel'] = rels_sorted[len(rels_sorted) // 2]
        rec['q10_rel'] = rels_sorted[max(0, len(rels_sorted) // 10)]
        rec['hist'] = hist(rels)
        rec['mean_iters'] = float(np.mean([c.get('iters', 0) for c in cands]))
    # candidates
    nc = 0
    for c in cands:
        if c.get('spurious', False) or c['rel'] >= args.cand_rel:
            continue
        nc += 1
        a1, r1 = (obj.polish_lm(c['a']) if spec != 'hess' else (c['a'], c['rel']))
        if spec != 'hess' and r1 > c['rel']:
            a1, r1 = obj.polish(c['a'])
        Dfl = obj.D_float(a1)
        tag = f'{args.job}:{spec}:d{d}:{lab}:b{batch}'
        crec = {'kind': 'cand', 'job': args.job, 'family': spec, 'd': d, 'sector': list(lab), 'batch': batch, 'seed': seed,
                'rel': c['rel'], 'polished_rel': r1, 'D_float': {mon_str(m): v for m, v in sorted(Dfl.items(), reverse=True)}, 'engine': ENGINE}
        if spec == 'hess':
            crec['note'] = 'Hessian relaxation candidate (N = A w); exact re-solve uses the full fixed-D system'
        att = reconstruct_and_confirm(inst, Dfl, log_path=CANDS, tag=tag)
        crec['confirm'] = [(x['max_den'], x['recon_err'], x['status']) for x in att]
        crec['solved'] = any(x['status'] == 'solved' for x in att)
        log_line(CANDS, crec)
        if crec['solved']:
            rec['WITNESS'] = True
            print('!!!!! WITNESS FOUND', tag, flush=True)
        if nc >= 3:
            break
    rec['n_cands'] = nc
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def done_set():
    s = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                o = json.loads(line)
            except Exception:
                continue
            if o.get('kind') == 'batch':
                s.add((o['family'], o['d'], tuple(o['sector']), o['batch']))
    return s


if __name__ == '__main__':
    done = done_set()
    tasks = []
    for b in range(args.batches):
        for d in DEGS:
            for lab in parse_sectors(args.sectors, d):
                if (args.family, d, lab, b) in done:
                    continue
                seed = args.seed0 + 7919 * b + 101 * d + int(hashlib.sha1(str(lab).encode()).hexdigest(), 16) % 1000
                tasks.append((args.family, d, lab, b, seed))
    print(f'[search] job {args.job}: family {args.family} degrees {DEGS} -> {len(tasks)} tasks ({len(done)} done before), procs {args.procs}', flush=True)
    t0 = time.time()
    n = 0
    with mp.get_context('fork').Pool(args.procs, maxtasksperchild=None) as pool:
        for rec in pool.imap_unordered(run_task, tasks):
            if rec is None:
                continue
            n += 1
            log_line(OUT, rec)
            br = rec.get('best_rel'); br = f'{br:.2e}' if br is not None else '-'
            print(f"[search] {rec['family']} d={rec['d']} sec={rec['sector']} b={rec['batch']} K={rec.get('K')} dimP={rec.get('dimP')} "
                  f"rows={rec.get('rows')} cols={rec.get('cols')} restarts={rec.get('restarts')} best={br} med={rec.get('median_rel', float('nan')):.1e} "
                  f"cands={rec.get('n_cands', 0)} {rec.get('seconds', 0)}s {rec.get('note', '')} {'WITNESS' if rec.get('WITNESS') else ''}", flush=True)
    print(f'[search] done: {n} batch records in {time.time()-t0:.0f}s', flush=True)
