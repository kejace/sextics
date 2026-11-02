#!/usr/bin/env python3
"""w2s_search.py -- random-restart joint search over a linear denominator family (wave 2, items 2-3, 5).

usage: python3 w2s_search.py --family J8 --restarts 200 --procs 8 --seed 100 [--S 1700] [--max-nfev 300]
                             [--method lm|lbfgs] [--tag name] [--stage1]
families: J<d>  (Jacobian-type D = sum w_i F_i, d = 5..9)
          K10   (numerical E1-kernel family from w2s_K10_family.npz; --sub nonJ|all|rand<k>)
Per restart one JSON line goes to w2s_runs_<tag>.jsonl immediately; a batch line goes to wave2-search.log.jsonl at the
end (and every --checkpoint restarts).  Candidates with rel < 1e-10 are rationally reconstructed and re-solved EXACTLY
(mod p reduced system); a consistent one is saved to witness-candidate-<tag>.json and flagged WITNESS-CANDIDATE.
--stage1 = Hessian relaxation (c = 0: minimise |p0(a)|^2 only, no Koszul unknowns).
"""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import numpy as np
import multiprocessing as mp

WD = '/tmp/hodge/results/quadratic-cubic-pullback'


def load_family(name, sub='all'):
    import w2s
    if name.startswith('J'):
        return w2s.jacobian_family(int(name[1:]))
    if name == 'K10':
        z = np.load(os.path.join(WD, 'w2s_K10_family.npz'))
        Mmat = z['Mmat']; Umat = [z[f'U{i}'] for i in range(6)]
        d = int(z['d'])
        if sub == 'nonJ':
            k = int(z['n_nonJ']); Mmat = Mmat[:, :k]; Umat = [U[:, :k] for U in Umat]
        elif sub.startswith('rand'):
            k = int(sub[4:]); rng = np.random.default_rng(12345)
            R = rng.standard_normal((Mmat.shape[1], k))
            Mmat = Mmat @ R; Umat = [U @ R for U in Umat]
        return w2s.Family(d, Mmat, Umat, name=f'K10-{sub}')
    raise ValueError(name)


class Stage1Wrapper:
    """c = 0 relaxation: residual p0(a) only (used with L-BFGS on the sphere)."""

    def __init__(self, pv):
        self.pv = pv

    def value_and_grad(self, uu):
        pv = self.pv
        a_full = pv.P @ uu; na = np.linalg.norm(a_full); a = a_full / na
        D, dD, u, divu = pv._state(a)
        p0 = D * divu - sum(u[i] * dD[i] for i in range(6))
        f = float(np.sum(np.abs(p0) ** 2))
        Pc = np.conj(p0)
        grad = (pv.TM @ (Pc * divu) + pv.Tdivu @ (Pc * D) - sum(pv.Tu[i] @ (Pc * dD[i]) for i in range(6))
                - sum(pv.TdM[i] @ (Pc * u[i]) for i in range(6)))
        grad = 2.0 * grad.real
        g_full = (grad - float(grad @ a) * a) / na
        return f, pv.P.T @ g_full


def worker(args):
    (famname, sub, S, seed, n_restarts, max_nfev, method, tag, stage1, wid) = args
    import w2s
    from scipy.optimize import minimize
    fam = load_family(famname, sub)
    pv = w2s.PointVarPro(fam, S=S, seed=1000 + seed)
    runfile = os.path.join(WD, f'w2s_runs_{tag}.jsonl')
    t0 = time.time()
    rels = []
    best = None
    cands_conf = []

    def record(c, extra=None):
        rec = {'family': fam.name, 'd': fam.d, 'K': fam.K, 'S': pv.S, 'nc': pv.nc_eff, 'seed': seed,
               'restart': c['restart'], 'rel': c['rel'], 'f': c['value'], 'spurious': c['spurious'],
               'nfev': c.get('nfev'), 'method': method, 'stage1': stage1, 'wid': wid,
               'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
        if extra:
            rec.update(extra)
        with open(runfile, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        return rec

    def confirm(c):
        """Rational reconstruction + exact reduced solve for a candidate (J-families only; K10 handled separately)."""
        out = {'confirm': []}
        if fam.exact is None:
            # reconstruct monomial coefficients of D directly
            Dv = fam.Mmat @ c['a']
            lst, _ = w2s.mons(fam.d, 0)
            fr = w2s.reconstruct_vector(Dv, max_dens=(1, 2, 3, 4, 6, 12, 24, 60, 120, 1000, 10000, 100000), tol=1e-6)
            Ds = []
            for f_ in fr:
                D = {lst[j]: f_[j] for j in range(len(lst)) if f_[j] != 0}
                Ds.append(D)
        else:
            Ds = [w2s.D_exact_from_family(fam, f_) for f_ in w2s.reconstruct_vector(c['a'])]
        for D in Ds[:3]:
            if not D:
                continue
            r = w2s.exact_reduced_solve(D, return_solution=True)
            rr = {k: v for k, v in r.items() if k not in ('solution', 'u')}
            rr['D'] = w2s.core.to_str(D)[:4000]
            out['confirm'].append(rr)
            if r.get('consistent'):
                divF, _ = w2s.jacobian.divisible_by(w2s.F, D)
                rr['F_divides_D'] = bool(divF)
                if not divF:
                    path = os.path.join(WD, f'witness-candidate-{tag}-{seed}-{c["restart"]}.json')
                    with open(path, 'w') as f:
                        json.dump({'d': fam.d, 'D': w2s.core.to_str(D), 'p': r['p'], 'note': 'mod-p consistent; needs rational lift',
                                   'rel': c['rel'], 'family': fam.name}, f)
                    rr['WITNESS-CANDIDATE'] = path
                    print('WITNESS-CANDIDATE', path, file=sys.stderr, flush=True)
                break
        return out

    def cb(c):
        nonlocal best
        rels.append(c['rel'])
        extra = None
        if c['rel'] < 1e-10 and not c['spurious']:
            extra = confirm(c)
            cands_conf.append(extra)
        record(c, extra)
        if best is None or (not c['spurious'] and c['rel'] < best['rel']):
            best = c

    if stage1:
        rng = np.random.default_rng(seed)
        w = Stage1Wrapper(pv)
        for it in range(n_restarts):
            u0 = rng.standard_normal(pv.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(w.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': max_nfev, 'gtol': 1e-14, 'ftol': 1e-18, 'maxcor': 30})
            a = pv.P @ res.x; a /= np.linalg.norm(a)
            f1, _ = w.value_and_grad(res.x)
            D, dD, u, divu = pv._state(a)
            scale_ = float(np.sum(np.abs(D * divu) ** 2) + sum(np.sum(np.abs(u[i] * dD[i]) ** 2) for i in range(6)))
            spur, _ = pv.is_spurious(a)
            c = {'value': f1, 'rel': f1 / max(scale_, 1e-300), 'a': a, 'spurious': spur, 'nfev': int(res.nfev),
                 'restart': it}
            cb(c)
    elif method == 'lm':
        pv.search_lm(n_restarts=n_restarts, seed=seed, max_nfev=max_nfev, callback=cb, stop_rel=-1)
    else:
        pv.search(n_restarts=n_restarts, seed=seed, maxiter=max_nfev, callback=cb, stop_rel=-1)
    rels = np.array(rels)
    summary = {'family': fam.name, 'd': fam.d, 'K': fam.K, 'S': pv.S, 'nc': pv.nc_eff, 'seed': seed, 'wid': wid,
               'restarts': int(len(rels)), 'best_rel': float(rels.min()) if len(rels) else None,
               'median_rel': float(np.median(rels)) if len(rels) else None,
               'quantiles_rel': [float(x) for x in np.quantile(rels, [0.0, 0.1, 0.5, 0.9, 1.0])] if len(rels) else None,
               'n_below_1e-6': int(np.sum(rels < 1e-6)), 'n_below_1e-10': int(np.sum(rels < 1e-10)),
               'confirmed': cands_conf, 'seconds': round(time.time() - t0, 1), 'method': method, 'stage1': stage1}
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--sub', default='all')
    ap.add_argument('--restarts', type=int, default=50, help='per worker')
    ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--S', type=int, default=0, help='points (0 = auto)')
    ap.add_argument('--max-nfev', type=int, default=300)
    ap.add_argument('--method', default='lm')
    ap.add_argument('--tag', default=None)
    ap.add_argument('--stage1', action='store_true')
    a = ap.parse_args()
    import w2s
    tag = a.tag or f'{a.family}{a.sub if a.sub != "all" else ""}-{a.method}{"-s1" if a.stage1 else ""}-{a.seed}'
    fam = load_family(a.family, a.sub)
    # auto S: 2S real equations >= 1.6 (nc + K)
    d = fam.d
    cpar = w2s.c_parities()
    nc = sum(len(w2s.mons(d - 4, par)[0]) for par in cpar.values())
    S = a.S or int(0.8 * (nc + fam.K)) + 50
    print(f'[search] family {fam.name} d={d} K={fam.K} nc={nc} S={S} procs={a.procs} restarts/worker={a.restarts} '
          f'method={a.method} stage1={a.stage1} tag={tag}', file=sys.stderr, flush=True)
    jobs = [(a.family, a.sub, S, a.seed + w, a.restarts, a.max_nfev, a.method, tag, a.stage1, w) for w in range(a.procs)]
    t0 = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(a.procs) as pool:
        sums = pool.map(worker, jobs, chunksize=1)
    rels_best = [s['best_rel'] for s in sums if s['best_rel'] is not None]
    batch = {'item': 'search', 'family': fam.name, 'sub': a.sub, 'd': d, 'K': fam.K, 'nc': nc, 'S': S, 'method': a.method,
             'stage1': a.stage1, 'tag': tag, 'seed0': a.seed, 'procs': a.procs,
             'restarts_total': sum(s['restarts'] for s in sums), 'best_rel': min(rels_best) if rels_best else None,
             'median_of_medians': float(np.median([s['median_rel'] for s in sums if s['median_rel'] is not None])) if sums else None,
             'n_below_1e-6': sum(s['n_below_1e-6'] for s in sums), 'n_below_1e-10': sum(s['n_below_1e-10'] for s in sums),
             'confirmed': [c for s in sums for c in s['confirmed']], 'seconds': round(time.time() - t0, 1),
             'workers': [{k: v for k, v in s.items() if k != 'confirmed'} for s in sums]}
    w2s.log_batch(batch)
    print(json.dumps({k: v for k, v in batch.items() if k != 'workers'}), flush=True)


if __name__ == '__main__':
    main()
