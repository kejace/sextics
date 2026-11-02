"""Structured families suggested by the geometry: cones with the G-stable vertex planes/lines
(ternary D(X), D(Y), D(X1-X0,X2-X0,e1Y), D(e1X,Y1-Y0,Y2-Y0); line vertices D(u1,u2,v1,v2), D(e1X,Y), D(X,e1Y)),
chi-symmetrised, degrees 7..12: restricted E1 kernel, single solves of random kernel points, varpro restarts.
python3 w2_struct.py --job struct --families ternX,ternY,uuv,evv,line_tt,line_sX,line_sY --chis triv,sgnX,sgnY,sgnsgn
   --dmin 7 --dmax 12 --restarts 150 --workers 8"""
import os, sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
import multiprocessing as mp
ap = argparse.ArgumentParser(); ap.add_argument('--job', default='struct'); ap.add_argument('--families', default='ternX,ternY,uuv,evv,line_tt,line_sX,line_sY')
ap.add_argument('--chis', default='triv,sgnX,sgnY,sgnsgn'); ap.add_argument('--dmin', type=int, default=7); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--restarts', type=int, default=150); ap.add_argument('--workers', type=int, default=8); ap.add_argument('--maxiter', type=int, default=400)
args = ap.parse_args()
JOBLOG = f'{RES}/w2b/{args.job}.log.jsonl'; os.makedirs(f'{RES}/w2b', exist_ok=True)
F, A, q = load(); p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
u1, u2 = sub(var(1), var(0)), sub(var(2), var(0))
v1, v2 = sub(var(4), var(3)), sub(var(5), var(3))
e1X, e1Y = e_k(X, 1), e_k(Y, 1)
COORDS = {'ternX': [var(0), var(1), var(2)], 'ternY': [var(3), var(4), var(5)], 'uuv': [u1, u2, e1Y], 'evv': [e1X, v1, v2],
          'line_tt': [u1, u2, v1, v2], 'line_sX': [e1X, var(3), var(4), var(5)], 'line_sY': [var(0), var(1), var(2), e1Y]}

def monomials_in(coords, d):
    n = len(coords); out = []
    def rec(i, left, cur):
        if i == n - 1:
            out.append(cur + [left]); return
        for k in range(left + 1):
            rec(i + 1, left - k, cur + [k])
    rec(0, d, [])
    polys = []
    for ex in out:
        P = {(0,) * 6: Fraction(1)}
        for c, k in zip(coords, ex):
            if k: P = mul(P, pw(c, k))
        polys.append(P)
    return polys

def worker(fam, seeds, out_path, maxiter):
    with open(out_path, 'a') as f:
        for s in seeds:
            def cb(c):
                rec = {'seed': c['seed'], 'rel': c['rel'], 'iters': c['iters'], 'spurious': bool(c['spurious'])}
                if c['rel'] < 1e-6: rec['a'] = [float(x) for x in c['a']]
                f.write(json.dumps(rec) + '\n'); f.flush()
            fam.search(n_restarts=1, seed=s, maxiter=maxiter, callback=cb)

for famname in args.families.split(','):
    coords = COORDS[famname]
    for d in range(args.dmin, args.dmax + 1):
        mons = monomials_in(coords, d)
        for chi in args.chis.split(','):
            G = GROUPS[chi]()
            for pD in (0, 1):
                t0 = time.time()
                sub_ = [P for P in mons if core.parity_of(P) == pD]
                basis = restricted_family(sub_, G, p) if sub_ else []
                label = f'{famname}/{chi}/d{d}/p{pD}'
                rec = {'job': args.job, 'family': 'structured-family', 'fam': famname, 'chi': chi, 'd': d, 'pD': pD, 'n_monomials': len(sub_), 'nD': len(basis)}
                if not basis:
                    rec['note'] = 'empty after symmetrisation'; wlog(rec, JOBLOG); continue
                try:
                    fam, info = build_family(F, A, q, G, basis, label, p=p)
                except Exception as ex:
                    rec['error'] = repr(ex)[:300]; wlog(rec, JOBLOG); print(json.dumps(rec), flush=True); continue
                rec.update({k: info.get(k) for k in ('dimK', 'dimF', 'dimP', 'rows', 'nN', 'nQ', 'build_s')})
                if fam is None or fam.dimP == 0:
                    rec['note'] = 'no admissible D (E1 kernel trivial or only F-multiples)'
                    wlog(rec, JOBLOG); print(json.dumps(rec), flush=True); continue
                # single mod-p solves at 3 random admissible points (float kernel -> rationalised small ints is not exact;
                # use the mod-p kernel vectors instead)
                dimK, kvecs, _, _ = e1_family_modp(F, A, q, basis, G, pD, p)
                rng = np.random.default_rng(d)
                singles = []
                for trial in range(min(3, len(kvecs))):
                    coeff = [int(rng.integers(-2, 3)) for _ in kvecs]
                    D = {}
                    for c, kv in zip(coeff, kvecs):
                        if c:
                            for j, Bj in enumerate(basis):
                                if kv[j]:
                                    D = core.add(D, core.scale(core.reduce_mod(Bj, p), (c * kv[j]) % p, p), p)
                    if not D: continue
                    r = solve_single(F, A, q, D, G, p)
                    singles.append({k: r.get(k) for k in ('e1', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')})
                rec['single_solves'] = singles
                # varpro
                n = min(args.restarts, 40 * fam.dimP); W = min(args.workers, n)
                seeds = list(range(5000, 5000 + n))
                files = [f'{RES}/w2b/{args.job}_{famname}-{chi}-{d}-{pD}_w{i}.jsonl' for i in range(W)]
                ctx = mp.get_context('fork'); procs = []
                for i in range(W):
                    pr = ctx.Process(target=worker, args=(fam, seeds[i::W], files[i], args.maxiter)); pr.start(); procs.append(pr)
                for pr in procs: pr.join()
                recs = []
                for fn in files:
                    if os.path.exists(fn):
                        for line in open(fn):
                            try: recs.append(json.loads(line))
                            except Exception: pass
                good = [r for r in recs if not r['spurious']] or recs
                rels = np.array([r['rel'] for r in good]) if good else np.array([1.0])
                rec.update(restarts=len(recs), best_rel=float(rels.min()), median_rel=float(np.median(rels)),
                           n_below_1e_6=int((rels < 1e-6).sum()), n_below_1e_10=int((rels < 1e-10).sum()), seconds=round(time.time() - t0, 1))
                cands = [r for r in good if r['rel'] < 1e-10 and 'a' in r]
                conf = []
                for c in cands[:10]:
                    ok, res = try_confirm(fam, {'a': np.array(c['a']), 'rel': c['rel']}, dict(rec, event='confirm'))
                    conf.append({'rel': c['rel'], 'ok': ok, 'res': res})
                    if ok: rec['WITNESS'] = True
                rec['confirmations'] = conf
                wlog(rec, JOBLOG); print(json.dumps({k: v for k, v in rec.items() if k != 'confirmations'}, default=str), flush=True)
print('ALL DONE', flush=True)
