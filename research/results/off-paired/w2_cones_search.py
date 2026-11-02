"""Step 4b: variable projection restricted to the E1-admissible cone subfamilies found by w2_cones.py.
usage: python3 -u w2_cones_search.py --dmin 6 --dmax 12 --restarts 200 [--minvars 3] [--name cones-A]"""
import sys, os, json, time, argparse, itertools
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import numpy as np
from fractions import Fraction
import flint
import w2lib as W, w2fast as WF
from hodge import core
ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=6); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--restarts', type=int, default=200); ap.add_argument('--minvars', type=int, default=3)
ap.add_argument('--name', type=str, default='cones'); ap.add_argument('--maxiter', type=int, default=400)
ap.add_argument('--only_nonjac', type=int, default=0)
args = ap.parse_args()
F, A, q = W.load()
names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
subs = [json.loads(l) for l in open('/tmp/hodge/results/off-paired/wave2-cones.jsonl') if l.strip()]
seen = set(); todo = []
for s in subs:
    key = (s['d'], tuple(s['chi']), tuple(s['vars']))
    if key in seen or s['d'] < args.dmin or s['d'] > args.dmax or len(s['vars']) < args.minvars:
        continue
    if args.only_nonjac and s['non_jacobian'] <= 0:
        continue
    seen.add(key); todo.append(s)
# drop subfamilies contained in a larger-variable-set subfamily with the same intersection dimension? keep all (cheap)
print(f'{len(todo)} cone subfamilies', flush=True)
cand_path = '/tmp/hodge/results/off-paired/wave2-candidates.jsonl'
for s in todo:
    d, c, V = s['d'], tuple(s['chi']), [names.index(v) for v in s['vars']]
    t0 = time.time()
    S = W.SectorSpaces(F, A, q, d, c)
    Dm = S.Dmons; K = len(Dm)
    cols = [j for j, m in enumerate(Dm) if all(m[i] == 0 for i in range(6) if i not in V)]
    # exact basis of V ∩ E1 modulo gauge: nullspace over Q of [E1 rows ; -unit rows] -> combinations of E1 vectors
    # supported on V.  Solve: x^T E1 = y^T U  <=>  [E1^T | -U^T] (x; y) = 0.
    E1 = S.E1
    if not E1:
        continue
    rowsM = []
    for jm in range(K):
        rowsM.append([int(v[jm]) for v in E1] + [(-1 if jm == j else 0) for j in cols])
    M = flint.fmpz_mat(rowsM)
    X, nul = M.nullspace()
    polys = []
    for j in range(nul):
        xv = [int(X[i, j]) for i in range(len(E1))]
        if not any(xv):
            continue
        vec = [sum(xv[t] * E1[t][jm] for t in range(len(E1))) for jm in range(K)]
        P = {Dm[jm]: Fraction(vec[jm]) for jm in range(K) if vec[jm]}
        if P:
            polys.append(P)
    # independent modulo gauge
    gvecs = [S._vec(P, Dm) for P in S.gauge]
    sel = S._extend(gvecs, [S._vec(P, Dm) for P in polys], K)
    basis = [polys[j] for j in sel]
    if not basis:
        continue
    L = W.E1Lifter(F, A, q, d, c)
    lifts = []
    for B in basis:
        N, res = L.lift(B); assert res < 1e-8, res; lifts.append(N)
    kos = [N for _, N in W.koszul_generators(F, d, c, W.poly_chi(A))]
    fam = WF.FastReducedVarPro(F, A, q, basis, lifts, kos, c, label=f'cone{s["vars"]}')
    res = fam.search(n_restarts=args.restarts if fam.K > 1 else 1, seed=7, maxiter=args.maxiter)
    rels = np.array([r['rel'] for r in res])
    best = res[0]
    rec = {'step': 'varpro', 'name': f'{args.name}-d{d}-{c[0]}{c[1]}-{"".join(s["vars"])}', 'd': d, 'chi': list(c),
           'kind': 'cone:' + ''.join(s['vars']), 'hess': 0, 'restarts': len(res), 'maxiter': args.maxiter, 'seed': 7,
           'workers': 1, 'best_rel': float(rels.min()), 'median_rel': float(np.median(rels)),
           'q10_log10': float(np.quantile(np.log10(np.maximum(rels, 1e-300)), 0.1)),
           'q90_log10': float(np.quantile(np.log10(np.maximum(rels, 1e-300)), 0.9)),
           'n_below_1e-4': int((rels < 1e-4).sum()), 'n_below_1e-7': int((rels < 1e-7).sum()),
           'n_below_1e-10': int((rels < 1e-10).sum()), 'n_spurious': int(sum(r['spurious'] for r in res)),
           'n_candidates': int((rels < 1e-7).sum()), 'seconds_wall': round(time.time() - t0, 1),
           'family': fam.info(), 'best_a': [float(v) for v in best['a']]}
    W.log_batch(rec)
    print(f"cone d={d} chi={c} V={s['vars']} K={fam.K} best rel {rels.min():.2e} median {np.median(rels):.2e} ({time.time()-t0:.0f}s)", flush=True)
    if best['rel'] < 1e-7 and not best['spurious']:
        with open(cand_path, 'a') as fc:
            fc.write(json.dumps({'name': rec['name'], 'd': d, 'chi': list(c), 'rel': float(best['rel']),
                                 'a': [float(v) for v in best['a']], 'basis': [W.D_str(B) for B in basis]}) + '\n')
        if best['rel'] < 1e-10:
            for md in (10, 100, 1000, 10000):
                ok = False
                for fr in fam.reconstruct(best['a'], max_den=md, tol=1e-8):
                    st, data = fam.confirm(fr)
                    print('  exact', md, st, W.D_str(fam.D_exact(fr))[:120], flush=True)
                    if st == 'solved':
                        wit = {'d': d, 'D': W.D_str(data['D']), 'N': [W.D_str(n) for n in data['N']], 'Q': W.D_str(data['Q'])}
                        json.dump(wit, open('/tmp/hodge/results/off-paired/witness-candidate.json', 'w'), indent=1)
                        W.log_batch({'step': 'WITNESS', 'name': rec['name'], 'd': d, 'chi': list(c), 'D': wit['D']})
                        print('WITNESS FOUND', wit['D'], flush=True); ok = True; break
                if ok:
                    break
os.system('python3 /tmp/hodge/results/off-paired/w2_report.py > /dev/null 2>&1')
print('done', flush=True)
