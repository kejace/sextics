"""Cross-check of the character-reduced varpro objective against the engine's unreduced varpro on the same family:
the invariant J^G_6 family (4 basis polys, dimP = 3 after removing F) and the sgnX d=7 family (D7 + J^{sgnX}_7).
Both objectives should have their minimisers at the same D (up to scale); the rel values differ by row weighting."""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
from hodge import varpro
JOBLOG = f'{RES}/w2b/xcheck.log.jsonl'
F, A, q = load()
for chi, d, pD in [('triv', 6, 0), ('sgnX', 7, 0)]:
    t0 = time.time()
    G = GROUPS[chi]()
    L = live_reps(d, pD, G); reps = [r for r, _ in L]; sos = [P for _, P in L]
    dimK, kvecs, _, _ = e1_family_modp(F, A, q, sos, G, pD)
    # exact rational kernel basis: lift the mod-p rref rows (small rationals expected) -> use float kernel instead
    Kb, finfo = e1_family_float(F, A, q, sos, G, pD, expect=dimK)
    # family basis for both engines: float combos of orbit sums (rational approximations with small denominators)
    basis = []
    for j in range(Kb.shape[1]):
        v = Kb[:, j]; v = v / np.max(np.abs(v))
        P = {}
        for k, B in enumerate(sos):
            c = Fraction(float(v[k])).limit_denominator(10**6)
            if c: P = core.add(P, core.scale(B, c))
        basis.append(P)
    mine = CharFamily(F, A, q, G, basis, verbose=True, label=f'xcheck {chi}/{d}/{pD}')
    Fv = generic_f_multiples(F, basis)
    if Fv.shape[1]:
        U, s, _ = np.linalg.svd(Fv, full_matrices=True); mine.set_subspace(U[:, Fv.shape[1]:])
    eng = varpro.VarProFamily(F, A, q, basis, exclude_F_multiples=True, verbose=True)
    print('dims: mine', mine.dimP, 'engine', eng.dimP, flush=True)
    cm = mine.search(n_restarts=20, seed=1, maxiter=400)
    ce = eng.search(n_restarts=4, seed=1, maxiter=200, verbose=True)
    am = cm[0]['a']; ae = ce[0]['a']
    # evaluate each engine at the other's minimiser
    rel_e_at_m = eng.value(am)[0] / max(eng.bnorm(am) ** 2, 1e-300)
    rel_m_at_e = mine.rel(ae)[0]
    cos = float(abs(am @ ae) / (np.linalg.norm(am) * np.linalg.norm(ae)))
    rec = {'job': 'xcheck', 'family': 'xcheck-reduced-vs-engine', 'chi': chi, 'd': d, 'pD': pD, 'K': len(basis), 'dimP': mine.dimP,
           'mine_best_rel': cm[0]['rel'], 'mine_minima': sorted(set(f"{c['rel']:.3g}" for c in cm)), 'engine_best_rel': ce[0]['rel'],
           'engine_minima': sorted(set(f"{c['rel']:.3g}" for c in ce)), 'engine_rel_at_my_min': rel_e_at_m, 'my_rel_at_engine_min': rel_m_at_e,
           'cos_angle_between_minimisers': cos, 'engine_rows_cols': (eng.rows, eng.ncols), 'mine_rows_cols': (mine.rows, mine.ncols), 'seconds': round(time.time() - t0, 1)}
    wlog(rec, JOBLOG); print(json.dumps(rec, default=str), flush=True)
print('ALL DONE', flush=True)

# ---------------------------------------------------------------- polish checks of the best restarts of running/finished sectors
import glob
for chi, d, pD, job in [('triv', 10, 0, 'vp2b'), ('triv', 9, 0, 'vp1b'), ('triv', 8, 0, 'vp1b'), ('triv', 8, 1, 'vp1b'), ('triv', 10, 1, 'vp2b'), ('triv', 9, 1, 'vp1b')]:
    t0 = time.time()
    recs = []
    for fn in glob.glob(f'{RES}/w2b/{job}_{chi}-{d}-{pD}_w*.jsonl'):
        for l in open(fn):
            try: recs.append(json.loads(l))
            except Exception: pass
    if not recs:
        continue
    recs.sort(key=lambda r: r['rel']); best = recs[0]
    try:
        fam, info = build_sector_family(F, A, q, chi, d, pD, label=f'{chi}/d{d}/p{pD}')
        if best.get('sparse_init'):
            c = fam.search(n_restarts=1, seed=best['seed'], maxiter=600)[0]   # sparse init not reproducible; plain init from the seed
        else:
            c = fam.search(n_restarts=1, seed=best['seed'], maxiter=600)[0]
        a, x, rel, hist = polish(fam, c['a'], iters=10)
        c2 = fam.search(n_restarts=1, seed=0, maxiter=2000, inits=(a,))[0]
        a2, x2, rel2, hist2 = polish(fam, c2['a'], iters=10)
        j = int(np.argmax(np.abs(a2))); an = a2 / a2[j]
        rec = {'job': 'xcheck', 'family': 'polish-check', 'chi': chi, 'd': d, 'pD': pD, 'seed': best['seed'], 'remote_rel': best['rel'],
               'remote_iters': best.get('iters'), 'rerun_rel': c['rel'], 'rerun_iters': c['iters'], 'polished_rel': rel, 'polish_hist': [f'{h:.2e}' for h in hist],
               'lbfgs2000_rel': c2['rel'], 'lbfgs2000_iters': c2['iters'], 'polished2_rel': rel2, 'dimP': fam.dimP,
               'dominant': [(list(fam.reps[k]), round(float(an[k]), 5)) for k in np.argsort(-np.abs(an))[:8]], 'seconds': round(time.time() - t0, 1)}
        if rel2 < 1e-10:
            ok, res = try_confirm(fam, {'a': a2, 'rel': rel2}, dict(rec, event='confirm'))
            rec['confirm'] = {'ok': ok, 'res': res}
    except Exception as ex:
        rec = {'job': 'xcheck', 'family': 'polish-check', 'chi': chi, 'd': d, 'pD': pD, 'error': repr(ex)[:300]}
    wlog(rec, JOBLOG); print(json.dumps(rec, default=str), flush=True)
print('POLISH CHECKS DONE', flush=True)
