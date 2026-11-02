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
    ce = eng.search(n_restarts=2, seed=1, maxiter=80, verbose=True)
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

print('ALL DONE', flush=True)
