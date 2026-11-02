"""Deep probe of the mod-J relaxation in one sector: multistart L-BFGS, LM polish of the best minima, full-E2 LM from
relaxation zeros.  python3 -u w2_relax_probe.py d parD k restarts [family J|KJ] [maxiter]"""
import os, sys, time, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import w2lib as L
from scipy.optimize import least_squares
F, A, q = L.load()
d, parD, k, nres = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
family = sys.argv[5] if len(sys.argv) > 5 else 'J'
maxiter = int(sys.argv[6]) if len(sys.argv) > 6 else 1500
t0 = time.time()
if family == 'J':
    fam = L.jacobian_family(F, A, d, parD, k)
else:
    from w2_varpro import build_family
    fam = build_family(d, parD, k, family, L.LANE + '/w2_sectors.pkl')
kv = L.KoszulVarPro(F, A, fam, restrict=True, k=k)
rx = L.Relax(kv, verbose=True)
out = rx.search(nres, 12345 + d * 100 + k, maxiter=maxiter)
rels = np.array([c['rel'] for c in out])
print(f'relax ({d},{parD},{k}) family {family}: K={kv.K} dimP={kv.dimP} nl={rx.nl} restarts={nres} best={rels.min():.3e} p10={np.percentile(rels,10):.3e} median={np.median(rels):.3e} [{time.time()-t0:.0f}s]', flush=True)

def fun(u):
    a = kv.P @ u; a /= np.linalg.norm(a); return rx.q_and_grad(a)[0]
def jac(u):
    a_full = kv.P @ u; na = np.linalg.norm(a_full); a = a_full / na; qv, dq = rx.q_and_grad(a); return ((dq - np.outer(dq @ a, a)) / na) @ kv.P
polished = []
seen = []
for c in out[:8]:
    u0 = kv.P.T @ c['a']
    res = least_squares(fun, u0, jac=jac, method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    a = kv.P @ res.x; a /= np.linalg.norm(a)
    rr = rx.rel(a)[0]; e2 = kv.rel(a)[0]
    dup = any(abs(abs(float(a @ s)) - 1) < 1e-6 for s in seen)
    seen.append(a)
    polished.append({'relax_rel': rr, 'e2_rel': e2, 'a': [float(x) for x in a], 'dup': dup})
    print(f'  minimum: relax rel {c["rel"]:.3e} -> polished {rr:.3e}; full-E2 rel at this point {e2:.3e}; {"(dup)" if dup else ""} top|a| {np.round(np.sort(np.abs(a))[::-1][:5],3)}', flush=True)
zeros = [p for p in polished if p['relax_rel'] < 1e-14 and not p['dup']]
e2_from_zero = []
for z in zeros[:5]:
    cand = kv.search_lm(n_restarts=1, init=np.array(z['a']), max_nfev=400)[0]
    e2_from_zero.append({'e2_rel': cand['rel'], 'relax_rel_end': rx.rel(cand['a'])[0], 'spurious': cand['spurious'], 'a': [float(x) for x in cand['a']]})
    print(f'  E2 LM from relaxation zero: rel {cand["rel"]:.3e} spurious={cand["spurious"]} relax at end {rx.rel(cand["a"])[0]:.2e}', flush=True)
    if cand['rel'] < 1e-10 and not cand['spurious']:
        from w2_varpro import try_candidate
        print('  CANDIDATE:', try_candidate(kv, cand, d, parD, k, 'relax-probe', None), flush=True)
L.wlog({'family': f'hess-relax-probe-{family}', 'd': d, 'parD': parD, 'k': k, 'K': kv.K, 'dimP': kv.dimP, 'n_equations': rx.nl,
        'restarts': nres, 'maxiter': maxiter, 'best_rel': float(rels.min()), 'median_rel': float(np.median(rels)),
        'polished': [{'relax_rel': p['relax_rel'], 'e2_rel': p['e2_rel'], 'dup': p['dup']} for p in polished],
        'n_relax_zeros': len(zeros), 'e2_from_zero': [{kk: v for kk, v in z.items() if kk != 'a'} for z in e2_from_zero],
        'seconds': round(time.time() - t0, 1), 'note': 'deep probe: LM-polished relaxation minima and full-E2 LM started from exact relaxation zeros'})
json.dump({'polished': polished, 'e2_from_zero': e2_from_zero}, open(f'{L.LANE}/w2-relaxprobe-d{d}-p{parD}-k{k}-{family}.json', 'w'))
print('PROBE_DONE', flush=True)
