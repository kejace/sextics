import sys, time
sys.path.insert(0, '.')
from hodge import core, fixed_d
for iid in ['collision-h3', 'hesse-triangle', 'xtype-sparse-rational-1m142', 'contact46006']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True)
    print(f'{iid}: d={d} consistent={out["consistent"]} rank={out["rank"]} nullity={out["nullity"]} exact={out.get("exact")} ({time.time()-t0:.1f}s)')
    D2 = dict(D); e = next(iter(D2)); D2[e] = D2[e] + 1
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, D2, want_exact=False)
    print(f'{iid}: perturbed consistent={out["consistent"]} rank={out["rank"]} ({time.time()-t0:.1f}s)')
