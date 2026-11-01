import sys, time
sys.path.insert(0, '.')
from hodge import core, fixed_d
for iid, Dexpr in [('dense-rational', 'X0**6 + X1**6 + X2**6 + 2*Y2**6'), ('rank-one-biquartic', 'X1**6 + Y0**6 + Y1**6 + Y2**6'), ('eight-plane-configuration', '(X0**2 + Y0**2)**3')]:
    F, A, q, rec = core.load_instance(iid)
    D = core.parse(Dexpr)
    t0 = time.time()
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True, want_exact=False)
    print(f'{iid}: D={Dexpr} consistent={out["consistent"]} rank={out["rank"]} nullity={out["nullity"]} rows={out["rows"]} cols={out["cols"]} ({time.time()-t0:.1f}s)', flush=True)
