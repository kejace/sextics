import os, sys, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
import numpy as np
import fptlib as new
import old_w2lib_remote as old
from fptlib import core, symmetry, label_of_poly
iid = 'family-hesse-web'
for name, lib in [('new', new), ('old', old)]:
    t0 = time.time()
    inst0 = lib.Instance(iid); d, Dw, Nw, Qw = core.load_witness(iid)
    sub = symmetry.stabilizer(inst0.G, [inst0.A, Dw]); inst = lib.Instance(iid, gens=sub)
    lab = lib.label_of_poly(inst, Dw); fam = lib.e1_family(inst, d, lab); svp = lib.SectorVarPro(inst, d, lab, fam['basis'])
    print(name, 'built', round(time.time() - t0, 1), 's; K', svp.K, 'dimP', svp.dimP, flush=True)
    t0 = time.time()
    res = svp.search(n_restarts=8, seed=9, maxiter=400)
    print(name, 'search 8 restarts:', round(time.time() - t0, 1), 's; rels', [f"{r['rel']:.1e}" for r in res], 'iters', [r['iters'] for r in res], flush=True)
    # confirmation pipeline on the best
    best = res[0]
    if best['rel'] < 1e-8 and name == 'new':
        a1, r1 = svp.polish(best['a']); Dfl = svp.D_float(a1); print('polished', r1, flush=True)
        for k, (md, err, Dfrac) in enumerate(new.rational_candidates(Dfl)):
            if k >= 6: break
            st, data = new.confirm_exact(inst, Dfrac)
            print(f'  max_den={md} err={err:.1e} nterms={len(Dfrac)} -> {st} rank={data.get("rank")} nullity={data.get("nullity")} note={data.get("note")} verified={data.get("verified")}', flush=True)
            if st == 'solved': break
