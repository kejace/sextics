import os, sys, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
import numpy as np
from fptlib import *
iid = 'family-hesse-web'
inst0 = Instance(iid); d, Dw, Nw, Qw = core.load_witness(iid)
sub = symmetry.stabilizer(inst0.G, [inst0.A, Dw]); inst = Instance(iid, gens=sub)
lab = label_of_poly(inst, Dw); fam = e1_family(inst, d, lab); svp = SectorVarPro(inst, d, lab, fam['basis'])
res = svp.search(n_restarts=8, seed=9, maxiter=400)
best = res[0]; print('search best rel', best['rel'], flush=True)
t0 = time.time(); a1, r1 = svp.polish_lm(best['a']); print('LM polished rel', r1, round(time.time()-t0, 1), 's', flush=True)
Dfl = svp.D_float(a1)
att = reconstruct_and_confirm(inst, Dfl, tag='local-test')
for x in att: print(x['max_den'], f"{x['recon_err']:.1e}", x['status'], x['nterms'], x['D'][:100])
