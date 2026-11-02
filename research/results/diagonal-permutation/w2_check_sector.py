"""Sanity check of the (01)-sector reduction: for random (01)-equivariant w (D invariant), compare min_c residual with
full c, with c restricted to the anti-equivariant sector (expected ~equal to full) and to the equivariant sector
(expected much larger).  Also the same for anti-equivariant w.  usage: python3 w2_check_sector.py d"""
import os, sys, time, json
os.environ['OPENBLAS_NUM_THREADS'] = '4'
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2pt import *
from w2sym import *
d = int(sys.argv[1])
F, A, q, rec = core.load_instance('diagonal-permutation')
K0 = sum(len(monomials(d - 5, PAR[i])) for i in range(NV))
pA = core.parity_of(A)
nK = sum(len(monomials(d - 4, ((pA - PAR[i]) % 2 - PAR[j]) % 2)) for i in range(NV) for j in range(i + 1, NV))
pts = Points(F, K0 + nK + 60, seed=5)
fam = jacobian_family(F, A, d, pts); kb = KoszulBlock(F, A, d, pts)
wb = [m['w'] for m in fam.members]
rng = np.random.default_rng(1)
out = {'d': d}
for wsign in (+1, -1):
    W = w_sector_matrix(wb, wsign)
    a = W @ rng.standard_normal(W.shape[1]); a /= np.linalg.norm(a)
    vals = {}
    for name, cs in (('full', None), ('c-anti', c_sector_matrix(kb.cols, -1)), ('c-equi', c_sector_matrix(kb.cols, +1))):
        vp = VarPro(fam, kb, pts, csector=cs)
        f, rel, c, rho = vp.value(a)
        vals[name] = rel
    out[f'w{wsign:+d}'] = vals
    print(f'w-sector {wsign:+d}: ', {k: f'{v:.4e}' for k, v in vals.items()}, flush=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'w2-check-sector.jsonl'), 'a') as fh:
    fh.write(json.dumps(out) + '\n')
