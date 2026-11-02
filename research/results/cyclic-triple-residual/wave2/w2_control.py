"""Controls: (1) Fermat F0 with A = X0^3 Y0^3, family sum c_i Z_i^6 (d=6, trivial character) -> join-type zeros
exist, the search must find rel ~ 0 and the exact confirm must succeed.  (2) coupled F: f at D8 > 0."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2')
from w2lib import *
import numpy as np
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
A = core.parse('X0**3*Y0**3')
inst0 = Inst(F=F0, A=A)
print('|G| Fermat', len(inst0.G), 'chiA', inst0.chiA, flush=True)
c0 = inst0.chi_of(core.parse('X0**6+X1**6+X2**6'))
t0 = time.time()
meta, Jb, reps = inst0.e1_space(6, c0, nprimes=2)
print('fermat d6 sector', c0, {k: v for k, v in meta.items() if k not in ('J_indep', 'reps')}, round(time.time()-t0, 1), 's', flush=True)
basis, w_of = basis_from_meta(inst0, meta)
for mode in (False, True):
    V = SectorVP(inst0, basis, c0, relax=mode, w_of_basis=w_of if mode else None, verbose=True)
    a = np.array([1.0 if core.to_str(Bk).count('X') else 0.0 for Bk in basis])
    print('mode relax=%s rel at D0 direction' % mode, V.rel(a), flush=True)
    res = V.search(n_restarts=6, seed=3, maxiter=300, verbose=True)
    best = res[0]
    print('best rel', best['rel'], 'a', np.round(best['a'] / np.abs(best['a']).max(), 6), flush=True)
    for fr in V.reconstruct(best['a']):
        st, data = V.confirm(fr)
        print('confirm', fr, '->', st, {k: data.get(k) for k in ('rank', 'rank_aug', 'verified')}, flush=True)
        break
