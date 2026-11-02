import sys, time, json
sys.path.insert(0,'/tmp/hodge/results/cyclic-triple-residual/wave2')
from w2lib import *
import numpy as np
T=time.time()
def P(*a): print(f'[{time.time()-T:7.1f}s]', *a, flush=True)
inst = Inst()
P('chars', inst.chars, 'chiA', inst.chiA)
meta, Jb, reps = inst.e1_space(8, (0,4))
P('d8 (0,4):', {k:v for k,v in meta.items() if k not in ('J_indep','reps')}, 'nJ', len(Jb), 'nreps', len(reps))
D8 = core.parse(json.load(open('/tmp/hodge/results/cyclic-triple-residual/kernel_d8_0_4_0.json'))['D'])
diff = core.sub(reps[0], D8)
P('rep == D8 ?', diff=={}, 'terms rep', len(reps[0]))
basis, w_of = basis_from_meta(inst, meta)
V = SectorVP(inst, basis, (0,4), verbose=True)
P('built')
a = np.zeros(V.K); a[-1]=1.0
P('rel at D8:', V.rel(a))
rng=np.random.default_rng(0)
for _ in range(3):
    a=rng.standard_normal(V.K); P('rel random:', V.rel(a)[0])
res = V.search(n_restarts=3, seed=1, maxiter=200, verbose=True)
P('search 3 restarts; best rel', res[0]['rel'], 'nfev', V.nfev)
# relax mode
Vr = SectorVP(inst, basis, (0,4), relax=True, w_of_basis=w_of, verbose=True)
res = Vr.search(n_restarts=3, seed=1, maxiter=200, verbose=True)
P('relax search 3 restarts; best rel', res[0]['rel'], 'nfev', Vr.nfev)
