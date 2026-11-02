import sys, time
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import numpy as np
import w2s
from w2s import core
F0 = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
w2s.configure(F0, w2s.A)
famF = w2s.jacobian_family(6)
pvF = w2s.PointVarPro(famF, S=250, seed=7, verbose=True)
t0 = time.time()
res = pvF.search_lm(n_restarts=40, seed=int(sys.argv[1]) if len(sys.argv) > 1 else 3, max_nfev=300, verbose=False, stop_rel=-1)
print('LM Fermat control: %d restarts in %.1fs' % (len(res), time.time() - t0))
print(' rels:', ['%.1e' % r['rel'] for r in res])
print(' nfev:', [r['nfev'] for r in res][:20])
good = [r for r in res if r['rel'] < 1e-12]
print(' zeros found:', len(good))
for r in good[:5]:
    a = r['a'] / np.max(np.abs(r['a']))
    print('  a (rounded):', [(famF.labels[k], round(float(a[k]), 4)) for k in range(famF.K) if abs(a[k]) > 1e-6])
