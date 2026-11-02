import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq
from hodge import core
F, A, q, rec = core.load_instance('cycle-q2')
D = core.parse('X0**6 + X1**6 + X2**6')
D2 = core.mul(D, D)
t0 = time.time()
r = eq.test_D_fast(F, A, q, D2, 'X-block', note='validation of fast path vs flint rank 2847', verbose=True, Dstr='(X0**6 + X1**6 + X2**6)**2 [fastval]')
print({k: v for k, v in r.items() if not k.startswith('_')}, time.time() - t0, flush=True)
