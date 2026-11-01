import sys, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
F, A, q, rec = load()
# Fermat denominator D0 = eX on this F: full solve + served space
t0 = time.time()
status, data = run(F, A, q, eX, 'cone-X: eX = X0^6+X1^6+X2^6 (Fermat subset-Euler denominator)', note='D_e with e=0')
print('eX:', status, {k: v for k, v in data.items() if k not in ('N', 'Q')}, f'{time.time()-t0:.1f}s', flush=True)
t0 = time.time()
dim, basis, S, dt = served(F, q, eX, verbose=True)
print(f'served dim for eX: {dim} ({dt:.1f}s)  [expect 18 Jacobian + geometric part]', flush=True)
