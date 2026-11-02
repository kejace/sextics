import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
F, A, q = load()
t0 = time.time()
fam, info = build_sector_family(F, A, q, 'triv', 6, 0, verbose=True)
print('real F triv/6/0 info:', info, f'{time.time()-t0:.1f}s')
if fam is not None:
    cands = fam.search(n_restarts=12, seed=1, verbose=False)
    print('real F best rel:', [f"{c['rel']:.2e}" for c in cands[:6]], 'nfev', fam._nfev)
# Fermat sanity: must find p6X
F0 = core.add(p_k(X, 6), p_k(Y, 6))
t0 = time.time()
fam0, info0 = build_sector_family(F0, A, q, 'triv', 6, 0, verbose=True)
print('Fermat triv/6/0 info:', info0, f'{time.time()-t0:.1f}s')
cands = fam0.search(n_restarts=8, seed=3, verbose=True)
best = cands[0]
print('Fermat best rel', best['rel'], 'spurious', best['spurious'])
a = best['a']; j = int(np.argmax(np.abs(a))); print('normalized a:', np.round(a / a[j], 6))
print('reps:', fam0.reps)
ok, res = try_confirm(fam0, best, {'family': 'test-fermat'}, verbose=False)
print('confirm:', ok, res)
