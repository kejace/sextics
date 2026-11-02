import time, sys, json
import flint
from flint import nmod_mat
p = 2147483647
for n in (2000, 4000):
    for th in (1, 8):
        flint.ctx.threads = th
        M = nmod_mat.randtest(n, n + 1, p) if hasattr(nmod_mat, 'randtest') else None
        if M is None:
            import random
            M = nmod_mat(n, n + 1, [random.randrange(p) for _ in range(n * (n + 1))], p)
        t0 = time.time(); Mr, rk = M.rref(); dt = time.time() - t0
        rec = {'bench': 'rref', 'n': n, 'threads': th, 'seconds': round(dt, 2), 'rank': rk}
        print(json.dumps(rec), flush=True)
        t0 = time.time(); X, nul = M.nullspace(); dt = time.time() - t0
        print(json.dumps({'bench': 'nullspace', 'n': n, 'threads': th, 'seconds': round(dt, 2)}), flush=True)
