"""Full served-target space fixed_d.target_image for random Jacobian denominators w(F) (d = 6, 7): is served(w(F)) = J_6^- always?"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from w2x_hess import field_basis
d = int(sys.argv[1]); n = int(sys.argv[2]); tag = sys.argv[3]; p = P_NEW
random.seed(d * 100 + n)
Fp, Fi = partials(p)
gens = J_generators(Fi, 6, 1, p)
for r in range(n):
    for par in (0, 1):
        B = field_basis(d, par)
        D = {}
        for (i, m) in B:
            c = random.randint(-3, 3)
            if c: D = add(D, scale(shift(Fi[i], m, 1, p), c % p, p), p)
        t0 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D, p=p)
        mons = {}
        for P_ in basis + gens:
            for ee in P_: mons.setdefault(ee, len(mons))
        def rk(L):
            flat = [0] * (len(L) * len(mons))
            for rr, P_ in enumerate(L):
                for ee, v in P_.items(): flat[rr * len(mons) + mons[ee]] = v % p
            return nmod_mat(len(L), len(mons), flat, p).rank() if L else 0
        rJ = rk(gens); rBJ = rk(basis + gens)
        e = {'exp': 'c-served-random-J', 'D': f'random w(F) d={d} par={par} r={r}', 'family': 'random Jacobian denominator', 'd': d, 'prime': p,
             'result': f'served dim {dim}, modulo J_6^-: {rBJ - rJ}', 'served_dim': int(dim), 'served_mod_J': int(rBJ - rJ), 'rank': None,
             'seconds': round(time.time() - t0, 1), 'job': tag, 'note': 'fixed_d.target_image (all odd degree-6 targets served by this D)'}
        print(json.dumps(e), flush=True)
        wlog(e, fn=f'{PK}/wave2-excl-remote-{tag}.jsonl' if tag != 'local' else None)
