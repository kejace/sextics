"""Served-target space (target_image) of Jacobian-type denominators D = w(F): is it only J_6^- (dim 18)?"""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
d = int(sys.argv[1]); reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
random.seed(d)
for par in (0, 1):
    for r in range(reps):
        D = {}
        for i in range(NV):
            pw = (par - (0 if i < 3 else 1)) % 2
            for m in monomials(d - 5, pw):
                c = random.randint(-2, 2)
                if c: D = add(D, scale(shift(Fi[i], m), c))
        t0 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D, p=PRIME)
        # how much of the image lies in J_6^-?
        from hodge import jacobian
        inJ = sum(1 for b in basis if jacobian.in_jacobian(F, b, p=PRIME)[0])
        rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'target_image of random Jacobian-type D (mod p)', 'd': d, 'par': par, 'rep': r,
               'target_image_dim': int(dim), 'basis_elements_in_J6': inJ, 'dim_J6_odd': 18, 'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 1)}
        print(json.dumps(rec), flush=True)
        jlog(f'{WD}/wave2-search.log.jsonl', rec)
