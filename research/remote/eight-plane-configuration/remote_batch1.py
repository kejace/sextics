import sys, time, json, random
sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration')
from lib import *
import lib
_log = lib.log
def log2(e): e['engine'] = 2; e['host'] = 'remote'; _log(e)
lib.log = log2
from hodge.core import monomials
Fi = [core.pderiv(F, i) for i in range(6)]
which = sys.argv[1]
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
if which == 'target':
    for name, D in [('dF/dX0', Fi[0]), ('(X0 dX0 + Y0 dY0)F', core.add(core.mul(X0, Fi[0]), core.mul(Y0, Fi[3])))]:
        t = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D, verbose=True)
        print(json.dumps({'D': name, 'd': core.degree(D), 'target_image_dim': int(dim), 'seconds': round(time.time()-t,1)}), flush=True)
        with open(f'{WD}/target_image.jsonl','a') as f: f.write(json.dumps({'D': name, 'd': core.degree(D), 'target_image_dim': int(dim), 'seconds': round(time.time()-t,1), 'host':'remote'})+'\n')
elif which == 'randomJ':
    for d, reps in [(5, 2), (6, 2), (7, 2), (8, 1)]:
        for r in range(reps):
            for par in (0, 1):
                D = {}
                for i in range(6):
                    parM = (par - (0 if i < 3 else 1)) % 2
                    for m in monomials(d - 5, parM):
                        c = random.randint(-3, 3)
                        if c: D = core.add(D, core.scale(core.shift(Fi[i], m), c))
                if not D: continue
                name = f'randomJ_d{d}_par{par}_r{r}'
                e = test_D(D, name, 'random-J', note=f'random element of J(F)_{d} (coeffs in [-3,3]), parity {par}', fast=True, want_exact=True)
                print(name, {k: e.get(k) for k in ('d','e1','consistent','rank','nullity','rows','cols','seconds')}, flush=True)
