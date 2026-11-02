import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_e6 import residual_space
from collections import Counter
F, A, q, rec = load()
out = {}
for lname, d, par in [('X0', 8, 0), ('Y0', 8, 0), ('X1', 8, 0), ('Y1', 8, 0)]:
    lists = {}
    rat = None
    for p in PGB + core.PRIMES[:3]:
        reps, jref = residual_space(lname, d, par, p); lists[p] = reps
        if len(lists) >= 2:
            rat, info = reconstruct_polys(lists)
            if rat: break
    print(f'Ann({lname})_{d}^+ : dim {len(lists[p])}, reconstruction {"ok" if rat else "failed"} with {len(lists)} primes', flush=True)
    if rat:
        E = rat[0]
        xd = sorted(Counter(sum(e[:3]) for e in E).items()); yd = sorted(Counter(sum(e[3:]) for e in E).items())
        print(f'  NF rep: terms {len(E)}, X-degrees {xd}, Y-degrees {yd}, min exp of {lname}:', min(e[0] if lname == "X0" else (e[3] if lname == "Y0" else (e[1] if lname == "X1" else e[4])) for e in E), flush=True)
        dY = [core.pderiv(E, i) for i in range(3, 6)]; dX = [core.pderiv(E, i) for i in range(3)]
        print('  Y only via M:', core.sub(dY[1], core.scale(dY[0], 3)) == {} and core.sub(dY[2], core.scale(dY[0], 5)) == {},
              ' X only via L:', core.sub(dX[1], core.scale(dX[0], 2)) == {} and core.sub(dX[2], core.scale(dX[0], 3)) == {}, flush=True)
        out[lname] = core.to_str(E)
json.dump(out, open('Ann8_reps.json', 'w')); print('saved Ann8_reps.json', flush=True)
