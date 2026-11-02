"""E3: E1 annihilator dimensions Ann_d = {[D] in R_d : [A D] = 0}, sector-wise, d = 0..12, at two primes.
Zero kernels are rigorous over Q.  Also records per-sector kernels (dim incl. J-part) for later use."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
for p in (P2, P1):
    JS = JacobianSectors(F2, p)
    for d in range(0, dmax + 1):
        t0 = time.time()
        tot = 0; totk = 0; totj = 0; per = {}
        for key in all_sector_keys(JS.gens, d):
            r = JS.e1_kernel(A, d, key)
            tot += r['dim']; totk += r['kernel_dim']; totj += r['jdim']
            if r['kernel_dim']:
                par = (sum(m[3:]) % 2 for m in r['basis'][0].keys()) if r['basis'] else None
                per[','.join(map(str, key))] = dict(dim=r['dim'], kernel_dim=r['kernel_dim'], jdim=r['jdim'], ring_dim=r['ring_dim'], nmonos=r['nmonos'])
        rec = dict(exp='E3-E1-annihilator', engine=ENGINE_TAG, p=p, d=d, Ann_dim=tot, E1_space_dim=totk, J_dim=totj, R_dim=len(monos(d)) - totj,
                   sectors_nonzero=len(per), per_sector=per if d <= 8 else {k: v for k, v in list(per.items())[:0]}, seconds=round(time.time() - t0, 1))
        print(p, d, 'Ann', tot, 'E1space', totk, 'J', totj, 'R', rec['R_dim'], f'{rec["seconds"]}s', flush=True)
        xlib.log(rec)
        if p == P1 and d >= 8:
            break
