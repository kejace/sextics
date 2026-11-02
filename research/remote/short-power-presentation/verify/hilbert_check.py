import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
H = hilb()
print('e  S+  S-  R+  R-  J+  J-')
for e in range(0, 25):
    print(e, H[e]['S_plus'], H[e]['S_minus'], H[e]['R_plus'], H[e]['R_minus'], H[e]['J_plus'], H[e]['J_minus'])
# compare with logged ranks
logged = []
for fn in ['/tmp/hodge/results/short-power-presentation/e1kernel.jsonl']:
    for line in open(fn):
        r = json.loads(line)
        d = r['d']
        if 'parD' in r: par = r['parD']
        elif 'parity_D' in r: par = r['parity_D']
        else: par = 0
        rk = r.get('rank J_{d+6}', r.get('rank J_{d+6}^-'))
        tpar = (par + 1) % 2
        exp = H[d + 6]['J_minus' if tpar == 1 else 'J_plus']
        expD = H[d]['J_minus' if par == 1 else 'J_plus']
        rkD = r.get('dim J_d', r.get('dim J_d^+'))
        logged.append((d, par, rk, exp, rk == exp, rkD, expD, rkD == expD))
ok = all(x[4] and x[7] for x in logged)
for x in logged: print('d=%d parD=%d rankJ_{d+6}=%d expected=%d match=%s | dimJ_d=%d expected=%d match=%s' % x)
print('ALL MATCH:', ok)
logjson(claim='rank_p J_{d+6}^{par} and dim J_d^{par} logged in e1kernel.jsonl equal the Hilbert-series values (needed to make mod-p E1 results rigorous over Q)', verdict='CONFIRMED' if ok else 'REFUTED', detail=[dict(d=x[0], parD=x[1], rank=x[2], expected=x[3]) for x in logged])
