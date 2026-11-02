"""Job: excl-extras<d>.  Screen: full sector-restricted fixed-D solve (E1+E2 exact mod p) for every basis element of the
non-Jacobian part of the E1 kernel in every G_aff sector at degree d (basis = rref-reduced kernel vectors modulo J_d[c]).
A basis-element screen is not a family exclusion (combinations are untested); it is logged as 'screen'.
Usage: extras_screen.py d [workers]"""
import sys, time, json, os
from multiprocessing import Pool
from w2lib import *
from family import e1_profile_sector
d = int(sys.argv[1]); workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
OUT = PK + f'/extras_d{d}'
os.makedirs(OUT, exist_ok=True)
chars = sorted({chi_of({e: 1}) for e in monomials(d)})

def one(args):
    c, i, P = args
    t0 = time.time()
    S = SectorSystem(P, verbose=False)
    res = S.solve()
    ok = None
    if res['consistent']:
        N, Q = S.unpack(res['x'])
        E = core.defect(Fp, Ap, P, 1, N, p)
        ok = (sub(E, mul(power(Fp, 2, p), Q, p), p) == {})
        json.dump({'D': core.to_str(P), 'N': [core.to_str(x) for x in N], 'Q': core.to_str(Q), 'p': p}, open(OUT + f'/CONSISTENT_{"".join(map(str,c))}_{i}.json', 'w'))
    return {'chi': str(c), 'i': i, 'consistent': res['consistent'], 'verified_mod_p': ok, 'rank': res['rank'], 'nullity': res['nullity'], 'rows': res['rows'], 'cols': res['cols'], 'seconds': round(time.time() - t0, 1), 'nterms': len(P)}

tasks = []
prof_all = {}
for c in chars:
    prof, extras = e1_profile_sector(d, c)
    prof_all[str(c)] = prof
    for i, P in enumerate(extras):
        tasks.append((c, i, P))
    print(c, prof, file=sys.stderr, flush=True)
log('extras-profile', d=d, n_extras=len(tasks), per_sector={k: v['extra'] for k, v in prof_all.items()})
json.dump({'profile': prof_all, 'extras': {str(c): core.to_str(P) for (c, i, P) in tasks}}, open(OUT + '/extras.json', 'w'))
n_cons = 0
with Pool(workers) as pool:
    for r in pool.imap_unordered(one, tasks):
        if r['consistent']: n_cons += 1
        log('extras-screen', d=d, **r, result='CONSISTENT' if r['consistent'] else 'inconsistent')
log('extras-summary', d=d, n_tested=len(tasks), n_consistent=n_cons, note='sector-restricted full fixed-D solves of basis elements of the non-Jacobian E1 kernel (screen, not a family exclusion)')
print('DONE', file=sys.stderr)
