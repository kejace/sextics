import sys, pickle, time, os, json
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from candidates import cands
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
fams = sys.argv[3].split(',') if len(sys.argv) > 3 else None
C = [(n, Pp, fam) for n, Pp, fam in cands() if dmin <= core.degree(Pp) <= dmax and (fams is None or any(fam.startswith(x) for x in fams))]
C.sort(key=lambda x: core.degree(x[1]))
os.makedirs('served', exist_ok=True)
done = set(os.listdir('served'))
print(f'{len(C)} candidates, degrees {dmin}..{dmax}', flush=True)
for name, Pp, fam in C:
    fn = name.replace('/', '_').replace('*', 'x').replace('|', 'I') + '.pkl'
    if fn in done: 
        print('skip', name); continue
    try:
        r = test5(name, fam, note=f'5var candidate {fam}', Dpoly=Pp)
    except Exception as ex:
        print('ERROR', name, ex, flush=True); continue
    with open('served/' + fn, 'wb') as fh:
        pickle.dump({'name': name, 'family': fam, 'D': {e: str(c) for e, c in Pp.items()}, 'p': r['p'], 'd': r['d'],
                     'served_basis': r['served_basis'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                     'serves_A': r['serves_A'], 'serves_A+': r['serves_A+'], 'serves_A-': r['serves_A-'],
                     'rank_N': r['rank_N'], 'nullity_N': r['nullity_N'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds']}, fh)
print('DONE', flush=True)
