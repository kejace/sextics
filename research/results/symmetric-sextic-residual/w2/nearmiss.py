"""Examine the best minimum of a sector: re-minimise with tighter lstsq cutoff, print the dominant a-coordinates
(in the orbit-sum basis), rationalise and solve exactly mod p."""
import sys, json, glob
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
F, A, q = load()
chi, d, pD, job = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
recs = []
for fn in glob.glob(f'{RES}/w2/{job}_{chi}-{d}-{pD}_w*.jsonl'):
    for line in open(fn):
        r = json.loads(line)
        if 'a' in r: recs.append(r)
recs.sort(key=lambda r: r['rel'])
best = recs[0]; a = np.array(best['a'])
fam, info = build_sector_family(F, A, q, chi, d, pD, label=f'{chi}/d{d}/p{pD}')
print('info', {k: info[k] for k in ('nD', 'dimK', 'dimJ', 'ker_beyond_J', 'dimF', 'dimP')})
print('best rel from log', best['rel'])
for cond in (1e-9, 1e-12, 1e-14):
    fam.cond = cond
    rel, f, x = fam.rel(a); print(f'cond={cond}: rel={rel:.3e} f={f:.3e}')
fam.cond = 1e-12
c = fam.search(n_restarts=1, seed=0, maxiter=2000, inits=(a,), verbose=True)[0]
print('re-minimised rel', c['rel'], 'iters', c['iters'])
a = c['a']; j = int(np.argmax(np.abs(a))); an = a / a[j]
order = np.argsort(-np.abs(an))
print('dominant coordinates (rep: value):')
for k in order[:12]:
    print('   ', fam.reps[k], f'{an[k]:+.6f}')
for fr in fam.reconstruct(a, max_dens=(12, 100, 1000)):
    D = fam.D_exact(fr)
    r = solve_single(F, A, q, D, fam.G, P0)
    print('rationalised den<=', max(x.denominator for x in fr), '->', {k: r[k] for k in ('e1', 'consistent', 'rank', 'nullity', 'rows', 'cols')})
