"""d=8 tests: (1) the unique non-Jacobian E1-kernel element D8; (2) S3-invariant / sign / random members of the
quadruple-along-curve families (Y)^4∩J_8 and (X)^4∩J_8."""
import sys, json, random, os, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
from hodge.core import add, scale, parse
tests = []
p8 = '/tmp/hodge/results/diagonal-permutation/D8_e1kernel.json'
if os.path.exists(p8):
    D8 = parse(json.load(open(p8))['D'])
    tests.append(('D8 = unique non-Jacobian E1-kernel element (deg 8)', D8, 'E1kernel-d8'))
fam = json.load(open('/tmp/hodge/results/diagonal-permutation/families_sym.json'))
random.seed(7)
def rnd(strs, k):
    D = {}
    for s in strs:
        D = add(D, scale(parse(s), Fraction(random.randint(-k, k))))
    return D
for key, label in [('idx[3, 4, 5]_k4_d8', '(Y)^4∩J_8'), ('idx[0, 1, 2]_k4_d8', '(X)^4∩J_8')]:
    tests.append((f'{label} S3-invariant random combo #1', rnd(fam[key]['inv'], 3), f'J8-{label}-S3inv'))
    tests.append((f'{label} S3-invariant random combo #2', rnd(fam[key]['inv'], 3), f'J8-{label}-S3inv'))
    tests.append((f'{label} S3-sign random combo', rnd(fam[key]['sign'], 3), f'J8-{label}-S3sign'))
    tests.append((f'{label} (01)-anti random combo', rnd(fam[key]['t01anti'], 3), f'J8-{label}-t01anti'))
for name, D, fam_ in tests:
    e = test_D(name, fam_, D=D, verbose=False)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'rank', e['rank'], 'null', e['nullity'], e['seconds'], 's', flush=True)
    if e['e1'] and e.get('consistent') is False:
        t0 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D)
        print(f'   target_image dim = {dim} (baseline 18)  [{time.time()-t0:.0f}s]', flush=True)
        log({'D': name, 'd': e['d'], 'family': fam_, 'e1': True, 'consistent': None, 'rank': -1, 'nullity': -1, 'seconds': round(time.time()-t0,1), 'note': f'target_image dim={dim}'})
