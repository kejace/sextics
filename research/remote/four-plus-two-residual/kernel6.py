"""Full fixed-D system (mod p) for each of the 20 d=6 E1-kernel basis elements (mod J) and random combinations."""
from h import *
import json, random
random.seed(11)
p = core.PRIMES[0]
basis = json.load(open(f'{RES}/e1kernel-basis-d6-par0.json'))
basis = [{tuple(int(x) for x in k.strip('()').split(',')): v % p for k, v in b.items()} for b in basis]
tests = []
for r in range(3):
    D = {}
    for b in basis:
        D = core.add(D, core.scale(b, random.randrange(1, p), p), p)
    tests.append((D, 'e1kernel-d6-random', f'random combination #{r} of the 20 d=6 E1-kernel basis elements'))
for i, b in enumerate(basis):
    lead = min(b)  # some monomial
    tests.append((b, 'e1kernel-d6-basis', f'd=6 E1-kernel (mod J) rref basis element #{i}, 17 terms, leading {lead}'))
for D, fam, note in tests:
    r = run_modp(D, p, fam, expr=note, note=note)
    print({k: r[k] for k in ('D','d','consistent','rank','nullity','rows','cols','seconds')}, flush=True)
