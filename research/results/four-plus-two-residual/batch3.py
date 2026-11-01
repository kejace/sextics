"""Batch 3: degree-5 Jacobian elements (the only degree-5 E1 survivors), then Y-odd d=6 E1-kernel."""
from h import *
import json, random, subprocess
random.seed(3)
Fs = core.to_str
tests = []
for i, v in enumerate(core.VARS):
    tests.append((Fs(core.pderiv(F, i)), 'partials', f'dF/d{v} (in J_5; E1 trivial)'))
for (i,j,name) in [(0,3,'X0,Y0'),(1,4,'X1,Y1'),(2,5,'X2,Y2'),(0,1,'X0,X1'),(3,4,'Y0,Y1')]:
    R = core.add(core.mul(P(core.VARS[i]), core.pderiv(F,i)), core.mul(P(core.VARS[j]), core.pderiv(F,j)))
    tests.append((Fs(R), 'radial', f'({name}) radial derivative of F'))
cX = [random.choice([1,2,3,5,7,-1,-2,-3]) for _ in range(3)]
cY = [random.choice([1,2,3,5,7,-1,-2,-3]) for _ in range(3)]
DX = {}
for i in range(3): DX = core.add(DX, core.scale(core.pderiv(F, i), cX[i]))
DY = {}
for j in range(3): DY = core.add(DY, core.scale(core.pderiv(F, 3 + j), cY[j]))
tests.append((Fs(DX), 'J5-combo', f'generic Y-even element of J_5: {cX}.(F_X0,F_X1,F_X2)'))
tests.append((Fs(DY), 'J5-combo', f'generic Y-odd element of J_5: {cY}.(F_Y0,F_Y1,F_Y2)'))
seen = set()
for expr, fam, note in tests:
    if expr in seen: continue
    seen.add(expr)
    r = run(expr, fam, note=note, fast=False)
    print({k: r[k] for k in ('D','d','e1','consistent','rank','nullity','seconds','note')}, flush=True)
# Y-odd d=6 kernel
subprocess.run(['python3', f'{RES}/e1kernel.py', '6', '6', '1'], check=True)
p = core.PRIMES[0]
basis = json.load(open(f'{RES}/e1kernel-basis-d6-par1.json'))
basis = [{tuple(int(x) for x in k.strip('()').split(',')): v % p for k, v in b.items()} for b in basis]
tests = []
for r in range(2):
    D = {}
    for b in basis:
        D = core.add(D, core.scale(b, random.randrange(1, p), p), p)
    tests.append((D, 'e1kernel-d6odd-random', f'random combination #{r} of the {len(basis)} Y-odd d=6 E1-kernel basis elements'))
for i, b in enumerate(basis):
    tests.append((b, 'e1kernel-d6odd-basis', f'Y-odd d=6 E1-kernel (mod J) rref basis element #{i}, {len(b)} terms'))
for D, fam, note in tests:
    r = run_modp(D, p, fam, expr=note, note=note)
    print({k: r[k] for k in ('D','d','consistent','rank','nullity','rows','cols','seconds')}, flush=True)
