"""Batch 2: cubic cones in all 3-variable subsets (vertex a plane), pure parity, plus a few random cubics."""
from h import *
import itertools, random
random.seed(2)
V = core.VARS
def mon(e): return core.parse('*'.join(f'{V[i]}**{k}' for i,k in enumerate(e) if k) or '1')
def is_smooth_plane_curve(D, idx):
    # quick check mod p: random points not needed; skip (we just run)
    return True
tests = []
for S in itertools.combinations(range(6), 3):
    names = [V[i] for i in S]
    ny = sum(1 for i in S if i >= 3)
    for par in (0, 1):
        # monomials of degree 3 in the 3 variables with Y-parity par
        mons = []
        for e in core.monomials(3, par):
            if all(e[i] == 0 for i in range(6) if i not in S):
                mons.append(e)
        if not mons: continue
        # canonical: sum of cubes if pure parity allows
        cubes = [e for e in mons if max(e) == 3]
        if len(cubes) == 3:
            tests.append(('+'.join(core.to_str(mon(e)) for e in cubes), 'cone3-cubic', f'Fermat cubic in {names} par{par}'))
        elif len(cubes) == 2:
            # add a square-times-linear term with the third variable
            others = [e for e in mons if max(e) == 2 and e not in cubes]
            for o in others[:2]:
                tests.append(('+'.join(core.to_str(mon(e)) for e in cubes + [o]), 'cone3-cubic', f'cubic cone in {names} par{par}'))
        elif len(cubes) == 1:
            others = [e for e in mons if e not in cubes]
            for o in others[:2]:
                tests.append(('+'.join(core.to_str(mon(e)) for e in cubes + [o]), 'cone3-cubic', f'cubic cone in {names} par{par}'))
        # random cubics (2 per subset/parity)
        for t in range(2):
            coeffs = [random.choice([-2,-1,1,2,3]) for _ in mons]
            expr = ' + '.join(f'{c}*{core.to_str(mon(e))}' for c, e in zip(coeffs, mons))
            tests.append((expr, 'cone3-cubic-random', f'random cubic in {names} par{par}'))
print(len(tests), 'tests')
seen=set()
for expr, fam, note in tests:
    if expr in seen: continue
    seen.add(expr)
    r = run(expr, fam, note=note)
    print({k: r[k] for k in ('D','d','e1','consistent','rank','nullity','seconds','note')}, flush=True)
