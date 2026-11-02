"""Structural family tests: D_alpha (d=7, triple along Z_Y = X ∩ {X=0}) and members of the (X)^3 ∩ J_7 family."""
import sys, json, random; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
from hodge.core import pderiv, mul, add, sub, scale, to_str, parse
Fi = [pderiv(F, i) for i in range(6)]
X = [parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
def D_alpha(al):
    D = {}
    for j in range(3):
        a, b = [t for t in range(3) if t != j]
        # G_{Xj}^{(0)} = alpha_j (Y_a^2 - Y_b^2);  G_{Y_a} += -(5/39) alpha_j X_j Y_b ; G_{Y_b} += (5/39) alpha_j X_j Y_a
        if al[j] == 0: continue
        D = add(D, scale(mul(sub(mul(X[3+a],X[3+a]), mul(X[3+b],X[3+b])), Fi[j]), Fraction(al[j])))
        D = add(D, scale(mul(mul(X[j], X[3+b]), Fi[3+a]), Fraction(-5, 39) * al[j]))
        D = add(D, scale(mul(mul(X[j], X[3+a]), Fi[3+b]), Fraction(5, 39) * al[j]))
    return D
def xdeg_min(D): return min(e[0]+e[1]+e[2] for e in D)
tests = []
for al, nm in [((1,0,0),'alpha=(1,0,0)'), ((1,1,1),'alpha=(1,1,1)'), ((1,-1,0),'alpha=(1,-1,0)'), ((1,2,3),'alpha=(1,2,3)')]:
    D = D_alpha(al)
    assert xdeg_min(D) >= 3, (nm, xdeg_min(D))
    tests.append((f'D_alpha {nm}', D, 'J7-(X)^3-alpha'))
fam = json.load(open('/tmp/hodge/results/diagonal-permutation/families.json'))
reps = [parse(s) for s in fam['idx[0, 1, 2]_k3_d7']]
random.seed(1)
for t in range(2):
    D = {}
    for P in reps:
        D = add(D, scale(P, Fraction(random.randint(-5, 5))))
    tests.append((f'random member #{t} of (X)^3∩J_7 (21-dim)', D, 'J7-(X)^3-random'))
# alpha family plus random Q-part
D = add(D_alpha((1,1,1)), add(mul(mul(X[0],X[1]), Fi[2]), mul(mul(X[2],X[2]), Fi[0])))
tests.append(('D_alpha(1,1,1) + X0X1 F_X2 + X2^2 F_X0', D, 'J7-(X)^3-mixed'))
for name, D, fam_ in tests:
    e = test_D(name, fam_, D=D, verbose=False)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'rank', e['rank'], 'null', e['nullity'], e['seconds'], 's', flush=True)
    if e['e1'] and e.get('consistent') is False:
        # served target space dimension
        import time; t0 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D)
        print(f'   target_image dim = {dim} (J_6^- = 18 expected baseline)  [{time.time()-t0:.0f}s]', flush=True)
        log({'D': name, 'd': e['d'], 'family': fam_, 'e1': True, 'consistent': None, 'rank': -1, 'nullity': -1, 'seconds': round(time.time()-t0,1), 'note': f'target_image dim={dim}'})
