"""Comparison on the SOLVED neighbour one-reflection (F0 = f0 + Y2^6): activity of the same low-degree families."""
import sys, os, pickle, itertools
os.environ['HODGE_INST'] = 'one-reflection'
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import five_generic as fg
from five_generic import *
from fractions import Fraction as Fr
P = core.parse
f0 = fg.f
inv = {'G0': P('X0**6+X1**6+X2**6'), 'c': P('X0**3*X1+X1**3*X2'), 'e': P('X0*X1'), 'p': P('Y0*Y1'), 't': P('Y0**2-Y1**2'),
       'm': P('Y0**2+Y1**2'), 's': P('Y0**3-Y1**3'), 'S': P('Y0**3+Y1**3'), 'Y6': P('Y0**6+Y1**6'), 'f': f0,
       'X0': X[0], 'X1': X[1], 'X2': X[2], 'Y0': X[3], 'Y1': X[4], 'cp': P('(X0**3*X1+X1**3*X2)*Y0*Y1')}
G0, c, p, s = inv['G0'], inv['c'], inv['p'], inv['s']
fc = add(add(G0, mul(c, p)), scale(mul(p, mul(p, p)), 2))          # companion's f = G0 + c p + 2 p^3  (f0 = fc + s^2)
assert sub(f0, add(fc, mul(s, s))) == {}, 'f0 = fc + s^2 check'
delta = add(scale(mul(G0, G0), 324), scale(mul(c, mul(c, c)), 24))  # deg 12
h = add(add(scale(mul(c, c), -1), scale(mul(G0, p), 18)), scale(mul(c, mul(p, p)), 6))  # deg 8
cands = []
for lam in [Fr(1,3), Fr(-1,3), Fr(1), Fr(-1), Fr(3), Fr(-3), Fr(2), Fr(1,2)]:
    cands.append((f'f0+({lam})*s^2', add(f0, scale(mul(s, s), lam)), 'one-reflection:cone-norm'))
    cands.append((f'f0+({lam})*S^2', add(f0, scale(mul(inv['S'], inv['S']), lam)), 'one-reflection:cone-norm'))
cands.append(('h', h, 'one-reflection:companion-piece'))
cands.append(('fc', fc, 'one-reflection:companion-piece'))
cands.append(('fc*X0', mul(fc, X[0]), 'one-reflection:companion-piece'))
cands.append(('h*X0', mul(h, X[0]), 'one-reflection:companion-piece'))
cands.append(('h*p', mul(h, p), 'one-reflection:companion-piece'))
cands.append(('h*e', mul(h, inv['e']), 'one-reflection:companion-piece'))
cands.append(('s^2-9fc', sub(mul(s, s), scale(fc, 9)), 'one-reflection:companion-piece'))
cands.append(('s^2*(s^2-9fc)', mul(mul(s, s), sub(mul(s, s), scale(fc, 9))), 'one-reflection:companion-piece'))
cands.append(('delta', delta, 'one-reflection:companion-piece'))
names = ['G0', 'c', 'e', 'p', 't', 'm', 's', 'S', 'Y6', 'f', 'X0', 'X1', 'X2', 'Y0', 'Y1', 'cp']
seen = set()
for k in (2, 3):
    for combo in itertools.combinations_with_replacement(names, k):
        Pp = {(0,)*6: Fr(1)}
        for nm in combo: Pp = mul(Pp, inv[nm])
        d = core.degree(Pp)
        if d < 7 or d > 8 or rparity_of(Pp) != 0: continue
        key = tuple(sorted(Pp.items()))
        if key in seen: continue
        seen.add(key); cands.append(('*'.join(combo), Pp, f'one-reflection:invariant-product:k{k}'))
for i in range(5):
    for j in range(i + 1, 5):
        H = [[pderiv(fi(a), b) for b in range(5)] for a in range(5)]
        cands.append((f'Hess0[{i},{j}]', sub(mul(H[i][i], H[j][j]), mul(H[i][j], H[i][j])), 'one-reflection:hessian-2x2'))
cands = [(n, D, fam) for n, D, fam in cands if D and rparity_of(D) == 0]
cands.sort(key=lambda x: core.degree(x[1]))
os.makedirs('served_onerefl', exist_ok=True)
print(len(cands), 'candidates on one-reflection', flush=True)
for name, D, fam in cands:
    fn = name.replace('/', '_').replace('*', 'x') + '.pkl'
    if os.path.exists('served_onerefl/' + fn): continue
    try:
        r = test5(name, fam, note='comparison on solved neighbour one-reflection', Dpoly=D, do_exact=False)
    except Exception as ex:
        print('ERR', name, ex, flush=True); continue
    with open('served_onerefl/' + fn, 'wb') as fh:
        pickle.dump({'name': name, 'family': fam, 'd': r['d'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                     'serves_A': r['serves_A'], 'served_basis': r['served_basis'], 'p': r['p']}, fh)
print('DONE')
