import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from collections import Counter
F, A, q, rec = load()
J = json.load(open('families_d7_p1_w2.json'))
fam = J['families']['Y0+Y1+Y2+L']
print('Y0+Y1+Y2+L: admissible', fam['admissible_dim'], 'kernel part', fam['kernel_part_rank'], 'rational:', fam.get('rational') is not None, flush=True)
els = [el for el in (fam.get('rational') or []) if el.get('nonjacobian')]
print('non-Jacobian reps:', len(els), flush=True)
for el in els:
    D = core.parse(el['D'])
    okE1, info = fixed_d.e1_solvable(F, A, 1, D)
    print('  E1:', okE1, 'terms', len(D), 'family coords', len(el['coords']), flush=True)
    print('  coords (exponents of Y0,Y1,Y2,L -> coeff):', el['coords'], flush=True)
    dX = [core.pderiv(D, i) for i in range(3)]
    print('  depends on X only via L:', core.sub(dX[1], core.scale(dX[0], 2)) == {} and core.sub(dX[2], core.scale(dX[0], 3)) == {}, flush=True)
    json.dump({'d': 7, 'par': 1, 'D_cone': core.to_str(D), 'coords': el['coords']}, open('D7odd_cone.json', 'w'))
