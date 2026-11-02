import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
F, A, q, rec = load()
J = json.load(open('families_d7_p0_localtest.json'))
print('keys', list(J['families'])[:5], flush=True)
for key in ['X0+X1', 'X0+Y0+Y1', 'X0+X1+X2+M']:
    fam = J['families'][key]
    print(key, 'admissible_dim', fam['admissible_dim'], 'kernel_part', fam['kernel_part_rank'], 'rational?', fam.get('rational') is not None, fam.get('recon_info'), flush=True)
    if fam.get('rational'):
        for el in fam['rational'][:3]:
            D = core.parse(el['D'])
            ok, info = jacobian.in_jacobian(F, D)
            okE1, infoE1 = fixed_d.e1_solvable(F, A, 1, D)
            print('   coords', el['coords'], '| in J_7:', ok, '| E1:', okE1, '| terms', len(D), flush=True)
