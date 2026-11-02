import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
F, A, q, rec = load()
J = json.load(open('families_d7_p1_w2.json'))
fam = J['families']['Y0+Y1+Y2+L']
print('Y0+Y1+Y2+L: admissible', fam['admissible_dim'], 'kernel part', fam['kernel_part_rank'], 'rational:', fam.get('rational') is not None, fam.get('recon_info'))
if fam.get('rational'):
    # find the element(s) with nonzero kernel part: test E1 and J-membership for each rational element
    found = []
    for el in fam['rational']:
        D = core.parse(el['D'])
        inJ, _ = jacobian.in_jacobian(F, D)
        if inJ: continue
        okE1, info = fixed_d.e1_solvable(F, A, 1, D)
        print('  element not in J: E1 =', okE1, 'terms', len(D), 'coords', len(el['coords']))
        if okE1: found.append(D)
    if found:
        D = found[0]
        # structure: L-degrees present, Y-degrees
        from collections import Counter
        print('  Y-degrees:', sorted(Counter(sum(e[3:]) for e in D).items()))
        dX = [core.pderiv(D, i) for i in range(3)]
        print('  dX1 == 2 dX0:', core.sub(dX[1], core.scale(dX[0], 2)) == {}, ' dX2 == 3 dX0:', core.sub(dX[2], core.scale(dX[0], 3)) == {})
        out = {'d': 7, 'par': 1, 'primes': [], 'info': {}, 'basis_modp': {}, 'basis_rational': [core.to_str(D)],
               'provisional': 'from family screen Y0+Y1+Y2+L (2-prime reconstruction in family coordinates), E1 verified mod p'}
        json.dump(out, open('kernel_d7_p1_fam.json', 'w'))
        print('  saved kernel_d7_p1_fam.json; coefficient heights: max num bits', max(abs(c.numerator).bit_length() for c in D.values()), 'max den bits', max(c.denominator.bit_length() for c in D.values()))
