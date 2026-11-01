"""Extra family sets (resume after crash). Usage: python3 families2.py <setname>"""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
import families as fam
from families import *
sets = {}
sets['deg6b'] = [('E_X(F)=6a+4b+2c', EX(F)), ('E_Y(F)=2b+4c+6d', EY(F)),
                 ('a - d', sub(B[0], B[6])), ('X0^6+Y0^6 pair', core.parse('X0**6+Y0**6')),
                 ('(X0^2+X1^2+X2^2)^3', power(core.parse('X0**2+X1**2+X2**2'), 3)),
                 ('(X0^2+Y0^2)^3', power(core.parse('X0**2+Y0**2'), 3)),
                 ('E_X(F) - E_Y(F)', sub(EX(F), EY(F))),
                 ('X0*F_X0 (block radial, 1 var)', mul(core.parse('X0'), Fi[0])),
                 ('X0*F_X0 + Y0*F_Y0 (pair radial)', add(mul(core.parse('X0'), Fi[0]), mul(core.parse('Y0'), Fi[3])))]
sets['deg7x'] = [('X0^2*F_X0', mul(core.parse('X0**2'), Fi[0])),
                 ('X0*F_X0 + Y0*F_Y0 ... times X0? no: X0*E_X(F)', mul(core.parse('X0'), EX(F)))]
sets['deg8x'] = [('X0^3*F_X0', mul(core.parse('X0**3'), Fi[0])),
                 ('X0^3 F_X0 + Y0^3 F_Y0', add(mul(core.parse('X0**3'), Fi[0]), mul(core.parse('Y0**3'), Fi[3]))),
                 ('Hess minor F_X0X0*F_Y0Y0 - F_X0Y0^2', fam.sets['deg8'][0][1]),
                 ('Hess minor F_X0X0*F_X1X1 - F_X0X1^2', fam.sets['deg8'][1][1])]
if __name__ == '__main__':
    name = sys.argv[1]
    for nm, D in sets[name]:
        if core.parity_of(D) != 0:
            print('skip (not Y-even):', nm); continue
        t0 = time.time(); d = core.degree(D)
        if d <= 9 and not in_J(D):
            log(nm, d, name, False, False, None, None, time.time()-t0, note='D not in J_d; E1 fails since ker([A] on R_d^+)=0 for d<=9 (E1-kernel sweep)')
            print(f"{nm}: d={d} NOT in J_d -> E1 fails (kernel mod J = 0) ({time.time()-t0:.1f}s)", flush=True); continue
        r = test_D(F, A, q, D, family=name, name=nm, full=True)
        print(f"{nm}: d={r['d']} e1={r['e1']} consistent={r.get('consistent')} rank={r.get('rank')} nullity={r.get('nullity')} ({time.time()-t0:.1f}s)", flush=True)
