"""JOB excl-dp-fam: structured candidate denominators of degree 8..12 (dossier §4.2 list and products of structured
factors): E1 test (A*D in J, mod p -- a failure is a proof), and for E1-passing candidates the reduced E2' system at
P1 and P2 (consistency; ranks), optionally restricted to the (01)-sector when D is a (01)-eigenvector.
Also: (01)-sector E1-kernel basis elements at d = 9..12 loaded from e1k_e1k.json (if present and small) are tested.
Output: fam_<job>.json; log lines exp in {'fam_e1','fam_e2'}."""
import sys, os, time, json
import numpy as np
import w2lib as W
from hodge import core, fixed_d
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree, parse

job = os.environ.get('W2JOB', 'fam')
OUT = f'{W.PKT}/fam_{job}.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
F, A = W.F, W.A


def save():
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


O1 = []; O2 = []; O3 = []
for (i, j, k) in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
    O1.append(f'X{i}**3*X{j}*Y{i}*Y{j}'); O2.append(f'X{i}**3*X{j}*Y{k}**2'); O3.append(f'Y{i}**3*Y{j}*X{k}**2')
O1 = '(' + '+'.join(O1) + ')'; O2 = '(' + '+'.join(O2) + ')'; O3 = '(' + '+'.join(O3) + ')'
Fc = f'(5*{O1}/7+5*{O2}/11+5*{O3}/13)'
p2X = '(X0**2+X1**2+X2**2)'; p2Y = '(Y0**2+Y1**2+Y2**2)'; pXY = '(X0*Y0+X1*Y1+X2*Y2)'
VX2 = '((X0-X1)*(X1-X2)*(X2-X0))**2'; VY2 = '((Y0-Y1)*(Y1-Y2)*(Y2-Y0))**2'
VXVY = '((X0-X1)*(X1-X2)*(X2-X0)*(Y0-Y1)*(Y1-Y2)*(Y2-Y0))'
p6X = '(X0**6+X1**6+X2**6)'; p6Y = '(Y0**6+Y1**6+Y2**6)'
qp = '(((X0-X1)**2+(Y0-Y1)**2)*((X1-X2)**2+(Y1-Y2)**2)*((X2-X0)**2+(Y2-Y0)**2))'
qm = '(((X0-X1)**2-(Y0-Y1)**2)*((X1-X2)**2-(Y1-Y2)**2)*((X2-X0)**2-(Y2-Y0)**2))'
sumT = '(X0**3*Y0**3+X1**3*Y1**3+X2**3*Y2**3)'; B2 = '(X0**3*Y0**3+X1**3*Y1**3-2*X2**3*Y2**3)'; Aexp = '(X0**3*Y0**3-X1**3*Y1**3)'
sumT2 = '(X0**6*Y0**6+X1**6*Y1**6+X2**6*Y2**6)'

cands = {
    # degree 8
    'VX2*p2X': f'{VX2}*{p2X}', 'VX2*p2Y': f'{VX2}*{p2Y}', 'VX2*pXY': f'{VX2}*{pXY}', 'O1*p2X': f'{O1}*{p2X}', 'O2*p2Y': f'{O2}*{p2Y}',
    'O3*p2X': f'{O3}*{p2X}', 'p6X*p2X': f'{p6X}*{p2X}', 'qp*p2X': f'{qp}*{p2X}', 'p2X^4': f'{p2X}**4', 'Fc*p2X': f'{Fc}*{p2X}',
    'Fc*pXY': f'{Fc}*{pXY}', 'VXVY*pXY': f'{VXVY}*{pXY}', 'sumT*pXY': f'{sumT}*{pXY}',
    # degree 10
    'O1*p2X^2': f'{O1}*{p2X}**2', 'O1*p2Xp2Y': f'{O1}*{p2X}*{p2Y}', 'O1*pXY^2': f'{O1}*{pXY}**2', 'VX2*p2X^2': f'{VX2}*{p2X}**2',
    'VX2*p2Y^2': f'{VX2}*{p2Y}**2', 'VX2*pXY^2': f'{VX2}*{pXY}**2', 'VXVY*p2Xp2Y': f'{VXVY}*{p2X}*{p2Y}', 'VXVY*pXY*p2X': f'{VXVY}*{pXY}*{p2X}',
    'qp*p2X^2': f'{qp}*{p2X}**2', 'qm*p2X^2': f'{qm}*{p2X}**2', 'p6X*p2X^2': f'{p6X}*{p2X}**2', 'p6X*pXY^2': f'{p6X}*{pXY}**2',
    'O2*p2Y^2': f'{O2}*{p2Y}**2', 'O3*p2X^2': f'{O3}*{p2X}**2', 'Fc*p2X^2': f'{Fc}*{p2X}**2', 'Fc*pXY^2': f'{Fc}*{pXY}**2',
    '(p6X+p6Y)*p2X^2': f'({p6X}+{p6Y})*{p2X}**2', 'p2X^5': f'{p2X}**5', 'sumT*p2X^2': f'{sumT}*{p2X}**2', 'B2*p2Xp2Y': f'{B2}*{p2X}*{p2Y}',
    # degree 12
    'VX2*VY2': f'{VX2}*{VY2}', 'O1*O2': f'{O1}*{O2}', 'O1*O3': f'{O1}*{O3}', 'O2*O3': f'{O2}*{O3}', 'O1^2': f'{O1}**2', 'O2^2': f'{O2}**2',
    'O3^2': f'{O3}**2', 'O1*p6X': f'{O1}*{p6X}', 'O1*p6Y': f'{O1}*{p6Y}', 'O2*p6X': f'{O2}*{p6X}', 'O3*p6Y': f'{O3}*{p6Y}',
    'p6X*p6Y': f'{p6X}*{p6Y}', 'p6X^2': f'{p6X}**2', '(p6X-p6Y)^2': f'({p6X}-{p6Y})**2', 'qp*qm': f'{qp}*{qm}', 'qp^2': f'{qp}**2',
    'qm^2': f'{qm}**2', 'VX2*p6X': f'{VX2}*{p6X}', 'VX2*p6Y': f'{VX2}*{p6Y}', 'VY2*p6X': f'{VY2}*{p6X}', 'VXVY*O1': f'{VXVY}*{O1}',
    'VXVY*p6X': f'{VXVY}*{p6X}', 'p2X^3p2Y^3': f'{p2X}**3*{p2Y}**3', 'Fc^2': f'{Fc}**2', 'Fc*p6X': f'{Fc}*{p6X}', 'Fc*p6Y': f'{Fc}*{p6Y}',
    'A^2': f'{Aexp}**2', 'sumT^2': f'{sumT}**2', 'sumT2': sumT2, 'B2^2': f'{B2}**2', 'sumBi^2': f'({B2}**2+(X0**3*Y0**3-2*X1**3*Y1**3+X2**3*Y2**3)**2+(-2*X0**3*Y0**3+X1**3*Y1**3+X2**3*Y2**3)**2)',
    'p2X^2p2Y^2pXY^2': f'{p2X}**2*{p2Y}**2*{pXY}**2', 'VX2VY2-qpqm': f'{VX2}*{VY2}-{qp}*{qm}', 'A*B2': f'{Aexp}*{B2}', 'A*sumT': f'{Aexp}*{sumT}',
}

mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
only_e1 = '--e1-only' in sys.argv

for label, expr in cands.items():
    if label in res and (only_e1 or 'e2' in res[label] or not res[label].get('e1')):
        continue
    D = parse(expr)
    d = degree(D)
    if not core.is_homogeneous(D) or core.parity_of(D) is None:
        res[label] = {'label': label, 'd': d, 'note': 'skipped (mixed parity / inhomogeneous)'}; save(); continue
    t0 = time.time()
    e1, info = fixed_d.e1_solvable(F, A, 1, D, p=W.P1)
    rec = res.get(label, {'label': label, 'd': d, 'expr': expr})
    rec['e1'] = bool(e1); rec['e1_info'] = {k: int(v) for k, v in info.items()}; rec['e1_prime'] = W.P1
    rec['e1_seconds'] = round(time.time() - t0, 1)
    res[label] = rec; save()
    W.log('fam_e1', D=label, d=d, prime=W.P1, result=('E1 passes' if e1 else 'E1 FAILS (excluded over Q)'),
          rank=int(info['rank']), seconds=rec['e1_seconds'], note=f'rank_aug {int(info["rank_aug"])}')
    print(label, 'd', d, 'E1', e1, flush=True)
    if not e1 or only_e1:
        continue
    # reduced E2' at P1 and P2
    rec['e2'] = {}
    for pr in (W.P1, W.P2):
        t1 = time.time()
        try:
            u = W.solve_E1(D, pr)
            S = W.ReducedE2(D, pr, u=u, verbose=True)
            r = S.solve(compressed=True, verbose=True)
        except Exception as ex:
            r = {'error': str(ex)}
        rec['e2'][str(pr)] = r
        save()
        W.log('fam_e2', D=label, d=d, prime=pr, result=('consistent' if r.get('consistent') else 'inconsistent' if 'consistent' in r else 'error'),
              rank=r.get('rank', -1), nullity=r.get('nullity', -1), seconds=round(time.time() - t1, 1), note=f'reduced E2 system {r.get("rows")}x{r.get("cols")}')
        print(label, 'E2 at', pr, r, flush=True)
        if r.get('consistent'):
            print('!!! CONSISTENT mod p -- run the exact lift with the engine', flush=True)
            W.log('fam_e2_ALERT', D=label, d=d, prime=pr, result='CONSISTENT mod p', note='needs exact lift')
print('DONE', flush=True)
