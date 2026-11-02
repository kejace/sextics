"""(a) Adversarial re-verification of the fixed-D exclusions of dossier A.5 (and new siblings) with the
independent solver, new primes, and exact Farkas certificates saved to farkas/."""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
Fi = [pderiv(F2, i) for i in range(6)]          # derivatives of 2F (scaling irrelevant)
X = [{tuple(int(j == i) for j in range(6)): 1} for i in range(6)]
cands = []
for i in range(3):
    cands.append((f'F_X{i}', 'deg5-J', Fi[i]))
for i in range(3):
    cands.append((f'F_Y{i}', 'deg5-J-Yodd', Fi[3 + i]))
B6 = {}
for i in range(3):
    B6 = padd(B6, pmul(X[i], Fi[i]))
cands.append(('B6=sum X_i F_Xi', 'deg6-J', B6))
for i in range(3):
    cands.append((f'Grad{i}=X{i}F_X{i}+Y{i}F_Y{i}', 'deg6-J', padd(pmul(X[i], Fi[i]), pmul(X[3 + i], Fi[3 + i]))))
for j in range(6):
    for i in range(6):
        cands.append((f'{VARS[j]}*F_{VARS[i]}', 'deg6-J-ZjFi', pmul(X[j], Fi[i])))
cands.append(('X0^2X1^2X2^2', 'deg6-monomial', {(2, 2, 2, 0, 0, 0): 1}))
cands.append(('Y0^2Y1^2Y2^2', 'deg6-monomial', {(0, 0, 0, 2, 2, 2): 1}))
K10 = padd(pmul({(0, 0, 0, 0, 0, 5): 1}, Fi[4]), pmul({(0, 0, 0, 0, 5, 0): 1}, Fi[5]), -1)
cands.append(('K10=Y2^5 F_Y1 - Y1^5 F_Y2', 'deg10', K10))
for i in range(6):
    for j in range(i, 6):
        cands.append((f'F_{VARS[i]}*F_{VARS[j]}', 'deg10-FiFj', pmul(Fi[i], Fi[j])))
done = set()
if os.path.exists(xlib.LOG):
    for l in open(xlib.LOG):
        r = json.loads(l)
        if r.get('exp') == 'A-reverify':
            done.add(r['D_hash'])
for name, fam, D in cands:
    h = pstr_hash(D)
    if h in done:
        continue
    t0 = time.time()
    r = solve_fixed_D(F2, A, q, D)
    r.update(exp='A-reverify', name=name, family=fam, yparity=ypar(D))
    if '_farkas_full' in r:
        fc = r['_farkas_full']
        json.dump(dict(name=name, D=to_str(D), q=q, F2=to_str(F2), A=to_str(A), H=int(r['H']), prime_used=r['primes'][-1], rows=fc['rows'], Y=fc['Y'], ell=fc['ell'],
                       note='Y^T M_rows = 0, Y^T b_rows = ell, exact integers; M,b from xlib.FixedDSystem (sector-restricted, ell=2)'),
                  open(f'/tmp/hodge/results/cycle-q2/wave2-excl/farkas/{fam}_{h}.json', 'w'))
    print(name, fam, 'd', r['d'], r['status'], r['ranks'], r.get('farkas'), f"{r['seconds']}s", flush=True)
    xlib.log(r)
