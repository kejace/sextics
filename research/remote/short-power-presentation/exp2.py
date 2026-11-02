import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
F, A, q, rec = load()
which = sys.argv[1]
L6 = core.power(L, 6); M6 = core.power(M, 6)
L4M2 = core.mul(core.power(L, 4), core.power(M, 2)); L2M4 = core.mul(core.power(L, 2), core.power(M, 4))
cands = {}
if which == 'a':
    cands['D_e=eX-L^6'] = core.sub(eX, L6)
    cands['D_e=eX+L^6'] = core.add(eX, L6)
    cands["D'_s=eY-M^6"] = core.sub(eY, M6)
    cands['eY'] = eY
    cands['eX+M^6'] = core.add(eX, M6)
    cands['eY+L^6'] = core.add(eY, L6)
else:
    EXF = {}
    for i in range(3):
        EXF = core.add(EXF, core.mul({tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}, core.pderiv(F, i)))
    EYF = {}
    for i in range(3, 6):
        EYF = core.add(EYF, core.mul({tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}, core.pderiv(F, i)))
    cands['E_X(F)'] = EXF
    cands['E_Y(F)'] = EYF
    cands['eX+L^4M^2'] = core.add(eX, L4M2)
    cands['eX+L^2M^4'] = core.add(eX, L2M4)
    cands['eX-2*L^6+M^6'] = core.add(core.sub(eX, core.scale(L6, 2)), M6)
    cands['eX+3*L^4M^2-M^6'] = core.sub(core.add(eX, core.scale(L4M2, 3)), M6)
for name, D in cands.items():
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    dim, basis, S, dt = served(F, q, D)
    sec = time.time() - t0
    log(core.to_str(D), core.degree(D), 'base-pullback: ' + name, ok, None, info['rank'], info['cols'] - info['rank'], sec,
        f'served-dim={dim} (18 = Jacobian only); E1 info {info}')
    print(f'{name}: E1={ok} served_dim={dim} ({sec:.0f}s)', flush=True)
