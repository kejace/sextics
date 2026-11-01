import sys, random
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc import *
F, A, q = load(); p = P0
random.seed(1)
D = core.parse('X0**5*X1 + X1**5*X2 + X0*X2**5 + 3*X0*X1**3*Y0*Y1')
Fp, Ap, Dp = core.reduce_mod(F,p), core.reduce_mod(A,p), core.reduce_mod(D,p)
S = WSystem(F, A, q, D, p, weights=True)
N = [{} for _ in range(NV)]
for (kind, i, m) in S.cols:
    if kind == 'N' and random.random() < 0.3: N[i][m] = random.randint(1, 50)
E = core.defect(Fp, Ap, Dp, q, N, p)
lhs = core.sub(E, core.mul(Ap, core.mul(Dp, Dp, p), p), p)
acc = {}
for j, (kind, i, m) in enumerate(S.cols):
    v = N[i].get(m, 0) if kind == 'N' else 0
    if v:
        for e, c in S.colpolys[j].items(): acc[e] = (acc.get(e, 0) + c*v) % p
acc = {e: c for e, c in acc.items() if c}
print('column construction matches core.defect (reduced inputs):', not core.sub(acc, lhs, p))
# klein-hesse D_15 and friends: E1 block test on OUR F
import json
w = json.load(open('/home/kejace/devel/sextics/data/klein-hesse/witness.json'))
D15 = core.parse(w['D'])
print('klein-hesse D: deg', core.degree(D15), 'sigma-weight', sigma_weight(D15), 'parity', core.parity_of(D15))
import time
for name, Dx in [('klein-hesse D_15', D15)]:
    t = time.time(); ok, info = e1_block(F, A, q, Dx, p); print(name, 'E1:', ok, info, f'{time.time()-t:.1f}s')
    log({'D': name, 'd': core.degree(Dx), 'family': 'companion', 'engine': 2, 'k': sigma_weight(Dx), 'e1': ok, 'consistent': None, 'rank': info['rank'], 'nullity': None, 'seconds': round(time.time()-t,2), 'note': 'E1 block only; ' + str(info)})
