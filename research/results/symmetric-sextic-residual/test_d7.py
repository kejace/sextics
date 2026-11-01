import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
F, A, q = load()
def dlt(vs):
    return core.mul(core.mul(core.sub(vs[0], vs[1]), core.sub(vs[0], vs[2])), core.sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
e2sqX = core.sub(core.power(e_k(X,2),2), core.scale(core.mul(e_k(X,1), e_k(X,3)), 2))
e2sqY = core.sub(core.power(e_k(Y,2),2), core.scale(core.mul(e_k(Y,1), e_k(Y,3)), 2))
cands = [('D7=DX*e2(X^2)', core.mul(DX, e2sqX)), ('D7Y=DY*e2(Y^2)', core.mul(DY, e2sqY))]
for lab, D in cands:
    t0 = time.time()
    r = test_D(F, A, q, D, 'e1-kernel-deg7', note=lab + ' (unique E1-kernel element of degree 7, this parity)', fast=True, verbose=True, label=lab)
    print(lab, {k: v for k, v in r.items() if k not in ('N', 'Q')}, f'{time.time()-t0:.1f}s', flush=True)
    if r.get('consistent'):
        print('!!!!! CONSISTENT', lab, flush=True)
        if 'N' in r:
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in r['N']], 'Q': core.to_str(r['Q'])}, open(f'{RES}/witness-candidate.json', 'w'), indent=1)
