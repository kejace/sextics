"""Batch: structured semi-invariant candidates (products of symmetric atoms, discriminants, E1 classes, Hessian blocks),
degrees 6..12: twisted E1 test, then exact reduced solve for every E1-passing D."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import w2lib as W, charsolve as cs
from h import var, X, Y, e_k, p_k
from hodge import core
from fractions import Fraction
F, A, q = W.load()
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
def radial(P, idx):
    out = {}
    for i in idx: out = add(out, mul(var(i), core.pderiv(P, i)))
    return out
def det3(M):
    t = {}
    for perm, s in [((0,1,2),1),((1,2,0),1),((2,0,1),1),((0,2,1),-1),((2,1,0),-1),((1,0,2),-1)]:
        t = add(t, sc(mul(mul(M[0][perm[0]], M[1][perm[1]]), M[2][perm[2]]), s))
    return t
DX, DY = dlt(X), dlt(Y)
e1X, e2X, e3X = e_k(X,1), e_k(X,2), e_k(X,3); e1Y, e2Y, e3Y = e_k(Y,1), e_k(Y,2), e_k(Y,3)
D7 = mul(DX, sub(pw(e2X,2), sc(mul(e1X,e3X),2))); D7Y = mul(DY, sub(pw(e2Y,2), sc(mul(e1Y,e3Y),2)))
EXF, EYF = radial(F,(0,1,2)), radial(F,(3,4,5))
H = [[core.pderiv(core.pderiv(F,i),j) for j in range(6)] for i in range(6)]
HXX = det3([[H[i][j] for j in range(3)] for i in range(3)]); HYY = det3([[H[i+3][j+3] for j in range(3)] for i in range(3)])
HXY = det3([[H[i][j+3] for j in range(3)] for i in range(3)])
s2X = sub(sc(p_k(X,2),3), pw(e1X,2)); s2Y = sub(sc(p_k(Y,2),3), pw(e1Y,2))
atoms = {'e1X':(e1X,1),'e2X':(e2X,2),'e3X':(e3X,3),'e1Y':(e1Y,1),'e2Y':(e2Y,2),'e3Y':(e3Y,3),'DX':(DX,3),'DY':(DY,3),
         'D7':(D7,7),'D7Y':(D7Y,7),'EXF':(EXF,6),'EYF':(EYF,6),'HXX':(HXX,12),'HYY':(HYY,12),'HXY':(HXY,12),'s2X':(s2X,2),'s2Y':(s2Y,2),
         'p6X':(p_k(X,6),6),'p6Y':(p_k(Y,6),6),'p3X':(p_k(X,3),3),'p3Y':(p_k(Y,3),3)}
names = list(atoms)
# enumerate exponent vectors with total degree in [6,12] (cap exponents of big atoms)
cands = {}
def rec(i, deg, cur):
    if deg > 12: return
    if i == len(names):
        if 6 <= deg: cands[tuple(cur)] = deg
        return
    nm = names[i]; a = atoms[nm][1]
    emax = 12 // a
    if nm in ('EXF','EYF','p6X','p6Y','HXX','HYY','HXY','D7','D7Y'): emax = min(emax, 1)
    if nm in ('p3X','p3Y','s2X','s2Y'): emax = min(emax, 2)
    for e in range(emax + 1):
        rec(i + 1, deg + e * a, cur + [e])
rec(0, 0, [])
print('candidates', len(cands), flush=True)
groups = {c: W.group(c) for c in W.CHARS}
t_all = time.time(); n_e1 = 0; n_solve = 0; passed = []
log = open('/tmp/hodge/results/symmetric-sextic-residual/w2struct.jsonl', 'a')
done = set()
try:
    for l in open('/tmp/hodge/results/symmetric-sextic-residual/w2struct.jsonl'):
        done.add(json.loads(l)['label'])
except Exception: pass
for ev, deg in sorted(cands.items(), key=lambda kv: kv[1]):
    label = '*'.join(f'{n}^{e}' if e > 1 else n for n, e in zip(names, ev) if e)
    if label in done: continue
    D = {(0,)*6: Fraction(1)}
    for n, e in zip(names, ev):
        if e: D = mul(D, pw(atoms[n][0], e))
    if not D: continue
    pD = core.parity_of(D)
    # character: sign atoms DX, D7, HXY (X-sign), DY, D7Y, HXY (Y-sign)
    sx = (-1) ** sum(e for n, e in zip(names, ev) if n in ('DX','D7','HXY'))
    sy = (-1) ** sum(e for n, e in zip(names, ev) if n in ('DY','D7Y','HXY'))
    chi = {(1,1):'triv',(-1,1):'sgnX',(1,-1):'sgnY',(-1,-1):'sgnsgn'}[(sx,sy)]
    G = groups[chi]
    t0 = time.time()
    ok, info = cs.e1_char(F, A, q, D, G)
    n_e1 += 1
    r = {'label': label, 'd': deg, 'pD': pD, 'chi': chi, 'e1': ok, 'seconds': round(time.time()-t0, 2)}
    if ok:
        divF, _ = __import__('hodge.jacobian', fromlist=['x']).divisible_by(F, D)
        r['F_divides_D'] = divF
        if not divF:
            st, data = cs.solve_char(F, A, q, D, G)
            n_solve += 1
            r.update(status=st, rank=data.get('rank'), nullity=data.get('nullity'), rows=data.get('rows'), cols=data.get('cols'))
            r['seconds'] = round(time.time()-t0, 2)
            if st == 'solved':
                okv, _ = core.verify_witness(F, A, q, D, data['N'], data['Q'])
                r['verified'] = okv
                if okv:
                    W.save_witness_candidate(D, data['N'], data['Q'], extra={'label': label, 'batch': 'struct'})
                    open('/tmp/hodge/results/symmetric-sextic-residual/WITNESS_FOUND.txt','a').write(json.dumps(r)+'\n')
        passed.append(r)
    log.write(json.dumps(r) + '\n'); log.flush()
    print(json.dumps(r), flush=True)
W.wlog({'batch': 'struct', 'family': 'atom-products', 'degrees': '6-12', 'n_candidates': len(cands), 'n_e1_tested': n_e1,
        'n_e1_pass': len(passed), 'n_full_solves': n_solve, 'passed': [p['label'] for p in passed],
        'solved': [p['label'] for p in passed if p.get('status') == 'solved'], 'seconds': round(time.time()-t_all, 1)})
