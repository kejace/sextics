"""Degree-7 E2 batch: all l*F_i (36), the 8 extra Y-even E1-kernel directions, and random combos."""
import sys, time, pickle, random
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from orslib import *
random.seed(1)
cands = []
for i in range(6):
    for j in range(6):
        D = core.mul(X[j], Fi(i))
        cands.append((f'{core.VARS[j]}*dF/d{core.VARS[i]}', 'J-multiple linear*partial', D))
for name, fam, D in cands:
    r = test(name, fam, note='deg7 J-multiple', Dpoly=D, verbose=False)
    print(name, r['d'], 'parity', core.parity_of(D), r['status'], r.get('rank'), r.get('nullity'), round(r['seconds'],1), flush=True)
# extra kernel directions (mod p basis, lift to small integers by symmetric remainder)
with open('/tmp/hodge/results/one-reflection-skew/e1ker_d7.pkl','rb') as f:
    K = pickle.load(f)
p = K['p']
def lift(poly):
    out = {}
    for e, v in poly.items():
        v = v if v <= p//2 else v - p
        out[e] = Fraction(v)
    return out
# reduce basis mod J to get representatives: basisD already rref'd; take them as is but note they contain J parts
for k, b in enumerate(K['basisD']):
    D = lift(b)
    # check it is genuinely in the kernel with small lift: test e1
    name = f'<E1ker d7 basis #{k}: {len(D)} terms>'
    r = test(name, 'E1-kernel-basis-d7', note='rref basis vector of E1 kernel (lifted symmetric residues)', Dpoly=D)
    print(name, r['e1'], r['status'], r.get('rank'), r.get('nullity'), round(r['seconds'],1), flush=True)
