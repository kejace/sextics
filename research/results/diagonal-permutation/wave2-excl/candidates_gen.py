"""Structured invariant candidate denominators (dossier §4.2 + attack1) as expressions -> candidates.jsonl"""
import itertools, json
from fractions import Fraction as Fr
from common import *
from hodge.core import mul, add, sub, scale, power, parse, to_str
F, A, q, rec = core.load_instance(IID)
def V(i): return {tuple(1 if j == i else 0 for j in range(6)): Fr(1)}
X = [V(0), V(1), V(2)]; Y = [V(3), V(4), V(5)]
def S(*ps):
    out = {}
    for p_ in ps: out = add(out, p_)
    return out
def prod(*ps):
    out = {(0,)*6: Fr(1)}
    for p_ in ps: out = mul(out, p_)
    return out
perms = list(itertools.permutations(range(3)))
O1 = S(*[prod(power(X[i], 3), X[j], Y[i], Y[j]) for (i, j, k) in perms])
O2 = S(*[prod(power(X[i], 3), X[j], power(Y[k], 2)) for (i, j, k) in perms])
O3 = S(*[prod(power(Y[i], 3), Y[j], power(X[k], 2)) for (i, j, k) in perms])
p2X = S(*[power(X[i], 2) for i in range(3)]); p2Y = S(*[power(Y[i], 2) for i in range(3)])
pXY = S(*[mul(X[i], Y[i]) for i in range(3)])
X6 = S(*[power(X[i], 6) for i in range(3)]); Y6 = S(*[power(Y[i], 6) for i in range(3)])
X4 = S(*[power(X[i], 4) for i in range(3)]); X2Y2 = S(*[mul(power(X[i], 2), power(Y[i], 2)) for i in range(3)])
DX = prod(sub(X[0], X[1]), sub(X[0], X[2]), sub(X[1], X[2])); DY = prod(sub(Y[0], Y[1]), sub(Y[0], Y[2]), sub(Y[1], Y[2]))
qp = prod(*[add(power(sub(X[i], X[j]), 2), power(sub(Y[i], Y[j]), 2)) for i, j in [(0, 1), (0, 2), (1, 2)]])
qm = prod(*[sub(power(sub(X[i], X[j]), 2), power(sub(Y[i], Y[j]), 2)) for i, j in [(0, 1), (0, 2), (1, 2)]])
Fc = sub(sub(F, X6), Y6)
T = [prod(power(X[i], 3), power(Y[i], 3)) for i in range(3)]
sumT = S(*T); B2 = sub(add(T[0], T[1]), scale(T[2], 2))
DX2, DY2 = power(DX, 2), power(DY, 2)
cands = {
 # degree 2
 'p2X': p2X, 'p2Y': p2Y, 'pXY': pXY,
 # degree 4
 'p2X^2': power(p2X, 2), 'p2X*p2Y': mul(p2X, p2Y), 'pXY^2': power(pXY, 2), 'X4': X4, 'X2Y2': X2Y2,
 # degree 6
 'DX^2': DX2, 'DY^2': DY2, 'DX*DY': mul(DX, DY), 'qp': qp, 'qm': qm, 'O1': O1, 'O2': O2, 'O3': O3, 'X6': X6, 'Y6': Y6,
 'X6-Y6': sub(X6, Y6), 'p2X^3': power(p2X, 3), 'p2X^2*p2Y': mul(power(p2X, 2), p2Y), 'p2X*pXY^2': mul(p2X, power(pXY, 2)),
 'Fc': Fc, 'sumT': sumT, 'pXY^3': power(pXY, 3), 'B2': B2,
 # degree 8
 'DX^2*p2X': mul(DX2, p2X), 'DX^2*p2Y': mul(DX2, p2Y), 'DX^2*pXY': mul(DX2, pXY), 'O1*p2X': mul(O1, p2X), 'O2*p2X': mul(O2, p2X),
 'O3*p2X': mul(O3, p2X), 'O1*p2Y': mul(O1, p2Y), 'X6*p2X': mul(X6, p2X), 'qp*p2X': mul(qp, p2X), 'p2X^4': power(p2X, 4),
 'sumT*pXY': mul(sumT, pXY), 'B2*pXY': mul(B2, pXY), 'A*pXY': mul(A, pXY),
 # degree 10
 'O1*p2X^2': mul(O1, power(p2X, 2)), 'O1*p2X*p2Y': mul(O1, mul(p2X, p2Y)), 'O1*pXY^2': mul(O1, power(pXY, 2)),
 'DX^2*p2X^2': mul(DX2, power(p2X, 2)), 'DX^2*p2Y^2': mul(DX2, power(p2Y, 2)), 'DX^2*pXY^2': mul(DX2, power(pXY, 2)),
 'DX*DY*p2X*p2Y': prod(DX, DY, p2X, p2Y), 'DX*DY*pXY*p2X': prod(DX, DY, pXY, p2X), 'qp*p2X^2': mul(qp, power(p2X, 2)),
 'qm*p2X^2': mul(qm, power(p2X, 2)), 'X6*p2X^2': mul(X6, power(p2X, 2)), 'X6*pXY^2': mul(X6, power(pXY, 2)),
 'O2*p2Y^2': mul(O2, power(p2Y, 2)), 'O3*p2X^2': mul(O3, power(p2X, 2)), 'Fc*p2X^2': mul(Fc, power(p2X, 2)),
 'Fc*pXY^2': mul(Fc, power(pXY, 2)), '(X6+Y6)*p2X^2': mul(add(X6, Y6), power(p2X, 2)), 'p2X^5': power(p2X, 5),
 'sumT*pXY*p2X': prod(sumT, pXY, p2X), 'B2*pXY*p2X': prod(B2, pXY, p2X),
 # degree 12
 'DX^2*DY^2': mul(DX2, DY2), 'O1*O2': mul(O1, O2), 'O1*O3': mul(O1, O3), 'O2*O3': mul(O2, O3), 'O1^2': power(O1, 2),
 'O2^2': power(O2, 2), 'O3^2': power(O3, 2), 'O1*X6': mul(O1, X6), 'O1*Y6': mul(O1, Y6), 'O2*X6': mul(O2, X6), 'O3*Y6': mul(O3, Y6),
 'X6*Y6': mul(X6, Y6), 'X6^2': power(X6, 2), '(X6-Y6)^2': power(sub(X6, Y6), 2), 'qp*qm': mul(qp, qm), 'qp^2': power(qp, 2),
 'qm^2': power(qm, 2), 'DX^2*X6': mul(DX2, X6), 'DX^2*Y6': mul(DX2, Y6), 'DY^2*X6': mul(DY2, X6), 'DX*DY*O1': prod(DX, DY, O1),
 'DX*DY*X6': prod(DX, DY, X6), 'p2X^3*p2Y^3': mul(power(p2X, 3), power(p2Y, 3)), 'Fc^2': power(Fc, 2), 'Fc*X6': mul(Fc, X6),
 'Fc*Y6': mul(Fc, Y6), 'A^2': power(A, 2), 'sumT^2': power(sumT, 2), 'sumT2': S(*[power(t, 2) for t in T]), 'B2^2': power(B2, 2),
 'p2X^2*p2Y^2*pXY^2': prod(power(p2X, 2), power(p2Y, 2), power(pXY, 2)), 'DX^2*DY^2-qp*qm': sub(mul(DX2, DY2), mul(qp, qm)),
 'O1*sumT': mul(O1, sumT), 'A*B2': mul(A, B2), 'T0*T1': mul(T[0], T[1]), 'sumT*(X6-Y6)': mul(sumT, sub(X6, Y6)),
}
with open('candidates.jsonl', 'w') as f:
    for k, v in cands.items():
        if core.parity_of(v) is None or not core.is_homogeneous(v): print('skip', k); continue
        f.write(json.dumps({'name': k, 'd': core.degree(v), 'parity': core.parity_of(v), 'D': to_str(v)}) + '\n')
print(len(cands), 'candidates written')
