import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
import numpy as np
from flint import nmod_mat
F, A, q = load()
p = P0
QX = e_k(X, 2)
A0 = {core.monomials(6, 1)[0]: 1}
S = InvSystem(F, A0, q, QX, p=p, targets=True)
M = np.stack(S.cols, axis=1)
Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
Xn, nullity = Mn.nullspace()
tidx = [j for j, c in enumerate(S.colinfo) if c[0] == 'A']
print('nullity', nullity, 'ncols', S.ncols, 'nN', S.nN, 'nQ', S.nQ, 'nT', S.nT)
# pick nullspace vector with nonzero A part
for k in range(nullity):
    x = [int(Xn[j, k]) for j in range(S.ncols)]
    if any(x[j] for j in tidx):
        break
else:
    print('no vector with A part'); sys.exit()
xs = [(v if v < p//2 else v - p) for v in x]
N, Q = S.unfold([Fraction(v) for v in xs[:S.nN + S.nQ]] + [Fraction(0)]*S.nT)
Ap = {}
for j in tidx:
    if xs[j]:
        Ap = core.add(Ap, core.scale(orbit_sum(S.colinfo[j][2]), Fraction(xs[j])))
print('A\' terms', len(Ap), 'N terms', [len(n) for n in N], 'Q', len(Q))
# the inv system equation: sum cols * x = 0 means: (N-cols)(x) + (Q-cols)(x) + D^2 A' = 0.
# N-col for (m, slot) encodes  -D F_s m - F D_s m + FD d_s m  ; so identity: -D N(F) - F N(D) + F D div N - F^2 Q + A' D^2 = 0
# which is exactly E = A' D^2 - q D N(F) + F (D div N - N(D)) - F^2 Q = 0 with the Q sign as built (Q cols = -F^2 m).
ok, info = core.verify_witness(F, Ap, q, QX, N, Q, p=p)
print('verify mod p:', ok, info)
ok2, info2 = core.verify_witness(F, Ap, q, QX, N, core.scale(Q, -1), p=p)
print('verify mod p with -Q:', ok2, info2)
# compute E directly mod p
D = QX
E = core.mul(Ap, core.mul(D, D, p), p)
NF = {}
for i in range(6): NF = core.add(NF, core.mul(N[i], core.pderiv(F, i, p), p), p)
E = core.sub(E, core.mul(D, NF, p), p)
divN = {}
for i in range(6): divN = core.add(divN, core.pderiv(N[i], i, p), p)
ND = {}
for i in range(6): ND = core.add(ND, core.mul(N[i], core.pderiv(D, i, p), p), p)
E = core.add(E, core.mul(F, core.sub(core.mul(D, divN, p), ND, p), p), p)
F2Q = core.mul(core.mul(F, F, p), Q, p)
r1 = core.sub(E, F2Q, p); r2 = core.add(E, F2Q, p)
print('E - F^2 Q terms:', len(core.reduce_mod(r1, p)), ' E + F^2 Q terms:', len(core.reduce_mod(r2, p)), ' E terms:', len(core.reduce_mod(E, p)))
# check parities of N
print('N parities:', [core.parity_of(n) if n else None for n in N], 'A\' parity', core.parity_of(Ap))
