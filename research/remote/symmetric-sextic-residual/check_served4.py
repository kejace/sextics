import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
import numpy as np
from flint import nmod_mat
F, A, q = load()
p = P0
def served_polys(D):
    A0 = {core.monomials(6, 1)[0]: 1}
    S = InvSystem(F, A0, q, D, p=p, targets=True)
    M = np.stack(S.cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    tidx = [j for j, c in enumerate(S.colinfo) if c[0] == 'A']
    out = []
    for k in range(nullity):
        x = [int(Xn[j, k]) for j in range(S.ncols)]
        if not any(x[j] for j in tidx): continue
        xs = [(v if v < p//2 else v - p) for v in x]
        N, Q = S.unfold([Fraction(v) for v in xs[:S.nN + S.nQ]] + [Fraction(0)]*S.nT)
        Ap = {}
        for j in tidx:
            if xs[j]:
                Ap = core.add(Ap, core.scale(sym({S.colinfo[j][2]: Fraction(1)}), Fraction(xs[j], 36)))
        out.append((Ap, N, Q))
    return out
for lab, D in [('QX', e_k(X, 2)), ('p3X+e3X', core.add(p_k(X,3), e_k(X,3)))]:
    sols = served_polys(D)
    print(lab, 'nullspace vectors with A part:', len(sols))
    for Ap, N, Q in sols:
        ok, info = core.verify_witness(F, Ap, q, D, N, Q, p=p)
        ok1, info1 = fixed_d.e1_solvable(F, Ap, q, core.parse('1'))
        print('   verify mod p:', ok, '| A\' in J (E1 with D=1):', ok1, '| A\' terms', len(Ap))
