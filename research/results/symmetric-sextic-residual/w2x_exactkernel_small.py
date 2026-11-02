"""Exact (rational) E1-kernel check for small degrees: {D in S_d^{pD}: A D in J_{d+6}} = J_d over Q, via fmpq rank."""
import sys, time, json, math
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import w2x_lib as L
from hodge import core
from fractions import Fraction
from flint import fmpq_mat, fmpq
F, A, q = L.load()
Fi = [core.pderiv(F, i) for i in range(6)]
for d in [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '1,2,3').split(',')]:
    for pD in (0, 1):
        t0 = time.time()
        tp = (1 + pD) % 2
        rowmons = core.monomials(d + 6, tp); idx = {m: i for i, m in enumerate(rowmons)}
        cols = []
        for i in range(6):
            pm = tp if i < 3 else (tp + 1) % 2
            for m in core.monomials(d + 1, pm): cols.append(core.shift(Fi[i], m, 1))
        nJ = len(cols)
        Dm = core.monomials(d, pD)
        for m in Dm: cols.append(core.shift(A, m, 1))
        M = fmpq_mat(len(rowmons), len(cols))
        for j, P in enumerate(cols):
            for e, c in P.items():
                c = Fraction(c); M[idx[e], j] = fmpq(c.numerator, c.denominator)
        rk_all = M.rank()
        MJ = fmpq_mat(len(rowmons), nJ)
        for j in range(nJ):
            for e, c in cols[j].items():
                c = Fraction(c); MJ[idx[e], j] = fmpq(c.numerator, c.denominator)
        rk_J = MJ.rank()
        dimK = len(Dm) - (rk_all - rk_J)     # dim {D : A D in J}
        rec = {'exp': 'e1-kernel-exactQ', 'd': d, 'pD': pD, 'prime': None, 'rows': len(rowmons), 'nJ': nJ, 'nA': len(Dm), 'rank_J_Q': rk_J, 'rank_JA_Q': rk_all,
               'dimK_Q': dimK, 'ci_rank_J': L.dimJ(d + 6, tp), 'result': 'kernel zero over Q' if dimK == 0 else f'kernel dim {dimK} over Q', 'seconds': round(time.time() - t0, 1), 'job': 'exactQ'}
        L.log(rec); print(json.dumps(rec), flush=True)
