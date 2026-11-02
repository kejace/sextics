"""Rigorous necessary test for structured polar candidates D = w.grad F:  A . d_w(D) in J(F)_{2d}  (T5).
A failure is rigorous when the degree-2d Jacobian block has its Hilbert rank mod p.  A pass means the candidate
survives this test (the full fixed-D system then decides).   usage: python3 wave2_excl_hesstest.py <d> [primeindex]
"""
import sys, os, json, time
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, Qs, var, mul, add, sub, scale, power, pderiv, logline, dim_J
from hodge import jacobian
from fractions import Fraction
D_DEG = int(sys.argv[1]); PI = int(sys.argv[2]) if len(sys.argv) > 2 else 1
P = core.PRIMES[PI]
LOG = os.path.join(w2x.LANE, 'wave2-excl.hesstest.jsonl')
X0, X1, X2, Y0, Y1, Y2 = [var(i) for i in range(6)]
Q0, Q1, Q2, Q3, Q4, Q5 = Qs
Z = [X0, X1, X2, Y0, Y1, Y2]


def grad(h):
    return [pderiv(h, i) for i in range(6)]


def cands(d):
    C = []
    def addw(name, fam, w, note=''):
        C.append((name, fam, w, note))
    if d == 6:
        addw('EX', 'polar6-blockEuler', [X0, X1, X2, {}, {}, {}])
        for k in range(6): addw(f'gradQ{k}', 'polar6-gradQ', grad(Qs[k]))
    if d == 7:
        addw('gradX0Q0', 'polar7-gradh', grad(mul(X0, Q0))); addw('gradX0Q3', 'polar7-gradh', grad(mul(X0, Q3)))
        addw('gradX1Q0', 'polar7-gradh', grad(mul(X1, Q0))); addw('gradY0Q0-odd', 'polar7-gradh-odd', grad(mul(Y0, Q0)))
        addw('gradX0cube', 'polar7-gradh', grad(power(X0, 3))); addw('gradX0Y0Y1', 'polar7-gradh', grad(mul(X0, mul(Y0, Y1))))
        addw('wQX', 'polar7-Qfield', [Q0, Q1, Q2, {}, {}, {}]); addw('wQY-odd', 'polar7-Qfield-odd', [{}, {}, {}, Q3, Q4, Q5])
        addw('wQshift', 'polar7-Qfield', [Q1, Q2, Q0, {}, {}, {}]); addw('wQcross', 'polar7-Qfield', [Q3, Q4, Q5, {}, {}, {}])
        addw('wQmix-odd', 'polar7-Qfield-odd', [{}, {}, {}, Q0, Q1, Q2])
        for k in range(6):
            addw(f'wQ{k}e{k}', 'polar7-Qdiag', [Qs[k] if i == k else {} for i in range(6)], 'w = Q_k e_k (D = Q_k F_k, reducible control)')
    if d == 8:
        addw('gradQ0Q1', 'polar8-gradh', grad(mul(Q0, Q1))); addw('gradQ0Q3', 'polar8-gradh', grad(mul(Q0, Q3)))
        addw('gradQ3Q4', 'polar8-gradh', grad(mul(Q3, Q4))); addw('gradQ0Q5', 'polar8-gradh', grad(mul(Q0, Q5)))
        sq = {}
        for k in range(6): sq = add(sq, mul(Qs[k], Qs[k]))
        addw('gradSumQsq', 'polar8-gradh', grad(sq)); addw('gradX0sqY0sq', 'polar8-gradh', grad(mul(mul(X0, X0), mul(Y0, Y0))))
        addw('gradX0quart', 'polar8-gradh', grad(power(X0, 4)))
        addw('wQEuler', 'polar8-Qfield', [mul(Qs[i], Z[i]) for i in range(6)]); addw('wQEulerX', 'polar8-Qfield', [mul(Qs[i], Z[i]) if i < 3 else {} for i in range(6)])
        addw('wQcyc', 'polar8-Qfield', [mul(Q0, X1), mul(Q1, X2), mul(Q2, X0), mul(Q3, Y1), mul(Q4, Y2), mul(Q5, Y0)])
    if d == 9:
        addw('gradQ0Q1Q2', 'polar9-gradh', grad(mul(mul(Q0, Q1), Q2))); addw('gradQ0Q1Q3', 'polar9-gradh', grad(mul(mul(Q0, Q1), Q3)))
        addw('gradQ3Q4Q5', 'polar9-gradh', grad(mul(mul(Q3, Q4), Q5))); addw('gradX0cubeY0sq', 'polar9-gradh', grad(mul(power(X0, 3), mul(Y0, Y0))))
        addw('wQsqX', 'polar9-Qfield', [mul(Qs[i], Qs[i]) if i < 3 else {} for i in range(6)])
        addw('wQsqCross', 'polar9-Qfield', [mul(Qs[i + 3], Qs[i + 3]) if i < 3 else {} for i in range(6)])
        addw('wQsqY-odd', 'polar9-Qfield-odd', [{} if i < 3 else mul(Qs[i], Qs[i]) for i in range(6)])
    return C


def main():
    t0 = time.time()
    d = D_DEG
    Fi = [pderiv(F, i) for i in range(6)]
    Fij = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    reducers = {}
    for name, fam, w, note in cands(d):
        t1 = time.time()
        D = w2x.polar(w)
        par = core.parity_of(D)
        # T = A sum_{k,l} w_k w_l F_kl
        T = {}
        for k in range(6):
            if not w[k]: continue
            for l in range(6):
                if not w[l]: continue
                T = add(T, mul(w[k], mul(w[l], Fij[k][l])))
        T = mul(A, T)
        tpar = core.parity_of(T)
        if tpar not in reducers:
            tt = time.time()
            reducers[tpar] = jacobian._Reducer(F, 2 * d, P, tpar)
            print(f'reducer deg {2*d} parity {tpar}: rank {reducers[tpar].rank} (Hilbert {dim_J(2*d, tpar)}) in {time.time()-tt:.0f}s', flush=True)
        red = reducers[tpar]
        nf = red.normal_form(core.reduce_mod(T, P))
        inJ = (len(nf) == 0)
        cert = (red.rank == dim_J(2 * d, tpar))
        rec = dict(exp='hesstest', name=name, family=fam, d=d, parity=par, prime=P, D_terms=len(D), hess_in_J=inJ, J_rank=red.rank,
                   hilbJ=dim_J(2 * d, tpar), J_rank_certified=cert, result=('survives-hessian-test' if inJ else ('EXCLUDED-rigorous' if cert else 'excluded-modular')),
                   note=note, seconds=round(time.time() - t1, 1))
        logline(LOG, rec); print(json.dumps(rec), flush=True)
    print('ALL DONE', round(time.time() - t0), 's', flush=True)


if __name__ == '__main__':
    main()
