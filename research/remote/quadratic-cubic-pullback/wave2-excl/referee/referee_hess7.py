#!/usr/bin/env python3
"""Referee replay, part 2: (a) three of the 17 structured degree-7 Jacobian-type candidates against the necessary condition
A (w^T Hess(F) w) in J_14 at the fresh prime p = 1000000009 (J_14^odd block rank must equal the Hilbert value for the
failure to transfer to characteristic 0); (b) E1 kernels K_8 (both parities) at the same prime.
Own matrix construction (hodge.core only for parsing and dict-polynomial arithmetic).  Appends to referee-replay.jsonl.
usage: python3 referee_hess7.py [hess7|kernel8 ...]
"""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hodge import core
import referee_replay as R

P = R.P
PARTS = sys.argv[1:] or ['hess7', 'kernel8']
X0, X1, X2, Y0, Y1, Y2 = [{tuple(1 if j == i else 0 for j in range(6)): 1} for i in range(6)]
Q0, Q1, Q2, Q3, Q4, Q5 = R.Qs


def grad(h):
    return [core.pderiv(h, i) for i in range(6)]


def hess_target(w):
    T = {}
    for k in range(6):
        if not w[k]: continue
        for l in range(6):
            if not w[l]: continue
            T = core.add(T, core.mul(core.mul(w[k], w[l]), R.Fij[k][l]))
    return core.mul(R.A, T)


def part_hess7():
    t0 = time.time()
    cands = [('gradX0Q0', grad(core.mul(X0, Q0)), 'w = grad(X0 Q0), D Y-even'),
             ('wQX', [Q0, Q1, Q2, {}, {}, {}], 'w = (Q0,Q1,Q2,0,0,0), D Y-even'),
             ('wQ3e3', [{}, {}, {}, Q3, {}, {}], 'w = Q3 e_3 (D = Q3 F_3, reducible control), D Y-odd')]
    cols = R.monos(14, 1)
    J = R.jrows(14, 1, P)
    MJ = R.tomat(J, cols, P); rJ = MJ.rank(); del MJ
    hJ = R.dimJ(14, 1)
    R.log(dict(check='J14_odd_block_rank', prime=P, rows=len(J), cols=len(cols), rank_modp=rJ, hilbert=hJ, certified=(rJ == hJ), seconds=round(time.time() - t0, 1)))
    for name, w, note in cands:
        t1 = time.time()
        D = {}
        for i in range(6):
            if w[i]: D = core.add(D, core.mul(w[i], R.Fi[i]))
        T = hess_target(w)
        Tp = core.reduce_mod(T, P)
        assert all((e[3] + e[4] + e[5]) % 2 == 1 for e in Tp), 'target parity'
        M = R.tomat(J + [Tp], cols, P); rM = M.rank(); del M
        inJ = (rM == rJ)
        R.log(dict(check='hess7_structured_candidate', name=name, note=note, d=7, prime=P, D_terms=len(D), D_parity=core.parity_of(D), T_terms=len(T),
                   rankJ=rJ, rank_stack=rM, hess_in_J_modp=inJ, J_rank_certified=(rJ == hJ),
                   verdict=('EXCLUDED (rigorous: A w^T Hess F w not in J_14, block at Hilbert rank)' if (not inJ and rJ == hJ) else ('survives' if inJ else 'excluded mod p only')),
                   seconds=round(time.time() - t1, 1)))


def part_kernel8():
    for par in (0, 1):
        R.log(R.e1_kernel(8, par, P))


if __name__ == '__main__':
    for part in PARTS:
        globals()['part_' + part]()
    print('REFEREE HESS7/KERNEL8 DONE', flush=True)
