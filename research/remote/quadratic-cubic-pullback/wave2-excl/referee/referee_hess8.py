#!/usr/bin/env python3
"""Referee replay, part 3 (remote): two of the 10 structured degree-8 Jacobian-type candidates against the necessary
condition A (w^T Hess(F) w) in J_16 at the fresh prime p = 1000000009, with the J_16^odd block rank certified against
the Hilbert value (9684).  Own matrix construction.  Appends to referee-replay-remote.jsonl next to this file.
"""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hodge import core
import referee_replay as R
R.OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'referee-replay-remote.jsonl')
P = R.P
Z = [{tuple(1 if j == i else 0 for j in range(6)): 1} for i in range(6)]
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


t0 = time.time()
cands = [('gradQ0Q1', grad(core.mul(Q0, Q1)), 'w = grad(Q0 Q1), D Y-even, degree 8'),
         ('wQEuler', [core.mul(R.Qs[i], Z[i]) for i in range(6)], 'w_i = Q_i z_i, D Y-even, degree 8')]
cols = R.monos(16, 1)
J = R.jrows(16, 1, P)
R.log(dict(check='J16_odd_block_build', rows=len(J), cols=len(cols), seconds=round(time.time() - t0, 1)))
MJ = R.tomat(J, cols, P); rJ = MJ.rank(); del MJ
hJ = R.dimJ(16, 1)
R.log(dict(check='J16_odd_block_rank', prime=P, rows=len(J), cols=len(cols), rank_modp=rJ, hilbert=hJ, certified=(rJ == hJ), seconds=round(time.time() - t0, 1)))
for name, w, note in cands:
    t1 = time.time()
    D = {}
    for i in range(6):
        if w[i]: D = core.add(D, core.mul(w[i], R.Fi[i]))
    T = hess_target(w)
    Tp = core.reduce_mod(T, P)
    M = R.tomat(J + [Tp], cols, P); rM = M.rank(); del M
    inJ = (rM == rJ)
    R.log(dict(check='hess8_structured_candidate', name=name, note=note, d=8, prime=P, D_terms=len(D), D_parity=core.parity_of(D), T_terms=len(T),
               rankJ=rJ, rank_stack=rM, hess_in_J_modp=inJ, J_rank_certified=(rJ == hJ),
               verdict=('EXCLUDED (rigorous: A w^T Hess F w not in J_16, block at Hilbert rank)' if (not inJ and rJ == hJ) else ('survives' if inJ else 'excluded mod p only')),
               seconds=round(time.time() - t1, 1)))
print('REFEREE HESS8 DONE', round(time.time() - t0), 's', flush=True)
