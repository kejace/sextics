#!/usr/bin/env python3
"""First-order survival of the WHOLE served space of a Fermat cone under F_t = F0 + t*Qdir:
for D0 = f_S * X0^(d-6), the (3,0)-sector classes served at t = 0 (join classes), find the subspace of target directions
alpha for which the t = 0 certificate deforms to first order:  L(dD,dN,dQ) + sum_t alpha_t T_t = 0 solvable.
  python3 w2_firstorder_free.py <d> <dir>     (dir in Qx, chain, extra, link1, link2)"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
import json, time, itertools
from fractions import Fraction
import w2lib, oplib
from hodge import core
from hodge.core import NV, add, sub, scale, mul, pderiv, shift
from flint import nmod_mat
F, A, q = w2lib.load(); p = oplib.P0
d = int(sys.argv[1]); DIR = sys.argv[2] if len(sys.argv) > 2 else 'Qx'
f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6'); F0 = add(f6, g6)
Qx = {'Qx': sub(F, F0), 'chain': core.parse('X0**3*X1*Y0*Y1/2 + X1**3*X2*Y1*Y2/2'), 'extra': core.parse('X0**2*X1*X2*Y0*Y2/3'),
      'link1': core.parse('X0**3*X1*Y0*Y1/2'), 'link2': core.parse('X1**3*X2*Y1*Y2/2')}[DIR]
F0i = [pderiv(F0, i) for i in range(NV)]; Qxi = [pderiv(Qx, i) for i in range(NV)]
c = (0, 0); cA = (3, 0)
m = {tuple([d - 6, 0, 0, 0, 0, 0]): Fraction(1)}
Tm = w2lib.monos_char(6, cA)
out = []
CONES = {'f6': (0, 1, 2), 'g6': (3, 4, 5), 'X0X1Y1': (0, 1, 4), 'X0X2Y2': (0, 2, 5)}
for S in [CONES[sys.argv[3]] if len(sys.argv) > 3 else (0, 1, 2)]:
    fS = {tuple(6 if t == i else 0 for t in range(NV)): Fraction(1) for i in S}
    D0 = mul(fS, m)
    t0 = time.time()
    served = []
    for Ap_mon in Tm:
        Ap = {Ap_mon: Fraction(1)}
        st, data = oplib.lift_solution(F0, Ap, q, D0, max_primes=6)
        if st == 'solved':
            served.append((Ap_mon, data['N'], data['Q']))
    print(f'd={d} dir={DIR}: {len(served)} of {len(Tm)} sector monomials served at t=0 by {core.to_str(D0)[:40]}', flush=True)
    # T_t for each served target (exact)
    D0i = [pderiv(D0, i) for i in range(NV)]
    Ts = []
    for Ap_mon, N0, Q0 in served:
        divN0 = {}; NgradD = {}; sNQ = {}
        for i in range(NV):
            divN0 = add(divN0, pderiv(N0[i], i)); NgradD = add(NgradD, mul(N0[i], D0i[i])); sNQ = add(sNQ, mul(N0[i], Qxi[i]))
        T = scale(mul(D0, sNQ), -1); T = add(T, mul(Qx, sub(mul(D0, divN0), NgradD))); T = sub(T, scale(mul(mul(F0, Qx), Q0), 2))
        Ts.append(core.reduce_mod(T, p))
    # linearization columns (mod p) -- same for all targets (depends on D0, N0?  NO: L depends on N0 through sumNF, divN0.
    # L(dD) = 2 A D0 dD - dD*sumNF + F0(dD divN0 - N0.grad dD) depends on (A, N0): handle by making the system block-wise:
    # unknowns: for each served target t a coefficient alpha_t; dD, dN, dQ shared?  The first-order equation for the
    # combination sum_t alpha_t (certificate_t) is linear in alpha with certificate-dependent L_t(dD) terms:
    #   sum_t alpha_t [2 A_t D0 dD - dD sumNF_t + F0 (dD divN0_t - N0_t.grad dD)] + [-D0 sum dN_i F0_i + F0(D0 div dN - dN.grad D0)] - F0^2 dQ + sum_t alpha_t T_t = 0
    # which is bilinear in (alpha, dD).  Restrict dD to the gauge-free colon slice is not needed; instead test each target
    # direction separately (alpha = e_t) and then the joint question with dD = 0 (certificate corrections only in N, Q).
    degN = d + 1; degQ = 2 * d - 6
    F0p = core.reduce_mod(F0, p); D0p = core.reduce_mod(D0, p); F0ip = [core.reduce_mod(x, p) for x in F0i]
    D0ip = [pderiv(D0p, i, p) for i in range(NV)]; F0sq = mul(F0p, F0p, p)
    NQcols = []
    for i in range(NV):
        for mm in w2lib.monos_char(degN, oplib.cadd(oplib.cadd(c, oplib.VARCHAR[i]), cA)):
            mp_ = {mm: 1}
            col = scale(mul(mul(D0p, mp_, p), F0ip[i], p), p - 1, p)
            col = add(col, mul(F0p, sub(mul(D0p, pderiv(mp_, i, p), p), mul(mp_, D0ip[i], p), p), p), p)
            NQcols.append(col)
    for mm in w2lib.monos_char(degQ, oplib.cadd(cA, oplib.cadd(c, c))):
        NQcols.append(scale(shift(F0sq, mm, 1, p), p - 1, p))
    def rank_of(cols):
        rows = {}
        for col in cols:
            for e in col: rows.setdefault(e, len(rows))
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, col in enumerate(cols):
            for e, v in col.items(): flat[rows[e] * C + j] = v
        return nmod_mat(R, C, flat, p).rank()
    # (i) dD = 0: which combinations of the T_t lie in the span of the (N,Q)-columns?  dim = (#served) - (rank([NQ|Ts]) - rank(NQ))
    rNQ = rank_of(NQcols); rAll = rank_of(NQcols + Ts)
    surv_noD = len(Ts) - (rAll - rNQ)
    # (ii) per target with dD free (full sector) -- as in w2_firstorder
    per = []
    for (Ap_mon, N0, Q0), T in zip(served, Ts):
        Ap = core.reduce_mod({Ap_mon: 1}, p); N0p = [core.reduce_mod(n, p) for n in N0]
        sumNF = {}; divN0p = {}
        for i in range(NV):
            sumNF = add(sumNF, mul(N0p[i], F0ip[i], p), p); divN0p = add(divN0p, pderiv(N0p[i], i, p), p)
        Dcols = []
        for mu in w2lib.monos_char(d, c):
            mup = {mu: 1}
            col = add(scale(mul(mul(Ap, D0p, p), mup, p), 2, p), scale(mul(mup, sumNF, p), p - 1, p), p)
            inner = mul(mup, divN0p, p)
            for i in range(NV): inner = sub(inner, mul(N0p[i], pderiv(mup, i, p), p), p)
            Dcols.append(add(col, mul(F0p, inner, p), p))
        rL = rank_of(Dcols + NQcols); rA = rank_of(Dcols + NQcols + [T])
        per.append((core.to_str({Ap_mon: Fraction(1)}), rA == rL))
    rec = {'d': d, 'dir': DIR, 'D0': core.to_str(D0)[:40], 'served_at_t0': len(served), 'survive_first_order_with_dD_free': [t for t, ok in per if ok],
           'dim_surviving_combinations_dD0': surv_noD, 'time': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True); out.append(rec)
w2lib.append_json(os.path.join(os.path.dirname(HERE), 'wave2-search.log.jsonl'), {'ts': w2lib.now(), 'step': 'first-order-served-survival', 'results': out})
