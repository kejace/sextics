#!/usr/bin/env python3
"""Exact first-order deformation test along F_t = F0 + t*Qx at t = 0 for the Fermat solutions D0 = f_S * X0^(d-6)
(S a 3-subset of variables containing X0 but not Y0; the cone over the Fermat curve of S serves X0^3 Y0^3 on F0).
Solves the t=0 certificate exactly (sector-reduced), then asks (mod p) whether the t-derivative of the cleared identity
lies in the image of the linearization in (dD, dN, dQ) with dD free in the sector.  Light (seconds)."""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
import json, time, itertools
from fractions import Fraction
import w2lib, oplib
from hodge import core
from hodge.core import NV, add, sub, scale, mul, pderiv, shift, power
from flint import nmod_mat
F, A, q = w2lib.load()
DEGS = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [6, 7, 8, 9]
PRIMES = core.PRIMES[:2]
f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6'); F0 = add(f6, g6)
DIR = sys.argv[2] if len(sys.argv) > 2 else 'Qx'
Qx = {'Qx': sub(F, F0), 'chain': core.parse('X0**3*X1*Y0*Y1/2 + X1**3*X2*Y1*Y2/2'), 'extra': core.parse('X0**2*X1*X2*Y0*Y2/3'),
      'link1': core.parse('X0**3*X1*Y0*Y1/2'), 'link2': core.parse('X1**3*X2*Y1*Y2/2')}[DIR]
F0i = [pderiv(F0, i) for i in range(NV)]; Qxi = [pderiv(Qx, i) for i in range(NV)]
c = (0, 0); cA = (3, 0)
out = []
for d in DEGS:
    m = {tuple([d - 6, 0, 0, 0, 0, 0]): Fraction(1)}
    for S in itertools.combinations(range(NV), 3):
        if 0 not in S or 3 in S: continue
        fS = {tuple(6 if t == i else 0 for t in range(NV)): Fraction(1) for i in S}
        D0 = mul(fS, m)
        if oplib.char_of(D0) != c: continue
        t0 = time.time()
        st, data = oplib.lift_solution(F0, A, q, D0, max_primes=6)
        if st != 'solved':
            out.append({'d': d, 'S': [core.VARS[i] for i in S], 'status': f't0-{st}'}); print(out[-1]); continue
        N0, Q0 = data['N'], data['Q']
        # t-derivative term (exact, degree 2d+6): -D0 sum N0_i Qx_i + Qx (D0 div N0 - N0.grad D0) - 2 F0 Qx Q0
        divN0 = {}
        for i in range(NV): divN0 = add(divN0, pderiv(N0[i], i))
        NgradD = {}
        D0i = [pderiv(D0, i) for i in range(NV)]
        for i in range(NV): NgradD = add(NgradD, mul(N0[i], D0i[i]))
        sNQ = {}
        for i in range(NV): sNQ = add(sNQ, mul(N0[i], Qxi[i]))
        T = scale(mul(D0, sNQ), -1)
        T = add(T, mul(Qx, sub(mul(D0, divN0), NgradD)))
        T = sub(T, scale(mul(mul(F0, Qx), Q0), 2))
        ranks = []
        for p in PRIMES:
          degN = d + 1; degQ = 2 * d - 6
          F0p = core.reduce_mod(F0, p); Ap = core.reduce_mod(A, p); D0p = core.reduce_mod(D0, p)
          N0p = [core.reduce_mod(n, p) for n in N0]; Q0p = core.reduce_mod(Q0, p)
          F0ip = [core.reduce_mod(x, p) for x in F0i]; D0ip = [pderiv(D0p, i, p) for i in range(NV)]
          divN0p = core.reduce_mod(divN0, p)
          sumNF = {}
          for i in range(NV): sumNF = add(sumNF, mul(N0p[i], F0ip[i], p), p)
          F0D0 = mul(F0p, D0p, p); F0sq = mul(F0p, F0p, p)
          cols = []
          # dD = monomial mu: 2 A D0 mu - mu*sumNF + F0 (mu divN0 - sum_i N0_i d_i mu)
          for mu in w2lib.monos_char(d, c):
              mup = {mu: 1}
              col = add(scale(mul(mul(Ap, D0p, p), mup, p), 2, p), scale(mul(mup, sumNF, p), p - 1, p), p)
              inner = mul(mup, divN0p, p)
              for i in range(NV):
                  inner = sub(inner, mul(N0p[i], pderiv(mup, i, p), p), p)
              col = add(col, mul(F0p, inner, p), p)
              cols.append(col)
          # dN_i = monomial m: -D0 m F0_i + F0 (D0 d_i m - m D0_i)
          for i in range(NV):
              for mm in w2lib.monos_char(degN, oplib.cadd(oplib.cadd(c, oplib.VARCHAR[i]), cA)):
                  mp_ = {mm: 1}
                  col = scale(mul(mul(D0p, mp_, p), F0ip[i], p), p - 1, p)
                  col = add(col, mul(F0p, sub(mul(D0p, pderiv(mp_, i, p), p), mul(mp_, D0ip[i], p), p), p), p)
                  cols.append(col)
          # dQ = monomial: -F0^2 m
          for mm in w2lib.monos_char(degQ, oplib.cadd(cA, oplib.cadd(c, c))):
              cols.append(scale(shift(F0sq, mm, 1, p), p - 1, p))
          rhs = scale(core.reduce_mod(T, p), p - 1, p)   # L(delta) + T = 0  ->  L(delta) = -T
          rows = {}
          for col in cols:
              for e in col: rows.setdefault(e, len(rows))
          for e in rhs: rows.setdefault(e, len(rows))
          R, C = len(rows), len(cols)
          flat = [0] * (R * (C + 1))
          for j, col in enumerate(cols):
              for e, v in col.items(): flat[rows[e] * (C + 1) + j] = v
          for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
          M = nmod_mat(R, C + 1, flat, p); ra = M.rank()
          M0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p); r0 = M0.rank()
          ranks.append((p, r0, ra))
        rec = {'dir': DIR, 'd': d, 'S': [core.VARS[i] for i in S], 'D0': core.to_str(D0)[:60], 'ranks_per_prime': ranks, 'rows': R, 'cols': C,
               'first_order_deformable': bool(any(ra == r0 for _, r0, ra in ranks)), 'T_terms': len(T), 'time': round(time.time() - t0, 1)}
        out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(f'wave2-firstorder-{DIR}-{"_".join(map(str, DEGS))}.json', 'w'), indent=1)
w2lib.append_json(w2lib.REMOTE_LOG,
                  {'ts': w2lib.now(), 'step': 'first-order-deformation', 'dir': DIR, 'results': [(r['d'], ''.join(r['S']), r.get('first_order_deformable')) for r in out]})
