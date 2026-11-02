#!/usr/bin/env python3
"""referee_singletons_hessian.py -- FAST variant (exact kernel + exact Hessian class only; no Farkas) of referee_singletons.py: repair of the nine degree-6 singleton sectors (dim K = 1, dim J = 0), lane cycle-q2.

sectors.py certified a SYMMETRIC MOD-p LIFT of the mod-p E1-kernel vector, which is in general not the exact Q-kernel element.
Here, for each singleton sector: (1) exact E1 kernel over Q (sector_exact machinery, fmpq rref); if its dimension is 0 the
sector carries no E1-admissible denominator; (2) if it is 1, D_Q is computed exactly, its Hessian class [H(D_Q)] in R_18 is
computed exactly (nonzero => no certificate, by the relaxation), and (3) an exact Farkas certificate for the integralised D_Q
is produced with w2lib.farkas (saved under certs_regen/) for an independent replay with referee_farkas.py.
"""
import sys, json, time, math
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, sector_exact as se
from hess6 import key_of, poly_key
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials
from fractions import Fraction
from flint import fmpq_mat, fmpq
OUT = '/tmp/hodge/referee/cycle-q2/singletons_hessian.results.jsonl'
SECTORS = [(0,0,3),(0,1,3),(0,1,4),(0,2,0),(0,5,5),(3,0,0),(3,0,3),(3,3,3),(3,5,1)]
F, A, q, G = w2lib.load(); gens = G.generators
F2 = scale(F, 2); Fi = [pderiv(F2, i) for i in range(6)]
kA = poly_key(A, gens); d = 6

def solve_E1(AD, kAD):
    cols = []; colp = []
    for t in range(6):
        kk = key_of(next(iter(Fi[t])), gens)
        need = tuple((x - y) % 6 for x, y in zip(kAD, kk))
        for m in monomials(13):
            if key_of(m, gens) == need: cols.append((t, m)); colp.append(shift(Fi[t], m, 1))
    rows = {}
    for P in colp:
        for e in P: rows.setdefault(e, len(rows))
    for e in AD: rows.setdefault(e, len(rows))
    Rn, Cc = len(rows), len(cols)
    Mq = fmpq_mat(Rn, Cc + 1)
    for j, P in enumerate(colp):
        for e, v in P.items(): fr = Fraction(v); Mq[rows[e], j] = fmpq(fr.numerator, fr.denominator)
    for e, v in AD.items(): fr = Fraction(v); Mq[rows[e], Cc] = fmpq(fr.numerator, fr.denominator)
    Mqr, rk = Mq.rref()
    pv = []
    for i in range(rk):
        j = 0
        while j <= Cc and Mqr[i, j] == 0: j += 1
        pv.append(j)
    assert pv[-1] != Cc, 'E1 inconsistent over Q?!'
    N = [{} for _ in range(6)]
    for i, j in enumerate(pv):
        val = Mqr[i, Cc]
        if val != 0:
            t, m = cols[j]; N[t][m] = Fraction(int(val.p), int(val.q))
    return N

for kD in SECTORS:
    t0 = time.time()
    Dmons = [m for m in monomials(d) if key_of(m, gens) == kD]
    kAD = tuple((x + y) % 6 for x, y in zip(kA, kD))
    R18 = se.ExactReducer(Fi, 18, kAD, gens)
    imgs = [R18.normal_form(shift(A, m, 1)) for m in Dmons]
    std = R18.standard; sidx = {m: i for i, m in enumerate(std)}
    M = fmpq_mat(len(std), len(Dmons))
    for j, nf in enumerate(imgs):
        for e, v in nf.items(): M[sidx[e], j] = fmpq(v.numerator, v.denominator)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < len(Dmons) and Mr[i, j] == 0: j += 1
        piv.append(j)
    free = [j for j in range(len(Dmons)) if j not in set(piv)]
    rec = {'sector': list(kD), 'dimS': len(Dmons), 'dimK_Q': len(free)}
    if len(free) == 1:
        fj = free[0]; D = {Dmons[fj]: Fraction(1)}
        for i, pj in enumerate(piv):
            c = Mr[i, fj]
            if c != 0: D[Dmons[pj]] = -Fraction(int(c.p), int(c.q))
        Dint, den = w2lib.integral(D)
        rec['D_Q_integral'] = core.to_str(Dint); rec['D_maxbits'] = max(abs(v).bit_length() for v in Dint.values())
        # exact Hessian class
        N = solve_E1(mul(A, D), kAD)
        divN = {}
        for t in range(6): divN = add(divN, pderiv(N[t], t))
        H = mul(D, divN)
        for t in range(6): H = core.sub(H, mul(N[t], pderiv(D, t)))
        H = {e: Fraction(v) for e, v in H.items() if Fraction(v) != 0}
        if H:
            kH = poly_key(H, gens); RH = se.ExactReducer(Fi, 18, kH, gens); nf = RH.normal_form(H)
            rec['H_class_in_R18_nonzero'] = bool(nf); rec['H_normal_form_terms'] = len(nf)
        else:
            rec['H_class_in_R18_nonzero'] = False; rec['H_normal_form_terms'] = 0
    rec['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(rec), flush=True)
    with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
