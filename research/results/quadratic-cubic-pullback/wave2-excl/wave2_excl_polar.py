"""Structured Jacobian-type (polar) and other structured denominators: certified E1 + served-space tests.

usage: python3 wave2_excl_polar.py <jobname> <workers> <primeindex> [dmax] [only=name1,name2]
Output: wave2-excl.<jobname>.jsonl
"""
import sys, os, json, time, traceback
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, Qs, var, polar, gradpolar, euler, mul, add, scale, power, logline
from fractions import Fraction
import multiprocessing as mp

JOB = sys.argv[1]
WORKERS = int(sys.argv[2])
PI = int(sys.argv[3])
DMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 10
ONLY = None
for a in sys.argv[5:]:
    if a.startswith('only='):
        ONLY = set(a[5:].split(','))
P = core.PRIMES[PI]
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')

X0, X1, X2, Y0, Y1, Y2 = [var(i) for i in range(6)]
Q0, Q1, Q2, Q3, Q4, Q5 = Qs


def cands():
    C = []
    def addc(name, fam, D, note=''):
        C.append((name, fam, D, note))
    # ---- degree 6: J_6 (linear vector fields)
    addc('EX', 'polar6-blockEuler', euler([1, 1, 1, 0, 0, 0]), 'E_X F; Fermat limit is the join X0^6+X1^6+X2^6 (serves 100 there)')
    addc('E01', 'polar6-blockEuler', euler([1, 1, 0, 0, 0, 0]))
    addc('E03', 'polar6-blockEuler', euler([1, 0, 0, 1, 0, 0]))
    addc('E0123', 'polar6-blockEuler', euler([1, 1, 1, 1, 0, 0]))
    addc('E012-3', 'polar6-blockEuler', euler([1, 1, 1, -1, 0, 0]))
    addc('Elam123', 'polar6-weightedEuler', euler([1, 2, 3, 0, 0, 0]))
    addc('Elam123456', 'polar6-weightedEuler', euler([1, 2, 3, 4, 5, 6]))
    addc('Elam1m12m3', 'polar6-weightedEuler', euler([1, -1, 0, 2, 0, -3]))
    addc('Elam10m1010', 'polar6-weightedEuler', euler([1, 0, -1, 0, 1, 0]))
    addc('Erho', 'polar6-shapeShift', polar([X1, X2, X0, Y1, Y2, Y0]), 'cyclic shape symmetry rho')
    addc('ErhoInv', 'polar6-shapeShift', polar([X2, X0, X1, Y2, Y0, Y1]))
    addc('Esigma', 'polar6-shapeSwap-odd', polar([Y0, Y1, Y2, X0, X1, X2]), 'Y-odd: exchange shape symmetry sigma')
    addc('EsigmaSymp', 'polar6-shapeSwap-odd', polar([Y0, Y1, Y2, scale(X0, -1), scale(X1, -1), scale(X2, -1)]))
    addc('Eswap01', 'polar6-linear', polar([X1, X0, {}, Y1, Y0, {}]))
    addc('Enilp', 'polar6-linear', polar([X1, X2, {}, Y1, Y2, {}]))
    addc('ElistD1', 'polar6-linear', polar([add(X0, scale(X1, 2)), core.sub(X2, X0), {}, add(Y0, Y1), core.sub(scale(Y2, 3), Y0), {}]))
    addc('EXY-odd', 'polar6-linear-odd', polar([Y0, {}, {}, X0, {}, {}]), 'Y0 F_X0 + X0 F_Y0 (odd)')
    for k in range(6):
        addc(f'gradQ{k}', 'polar6-gradQ', gradpolar(Qs[k]), 'grad Q_k . grad F')
    addc('gradQ0pQ3', 'polar6-gradQ', gradpolar(add(Q0, Q3)))
    addc('gradQ0mQ3', 'polar6-gradQ', gradpolar(core.sub(Q0, Q3)))
    addc('gradSumQ', 'polar6-gradQ', gradpolar(add(add(add(Q0, Q1), add(Q2, Q3)), add(Q4, Q5))))
    addc('gradX0sq', 'polar6-gradmono', gradpolar(mul(X0, X0)), 'grad(X0^2).grad F = 2 X0 F_0 (reducible control)')
    addc('gradX0Y0-odd', 'polar6-gradmono-odd', gradpolar(mul(X0, Y0)))
    addc('gradX0X1', 'polar6-gradmono', gradpolar(mul(X0, X1)))
    # ---- degree 7: J_7 (quadratic vector fields)
    addc('gradX0Q0', 'polar7-gradh', gradpolar(mul(X0, Q0)))
    addc('gradX0Q3', 'polar7-gradh', gradpolar(mul(X0, Q3)))
    addc('gradX1Q0', 'polar7-gradh', gradpolar(mul(X1, Q0)))
    addc('gradY0Q0-odd', 'polar7-gradh-odd', gradpolar(mul(Y0, Q0)))
    addc('gradX0cube', 'polar7-gradh', gradpolar(power(X0, 3)))
    addc('gradX0Y0Y1', 'polar7-gradh', gradpolar(mul(X0, mul(Y0, Y1))))
    addc('wQX', 'polar7-Qfield', polar([Q0, Q1, Q2, {}, {}, {}]), 'sum_{i<3} Q_i F_i')
    addc('wQY-odd', 'polar7-Qfield-odd', polar([{}, {}, {}, Q3, Q4, Q5]))
    addc('wQshift', 'polar7-Qfield', polar([Q1, Q2, Q0, {}, {}, {}]))
    addc('wQcross', 'polar7-Qfield', polar([Q3, Q4, Q5, {}, {}, {}]), 'sum_{i<3} Q_{i+3} F_i')
    addc('wQmix-odd', 'polar7-Qfield-odd', polar([{}, {}, {}, Q0, Q1, Q2]))
    # ---- degree 8: J_8 (cubic vector fields)
    addc('gradQ0Q1', 'polar8-gradh', gradpolar(mul(Q0, Q1)))
    addc('gradQ0Q3', 'polar8-gradh', gradpolar(mul(Q0, Q3)))
    addc('gradQ3Q4', 'polar8-gradh', gradpolar(mul(Q3, Q4)))
    addc('gradQ0Q5', 'polar8-gradh', gradpolar(mul(Q0, Q5)))
    addc('gradSumQsq', 'polar8-gradh', gradpolar(add(add(add(mul(Q0, Q0), mul(Q1, Q1)), add(mul(Q2, Q2), mul(Q3, Q3))), add(mul(Q4, Q4), mul(Q5, Q5)))), 'grad(sum Q_k^2).grad F')
    addc('gradX0sqY0sq', 'polar8-gradh', gradpolar(mul(mul(X0, X0), mul(Y0, Y0))))
    addc('gradX0quart', 'polar8-gradh', gradpolar(power(X0, 4)))
    addc('wQEuler', 'polar8-Qfield', polar([mul(Q0, X0), mul(Q1, X1), mul(Q2, X2), mul(Q3, Y0), mul(Q4, Y1), mul(Q5, Y2)]), 'sum Q_i z_i F_i')
    addc('wQEulerX', 'polar8-Qfield', polar([mul(Q0, X0), mul(Q1, X1), mul(Q2, X2), {}, {}, {}]))
    addc('wQcyc', 'polar8-Qfield', polar([mul(Q0, X1), mul(Q1, X2), mul(Q2, X0), mul(Q3, Y1), mul(Q4, Y2), mul(Q5, Y0)]))
    # ---- degree 9: J_9 (quartic vector fields)
    addc('gradQ0Q1Q2', 'polar9-gradh', gradpolar(mul(mul(Q0, Q1), Q2)))
    addc('gradQ0Q1Q3', 'polar9-gradh', gradpolar(mul(mul(Q0, Q1), Q3)))
    addc('gradQ3Q4Q5', 'polar9-gradh', gradpolar(mul(mul(Q3, Q4), Q5)))
    addc('gradX0cubeY0sq', 'polar9-gradh', gradpolar(mul(power(X0, 3), mul(Y0, Y0))))
    addc('wQsqX', 'polar9-Qfield', polar([mul(Q0, Q0), mul(Q1, Q1), mul(Q2, Q2), {}, {}, {}]), 'sum_{i<3} Q_i^2 F_i')
    addc('wQsqCross', 'polar9-Qfield', polar([mul(Q3, Q3), mul(Q4, Q4), mul(Q5, Q5), {}, {}, {}]))
    addc('wQsqY-odd', 'polar9-Qfield-odd', polar([{}, {}, {}, mul(Q3, Q3), mul(Q4, Q4), mul(Q5, Q5)]))
    # ---- degree 10: J_10 (quintic vector fields)
    addc('gradF', 'polar10-gradh', gradpolar(F), '|grad F|^2')
    addc('gradFX', 'polar10-gradh', polar([core.pderiv(F, 0), core.pderiv(F, 1), core.pderiv(F, 2), {}, {}, {}]), 'sum_{i<3} F_i^2')
    addc('gradJoin012', 'polar10-gradh', gradpolar(add(add(power(Q0, 3), power(Q1, 3)), power(Q2, 3))), 'grad(Q0^3+Q1^3+Q2^3).grad F')
    addc('gradQ0Q1Q2sq', 'polar10-gradh', gradpolar(mul(mul(Q0, Q1), mul(Q2, Q2))))
    addc('graddelta-odd', 'polar10-gradh-odd', gradpolar(w2x.delta_poly()), 'grad delta . grad F (Y-odd)')
    # ---- random members of J_d (generic polar at each degree), fixed seeds, small integer coefficients
    import random as _r
    parz = [0, 0, 0, 1, 1, 1]
    for dd in (6, 7, 8, 9):
        for par in (0, 1):
            for s in range(2 if dd <= 7 else 1):
                _r.seed(1000 * dd + 10 * par + s)
                w = []
                for i in range(6):
                    pw = (par - parz[i]) % 2
                    poly = {}
                    for m in core.monomials(dd - 5, pw):
                        c = _r.randint(-3, 3)
                        if c:
                            poly[m] = Fraction(c)
                    w.append(poly)
                addc(f'randJ{dd}_par{par}_{s}', f'polar{dd}-random' + ('-odd' if par else ''), polar(w), f'random element of J_{dd} (seed {1000*dd+10*par+s})')
    return C


def run(item):
    name, fam, D, note = item
    try:
        rec = w2x.full_test(name, fam, D, P, LOG, note=note)
    except Exception:
        rec = logline(LOG, dict(exp='fixedD', name=name, family=fam, prime=P, result='ERROR', error=traceback.format_exc()[-800:], job=JOB))
    print(json.dumps({k: v for k, v in rec.items() if k not in ('D_head',)}), flush=True)
    return rec


if __name__ == '__main__':
    items = [c for c in cands() if core.degree(c[2]) <= DMAX and (ONLY is None or c[0] in ONLY)]
    items.sort(key=lambda c: -core.degree(c[2]))
    print(f'{len(items)} candidates, prime {P}, workers {WORKERS}', flush=True)
    with mp.get_context('fork').Pool(WORKERS) as pool:
        for _ in pool.imap_unordered(run, items):
            pass
    print('ALL DONE', flush=True)
