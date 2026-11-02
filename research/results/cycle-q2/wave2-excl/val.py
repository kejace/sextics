"""E0: validation of xlib against known accepted witnesses and direct polynomial evaluation."""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
from fractions import Fraction
from hodge import core
import numpy as np

def witness_check(inst):
    F, A, q, rec = core.load_instance(inst)
    w = json.load(open(f'/home/kejace/devel/sextics/data/{inst}/witness.json'))
    D = core.parse(w['D']); N = [core.parse(s) for s in w['N']]; Q = core.parse(w['Q'])
    # scale to integers: F2 = 2F; common denominator L for D,N,Q
    import math
    L = 1
    for P in [D] + N + [Q]:
        for c in P.values():
            L = L * Fraction(c).denominator // math.gcd(L, Fraction(c).denominator)
    Di = pint({e: Fraction(c) * L for e, c in D.items()})
    # the identity is homogeneous of degree 2 in (D,N,Q) jointly?  A D^2, D N, N D, Q: scaling D,N,Q by L scales E by L^2 -> fine
    Ni = [pint({e: Fraction(c) * L for e, c in P.items()}) for P in N]
    Qi = pint({e: Fraction(c) * L * L for e, c in Q.items()})
    # careful: E = A D^2 - qD N(F) + F(D divN - N(D)) - F^{q+1}Q : D->LD, N->LN gives L^2 on all but Q; so Q->L^2 Q
    ellF = 1
    for c in F.values(): ellF = ellF * Fraction(c).denominator // math.gcd(ellF, Fraction(c).denominator)
    ellA = 1
    for c in A.values(): ellA = ellA * Fraction(c).denominator // math.gcd(ellA, Fraction(c).denominator)
    F2 = pint({e: Fraction(c) * ellF for e, c in F.items()}); Ai = pint({e: Fraction(c) * ellA for e, c in A.items()})
    # certificate for A  <=> certificate for ellA*A (scale N,Q by ellA): scale Ni, Qi
    Ni = [pscale(P, ellA) for P in Ni]; Qi = pscale(Qi, ellA)
    S = FixedDSystem(F2, Ai, q, Di, ell=ellF)
    E = S.evaluate(Ni, Qi)
    print(inst, 'd=', S.d, 'identity holds exactly (my arithmetic):', len(E) == 0, '| sector H size', S.Hsize, 'cols', S.ncols, 'rows', S.nrows)
    r = solve_fixed_D(F2, Ai, q, Di, ell=ellF)
    print(inst, 'solver status:', r['status'], 'ranks', r['ranks'], 'primes', r['primes'], 'secs', r['seconds'])
    xlib.log(dict(exp='E0-validation', instance=inst, identity_exact=len(E) == 0, **{k: v for k, v in r.items() if not k.startswith('_')}))
    return r

def matrix_selftest():
    F2, A, q = load_cycle_q2()
    D = pderiv(F2, 2)  # F2_{X2}
    S = FixedDSystem(F2, A, q, D)
    M, b = S.int_matrix()
    random.seed(1)
    x = [random.randint(-5, 5) for _ in range(S.ncols)]
    N = [dict() for _ in range(6)]; Q = {}
    for j, (kind, i, m) in enumerate(S.cols):
        if x[j] == 0: continue
        if kind == 'N': N[i][m] = x[j]
        else: Q[m] = x[j]
    Mx = {}
    for i in range(S.nrows):
        v = sum(int(M[i, j]) * x[j] for j in range(S.ncols))
        if v: Mx[S.Emonos[i]] = v
    E = S.evaluate(N, Q)
    # E = 2^{q+1} A D^2 + (M x)  =>  M x - (E - 2^{q+1} A D^2) == 0
    AD2 = pscale(pmul(A, pmul(D, D)), 2 ** (q + 1))
    diff = padd(Mx, padd(E, AD2, -1), -1)
    print('matrix selftest (M x == E - 2^{q+1} A D^2):', len(diff) == 0, 'cols', S.ncols, 'rows', S.nrows, 'H', S.Hsize)
    xlib.log(dict(exp='E0-matrix-selftest', ok=len(diff) == 0, cols=S.ncols, rows=S.nrows, H=int(S.Hsize)))

if __name__ == '__main__':
    matrix_selftest()
    for inst in ['block-q2', 'fermat-q2']:
        witness_check(inst)
