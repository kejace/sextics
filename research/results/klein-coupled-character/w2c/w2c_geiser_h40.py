"""(4,0)-forms on the Geiser-twisted fibre product Z = X x_{W,gamma} X  (cyclic 14-cover w^14 = phi^5 x^3 y^9 of X),
decomposed by the characters (j,k) of <d,e> = (Z/14)^2 (d deck, e lifts sigma*tau with w fixed), with the sign of
r = (tau: Y -> -Y, w -> eps w).  Esnault-Viehweg count on X (caveat: branch divisor not assumed snc; the count is the
standard eigensheaf formula H^0(X, L^j(-floor(jD/14))) with L = O(3), D = 5[phi]+3[x]+9[y]).
Off-diagonal (j != k) r-invariant pieces give <r,swap>-invariant (4,0)-forms on Z, i.e. h^{4,0}(Z/<r,swap>) > 0.
"""
import itertools
from collections import defaultdict
def monomials(deg, n=6):
    if n == 1:
        yield (deg,); return
    for a in range(deg+1):
        for rest in monomials(deg-a, n-1):
            yield (a,)+rest
res = {}
for eps in (+1, -1):
    offdiag = defaultdict(int); diag = defaultdict(int)
    for j in range(14):
        a = (5*j)//14; b = (3*j)//14; c = (9*j)//14
        deg = 3*j - 6*a - b - c
        for m in monomials(deg):
            al, be, ga, d1, d2, d3 = m
            delta = d1+d2+d3
            E = (9*al + 11*be + 15*ga + 9*b + 11*c + 7) % 14
            k = (-E) % 14
            rsign = (-1)**(1+delta) * (eps**j)
            if deg == 6 and m == (5,1,0,0,0,0):   # remove the F-direction (same character as F's monomials: weight 0)
                pass
            if rsign == 1:
                if j != k: offdiag[(j,k)] += 1
                else: diag[j] += 1
        if deg == 6:
            # F spans one dimension in the weight-0 class of S_6; its (j,k) index:
            E = (9*b + 11*c + 7) % 14; k = (-E) % 14
            rsign = (-1)**(1+0) * (eps**j)
            if rsign == 1:
                if j != k: offdiag[(j,k)] -= 1
                else: diag[j] -= 1
    tot_off = sum(offdiag.values()); tot_diag = sum(diag.values())
    print('eps =', eps, ': r-invariant off-diagonal (j,k) dims:', dict(sorted(offdiag.items())), ' total', tot_off)
    print('          r-invariant diagonal (j,j) dims (swap action not resolved):', dict(sorted(diag.items())), ' total', tot_diag)
