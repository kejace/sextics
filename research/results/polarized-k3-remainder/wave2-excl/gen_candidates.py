"""Structured irreducible candidate denominators in the UV frame (expressions use X=U, Y=V names)."""
import sys, json
sys.path.insert(0, '.')
from k3lib import *
f = ternary_f()
def U(i): return var(i)
def V(i): return var(3 + i)
def fU(P3): return {(e[0], e[1], e[2], 0, 0, 0): Fr(c) for e, c in P3.items()}
def fV(P3): return {(0, 0, 0, e[0], e[1], e[2]): Fr(c) for e, c in P3.items()}
fi = [tderiv(f, i) for i in range(3)]
fij = [[tderiv(fi[i], j) for j in range(3)] for i in range(3)]
def pw(P, k): return power(P, k)
def ssum(polys):
    out = {}
    for P in polys: out = add(out, P)
    return out
def det3(M):
    import itertools
    out = {}
    for perm in itertools.permutations(range(3)):
        sgn = 1
        for a in range(3):
            for b in range(a + 1, 3):
                if perm[a] > perm[b]: sgn = -sgn
        term = {(0, 0, 0): Fr(sgn)}
        for i in range(3): term = tmul(term, M[i][perm[i]])
        for e, c in term.items(): out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}
hess = det3(fij)                      # degree 12 ternary
M01 = {}
for e, c in tmul(fij[0][0], fij[1][1]).items(): M01[e] = M01.get(e, 0) + c
for e, c in tmul(fij[0][1], fij[0][1]).items(): M01[e] = M01.get(e, 0) - c
M01 = {e: c for e, c in M01.items() if c}   # degree 8 Hessian minor
C = []
def addc(label, family, P, d):
    assert core.degree(P) == d, (label, core.degree(P))
    C.append(dict(label=label, family=family, d=d, expr=core.to_str(P), frame='uv'))
# degree 8
addc('hess-pairing (4,4)', 'd8-eigen-sym', ssum([mul(fU(fij[i][j]), fV(fij[i][j])) for i in range(3) for j in range(3)]), 8)
addc('sum U_i^4 V_i^4', 'd8-eigen-sym', ssum([mul(pw(U(i), 4), pw(V(i), 4)) for i in range(3)]), 8)
addc('sum U_i^7 V_i + swap', 'd8-eigen-sym', ssum([add(mul(pw(U(i), 7), V(i)), mul(U(i), pw(V(i), 7))) for i in range(3)]), 8)
addc('sum f_i(U) U_i^2 V_i + swap', 'd8-eigen-sym', ssum([add(mul(mul(fU(fi[i]), pw(U(i), 2)), V(i)), mul(mul(fV(fi[i]), pw(V(i), 2)), U(i))) for i in range(3)]), 8)
addc('sum U^8 + sum V^8', 'd8-sym', ssum([pw(U(i), 8) for i in range(3)] + [pw(V(i), 8) for i in range(3)]), 8)
addc('sum U^8 - sum V^8', 'd8-anti', sub(ssum([pw(U(i), 8) for i in range(3)]), ssum([pw(V(i), 8) for i in range(3)])), 8)
addc('sum f_i(U) V_i^3 + swap', 'd8-sym', ssum([add(mul(fU(fi[i]), pw(V(i), 3)), mul(fV(fi[i]), pw(U(i), 3))) for i in range(3)]), 8)
addc('f(U) U0^2 + f(V) V0^2', 'd8-sym-threadfamily', add(mul(fU(f), pw(U(0), 2)), mul(fV(f), pw(V(0), 2))), 8)
addc('sum U^8 (cone)', 'd8-eigen-cone', ssum([pw(U(i), 8) for i in range(3)]), 8)
addc('M01(f)(U) (Hessian-minor cone)', 'd8-eigen-cone', fU(M01), 8)
addc('M01(f)(U) + M01(f)(V)', 'd8-sym', add(fU(M01), fV(M01)), 8)
# degree 9
addc('sum U^9 (cone)', 'd9-eigen-cone', ssum([pw(U(i), 9) for i in range(3)]), 9)
addc('sum f_i(U) V_i^4', 'd9-eigen', ssum([mul(fU(fi[i]), pw(V(i), 4)) for i in range(3)]), 9)
addc('sum f_i(U) U_i^3 V_i + swap-free', 'd9-eigen', ssum([mul(mul(fU(fi[i]), pw(U(i), 3)), V(i)) for i in range(3)]), 9)
# degree 10
addc('polar pairing sum f_i(U) f_i(V)', 'd10-eigen-sym', ssum([mul(fU(fi[i]), fV(fi[i])) for i in range(3)]), 10)
addc('sum f_i(U) f_{i+1}(V) + swap', 'd10-eigen-sym', ssum([add(mul(fU(fi[i]), fV(fi[(i + 1) % 3])), mul(fV(fi[i]), fU(fi[(i + 1) % 3]))) for i in range(3)]), 10)
addc('sum U_i^5 V_i^5', 'd10-eigen-sym', ssum([mul(pw(U(i), 5), pw(V(i), 5)) for i in range(3)]), 10)
addc('sum U^10 + sum V^10', 'd10-sym', ssum([pw(U(i), 10) for i in range(3)] + [pw(V(i), 10) for i in range(3)]), 10)
addc('sum f_i(U) U_i^2 V_i^3 + swap', 'd10-sym', ssum([add(mul(mul(fU(fi[i]), pw(U(i), 2)), pw(V(i), 3)), mul(mul(fV(fi[i]), pw(V(i), 2)), pw(U(i), 3))) for i in range(3)]), 10)
# degree 12
addc('Hess(f)(U) (cone)', 'd12-eigen-cone', fU(hess), 12)
addc('Hess(f)(U) + Hess(f)(V)', 'd12-eigen-sym', add(fU(hess), fV(hess)), 12)
addc('Hess(f)(U) - Hess(f)(V)', 'd12-eigen-anti', sub(fU(hess), fV(hess)), 12)
addc('sum U_i^6 V_i^6', 'd12-eigen-sym', ssum([mul(pw(U(i), 6), pw(V(i), 6)) for i in range(3)]), 12)
addc('sum_ij f_ij(U) f_ij(V) (U.V)^2', 'd12-eigen-sym-reducible', mul(ssum([mul(fU(fij[i][j]), fV(fij[i][j])) for i in range(3) for j in range(3)]), pw(ssum([mul(U(i), V(i)) for i in range(3)]), 2)), 12)
with open('candidates_wave2.jsonl', 'w') as fh:
    for c in C: fh.write(json.dumps(c) + '\n')
for c in C: print(c['d'], c['family'], '|', c['label'], '| terms', len(core.parse(c['expr'])))
