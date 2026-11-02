"""Structured-family E1 screen by Klein character.

For each degree d and each structured subspace U of S_d (monomial-defined shapes: cones with vertex a line/plane,
multiples of a coordinate, binary/ternary forms, polynomials in squares, u = Y2^2 polynomials, bounded X-degree, ...),
compute per character c:   adm(U,c) = {D in U_c : A D in J}  (mod p, good reduction checked),  and  adm cap J_d.
Survivors outside J_d are the structured candidates for fixed-D tests.
Usage: python3 structscreen.py <dmin> <dmax> <prime>
"""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
from hodge.core import NV, monomials, shift
from flint import nmod_mat

dmin, dmax, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
Ap = core.reduce_mod(L.A, p)


def shapes():
    out = []
    for i, j in itertools.combinations(range(6), 2):
        out.append((f'cone-line:indep-of-{V[i]},{V[j]}', lambda m, i=i, j=j: m[i] == 0 and m[j] == 0))
    for t in itertools.combinations(range(6), 3):
        out.append((f'cone-plane:indep-of-{",".join(V[k] for k in t)}', lambda m, t=t: all(m[k] == 0 for k in t)))
    for i in range(6):
        out.append((f'multiple-of-{V[i]}', lambda m, i=i: m[i] >= 1))
        out.append((f'multiple-of-{V[i]}^2', lambda m, i=i: m[i] >= 2))
    out.append(('poly-in-squares', lambda m: all(x % 2 == 0 for x in m)))
    out.append(('X-only', lambda m: m[3] + m[4] + m[5] == 0))
    out.append(('Y-only', lambda m: m[0] + m[1] + m[2] == 0))
    out.append(('X-degree<=2', lambda m: m[0] + m[1] + m[2] <= 2))
    out.append(('Y-degree<=2', lambda m: m[3] + m[4] + m[5] <= 2))
    out.append(('X-degree<=4', lambda m: m[0] + m[1] + m[2] <= 4))
    out.append(('Y-degree<=4', lambda m: m[3] + m[4] + m[5] <= 4))
    out.append(('Y2-even-and-Y0Y1-even-exponents', lambda m: m[5] % 2 == 0 and m[3] % 2 == 0 and m[4] % 2 == 0))
    out.append(('in-C[X,Y0Y1,Y0^2,Y1^2,Y2^2] (rho,r-invariant ring)', lambda m: m[5] % 2 == 0 and (m[3] + m[4]) % 2 == 0))
    out.append(('multiple-of-X0*Y0', lambda m: m[0] >= 1 and m[3] >= 1))
    out.append(('multiple-of-X0^3', lambda m: m[0] >= 3))
    out.append(('multiple-of-Y0^3', lambda m: m[3] >= 3))
    out.append(('in-(X0,Y0)^3', lambda m: m[0] + m[3] >= 3))
    out.append(('in-(X0,X1,Y0,Y1)^4', lambda m: m[0] + m[1] + m[3] + m[4] >= 4))
    return out


for d in range(dmin, dmax + 1):
    RT = {}
    Rd = {}
    for c in L.CHARS:
        RT[c] = L.CharReducer(d + 6, L.cadd(c, L.CHI_A), p)
        Rd[c] = L.CharReducer(d, c, p, need_W=True)
    for name, pred in shapes():
        t0 = time.time()
        res = {}
        total_adm, total_nonJ = 0, 0
        for c in L.CHARS:
            Um = [m for m in L.mons_chi(d, c) if pred(m)]
            if not Um:
                res[str(c)] = [0, 0, 0]; continue
            K = RT[c].cls_many([shift(Ap, m, 1, p) for m in Um])
            X, nullity = K.transpose().nullspace()
            kers = [{Um[j]: int(X[j, k]) for j in range(len(Um)) if int(X[j, k])} for k in range(nullity)]
            nonJ = 0
            if kers:
                KK = Rd[c].cls_many(kers)      # classes mod J_d^c
                nonJ = KK.rank()
            res[str(c)] = [len(Um), nullity, nonJ]
            total_adm += nullity; total_nonJ += nonJ
            if nonJ and d >= 8:
                json.dump([L.sym_str(P, p) for P in kers], open(f'{L.OUT}/struct_d{d}_{name.split(":")[0].replace("/","-")}_chi{c[0]}{c[1]}.json', 'w'))
        L.receipt('structscreen', exp='structured-E1-screen', d=d, prime=p, shape=name,
                  per_char={k: {'dimU': v[0], 'e1_adm': v[1], 'e1_adm_mod_J': v[2]} for k, v in res.items()},
                  e1_adm_total=total_adm, nonJ_total=total_nonJ,
                  good_reduction=all(RT[c].good and Rd[c].good for c in L.CHARS), seconds=round(time.time() - t0, 1))
print('STRUCT DONE', flush=True)
