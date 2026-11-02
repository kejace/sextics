"""More cheap independent checks: point matroid of the 8 power points; membership of structured sextics in J_6
with a new prime (rigorous non-membership via Hilbert rank match); D=1 test; Fermat-P^3 lines count (numerical)."""
import sys, json, itertools, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat
import numpy as np
X0, X1, X2, Y0, Y1, Y2 = VARS
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs); A = poly_from_sympy(As)
L = X0 + 2*X1 + 3*X2; M = Y0 + 3*Y1 + 5*Y2
pts = sp.Matrix([[1 if k == i else 0 for k in range(6)] for i in range(6)] + [[1, 2, 3, 1, 3, 5], [1, 2, 3, -1, -3, -5]])
# every 4 independent? minimal dependent sets (circuits)
dep4 = [S for S in itertools.combinations(range(8), 4) if pts[list(S), :].rank() < 4]
circ5 = [S for S in itertools.combinations(range(8), 5) if pts[list(S), :].rank() < 5]
print('dependent 4-sets:', dep4, ' dependent 5-sets:', circ5)
logjson(claim='8 power points: every 4 independent; dependent 5-sets exactly {0,1,2,6,7},{3,4,5,6,7}', verdict='CONFIRMED' if (not dep4 and circ5 == [(0,1,2,6,7),(3,4,5,6,7)]) else 'REFUTED', dep4=str(dep4), circ5=str(circ5))
# membership tests in J_6^{par} with p = 1000000007 (rank match => rigorous)
p = 1000000007
Fp = reduce_mod(F, p); H = hilb()
def in_J(P, deg, par):
    Pp = reduce_mod(P, p)
    J = jac_gens(Fp, deg, par, p)
    rows = {}
    for Q in J + [Pp]:
        for e in Q: rows.setdefault(e, len(rows))
    rJ = nmod_from_polys(J, p, rows)[0].rank(); rJP = nmod_from_polys(J + [Pp], p, rows)[0].rank()
    exp = H[deg]['J_plus' if par == 0 else 'J_minus']
    return rJP == rJ, rJ, exp
B = poly_from_sympy((L+M)**6 + (L-M)**6)
Bodd = poly_from_sympy(sp.expand(((L+M)**6 - (L-M)**6)/2))
eXY = poly_from_sympy(sum(v**6 for v in VARS))
res = {}
res['B_in_J6plus'] = in_J(B, 6, 0); res['eXeY_in_J6plus'] = in_J(eXY, 6, 0); res['Bodd_in_J6minus'] = in_J(Bodd, 6, 1)
res['A_in_J6minus'] = in_J(A, 6, 1)
print('B=(L+M)^6+(L-M)^6 in J_6^+ ? (member, rankJ, expected):', res['B_in_J6plus'])
print('e(X)+e(Y) in J_6^+ ?:', res['eXeY_in_J6plus'])
print('odd part of (L+M)^6 in J_6^- ?:', res['Bodd_in_J6minus'])
print('A in J_6^- ? (D=1 E1):', res['A_in_J6minus'])
# (L+M)^6 + sum_{k in I} Z_k^6 for all subsets I: even part in J_6^+ ?
hits = []
Bhalf = scale(B, Fraction(1, 2))
for r in range(0, 7):
    for I in itertools.combinations(range(6), r):
        P = dict(Bhalf)
        for k in I:
            m = [0]*6; m[k] = 6; P = add(P, {tuple(m): Fraction(1)})
        for cc in (Fraction(5, 7**7), Fraction(1)):   # with the true coefficient c and with 1
            Pc = add(scale(Bhalf, cc), add(P, scale(Bhalf, -1)))
            ok, rJ, exp = in_J(Pc, 6, 0)
            if ok: hits.append((I, str(cc)))
print('subsets I with  c*evenpart((L+M)^6) + sum_I Z_k^6 in J_6^+  (c in {5/7^7, 1}):', hits)
logjson(claim='(L+-M)^6 and (L+-M)^6 + coordinate sixth powers (all 64 subsets, scale c and 1) excluded as denominators: even part not in J_6^+', verdict='CONFIRMED' if not hits else 'REFUTED', hits=str(hits), rankJ6plus=res['B_in_J6plus'][1], expected=res['B_in_J6plus'][2])
for k, v in res.items(): logjson(claim='membership:' + k, verdict='DATA', member=v[0], rankJ=v[1], expected=v[2], p=p)
# quartic / sextic power independence, exact
forms = [v for v in VARS] + [L + M, L - M]
for dd in (4, 5, 6):
    mons = monomials(dd)
    Mx = sp.Matrix([[sp.Poly(f**dd, *VARS).coeff_monomial(sp.prod([v**e for v, e in zip(VARS, m)])) for m in mons] for f in forms])
    print(f'rank of the 8 forms l_k^{dd}:', Mx.rank())
    logjson(claim=f'rank of l_k^{dd} (8 powers) exact over Q', verdict='DATA', rank=int(Mx.rank()))
# lines: Fermat P^3's in P^7 meeting Lambda; numerical with singular values
k = (5.0 / 7**7) ** (1/6)
Lv = np.array([1,2,3,0,0,0], float); Mv = np.array([0,0,0,1,3,5], float)
fv = [np.eye(6)[i] for i in range(6)] + [k*(Lv+Mv), k*(Lv-Mv)]
zetas = [np.exp(1j*np.pi*(2*j+1)/6) for j in range(6)]
def pairings(items):
    if not items: yield (); return
    a = items[0]
    for b in items[1:]:
        for r in pairings([i for i in items if i not in (a, b)]): yield ((a, b),) + r
t0 = time.time(); minsv = 1e9; cnt = 0; lines = []
for pr in pairings(list(range(8))):
    for zs in itertools.product(range(6), repeat=4):
        Mx = np.array([fv[a] - zetas[z]*fv[b] for (a, b), z in zip(pr, zs)])
        sv = np.linalg.svd(Mx, compute_uv=False)
        minsv = min(minsv, sv[3]); cnt += 1
        # line = 2-dim null space; store an orthonormal basis for dedupe
        _, _, Vh = np.linalg.svd(Mx); lines.append(Vh[4:6].conj().T)
print(f'{cnt} Fermat P^3s; min 4th singular value = {minsv:.3e} ({time.time()-t0:.0f}s)')
# dedupe lines: projector P = Q Q^*; hash rounded
proj = {}
for Qb in lines:
    P = Qb @ Qb.conj().T
    key = tuple(np.round(P.real, 6).ravel()) + tuple(np.round(P.imag, 6).ravel())
    proj[key] = proj.get(key, 0) + 1
print('distinct lines (projector dedupe, 6 decimals):', len(proj))
logjson(claim='136080 Fermat P^3 cut Lambda in exactly a line each; 136066 distinct lines', verdict='CONFIRMED' if (minsv > 1e-6 and len(proj) == 136066) else 'WEAKENED', min_sv4=minsv, distinct=len(proj), total=cnt)
