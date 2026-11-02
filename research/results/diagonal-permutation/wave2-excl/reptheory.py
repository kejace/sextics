"""Exact (sympy/Fraction) re-check of the elementary structural claims of the dossier §1 and attack1 §1."""
import itertools, time
from fractions import Fraction as Fr
from common import *
F, A, q, rec = core.load_instance(IID)
t0 = time.time()
res = {}
# (1) sparse presentation: F = Fermat + (5/7) O1 + (5/11) O2 + (5/13) O3
def mono(e): return {tuple(e): Fr(1)}
X = lambda i: [1 if j == i else 0 for j in range(6)]
def add_e(*es): return tuple(sum(c) for c in zip(*es))
perms = list(itertools.permutations(range(3)))
O1 = {}; O2 = {}; O3 = {}
for (i, j, k) in perms:
    # O1: X_i^3 X_j Y_i Y_j ; O2: X_i^3 X_j Y_k^2 ; O3: Y_i^3 Y_j X_k^2
    e1 = [0]*6; e1[i] += 3; e1[j] += 1; e1[3+i] += 1; e1[3+j] += 1
    e2 = [0]*6; e2[i] += 3; e2[j] += 1; e2[3+k] += 2
    e3 = [0]*6; e3[3+i] += 3; e3[3+j] += 1; e3[k] += 2
    for O, e in ((O1, e1), (O2, e2), (O3, e3)):
        O[tuple(e)] = O.get(tuple(e), 0) + 1
Fermat = {tuple(6 if j == i else 0 for j in range(6)): Fr(1) for i in range(6)}
Fre = dict(Fermat)
for O, c in ((O1, Fr(5, 7)), (O2, Fr(5, 11)), (O3, Fr(5, 13))):
    for e, v in O.items(): Fre[e] = Fre.get(e, 0) + c * v
res['F_equals_Fermat_plus_orbits'] = all(Fr(F.get(e, 0)) == Fre.get(e, 0) for e in set(F) | set(Fre))
# (2) F = Fermat mod 5 (all coupling coefficients have 5-adic valuation >= 1)
res['F_mod5_is_Fermat'] = all((Fr(c).numerator % 5 == 0) for e, c in F.items() if e not in Fermat) and all(Fr(F[e]) == 1 for e in Fermat)
# (3) signed permutation symmetries of F (up to scalar): enumerate all 720*64
def act(P, perm, signs):
    out = {}
    for e, c in P.items():
        ne = [0]*6
        s = 1
        for i in range(6):
            ne[perm[i]] = e[i]
            if signs[i] == -1 and e[i] % 2: s = -s
        out[tuple(ne)] = out.get(tuple(ne), 0) + s * c
    return out
def proportional(P, Q):
    if set(P) != set(Q): return None
    k = None
    for e in P:
        r = Fr(Q[e]) / Fr(P[e])
        if k is None: k = r
        elif r != k: return None
    return k
syms = []
for perm in itertools.permutations(range(6)):
    for signs in itertools.product((1, -1), repeat=6):
        if signs[0] == -1: continue  # projective: normalize first sign
        k = proportional(F, act(F, perm, signs))
        if k is not None: syms.append((perm, signs, k))
res['num_signed_perm_symmetries_projective'] = len(syms)
res['signed_perm_symmetries'] = [(p, s, str(k)) for p, s, k in syms]
# (4) A and its S3 translates (simultaneous permutations), det = +1
def simperm(sig): return tuple([sig[i] for i in range(3)] + [3 + sig[i] for i in range(3)])
def det_perm(perm):
    # sign of the permutation
    s = 1; seen = [False]*6
    for i in range(6):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = perm[j]; L += 1
            if L % 2 == 0: s = -s
    return s
res['S3_simultaneous_dets'] = {str(sig): det_perm(simperm(sig)) for sig in perms}
T = [{tuple([3 if j == i else 0 for j in range(3)] + [3 if j == i else 0 for j in range(3)]): Fr(1)} for i in range(3)]
Aimgs = {str(sig): core.to_str(act(A, simperm(sig), (1,)*6)) for sig in perms}
res['A_translates'] = Aimgs
avg = {}
for sig in perms:
    for e, c in act(A, simperm(sig), (1,)*6).items(): avg[e] = avg.get(e, 0) + c
res['S3_average_of_A_is_zero'] = all(v == 0 for v in avg.values())
# (5) e_std identity in Q[S3]: (2/3) sum_s (1+s)/2 - 2 e_S3 == (1/3)(2e - c - c^2)
def compose(a, b): return tuple(a[b[i]] for i in range(3))
e = (0, 1, 2); c = (1, 2, 0); c2 = compose(c, c)
transp = [(1, 0, 2), (2, 1, 0), (0, 2, 1)]
lhs = {}
for s in transp:
    lhs[e] = lhs.get(e, 0) + Fr(2, 3) * Fr(1, 2); lhs[s] = lhs.get(s, 0) + Fr(2, 3) * Fr(1, 2)
for g in perms: lhs[g] = lhs.get(g, 0) - Fr(2, 6)
rhs = {e: Fr(2, 3), c: Fr(-1, 3), c2: Fr(-1, 3)}
res['e_std_identity'] = all(lhs.get(g, 0) == rhs.get(g, 0) for g in perms)
# (6) B_2 = T0+T1-2T2 fixed by (01); A = (B1 - B0)/3
B = [None]*3
for i in range(3):
    Bi = {}
    for j in range(3):
        for e_, v in T[j].items(): Bi[e_] = Bi.get(e_, 0) + v * (1 - 3 * (j == i))
    B[i] = Bi
comb = {}
for e_, v in B[1].items(): comb[e_] = comb.get(e_, 0) + v / 3
for e_, v in B[0].items(): comb[e_] = comb.get(e_, 0) - v / 3
res['A_eq_(B1-B0)/3'] = all(Fr(A.get(e_, 0)) == comb.get(e_, 0) for e_ in set(A) | set(comb))
res['B2_fixed_by_(01)'] = act(B[2], simperm((1, 0, 2)), (1,)*6) == B[2]
# (7) tau character of A: odd
res['A_Y_parity'] = core.parity_of(A)
# (8) Parity-split Koszul Hilbert function of R_F (dossier table), from ((1-t^5)/(1-t))^3 ((1-(ty)^5)/(1-ty))^3
import sympy as sp
t, y = sp.symbols('t y')
H = sp.expand(((1 - t**5) / (1 - t))**3 * ((1 - (t*y)**5) / (1 - t*y))**3)
H = sp.Poly(sp.cancel(H), t, y)
hp = {}; hm = {}
for (i, j), cval in H.terms():
    (hm if j % 2 else hp)[i] = (hm if j % 2 else hp).get(i, 0) + int(cval)
res['R_plus'] = [hp.get(i, 0) for i in range(25)]
res['R_minus'] = [hm.get(i, 0) for i in range(25)]
res['seconds'] = round(time.time() - t0, 1)
receipt('reptheory.json', res)
log({'exp': 'reptheory-exact-checks', **{k: v for k, v in res.items() if k not in ('signed_perm_symmetries', 'A_translates')}})
