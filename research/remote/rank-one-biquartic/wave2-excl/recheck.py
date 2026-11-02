"""Adversarial re-verification of the structural / elementary claims of the wave-1 report and dossier
(new prime p = 1000000007 where modular, exact rational arithmetic where possible)."""
import sys, json, time, itertools
sys.path.insert(0, '.')
import w2
from hodge import core, symmetry
from fractions import Fraction
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
R, Fp, Ap, Fi, F, A = w2.load(p)
res = {}
P = core.parse
f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
L = P('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2'); M = P('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
t = Fraction(5, 7)
L2M = core.mul(core.mul(L, L), M)
res['R1_decomposition_F=f+g+tL2M'] = core.sub(F, core.add(core.add(f, g), core.scale(L2M, t))) == {}
# M = L(Y); g = f(Y2,Y0,Y1) + t Y0^3 Y1 Y2^2
def subst_perm(Pd, perm):  # perm maps old var index -> new var index
    out = {}
    for e, c in Pd.items():
        e2 = [0] * 6
        for i, x in enumerate(e):
            e2[perm[i]] += x
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}
res['R1_M_equals_L_of_Y'] = subst_perm(L, [3, 4, 5, 0, 1, 2]) == M
res['R1_g_relation'] = core.sub(g, core.add(subst_perm(f, [5, 3, 4, 0, 1, 2]), core.scale(P('Y0**3*Y1*Y2**2'), t))) == {}
# R2 smoothness via mod 5: all non-pure coefficients are multiples of 5/7
pure = {tuple(6 if j == i else 0 for j in range(6)) for i in range(6)}
res['R2_nonpure_coeffs_in_5Z/7'] = all((c * 7) % 5 == 0 and (c * 7).denominator == 1 for e, c in F.items() if e not in pure) and all(F[e] == 1 for e in pure)
# R3 symmetry: diagonal group (engine) and permutation screen
G = symmetry.diagonal_group(F)
res['R3_diagonal_group_order'] = G.order
perms_ok = []
supp = set(F)
for perm in itertools.permutations(range(6)):
    if all(tuple(e[perm.index(j)] for j in range(6)) in supp for e in supp):
        perms_ok.append(perm)
res['R3_support_preserving_permutations'] = perms_ok
# R7 Euler residual: B_{R_X}(A X_i e_i) = 8 t A L^2 M exactly
def eps(Pp, idxs):
    out = {}
    for i in idxs:
        out = core.add(out, core.mul({tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}, core.pderiv(Pp, i)))
    return out
RX = eps(F, [0, 1, 2]); RY = eps(F, [3, 4, 5])
res['R7_RX=6f+4tL2M'] = core.sub(RX, core.add(core.scale(f, 6), core.scale(L2M, 4 * t))) == {}
res['R7_RY=6g+2tL2M'] = core.sub(RY, core.add(core.scale(g, 6), core.scale(L2M, 2 * t))) == {}
def B_exact(D, N):
    divN = {}
    ND = {}
    for i in range(6):
        divN = core.add(divN, core.pderiv(N[i], i))
        ND = core.add(ND, core.mul(N[i], core.pderiv(D, i)))
    return core.sub(core.mul(D, divN), ND)
NX = [core.mul(A, {tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}) if i < 3 else {} for i in range(6)]
res['R7_E1_RX'] = core.sub(sum((core.mul(NX[i], core.pderiv(F, i)) for i in range(6)), {}) if False else core.add(core.add(core.mul(NX[0], core.pderiv(F, 0)), core.mul(NX[1], core.pderiv(F, 1))), core.mul(NX[2], core.pderiv(F, 2))), core.mul(A, RX)) == {}
res['R7_B_RX(AX)=8tAL2M'] = core.sub(B_exact(RX, NX), core.scale(core.mul(A, L2M), 8 * t)) == {}
NY = [{} if i < 3 else core.mul(A, {tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}) for i in range(6)]
res['R7_B_RY(AY)=8tAL2M'] = core.sub(B_exact(RY, NY), core.scale(core.mul(A, L2M), 8 * t)) == {}
# R8 Prop 4.1: F_i - f_i = 2t L L_i M, F_{3+j} - g_j = t L^2 M_j (so J(F)+(L) = J(f+g)+(L))
ok = True
for i in range(3):
    ok &= core.sub(core.sub(core.pderiv(F, i), core.pderiv(f, i)), core.scale(core.mul(core.mul(L, core.pderiv(L, i)), M), 2 * t)) == {}
for j in range(3):
    ok &= core.sub(core.sub(core.pderiv(F, 3 + j), core.pderiv(g, 3 + j)), core.scale(core.mul(core.mul(L, L), core.pderiv(M, 3 + j)), t)) == {}
res['R8_Fi-F0i_in_(L)'] = ok
print(json.dumps(res, default=str, indent=1), flush=True)
w2.log(dict(exp='recheck-elementary', p=p, result=res))
# --- modular structure claims on K_7 (p new) ---
t0 = time.time()
fam7 = w2.e1_family(R, Fp, Ap, Fi, 7, 0, verbose=False)
rJ7, _ = w2.jacobian_dim(R, Fi, 7, 0)
res7 = dict(dimK7=fam7['dimK'], dimJ7=rJ7)
# non-Jacobian directions: complement of J_7 in K_7 (greedy)
from flint import nmod_mat
idx7 = w2.Index(7, 0)
gens7 = []
for i in range(6):
    pi = (0 - (0 if i < 3 else 1)) % 2
    for e in w2.monomials(2, pi):
        gens7.append(Fi[i] * R.mono(e))
rows = [idx7.vec(Pg) for Pg in gens7]
cur = rows[:]
currank = nmod_mat(len(cur), idx7.n, [x for r_ in cur for x in r_], p).rank()
nonJ = []
for D, n in fam7['basis']:
    v = idx7.vec(D)
    rk = nmod_mat(len(cur) + 1, idx7.n, [x for r_ in cur for x in r_] + v, p).rank()
    if rk > currank:
        nonJ.append((D, n)); cur.append(v); currank = rk
res7['n_nonJ_directions'] = len(nonJ)
# multiplier tests: m * K in J?
mults = {'L': L, 'M': M, 'X0^3': P('X0**3'), 'X0*L': core.mul(P('X0'), L), 'X0': P('X0'), 'X1': P('X1'), 'Y0': P('Y0'),
         'X0^2': P('X0**2'), 'L^2': core.mul(L, L), 'L*M': core.mul(L, M), 'f': f, 'g': g, 'Y0^3': P('Y0**3'), 'X0^3*Y0': P('X0**3*Y0')}
for name, mpoly in mults.items():
    mp = R.from_dict(mpoly)
    cnt = 0
    for D, n in nonJ:
        ok_, info = w2.in_J(R, Fi, mp * D)
        cnt += ok_
    res7[f'mult_{name}_kills'] = f'{cnt}/{len(nonJ)} basis directions'
# base-type / shape subspaces
def xmons(k):
    return [R.mono(e) for e in w2.monomials(k) if e[3] == e[4] == e[5] == 0]
def ymons(k):
    return [R.mono(e) for e in w2.monomials(k) if e[0] == e[1] == e[2] == 0]
Mp, Lp, fp, gp = R.from_dict(M), R.from_dict(L), R.from_dict(f), R.from_dict(g)
subspaces = {
    'S7(X)+M*S5(X)': xmons(7) + [Mp * x for x in xmons(5)],
    'S7(X)+L*M*S3(X)': xmons(7) + [Lp * Mp * x for x in xmons(3)],
    'base-type P(X)M^a g^b (deg7)': xmons(7) + [Mp * x for x in xmons(5)] + [Mp * Mp * x for x in xmons(3)] + [Mp ** 3 * x for x in xmons(1)] + [gp * x for x in xmons(1)],
    'X-shape Q(Y)L^a f^b (deg7)': ymons(7) + [Lp * y for y in ymons(5)] + [Lp * Lp * y for y in ymons(3)] + [Lp ** 3 * y for y in ymons(1)] + [fp * y for y in ymons(1)],
}
for name, V in subspaces.items():
    famV = w2.e1_family(R, Fp, Ap, Fi, 7, 0, verbose=False, dbasis=V)
    rJ, rJV = w2.jacobian_dim(R, Fi, 7, 0, dbasis=V)
    res7[f'sub[{name}]'] = dict(dimV=len(V), dim_K_cap_V=famV['dimK'], dim_J_cap_V=rJV, kernel_mod_J=famV['dimK'] - rJV)
res7['seconds'] = round(time.time() - t0, 1)
print(json.dumps(res7, default=str, indent=1), flush=True)
w2.log(dict(exp='recheck-K7-structure', p=p, result=res7))
# d = 8 shapes
t0 = time.time()
V8 = xmons(8) + [Mp * x for x in xmons(6)] + [Mp * Mp * x for x in xmons(4)] + [Mp ** 3 * x for x in xmons(2)] + [Mp ** 4] + [gp * x for x in xmons(2)] + [gp * Mp]
V8x = ymons(8) + [Lp * y for y in ymons(6)] + [Lp * Lp * y for y in ymons(4)] + [Lp ** 3 * y for y in ymons(2)] + [Lp ** 4] + [fp * y for y in ymons(2)] + [fp * Lp]
res8 = {}
for name, V in {'base-type (deg8)': V8, 'X-shape (deg8)': V8x}.items():
    famV = w2.e1_family(R, Fp, Ap, Fi, 8, 0, verbose=False, dbasis=V)
    rJ, rJV = w2.jacobian_dim(R, Fi, 8, 0, dbasis=V)
    res8[f'sub[{name}]'] = dict(dimV=len(V), dim_K_cap_V=famV['dimK'], dim_J_cap_V=rJV, kernel_mod_J=famV['dimK'] - rJV)
res8['seconds'] = round(time.time() - t0, 1)
print(json.dumps(res8, default=str, indent=1), flush=True)
w2.log(dict(exp='recheck-K8-structure', p=p, result=res8))
