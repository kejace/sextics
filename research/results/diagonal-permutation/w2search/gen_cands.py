"""Structured irreducible-type candidate denominators for diagonal-permutation (degrees 6-12), written to cands-structured.jsonl."""
import sys, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
F, A, q, rec = core.load_instance('diagonal-permutation')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
def S(P): return core.to_str(P)
def pd(P, i): return core.pderiv(P, i)
def restrict(P, zero):
    return {e: v for e, v in P.items() if all(e[i] == 0 for i in zero)}
cands = []
def add(D, fam):
    if isinstance(D, dict):
        if not D: return
        D = S(D)
    cands.append({'D': D, 'family': fam})
# A. cones with vertex a plane (3-variable polys) and vertex a line (4-variable restrictions of F)
add(restrict(F, [3, 4, 5]), 'cone-plane F|Y=0 (Fermat X)')
add(restrict(F, [0, 1, 2]), 'cone-plane F|X=0 (Fermat Y)')
for zero, tag in [([4, 5], 'F|Y1=Y2=0'), ([3, 5], 'F|Y0=Y2=0'), ([3, 4], 'F|Y0=Y1=0'), ([1, 2], 'F|X1=X2=0'), ([0, 2], 'F|X0=X2=0'), ([0, 1], 'F|X0=X1=0'),
                  ([2, 5], 'F|X2=Y2=0 (01-block K3 cone)'), ([0, 3], 'F|X0=Y0=0'), ([1, 4], 'F|X1=Y1=0'),
                  ([2], 'F|X2=0'), ([5], 'F|Y2=0'), ([0], 'F|X0=0'), ([3], 'F|Y0=0'), ([2, 4], 'F|X2=Y1=0'), ([1, 5], 'F|X1=Y2=0')]:
    add(restrict(F, zero), 'cone-line ' + tag)
# Fermat/Hesse-type cones over plane sextics (X-block, Y-block, mixed pair blocks)
for a, b, c in [('X0', 'X1', 'X2'), ('Y0', 'Y1', 'Y2'), ('X0', 'X1', 'Y2'), ('X0', 'Y0', 'X1'), ('X0', 'Y0', 'Y1'), ('X2', 'Y2', 'X0')]:
    add(f'{a}**3*{b}**3+{b}**3*{c}**3+{c}**3*{a}**3', f'cone-plane cubic-pairs ({a},{b},{c})')
    add(f'{a}**6+{b}**6+{c}**6-10*{a}**3*{b}**3-10*{b}**3*{c}**3-10*{c}**3*{a}**3', f'cone-plane sextic ({a},{b},{c})')
    add(f'({a}**3+{b}**3+{c}**3)**2-27*{a}**2*{b}**2*{c}**2', f'cone-plane Hesse-disc ({a},{b},{c})')
    add(f'{a}**4*{b}**2+{b}**4*{c}**2+{c}**4*{a}**2', f'cone-plane Klein-type ({a},{b},{c})')
    add(f'{a}**5*{b}+{b}**5*{c}+{c}**5*{a}', f'cone-plane Klein-sextic ({a},{b},{c})')
# B. (01)-adapted: r=X0-X1, s=Y0-Y1, p=X0+X1, qq=Y0+Y1
r, s, p_, qq, a_, b_ = '(X0-X1)', '(Y0-Y1)', '(X0+X1)', '(Y0+Y1)', 'X2', 'Y2'
add(f'{r}**6+{s}**6+{a_}**6', 'cone (r,s,a) Fermat'); add(f'{r}**6+{s}**6+{b_}**6', 'cone (r,s,b) Fermat')
add(f'{p_}**6+{qq}**6+{a_}**6', 'cone (p,q,a) Fermat'); add(f'{r}**3*{s}**3+{s}**3*{a_}**3+{a_}**3*{r}**3', 'cone (r,s,a) cubic pairs')
add(f'{r}**2*{s}**2*({r}**2+{s}**2)+{a_}**6', 'cone (r,s,a) mixed')
add(f'{r}**4*{s}**2+{s}**4*{a_}**2+{a_}**4*{r}**2', 'cone (r,s,a) Klein-type')
add(f'{r}**6+{s}**6+{a_}**6+{b_}**6', 'cone-line (r,s,a,b) Fermat')
add(f'{r}**6+{s}**6+{a_}**3*{b_}**3', 'cone-line (r,s,a,b) mixed')
# C. Hessian 2x2 minors (degree 8) and 3x3 block determinants (degree 12), bordered forms
Fij = [[pd(pd(F, i), j) for j in range(6)] for i in range(6)]
for (i, j) in [(0, 1), (0, 3), (3, 4), (0, 4), (2, 5), (0, 2), (3, 5)]:
    M = core.sub(core.mul(Fij[i][i], Fij[j][j]), core.mul(Fij[i][j], Fij[i][j]))
    add(M, f'Hessian 2x2 minor ({V[i]},{V[j]})')
def det3(idx_r, idx_c):
    tot = {}
    for perm in itertools.permutations(range(3)):
        sign = 1
        for x in range(3):
            for y in range(x + 1, 3):
                if perm[x] > perm[y]: sign = -sign
        term = {(0,) * 6: 1}
        for x in range(3):
            term = core.mul(term, Fij[idx_r[x]][idx_c[perm[x]]])
        tot = core.add(tot, core.scale(term, sign))
    return tot
add(det3([0, 1, 2], [0, 1, 2]), 'Hessian 3x3 det XX (d=12)')
add(det3([3, 4, 5], [3, 4, 5]), 'Hessian 3x3 det YY (d=12)')
add(det3([0, 1, 2], [3, 4, 5]), 'Hessian 3x3 det XY (d=12)')
add(det3([0, 3, 2], [0, 3, 2]), 'Hessian 3x3 det (X0,Y0,X2) (d=12)')
# bordered 2x2 determinants F_ii F_j^2 - 2 F_ij F_i F_j + F_jj F_i^2 (degree 14) -- too big; skip
# D. block Fermat-annihilator type at d=8,10 (non-products, S3/tau adapted)
add('X0**2*X1**2*Y0**2*Y1**2+X1**2*X2**2*Y1**2*Y2**2+X2**2*X0**2*Y2**2*Y0**2', 'orbit sum X_i^2X_j^2Y_i^2Y_j^2 (d=8)')
add('X0**4*Y0**4+X1**4*Y1**4+X2**4*Y2**4', 'sum T_i^(4/3) (d=8)')
add('X0**4*Y0**4+X1**4*Y1**4-2*X2**4*Y2**4', 'B-type d=8')
add('X0**3*Y0**3*X1*Y1+X1**3*Y1**3*X0*Y0', 'T0 X1Y1 + T1 X0Y0 (d=8)')
add('X0**5*Y0**5+X1**5*Y1**5+X2**5*Y2**5', 'sum (X_iY_i)^5 (d=10)')
add('X0**3*Y0**3*X1**2*Y1**2+X1**3*Y1**3*X0**2*Y0**2', 'T0 (X1Y1)^2 + T1 (X0Y0)^2 (d=10)')
# E. polar-type non-Jacobian: F restricted to the diagonal directions / S3-fixed data
add('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6', 'Fermat part P (d=6)')
add(core.sub(F, core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')), '5G = F - Fermat (d=6)')
for k, O in [('O1', 'X0**3*X1*Y0*Y1+X0*X1**3*Y0*Y1+X0**3*X2*Y0*Y2+X0*X2**3*Y0*Y2+X1**3*X2*Y1*Y2+X1*X2**3*Y1*Y2'),
             ('O2', 'X0**3*X1*Y2**2+X0*X1**3*Y2**2+X0**3*X2*Y1**2+X0*X2**3*Y1**2+X1**3*X2*Y0**2+X1*X2**3*Y0**2'),
             ('O3', 'X0**2*Y1**3*Y2+X0**2*Y1*Y2**3+X1**2*Y0**3*Y2+X1**2*Y0*Y2**3+X2**2*Y0**3*Y1+X2**2*Y0*Y1**3')]:
    add(O, f'orbit sum {k} (d=6)')
    add(f'X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6+({O})', f'Fermat + {k} (d=6)')
# F. Y-odd candidates d=6..9 (parity of D free)
add('X0**3*Y0**3+X1**3*Y1**3+X2**3*Y2**3', 'T0+T1+T2 (d=6, Y-odd)')
add('X0**3*Y0**3+X1**3*Y1**3-2*X2**3*Y2**3', 'B2 (d=6, Y-odd)')
add('X0**3*Y0**3', 'T0 (d=6, Y-odd)')
add('X0**2*Y0*X1**2*Y1*X2**2*Y2', 'X^2 Y product (d=9, Y-odd)')
add('X0**3*Y0**3*X2**2*Y2+X1**3*Y1**3*X2**2*Y2', '(T0+T1) X2^2 Y2 (d=9 Y-odd)')
add('X0**4*Y0**4*X2*Y2+X1**4*Y1**4*X2*Y2', '(d=10 Y-even) (T0'+'^(4/3)+T1^(4/3)) X2Y2')
with open('cands-structured.jsonl', 'w') as f:
    for c in cands: f.write(json.dumps(c) + '\n')
import collections
print(len(cands), 'candidates; degrees', collections.Counter(core.degree(core.parse(c['D'])) for c in cands))
