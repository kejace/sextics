"""Run fixed-D tests for named families.  Usage: python3 families.py <setname>"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
F, A, q = load()
X = ['X0', 'X1', 'X2']; Y = ['Y0', 'Y1', 'Y2']; Z = X + Y
Fi = [pderiv(F, i) for i in range(NV)]
def blocks(F):
    """F_{6-2j,2j} pieces by Y-degree 2j."""
    out = {0: {}, 2: {}, 4: {}, 6: {}}
    for e, c in F.items():
        out[e[3] + e[4] + e[5]][e] = c
    return out
B = blocks(F)
def EX(P_): return add(add(mul(core.parse('X0'), pderiv(P_, 0)), mul(core.parse('X1'), pderiv(P_, 1))), mul(core.parse('X2'), pderiv(P_, 2)))
def EY(P_): return add(add(mul(core.parse('Y0'), pderiv(P_, 3)), mul(core.parse('Y1'), pderiv(P_, 4))), mul(core.parse('Y2'), pderiv(P_, 5)))
sets = {}
# --- degree 6 structural family ---
s6 = []
subsets = []
for I in itertools.combinations(range(6), 3):
    if (0 in I) != (3 in I):   # exactly one of X0, Y0
        subsets.append(I)
for I in subsets:
    s6.append(('fermat-subset ' + '+'.join(Z[i] for i in I), {tuple(6 if k == i else 0 for k in range(6)): Fraction(1) for i in I}))
s6.append(('pure-powers sum Z_i^6', core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')))
s6.append(('F_{6,0}=a(X)', B[0])); s6.append(('F_{0,6}=d(Y)', B[6]))
s6.append(('F_{4,2}=b', B[2])); s6.append(('F_{2,4}=c', B[4]))
s6.append(('a+d (F minus coupling)', add(B[0], B[6])))
s6.append(('E_X(F)=6a+4b+2c', EX(F))); s6.append(('E_Y(F)=2b+4c+6d', EY(F)))
s6.append(('b+2c+3d (dt-derivative)', add(add(B[2], scale(B[4], 2)), scale(B[6], 3))))
s6.append(('a - d', sub(B[0], B[6])))
s6.append(('X0^6+Y0^6 pair', core.parse('X0**6+Y0**6')))
s6.append(('(X0^2+X1^2+X2^2)^3', power(core.parse('X0**2+X1**2+X2**2'), 3)))
s6.append(('(X0^2+Y0^2)^3', power(core.parse('X0**2+Y0**2'), 3)))
sets['deg6'] = s6
# --- degree 5: partials (in J, pass E1) ---
s5 = [(f'F_{Z[i]}', Fi[i]) for i in range(3)]   # even ones: F_X0, F_X1, F_X2 ; F_Yj are odd -> D would be odd: skip
s5.append(('X0^5+X1^5+X2^5', core.parse('X0**5+X1**5+X2**5')))
s5.append(('X0^5+X0*Y0^4', core.parse('X0**5+X0*Y0**4')))
sets['deg5'] = s5
# --- degree <=4 cones ---
s4 = [('X0^4+X1^4+X2^4', core.parse('X0**4+X1**4+X2**4')), ('X0^3+X1^3+X2^3', core.parse('X0**3+X1**3+X2**3')),
      ('X0^4+Y0^4', core.parse('X0**4+Y0**4')), ('X0^2+X1^2+X2^2', core.parse('X0**2+X1**2+X2**2')),
      ('X0^4+X1^4+X2^4+Y0^4+Y1^4+Y2^4', core.parse('X0**4+X1**4+X2**4+Y0**4+Y1**4+Y2**4')),
      ('X0*X1*X2*Y0*Y1*Y2? no (odd) -> X0*X1*X2*(X0+X1+X2)', core.parse('X0*X1*X2*(X0+X1+X2)'))]
sets['deg4'] = s4
# --- degree 7: J-elements w(F) for simple quadratic w, and l*D6 ---
s7 = []
s7.append(('X0^2*F_X0', mul(core.parse('X0**2'), Fi[0])))
s7.append(('X0*E_X(F)', mul(core.parse('X0'), EX(F))))
s7.append(('X0*(X0^6+X1^6+X2^6)', core.parse('X0*(X0**6+X1**6+X2**6)')))
s7.append(('X0*(Y0^6+Y1^6+Y2^6)', core.parse('X0*(Y0**6+Y1**6+Y2**6)')))
s7.append(('X0*F_X0 - Y0*F_Y0 (paired diff)', sub(mul(core.parse('X0'), Fi[0]), mul(core.parse('Y0'), Fi[3]))))
s7.append(('X0*F_X0 + Y0*F_Y0 (paired Euler)', add(mul(core.parse('X0'), Fi[0]), mul(core.parse('Y0'), Fi[3]))))
s7.append(('Y0^2*F_X0', mul(core.parse('Y0**2'), Fi[0])))
s7.append(('X0*Y0*F_Y0', mul(core.parse('X0*Y0'), Fi[3])))
s7.append(('X0^7+X1^7+X2^7', core.parse('X0**7+X1**7+X2**7')))
s7.append(('X0^7+X0*Y0^6', core.parse('X0**7+X0*Y0**6')))
s7.append(('X0*a(X)', mul(core.parse('X0'), B[0])))
s7.append(('X0*d(Y)', mul(core.parse('X0'), B[6])))
s7.append(('Y0^2*F_X0 + X0^2*F_X0', mul(core.parse('X0**2+Y0**2'), Fi[0])))
sets['deg7'] = s7
# --- degree 8: Hessian minors, products, w(F) cubic w ---
s8 = []
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
s8.append(('Hess minor F_X0X0*F_Y0Y0 - F_X0Y0^2', sub(mul(Fij[0][0], Fij[3][3]), mul(Fij[0][3], Fij[0][3]))))
s8.append(('Hess minor F_X0X0*F_X1X1 - F_X0X1^2', sub(mul(Fij[0][0], Fij[1][1]), mul(Fij[0][1], Fij[0][1]))))
s8.append(('X0^3*F_X0', mul(core.parse('X0**3'), Fi[0])))
s8.append(('X0^2*E_X(F)', mul(core.parse('X0**2'), EX(F))))
s8.append(('Y0^3*F_Y0', mul(core.parse('Y0**3'), Fi[3])))
s8.append(('X0^3 F_X0 + Y0^3 F_Y0', add(mul(core.parse('X0**3'), Fi[0]), mul(core.parse('Y0**3'), Fi[3]))))
s8.append(('X0^8+X1^8+X2^8', core.parse('X0**8+X1**8+X2**8')))
s8.append(('(X0^6+X1^6+X2^6)(X0^2+X1^2+X2^2)', core.parse('(X0**6+X1**6+X2**6)*(X0**2+X1**2+X2**2)')))
s8.append(('a(X)*(Y0^2+Y1^2+Y2^2)', mul(B[0], core.parse('Y0**2+Y1**2+Y2**2'))))
s8.append(('E_X(F)*(X0^2+Y0^2)', mul(EX(F), core.parse('X0**2+Y0**2'))))
sets['deg8'] = s8
def in_J(D):
    """Is D in J_d = span{m F_i} (mod p)?  Small rank test."""
    p = P
    d = core.degree(D)
    if d < 5: return False
    Dp = core.reduce_mod(D, p)
    Fp = core.reduce_mod(F, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    par = core.parity_of(D)
    cols = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d - 5, (par - pi) % 2):
            cols.append(shift(Fip[i], m, 1, p))
    rows = {}
    for c in cols + [Dp]:
        for e in c: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*(C+1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e]*(C+1)+j] = v
    for e, v in Dp.items(): flat[rows[e]*(C+1)+C] = v
    M = nmod_mat(R, C+1, flat, p)
    M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    return M.rank() == M0.rank()

if __name__ == '__main__':
    name = sys.argv[1]
    full = '--e1only' not in sys.argv
    for nm, D in sets[name]:
        if core.parity_of(D) != 0:
            print('skip (not Y-even):', nm); continue
        t0 = time.time()
        d = core.degree(D)
        if d <= 9 and not in_J(D):
            log(nm, d, name, False, False, None, None, time.time()-t0, note='D not in J_d; E1 fails since ker([A] on R_d^+)=0 for d<=9 (E1-kernel sweep)')
            print(f"{nm}: d={d} NOT in J_d -> E1 fails (kernel mod J = 0) ({time.time()-t0:.1f}s)", flush=True)
            continue
        r = test_D(F, A, q, D, family=name, name=nm, full=full)
        print(f"{nm}: d={r['d']} e1={r['e1']} consistent={r.get('consistent')} rank={r.get('rank')} nullity={r.get('nullity')} ({time.time()-t0:.1f}s)", flush=True)
