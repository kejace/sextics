"""Generate degree-8 (and a few other) denominator candidates in the (U,V) frame (X slots = U, Y slots = V)."""
import sys, json; sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction as Fr
uv = json.load(open('uv_forms.json'))
f = core.parse(uv['f'])                       # ternary sextic in slots 0..2 (U)
def swapP(P):  # U<->V
    return {(e[3], e[4], e[5], e[0], e[1], e[2]): c for e, c in P.items()}
def sym(P): return core.add(P, swapP(P))
def anti(P): return core.sub(P, swapP(P))
fV = swapP(f)
fi = [core.pderiv(f, i) for i in range(3)]
fij = [[core.pderiv(fi[i], j) for j in range(3)] for i in range(3)]
def var(i): e=[0]*6; e[i]=1; return {tuple(e): Fr(1)}
U = [var(i) for i in range(3)]; V = [var(3+i) for i in range(3)]
def pw(P, k): return core.power(P, k)
S = core.to_str
C = []  # (family, label, expr, frame)
# --- controls (reducible; expected inconsistent by the sum principle / degree<=7 exclusion)
C.append(('control-reducible', 'fU-U0sq-sym', S(sym(core.mul(f, pw(U[0], 2)))), 'uv'))
C.append(('control-reducible', 'fU-U0V0-sym', S(sym(core.mul(f, core.mul(U[0], V[0])))), 'uv'))
C.append(('control-reducible', 'fU-UdotV-sym', S(sym(core.mul(f, core.add(core.add(core.mul(U[0],V[0]),core.mul(U[1],V[1])),core.mul(U[2],V[2]))))), 'uv'))
# --- Fermat octics (irreducible)
C.append(('fermat-octic', 'sumU8+sumV8', S(sym(core.add(core.add(pw(U[0],8), pw(U[1],8)), pw(U[2],8)))), 'uv'))
C.append(('fermat-octic', 'sumU8-sumV8', S(anti(core.add(core.add(pw(U[0],8), pw(U[1],8)), pw(U[2],8)))), 'uv'))
C.append(('fermat-octic', 'U8w123+V8w123', S(sym(core.add(core.add(pw(U[0],8), core.scale(pw(U[1],8),2)), core.scale(pw(U[2],8),3)))), 'uv'))
# cone with vertex a line (4 linear forms): U0^8+U1^8+V0^8+V1^8
C.append(('cone-line-fermat', 'U0^8+U1^8+V0^8+V1^8', S(core.add(core.add(pw(U[0],8), pw(U[1],8)), core.add(pw(V[0],8), pw(V[1],8)))), 'uv'))
# cone with vertex the plane {Y=0} (XY frame): sum Y_i^8 ; and {X=0}: sum X_i^8
C.append(('cone-plane-fermat', 'Y0^8+Y1^8+Y2^8', 'Y0**8 + Y1**8 + Y2**8', 'xy'))
C.append(('cone-plane-fermat', 'X0^8+X1^8+X2^8', 'X0**8 + X1**8 + X2**8', 'xy'))
C.append(('cone-plane-fermat', 'Y0^8+2Y1^8+3Y2^8', 'Y0**8 + 2*Y1**8 + 3*Y2**8', 'xy'))
# cone with vertex plane {Y=0} over the octic Hessian minor curve of f and the polar octic
M01 = core.sub(core.mul(fij[0][0], fij[1][1]), core.mul(fij[0][1], fij[0][1]))
M02 = core.sub(core.mul(fij[0][0], fij[2][2]), core.mul(fij[0][2], fij[0][2]))
M12 = core.sub(core.mul(fij[1][1], fij[2][2]), core.mul(fij[1][2], fij[1][2]))
def toY(P): return {(0,0,0,e[0],e[1],e[2]): c for e, c in P.items()}
def toX(P): return P
C.append(('cone-plane-hessminor', 'M01(f)(Y)', S(toY(M01)), 'xy'))
C.append(('cone-plane-hessminor', 'M12(f)(Y)', S(toY(M12)), 'xy'))
C.append(('cone-plane-hessminor', 'M01(f)(X)', S(toX(M01)), 'xy'))
# --- Hessian minors of f symmetrized in the UV frame
C.append(('hessminor-sym', 'M01(U)+M01(V)', S(sym(M01)), 'uv'))
C.append(('hessminor-sym', 'M01(U)-M01(V)', S(anti(M01)), 'uv'))
C.append(('hessminor-sym', 'M12(U)+M12(V)', S(sym(M12)), 'uv'))
C.append(('hessminor-sym', 'M02(U)+M02(V)', S(sym(M02)), 'uv'))
C.append(('hessminor-sym', 'M01+M02+M12(U)+(V)', S(sym(core.add(core.add(M01, M02), M12))), 'uv'))
# --- polar-type octics: sum_i d_i f(U) V_i^3 + swap
pol3 = {}
for i in range(3):
    pol3 = core.add(pol3, core.mul(fi[i], pw(V[i], 3)))
C.append(('polar', 'sum_i f_i(U)V_i^3 + swap', S(sym(pol3)), 'uv'))
C.append(('polar', 'sum_i f_i(U)V_i^3 - swap', S(anti(pol3)), 'uv'))
# sum_i f_i(U) U_i^3  + swap  (pure block octic covariant)
pol3U = {}
for i in range(3):
    pol3U = core.add(pol3U, core.mul(fi[i], pw(U[i], 3)))
C.append(('polar', 'sum_i f_i(U)U_i^3 + swap', S(sym(pol3U)), 'uv'))
# second polar squared-type: (sum_ij f_ij(U) V_i V_j) is degree 6 (excluded); times nothing.
# --- perturbations of the reducible fU*U0^2: + t*(U0^8+V0^8)
base = sym(core.mul(f, pw(U[0], 2)))
for t in [Fr(1), Fr(-1), Fr(1,7)]:
    C.append(('perturb', f'fU*U0^2+swap + {t}*(U0^8+V0^8)', S(core.add(base, core.scale(core.add(pw(U[0],8), pw(V[0],8)), t))), 'uv'))
json.dump(C, open('candidates_d8.json', 'w'), indent=0)
print(len(C), 'candidates')
for c in C: print(c[0], '|', c[1], '| terms', len(core.parse(c[2])))
