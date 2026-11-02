import sys, json; sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction as Fr
uv = json.load(open('uv_forms.json'))
f = core.parse(uv['f'])
def swapP(P): return {(e[3], e[4], e[5], e[0], e[1], e[2]): c for e, c in P.items()}
def sym(P): return core.add(P, swapP(P))
def anti(P): return core.sub(P, swapP(P))
fi = [core.pderiv(f, i) for i in range(3)]
fiV = [swapP(g) for g in fi]
S = core.to_str
C = []
# polar pairing  sum_i f_i(U) f_i(V)  (swap-symmetric, bidegree (5,5))
pp = {}
for i in range(3): pp = core.add(pp, core.mul(fi[i], fiV[i]))
C.append(('polar-pairing-d10', 'sum_i f_i(U) f_i(V)', S(pp), 'uv'))
# twisted pairing with a permutation: f_0(U)f_1(V)+f_1(U)f_2(V)+f_2(U)f_0(V) + swap
tp = {}
for i in range(3): tp = core.add(tp, core.mul(fi[i], fiV[(i+1)%3]))
C.append(('polar-pairing-d10', 'sum_i f_i(U) f_{i+1}(V) + swap', S(sym(tp)), 'uv'))
json.dump(C, open('candidates_d10.json', 'w'), indent=0)
for c in C: print(c[0], '|', c[1], '| terms', len(core.parse(c[2])))
