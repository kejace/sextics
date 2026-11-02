"""J5: kernel of multiplication by cones over Z_X (= X cap {Y=0}) with a RANDOM vertex plane P (not tau-stable),
on the full R_6 (both parities) -- if zero, V(D_P) serves no primitive (3,1)-class of either tau-parity."""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.jacobian import jacobian_generators, _vectorize, reduce_mod
from hodge.core import monomials, mul
from flint import nmod_mat
import sympy as sp
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2'); V=[X0,X1,X2,Y0,Y1,Y2]
F, A, q, rec = core.load_instance('diagonal-permutation')
p = 2147483647
random.seed(1)
def kernel_full(F, D, d, p):
    Dp = reduce_mod(D, p)
    cols_D = []
    Dm = monomials(d, 0) + monomials(d, 1)
    for m in Dm:
        cols_D.append(mul(Dp, {m:1}, p))
    gens = jacobian_generators(F, d+core.degree(D), p, None)
    cols = cols_D + [g[2] for g in gens]
    rows = _vectorize(cols)
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j,col in enumerate(cols):
        for e,v in col.items(): flat[rows[e]*C+j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    K = len(Dm)
    if nullity == 0: return 0, 0
    proj = nmod_mat(nullity, K, [int(X[j,k]) for k in range(nullity) for j in range(K)], p)
    Pr, rk = proj.rref()
    # rk = dim {D' in S_d : D D' in J}; subtract dim J_d
    Jd = jacobian_generators(F, d, p, None)
    colsJ = [g[2] for g in Jd]; rowsJ = _vectorize(colsJ)
    flatJ = [0]*(len(rowsJ)*len(colsJ))
    for j,col in enumerate(colsJ):
        for e,v in col.items(): flatJ[rowsJ[e]*len(colsJ)+j] = v
    rJ = nmod_mat(len(rowsJ), len(colsJ), flatJ, p).rank()
    return rk - rJ, rk
log = open('wave2-construct.log.jsonl','a')
for trial in range(2):
    L = sp.Matrix(3,6, lambda i,j: random.randint(-3,3))
    lin = [sum(L[i,j]*V[j] for j in range(6)) for i in range(3)]
    D = sp.expand(sum(l**6 for l in lin))
    Dd = core.parse(str(D))
    t0=time.time(); k, rk = kernel_full(F, Dd, 6, p); dt=time.time()-t0
    r = {"exp":"J5-randcone-kernel","D":"cone over Z_X, vertex ker(L), L=%s"%(str(L.tolist())),"d":6,"prime":p,
         "result":"ker dim %d (both parities)"%k,"seconds":round(dt,1),"job":"local","note":"random vertex plane; 0 => V(D) serves no primitive (3,1) class"}
    print(r, flush=True); log.write(json.dumps(r)+"\n"); log.flush()
