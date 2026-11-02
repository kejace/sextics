"""Is the degree-8 primitive E1 kernel (dim 1 at (a,b,c)=(5/7,5/11,5/13)) special or generic in the family
F(a,b,c) = Fermat + a*O1 + b*O2 + c*O3  (O1: X_i^3 X_j Y_i Y_j, O2: X_i^3 X_j Y_k^2, O3: X_i^2 Y_j^3 Y_k)?  Also at the Fermat point."""
import sys, time, json; sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, parse, add, scale
from fractions import Fraction as Fr
from flint import nmod_mat
p = core.PRIMES[0]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fermat = parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
O1 = parse('X0**3*X1*Y0*Y1 + X0**3*X2*Y0*Y2 + X1**3*X0*Y0*Y1 + X1**3*X2*Y1*Y2 + X2**3*X0*Y0*Y2 + X2**3*X1*Y1*Y2')
O2 = parse('X0**3*X1*Y2**2 + X0**3*X2*Y1**2 + X1**3*X0*Y2**2 + X1**3*X2*Y0**2 + X2**3*X0*Y1**2 + X2**3*X1*Y0**2')
O3 = parse('X0**2*(Y1**3*Y2 + Y1*Y2**3) + X1**2*(Y0**3*Y2 + Y0*Y2**3) + X2**2*(Y0**3*Y1 + Y0*Y1**3)')
def Ffam(a, b, c): return add(add(add(Fermat, scale(O1, a)), scale(O2, b)), scale(O3, c))
assert core.sub(Ffam(Fr(5,7), Fr(5,11), Fr(5,13)), F) == {}, 'family does not reproduce F'
def e1_kernel_dims(Fq, d=8):
    Fp = core.reduce_mod(Fq, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    Dmons = monomials(d, 0)
    jcols = []
    for i in range(6):
        par_needed = (1 - (0 if i < 3 else 1)) % 2
        for m in monomials(d + 1, par_needed): jcols.append(shift(Fi[i], m, 1, p))
    acols = [shift(Ap, m, 1, p) for m in Dmons]
    rows = {}
    for col in jcols + acols:
        for e in col: rows.setdefault(e, len(rows))
    R = len(rows); CJ = len(jcols); CA = len(acols); C = CJ + CA
    flat = [0] * (R * C)
    for j, col in enumerate(jcols + acols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    X, nullity = nmod_mat(R, C, flat, p).nullspace()
    if nullity == 0: return 0, 0
    K = nmod_mat(CA, nullity, [int(X[CJ + j, k]) for j in range(CA) for k in range(nullity)], p).rank()
    jd = []
    for i in range(6):
        par_needed = (0 - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, par_needed): jd.append(shift(Fi[i], m, 1, p))
    idx = {m: j for j, m in enumerate(Dmons)}
    flat = [0] * (len(jd) * CA)
    for r, col in enumerate(jd):
        for e, v in col.items(): flat[r * CA + idx[e]] = v
    rkJ = nmod_mat(len(jd), CA, flat, p).rank()
    return K, rkJ
pts = [('instance (5/7,5/11,5/13)', (Fr(5,7), Fr(5,11), Fr(5,13))),
       ('Fermat (0,0,0)', (0, 0, 0)),
       ('O1 only (5/7,0,0)', (Fr(5,7), 0, 0)), ('O2 only (0,5/11,0)', (0, Fr(5,11), 0)), ('O3 only (0,0,5/13)', (0, 0, Fr(5,13))),
       ('a=c=0 line (0,5/11,0) same as O2', None),
       ('b=c=0 ... generic a (0,0) (2/3,0,0)', (Fr(2,3), 0, 0)),
       ('no O1: (0,5/11,5/13)', (0, Fr(5,11), Fr(5,13))), ('no O2: (5/7,0,5/13)', (Fr(5,7), 0, Fr(5,13))), ('no O3: (5/7,5/11,0)', (Fr(5,7), Fr(5,11), 0)),
       ('random (3/5,-7/9,2/7)', (Fr(3,5), Fr(-7,9), Fr(2,7))), ('random (-1/2,4/3,5/6)', (Fr(-1,2), Fr(4,3), Fr(5,6))),
       ('scaled t=1/2', (Fr(5,14), Fr(5,22), Fr(5,26))), ('scaled t=2', (Fr(10,7), Fr(10,11), Fr(10,13)))]
d = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for name, abc in pts:
    if abc is None: continue
    t0 = time.time()
    K, rkJ = e1_kernel_dims(Ffam(*abc), d)
    print(f'd={d} {name:40s} dim K_d={K} dim J_d^+={rkJ} primitive kernel={K-rkJ}  [{time.time()-t0:.0f}s]', flush=True)
