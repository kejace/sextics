"""Build the necessary quadratic system for polar denominators D = w(F) of degree d (w of degree e = d-5):
   w^T Hess(F) w  in  J_{d+1}^+   (valid when [A] is injective on R_{d+1}^+ -> R_{d+7}^-).
Unknowns: coefficients of w_i (parity: w_i Y-even for i<3, Y-odd for i>=3), modulo the gauge w -> w + l*Euler.
Output: msolve input file (char p), with a random affine slice L(c) = 1.  Usage: python3 polar_quadsys.py d p"""
import sys, time, random, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
d = int(sys.argv[1]); p = int(sys.argv[2])
e = d - 5
F, A, q = load()
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
# basis of w
basis = []
for i in range(NV):
    par = 0 if i < 3 else 1
    for m in monomials(e, par):
        basis.append((i, m))
nvar = len(basis)
# gauge vectors: l in S_{e-1}^+ (so that l*Z_i has the right parity), w = l * Z_i e_i
gauge = []
for l in monomials(e - 1, 0):
    v = [0]*nvar
    for k, (i, m) in enumerate(basis):
        ei = [0]*NV; ei[i] = 1
        lm = tuple(a + b for a, b in zip(l, ei))
        if lm == m: v[k] = 1
    gauge.append(v)
ng = len(gauge)
# J_{d+1}^+ generators and normal-form map
deg = d + 1
mons = monomials(deg, 0); idx = {m: j for j, m in enumerate(mons)}
gens = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(deg - 5, pi):
        gens.append(shift(Fi[i], m, 1, p))
G = nmod_mat(len(gens), len(mons), [0]*(len(gens)*len(mons)), p)
for r, g in enumerate(gens):
    for mm, v in g.items(): G[r, idx[mm]] = v
Gr, rk = G.rref()
pivots = []
for r in range(rk):
    j = 0
    while int(Gr[r, j]) == 0: j += 1
    pivots.append(j)
pivset = set(pivots)
nonpiv = [j for j in range(len(mons)) if j not in pivset]
print(f'd={d}: nvar={nvar} gauge={ng} dim J_{deg}^+ = {rk}  dim R_{deg}^+ = {len(nonpiv)}', flush=True)
# vectorized normal forms: V (npairs x nmons) ; NF = V - V[:,pivots] @ Gr  ; keep non-pivot columns
pairs = []
rowsV = []
for k in range(nvar):
    ik, mk = basis[k]
    for l in range(k, nvar):
        il, ml = basis[l]
        P = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1, p)
        if P:
            pairs.append((k, l)); rowsV.append(P)
npairs = len(pairs); nm = len(mons)
flat = [0]*(npairs*nm)
for r, P in enumerate(rowsV):
    for mm, c in P.items(): flat[r*nm + idx[mm]] = c % p
V = nmod_mat(npairs, nm, flat, p)
VP = nmod_mat(npairs, rk, [int(V[r, j]) for r in range(npairs) for j in pivots], p)
NF = V - VP * Gr
Qcoef = {}
for r, (k, l) in enumerate(pairs):
    nf = [int(NF[r, j]) for j in nonpiv]
    if any(nf): Qcoef[(k, l)] = nf
print(f'normal forms done in {time.time()-t0:.1f}s; {len(Qcoef)} nonzero pairs', flush=True)
# eliminate gauge: find ng coordinates with invertible gauge submatrix, set them to 0 (drop variables)
Gm = nmod_mat(ng, nvar, [x for v in gauge for x in v], p)
Gmr, grk = Gm.rref()
gpiv = []
for r in range(grk):
    j = 0
    while int(Gmr[r, j]) == 0: j += 1
    gpiv.append(j)
keep = [k for k in range(nvar) if k not in set(gpiv)]
varnames = [f'c{k}' for k in keep]
kpos = {k: t for t, k in enumerate(keep)}
# equations
eqs = []
neq = len(nonpiv)
for r in range(neq):
    terms = []
    for (k, l), nf in Qcoef.items():
        if k not in kpos or l not in kpos: continue
        c = nf[r]
        if k != l: c = (2 * c) % p
        if c: terms.append(f'{c}*c{k}*c{l}' if k != l else f'{c}*c{k}^2')
    if terms: eqs.append('+'.join(terms))
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
L = '+'.join(f'{random.randrange(1, p)}*{v}' for v in varnames) + f'-1'
eqs.append(L)
fn = f'/tmp/hodge/results/dense-rational/polar_d{d}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n')
    f.write(f'{p}\n')
    f.write(',\n'.join(eqs) + '\n')
print(f'wrote {fn}: {len(varnames)} vars, {len(eqs)} eqs (incl. affine slice)', flush=True)
