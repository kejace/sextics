"""Necessary quadratic system for polar denominators D = w(F) of degree d = 7 in a tau-pure sector.
   parity(D) = pd in {0,1};  w_i has parity (pd + p_i) mod 2, p_i = 0 (i<3), 1 (i>=3).
   Condition: w^T Hess(F) w in (J:A)_8^+ = J_8^+   (valid because [A] is injective on R_8^+, R_8^-).
   Gauge w -> w + l*Euler, l in S_1^{pd}, removed by dropping 3 pivot coordinates.
Usage: python3 polar_quadsys2.py d pd p [slice|homog] [seed]"""
import sys, time, random
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
t0 = time.time()
d = int(sys.argv[1]); pd = int(sys.argv[2]); p = int(sys.argv[3])
mode = sys.argv[4] if len(sys.argv) > 4 else 'slice'
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
assert d == 7, 'only d=7 uses J_{2d-6} = (J:A)_{2d-6}'
e = d - 5
F, A, q = load()
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
basis = []
for i in range(NV):
    par = (pd + (0 if i < 3 else 1)) % 2
    for m in monomials(e, par):
        basis.append((i, m))
nvar = len(basis)
gauge = []
for l in monomials(e - 1, pd):
    v = [0]*nvar
    for k, (i, m) in enumerate(basis):
        ei = [0]*NV; ei[i] = 1
        if tuple(a + b for a, b in zip(l, ei)) == m: v[k] = 1
    gauge.append(v)
ng = len(gauge)
deg = 2*e + 4
mons = monomials(deg, 0); idx = {m: j for j, m in enumerate(mons)}
gens = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(deg - 5, pi):
        gens.append(shift(Fi[i], m, 1, p))
G = nmod_mat(len(gens), len(mons), p)
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
print(f'd={d} pd={pd}: nvar={nvar} gauge={ng} dim J_{deg}^+ = {rk}  dim R_{deg}^+ = {len(nonpiv)}', flush=True)
pairs = []; rowsV = []
for k in range(nvar):
    ik, mk = basis[k]
    for l in range(k, nvar):
        il, ml = basis[l]
        Pm = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1, p)
        if Pm: pairs.append((k, l)); rowsV.append(Pm)
npairs = len(pairs); nm = len(mons)
V = nmod_mat(npairs, nm, p)
for r, Pm in enumerate(rowsV):
    for mm, c in Pm.items(): V[r, idx[mm]] = c % p
VP = nmod_mat(npairs, rk, [int(V[r, j]) for r in range(npairs) for j in pivots], p)
NF = V - VP * Gr
Qcoef = {}
for r, (k, l) in enumerate(pairs):
    nf = [int(NF[r, j]) for j in nonpiv]
    if any(nf): Qcoef[(k, l)] = nf
print(f'normal forms done in {time.time()-t0:.1f}s; {len(Qcoef)} nonzero pairs', flush=True)
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
eqs = []
for r in range(len(nonpiv)):
    terms = []
    for (k, l), nf in Qcoef.items():
        if k not in kpos or l not in kpos: continue
        c = nf[r]
        if k != l: c = (2 * c) % p
        if c: terms.append(f'{c}*c{k}*c{l}' if k != l else f'{c}*c{k}^2')
    if terms: eqs.append('+'.join(terms))
random.seed(seed)
if mode == 'slice':
    eqs.append('+'.join(f'{random.randrange(1, p)}*{v}' for v in varnames) + '-1')
fn = f'/tmp/hodge/results/dense-rational/polar_d{d}_pd{pd}_{mode}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n'); f.write(f'{p}\n'); f.write(',\n'.join(eqs) + '\n')
print(f'wrote {fn}: {len(varnames)} vars, {len(eqs)} eqs, mode={mode} ({time.time()-t0:.1f}s)', flush=True)
