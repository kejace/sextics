import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_families import form_polys, family_monomials
F, A, q, rec = load()
p = 536870909
D7 = core.parse(json.load(open('kernel_d7_p0.json'))['basis_rational'][0])
D7o = core.parse(json.load(open('D7odd_cone.json'))['D_cone'])
K8 = [modp_poly_from_json(o) for o in json.load(open('kernel_d8_p0.json'))['basis_modp'][str(p)]]
jref = JRef(F, 8, 0, p, verbose=False)
sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
def vec(nf): 
    v = [0] * ns
    for m, c in nf.items(): v[sidx[m]] = c
    return v
prods = []; labels = []
for i in range(3):
    prods.append(shift(reduce_mod(D7, p), tuple(1 if t == i else 0 for t in range(NV)), 1, p)); labels.append(f'X{i}*D7')
for i in range(3):
    prods.append(shift(reduce_mod(D7o, p), tuple(1 if t == 3 + i else 0 for t in range(NV)), 1, p)); labels.append(f"Y{i}*D7'")
nf_prods = jref.nf_polys(prods)
Mp = nmod_mat(6, ns, [x for nf in nf_prods for x in vec(nf)], p)
X, nul = Mp.transpose().nullspace()
print('relations among products mod J_8:', nul, flush=True)
for c in range(nul):
    col = [int(X[i, c]) for i in range(6)]
    # normalize by first nonzero and try small rational reconstruction
    j0 = next(i for i in range(6) if col[i]); inv = pow(col[j0], -1, p)
    norm = [(x * inv) % p for x in col]
    rats = [core.ratrecon(x, p) for x in norm]
    print('  relation', c, ':', {labels[i]: str(rats[i]) if rats[i] is not None else norm[i] for i in range(6) if norm[i]}, flush=True)
# new elements: K8 NFs modulo span of product NFs
nf_K = jref.nf_polys(K8)
allm = nmod_mat(6 + 6, ns, [x for nf in nf_prods for x in vec(nf)] + [x for nf in nf_K for x in vec(nf)], p)
print('rank products', Mp.rank(), ' rank products+K8', allm.rank(), flush=True)
# hypothesis: new elements in NF(C[X,M]_8) + NF(C[Y,L]_8)
FP = form_polys()
famXM = [reduce_mod(P, p) for _, P in family_monomials(['X0', 'X1', 'X2', 'M'], 8, 0, FP)]
famYL = [reduce_mod(P, p) for _, P in family_monomials(['Y0', 'Y1', 'Y2', 'L'], 8, 0, FP)]
nf_f = jref.nf_polys(famXM + famYL)
base = nmod_mat(len(nf_f), ns, [x for nf in nf_f for x in vec(nf)], p)
rb = base.rank()
tot = nmod_mat(len(nf_f) + 6, ns, [x for nf in nf_f for x in vec(nf)] + [x for nf in nf_K for x in vec(nf)], p)
print(f'sum-of-two-cones family: rank {rb}; with K8: {tot.rank()} -> K8 directions outside the family: {tot.rank() - rb}', flush=True)
# hypothesis 2: C[X,M]_8 + C[Y,L]_8 + C[L,M]*S_6 ... try C[X,M] + C[Y,L] + L*M*S_6? too big. Try: (L,M)-binary times products? skip.
# hypothesis 3: new elements in C[X0,X1,X2,M] * C[L] + C[Y0,Y1,Y2,L]*C[M]? i.e. C[X,M,L]=C[X,M]. skip.
# hypothesis 4: new elements in the span of {monomials in X,M of degree a} x {monomials in Y,L of degree 8-a}, a=1..7 (products of cones)
for a_ in range(1, 8):
    fa = [reduce_mod(P, p) for _, P in family_monomials(['X0', 'X1', 'X2', 'M'], a_, None, FP)] if False else None
