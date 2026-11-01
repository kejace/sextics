"""Cones with prescribed vertex v in P^3 that pass E1: intersect span{L^alpha : |alpha| = 6} (L1,L2,L3 vanishing at v)
with {D : a D in J(H)}; test every basis element and random combos in the reduced model."""
import sys, random
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
from threefold import *
from flint import nmod_poly
random.seed(31)
p = core.PRIMES[0]
Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
Hi = [pderiv(Hp, i, p) for i in CORE]
Jcols = []
for idx, i in enumerate(CORE):
    for m2 in cmon(7):
        Jcols.append(shift(Hi[idx], m2, 1, p))
MJ, rowsJ = build(Jcols, p); rJ = MJ.rank()
def lin(c): return {tuple(1 if j == i else 0 for j in range(6)): c % p for i, c in c.items() if c % p}
def cones_in_kernel(v, parity=None):
    """v: 4-tuple (x0,x1,y0,y1) mod p. Returns basis (mod p) of cones with vertex v satisfying a D in J(H)."""
    # three independent linear forms vanishing at v: pick pivot coordinate with v_k != 0
    k = max(range(4), key=lambda i: v[i] != 0)
    idxs = [CORE[i] for i in range(4)]
    inv = pow(v[k], -1, p)
    Ls = []
    for i in range(4):
        if i == k: continue
        Ls.append(lin({idxs[i]: 1, idxs[k]: -v[i] * inv}))
    mons = []
    for a_ in range(7):
        for b in range(7 - a_):
            c = 6 - a_ - b
            mons.append(mul(mul(core.power(Ls[0], a_, p), core.power(Ls[1], b, p), p), core.power(Ls[2], c, p), p))
    if parity is not None:
        mons = [m for m in mons if m and core.parity_of(m) == parity]
    Acols = [mul(ap, m, p) for m in mons]
    MAJ, _ = build(Acols + Jcols, p, dict(rowsJ))
    X, nullity = MAJ.nullspace()
    nA = len(mons)
    vecs = [[int(X[j, kk]) for j in range(nA)] for kk in range(nullity)]
    vecs = [vv for vv in vecs if any(vv)]
    if not vecs: return [], Ls
    Mv = nmod_mat(len(vecs), nA, [x for vv in vecs for x in vv], p)
    Mr, rk = Mv.rref()
    basis = []
    for i in range(rk):
        D = {}
        for j in range(nA):
            cij = int(Mr[i, j])
            if cij: D = add(D, scale(mons[j], cij, p), p)
        basis.append(D)
    return basis, Ls
def roots_binary(poly, i, j):
    coeffs = [0] * 7
    for e, c in poly.items(): coeffs[e[i]] = c % p
    return [int(r) for r, m in nmod_poly(coeffs, p).roots()]
f = reduce_mod(core.parse('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6'), p)
g = reduce_mod(core.parse('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6'), p)
tf = roots_binary(f, 0, 1); sg = roots_binary(g, 3, 4)
vertices = [((tf[0], 1, 0, 0), 'v=(t:1:0:0) on S_H, f(t,1)=0'), ((0, 0, sg[0], 1), 'v=(0:0:s:1) on S_H, g(s,1)=0'),
            ((1, 0, 0, 0), 'v=(1:0:0:0) off S_H'), ((0, 0, 1, 0), 'v=(0:0:1:0) off S_H'), ((1, 1, 1, 1), 'v=(1:1:1:1)')]
# general point of S_H
for _ in range(500):
    cand = [random.randrange(1, p) for _ in range(3)]
    coeffs = [0] * 7
    for e, c in Hp.items():
        coeffs[e[4]] = (coeffs[e[4]] + c * pow(cand[0], e[0], p) * pow(cand[1], e[1], p) * pow(cand[2], e[3], p)) % p
    rts = nmod_poly(coeffs, p).roots()
    if rts:
        vertices.append(((cand[0], cand[1], cand[2], int(rts[0][0])), 'v = general F_p-point of S_H')); break
vertices.append(((random.randrange(1, p), random.randrange(1, p), random.randrange(1, p), random.randrange(1, p)), 'v = random point off S_H'))
for v, name in vertices:
    for par in (0, 1, None):
        basis, Ls = cones_in_kernel(v, par)
        print(f'{name}, parity {par}: dim of E1-passing cones = {len(basis)}', flush=True)
        log({'D': f'cone space, {name}, parity {par}', 'd': 6, 'family': 'cone-E1-kernel-dim', 'e1': True, 'consistent': None, 'rank': len(basis), 'nullity': None, 'seconds': 0, 'note': 'dimension of {cones with vertex v} cap {a D in J(H)} (incl. J_6)', 'engine': 'threefold-v1', 'p': p})
        if par is None or not basis: continue
        tests = []
        for r in range(2):
            D = {}
            for b in basis: D = add(D, scale(b, random.randrange(1, p), p), p)
            tests.append((D, f'random E1-passing cone with vertex {name}, parity {par} (#{r})'))
        for i, b in enumerate(basis):
            tests.append((b, f'E1-passing cone basis #{i}, vertex {name}, parity {par}, {len(b)} terms'))
        for D, note in tests:
            run3(D, note, 'threefold-cones-E1kernel', p)
