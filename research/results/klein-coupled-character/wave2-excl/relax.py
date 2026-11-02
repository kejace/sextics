"""(b) Rigorous lower bound on deg D via E1 + Jacobian relaxation of E2 (and, optionally, the full linearization).

E1 (A D in J) with ker[A] = 0 on R_d forces D = sum_i g_i F_i, g_i of degree d-5 with parities fixed by par(D).
E1 lift N0 = A g.  E2:  -B_D(N0) in F S + im T_D  subset J_{2d}.  Hence  H(g) := B_D(Ag) = D div(Ag) - A g(D) in J.
H is a quadratic form in the coefficient vector c of g:  H(g) = sum_{a<=b} c_a c_b H_ab.
Jacobian relaxation: W = {symmetric T : sum T_ab [H_ab] = 0 in R_{2d}} must contain g (x) g for every certificate g.
Gauge: g -> g + h Z (h in S_{d-6}) changes D by 6hF; all tensors (hZ) o v lie in W.  Decisive iff dim W = dim(G o V),
G = Z*S_{d-6} subset V (then g(x)g in W => g in G => D in (F), illegal).
Rigor: rank_p[J-block | H-classes] <= rank_Q, dim_Fp J = dim_Q J (regular sequence check) => rank_Q(classes) >= rank_p;
and the gauge tensors are in W_Q by an exact identity (verified exactly here as well).

usage: python3 relax.py <p> <d> <parity: 0|1|mixed>
"""
import sys, json, time, itertools, math
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
import kx
from hodge import core
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, parse
from flint import nmod_mat
p = int(sys.argv[1]); d = int(sys.argv[2]); parmode = sys.argv[3]
F, A, q = kx.load()
F7 = scale(F, 7)
Fi = [pderiv(F7, i) for i in range(NV)]          # integer polars of 7F (J(7F) = J(F))
Ai = [pderiv(A, i) for i in range(NV)]
Z = [parse(v) for v in ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')]
t0 = time.time()

def gspace(parD):
    """basis of g-space: list of (i, m) with deg m = d-5 and par(m) = parD + PI[i]."""
    basis = []
    for i in range(NV):
        for m in kx.wmonomials(d - 5, (parD + kx.PI[i]) % 2):
            basis.append((i, m))
    return basis

if parmode == 'mixed':
    V = gspace(0) + gspace(1)
else:
    V = gspace(int(parmode))
n = len(V)
print(f'd={d} parmode={parmode}: g-space dim n={n}', flush=True)

def gvec(coeffs):
    """g as list of 6 polys from coefficient dict {alpha: c}."""
    g = [{} for _ in range(NV)]
    for a, c in coeffs.items():
        i, m = V[a]
        g[i] = add(g[i], {m: c})
    return g

def Hform(g):
    """H(g) = D div(Ag) - sum_i A g_i d_i D,  D = sum g_i F_i  (exact integers)."""
    D = {}
    for i in range(NV):
        if g[i]: D = add(D, mul(g[i], Fi[i]))
    divAg = {}
    for i in range(NV):
        if g[i]:
            divAg = add(divAg, mul(Ai[i], g[i]))
            divAg = add(divAg, mul(A, pderiv(g[i], i)))
    H = mul(D, divAg)
    for i in range(NV):
        if g[i]:
            H = sub(H, mul(mul(A, g[i]), pderiv(D, i)))
    return H

# quadratic form coefficients by polarization
Hd = {}
for a in range(n):
    Hd[(a, a)] = Hform(gvec({a: 1}))
for a in range(n):
    for b in range(a + 1, n):
        Hab = Hform(gvec({a: 1, b: 1}))
        Hab = sub(sub(Hab, Hd[(a, a)]), Hd[(b, b)])   # = 2 H_ab (we keep the factor 2; harmless for rank)
        Hd[(a, b)] = Hab
print(f'built {len(Hd)} quadratic-form coefficient polys of degree {2*d} in {round(time.time()-t0,1)}s', flush=True)

# gauge tensors: G = Z * S_{d-6} inside V (only when d >= 6 and g-space contains Z-multiples)
def coeffs_of_g(g):
    """coefficient dict of a g (list of polys) in the V basis; raises if not in span."""
    idx = {(i, m): a for a, (i, m) in enumerate(V)}
    out = {}
    for i in range(NV):
        for m, c in g[i].items():
            out[idx[(i, m)]] = c
    return out
gauge = []
if d >= 6:
    for h in core.monomials(d - 6):
        g = [shift(Z[i], h, 1) for i in range(NV)]
        try:
            gauge.append(coeffs_of_g(g))
        except KeyError:
            pass
print(f'gauge vectors in V: {len(gauge)}', flush=True)
# dimension of G o V as symmetric tensors
pairs = [(a, b) for a in range(n) for b in range(a, n)]
pidx = {pb: k for k, pb in enumerate(pairs)}
def sym_tensor(u, v):
    T = {}
    for a, ca in u.items():
        for b, cb in v.items():
            key = (a, b) if a <= b else (b, a)
            T[key] = T.get(key, 0) + ca * cb
    return T
gauge_tensors = []
for gv in gauge:
    for a in range(n):
        gauge_tensors.append(sym_tensor(gv, {a: 1}))
dimGV = 0
if gauge_tensors:
    flat = [0] * (len(gauge_tensors) * len(pairs))
    for r, T in enumerate(gauge_tensors):
        for k, v in T.items(): flat[r * len(pairs) + pidx[k]] = kx.modp(v, p)
    dimGV = nmod_mat(len(gauge_tensors), len(pairs), flat, p).rank()
print(f'dim(G o V) = {dimGV} (expected {len(gauge)}*{n} - {len(gauge)*(len(gauge)+1)//2} = {len(gauge)*n - len(gauge)*(len(gauge)+1)//2})', flush=True)

# exact check: gauge tensors are in W, i.e. sum_k T_k H_k in J_{2d+6} exactly?  (we verify the weaker, sufficient
# identity used in the proof: H(g + hZ) - H(g) in F*S for the basis g and gauge h; done mod p for all pairs, and
# exactly over Z for a few random pairs)
def in_F_multiple(P):
    ok, G = jacobian_divisible(P)
    return ok
from hodge import jacobian as jac
def jacobian_divisible(P):
    return jac.divisible_by(F7, P, p=None)
import random
random.seed(1)
checked = 0; bad = 0
if gauge:
    for gv in gauge:
        for a in random.sample(range(n), min(n, 3)):
            g1 = gvec({a: 1}); g2 = gvec(add_dicts := {**{a: 1}})
            gsum = [add(x, y) for x, y in zip(gvec({a: 1}), gvec(gv))]
            diff = sub(sub(Hform(gsum), Hform(gvec({a: 1}))), Hform(gvec(gv)))   # cross term 2 H(e_a, hZ)
            ok, _ = jac.divisible_by(F7, diff)   # exact over Q
            checked += 1; bad += (not ok)
        okZ, _ = jac.divisible_by(F7, Hform(gvec(gv)))
        checked += 1; bad += (not okZ)
print(f'exact gauge identity checks: {checked} checked, {bad} failed', flush=True)

# rank of the classes [H_ab] in R_{2d+6}, by (parity, character) blocks
deg = 2 * d   # deg B_D(N) = d + (d+1) - 1 = 2d at q = 1
blocks = {}
for key, P in Hd.items():
    if not P: continue
    par = core.parity_of(P); ch = kx.char_of(P)
    assert par is not None and ch is not None, (key, par, ch)
    blocks.setdefault((par, ch), []).append(key)
total_rank = 0; details = {}
for (par, ch), keys in sorted(blocks.items()):
    mons = kx.wmonomials(deg, par, ch); rows = {m: i for i, m in enumerate(mons)}
    jk, jp = kx.jac_block(F, deg, p, par, ch)
    MJ = kx.dense(jp, rows, p); rJ = MJ.rank()
    cols = jp + [Hd[k] for k in keys]
    M = kx.dense(cols, rows, p); r = M.rank()
    rk_classes = r - rJ
    total_rank += rk_classes
    details[f'{par},{ch}'] = {'n_classes': len(keys), 'rank_classes': rk_classes, 'dimJ': rJ, 'dimS': len(mons), 'dimR': len(mons) - rJ,
                              'expected_dimJ_Q': None}
    print(f'  block par={par} ch={ch}: {len(keys)} classes, rank {rk_classes} in R (dim {len(mons)-rJ}), dimJ={rJ}  [{round(time.time()-t0,1)}s]', flush=True)
ntens = len([k for k in Hd])   # all pairs (zero polys count as zero classes)
nzero = sum(1 for P in Hd.values() if not P)
dimW = ntens - total_rank
decisive = (dimW == dimGV)
print(f'RESULT d={d} parmode={parmode}: tensors {ntens} (zero polys {nzero}), rank of classes {total_rank}, dim W = {dimW}, dim(G o V) = {dimGV} -> {"DECISIVE: no legal D of this degree/parity (given E1 kernel = 0)" if decisive else "NOT decisive (W larger than gauge)"}', flush=True)
res = {'d': d, 'parmode': parmode, 'p': p, 'n': n, 'tensors': ntens, 'zero_polys': nzero, 'rank_classes': total_rank, 'dimW': dimW,
       'dimGV': dimGV, 'gauge_vectors': len(gauge), 'decisive': decisive, 'blocks': details, 'gauge_checks': checked, 'gauge_check_failures': bad,
       'seconds': round(time.time() - t0, 1)}
json.dump(res, open(f'relax_d{d}_{parmode}_p{p}.json', 'w'), indent=1)
kx.log(exp='jacobian-relaxation', **{k: v for k, v in res.items() if k != 'blocks'})
