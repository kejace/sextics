"""Hessian relaxation for Jacobian-type denominators (exclusion lens, diagonal-permutation, q = 1).

Setting.  Every E1-admissible D of degree d satisfies [A][D] = 0 in R = S/J(F).  When the primitive E1 kernel in
degree d is zero (d <= 7, and d = 8 Y-odd), D = w.grad F for a vector w of six polynomials of degree d-5; at d = 8 the
Y-even part may also contain a*D8 (A*D8 = sum U_i F_i).  With N = N0 + K + G z (N0 = A w [+ a U], K Koszul,
G z Euler gauge) one checks:  the Euler gauge G cancels exactly in E2 := D div N - N(D) - 6 G D, and the Koszul
terms D div K - K(D) lie in J(F).  Hence a NECESSARY condition for a certificate is

      E(u) := D(u) div N0(u) - N0(u)(D(u))  in  J(F)_{2d}        (u = coefficients of w [and a]),

a system of dim R_{2d} quadrics in u.  For D = w.grad F this is  A * w^T H_F w  in J(F).

This script builds the quadrics mod p (p < 2^26 so sparse int64 products are exact), reports the rank of the quadric
span inside Sym^2(U*), the dimension count of the trivial cone {M z : M in S_{d-6}} (D = 6 M F, not coprime), and
writes msolve / Singular input files.  Self-test: for random (w, c, G) the full E2 agrees with E(u) modulo J.

usage: python3 hess_relax.py d [--sector even|odd|all] [--p P] [--with-d8]"""
import sys, time, json, argparse, itertools, random
import numpy as np
import scipy.sparse as sps
from flint import nmod_mat
from common import *
from hodge import jacobian
from hodge.core import monomials, shift, pderiv, mul, add, sub, scale, reduce_mod, degree

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int)
ap.add_argument('--sector', default='all')       # parity of D: even | odd | all (mixed)
ap.add_argument('--p', type=int, default=67108859)  # prime < 2^26
ap.add_argument('--with-d8', action='store_true')
ap.add_argument('--selftest', action='store_true')
ap.add_argument('--tag', default='')
args = ap.parse_args()
d, p = args.d, args.p
F, A, q, rec = core.load_instance(IID)
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
t0 = time.time()

# ---- unknowns ------------------------------------------------------------------------------------
# w_i in S_{d-5}; parity restriction: D even <=> w_X even, w_Y odd ; D odd <=> w_X odd, w_Y even
basis = []   # entries: ('w', i, m) or ('d8',)
for i in range(6):
    for m in monomials(d - 5):
        par_m = core.yparity(m)
        par_D = (par_m + (0 if i < 3 else 1)) % 2       # parity of m*F_i
        if args.sector == 'even' and par_D != 0: continue
        if args.sector == 'odd' and par_D != 1: continue
        basis.append(('w', i, m))
extra = []
if args.with_d8:
    assert d == 8 and args.sector in ('even', 'all')
    D8 = reduce_mod(core.parse(json.load(open(os.path.join(LANE, 'D8_e1kernel.json')))['D']), p)
    # cofactors U mod p: A*D8 = sum U_i F_i  (solve dense)
    AD = mul(Ap, D8, p)
    cols = [(i, m) for i in range(6) for m in monomials(9, (1 - (0 if i < 3 else 1)) % 2)]
    rows = monomials(14, 1); ridx = {e: k for k, e in enumerate(rows)}
    R_, C_ = len(rows), len(cols)
    flat = [0] * (R_ * (C_ + 1))
    for j, (i, m) in enumerate(cols):
        for e, v in shift(Fi[i], m, 1, p).items(): flat[ridx[e] * (C_ + 1) + j] = v
    for e, v in AD.items(): flat[ridx[e] * (C_ + 1) + C_] = v
    Mb = nmod_mat(R_, C_ + 1, flat, p)
    Rr, rk = Mb.rref()
    U = [dict() for _ in range(6)]
    jj = 0
    for i_ in range(rk):
        while jj < C_ + 1 and int(Rr[i_, jj]) == 0: jj += 1
        assert jj < C_, 'A*D8 not in J mod p ?!'
        v = int(Rr[i_, C_])
        if v:
            i, m = cols[jj]; U[i][m] = v
        jj += 1
    chk = {}
    for i in range(6):
        for e, v in mul(U[i], Fi[i], p).items(): chk[e] = (chk.get(e, 0) + v) % p
    assert all(((chk.get(e, 0) - AD.get(e, 0)) % p) == 0 for e in set(chk) | set(AD))
    extra.append(('d8', D8, U))
    basis.append(('d8',))
n = len(basis)
print(f'd={d} sector={args.sector} unknowns n={n} (w-part {n - len(extra)})', flush=True)

def Dpoly(b):
    if b[0] == 'w': return shift(Fi[b[1]], b[2], 1, p)
    return extra[0][1]
def Npoly(b):
    if b[0] == 'w':
        N = [dict() for _ in range(6)]; N[b[1]] = shift(Ap, b[2], 1, p); return N
    return extra[0][2]
Ds = [Dpoly(b) for b in basis]
Ns = [Npoly(b) for b in basis]
dDs = [[pderiv(D_, i, p) for i in range(6)] for D_ in Ds]
divNs = []
for N in Ns:
    s = {}
    for i in range(6):
        for e, v in pderiv(N[i], i, p).items(): s[e] = (s.get(e, 0) + v) % p
    divNs.append({e: v for e, v in s.items() if v})

def Pkl(k, l):
    """D_k div(N_l) - N_l(D_k)"""
    out = mul(Ds[k], divNs[l], p)
    for i in range(6):
        if Ns[l][i]:
            out = sub(out, mul(Ns[l][i], dDs[k][i], p), p)
    return out

# ---- monomial indexing in degree 2d, both parities -----------------------------------------------
mons = monomials(2 * d); midx = {e: j for j, e in enumerate(mons)}
M = len(mons)
pairs = [(k, l) for k in range(n) for l in range(k, n)]
print('pairs', len(pairs), 'monomials', M, flush=True)
rows_i = []; cols_j = []; vals = []
for r, (k, l) in enumerate(pairs):
    Pq = Pkl(k, l) if k == l else add(Pkl(k, l), Pkl(l, k), p)
    for e, v in Pq.items():
        rows_i.append(r); cols_j.append(midx[e]); vals.append(v)
P = sps.csr_matrix((np.array(vals, dtype=np.int64), (rows_i, cols_j)), shape=(len(pairs), M))
print('quadric tensor built nnz', P.nnz, round(time.time() - t0, 1), 's', flush=True)

# ---- normal form modulo J_{2d}, per parity ------------------------------------------------------
quad_rows = []   # list of (parity, dense matrix pairs x standard-monomials)
std_all = []
Qmats = []
for par in (0, 1):
    red = jacobian._Reducer(Fp, 2 * d, p, par)
    pm = red.mons; pidx = {e: j for j, e in enumerate(pm)}
    glob = np.array([midx[e] for e in pm])          # columns of P for this parity
    Ppar = P[:, glob].tocsc()
    piv = red.pivots; r = len(piv)
    stdc = [j for j in range(len(pm)) if j not in set(piv)]
    # Jr[:, std] dense r x |std|
    Jr = np.zeros((r, len(stdc)), dtype=np.int64)
    sidx = {j: t for t, j in enumerate(stdc)}
    for i_, row in enumerate(red.rows):
        for e, v in row.items():
            j = pidx[e]
            if j in sidx: Jr[i_, sidx[j]] = v
    Ppiv = Ppar[:, piv]; Pstd = Ppar[:, stdc].toarray()
    # exact mod p: entries < 2^26, row nnz <= ~500 -> products < 2^52 * 500 < 2^61
    NF = (Pstd - (Ppiv @ Jr)) % p
    Qmats.append((par, NF, [pm[j] for j in stdc]))
    print(f'parity {par}: J rank {r}, R dim {len(stdc)}, NF block {NF.shape}', round(time.time() - t0, 1), 's', flush=True)
# the quadric coefficient matrix: rows = quadrics (standard monomials of both parities), cols = pairs
Qc = np.concatenate([NF.T for _, NF, _ in Qmats], axis=0)      # (dimR) x pairs
nq = Qc.shape[0]
if args.with_d8:
    # the relaxation is modulo J only (D8 is not in J); nothing to project.  (E(u) in J is the condition.)
    pass
Qn = nmod_mat(nq, len(pairs), [int(x) for x in Qc.ravel()], p)
rank = Qn.rank()
sym2 = n * (n + 1) // 2
# trivial cone: w = M z with M in S_{d-6} (parity-compatible with the sector): dimension
cone = 0
if d >= 6:
    for Mm in monomials(d - 6):
        parM = core.yparity(Mm)
        if args.sector == 'even' and parM != 0: continue
        if args.sector == 'odd' and parM != 1: continue
        cone += 1
# verify cone members are solutions: pick M = each monomial, w = M z
def eval_quadrics(u):
    uu = np.array([(u[k] * u[l]) % p for (k, l) in pairs], dtype=np.int64)
    return (Qc @ uu) % p if False else np.array([int(sum(int(Qc[i, j]) * int(uu[j]) for j in np.nonzero(uu)[0]) % p) for i in range(nq)])
cone_ok = True
bidx = {b: k for k, b in enumerate(basis)}
for Mm in monomials(d - 6) if d >= 6 else []:
    parM = core.yparity(Mm)
    if args.sector == 'even' and parM != 0: continue
    if args.sector == 'odd' and parM != 1: continue
    u = [0] * n
    for i in range(6):
        mi = list(Mm); mi[i] += 1
        key = ('w', i, tuple(mi))
        if key in bidx: u[bidx[key]] = 1
    vals_ = eval_quadrics(u)
    if any(vals_): cone_ok = False
# ---- self-test of the gauge cancellation (random w, c, G) ---------------------------------------
selftest = None
if args.selftest:
    random.seed(1)
    u = [random.randrange(p) for _ in range(n)]
    D = {}; N = [dict() for _ in range(6)]
    for k, b in enumerate(basis):
        D = add(D, scale(Ds[k], u[k], p), p)
        for i in range(6): N[i] = add(N[i], scale(Ns[k][i], u[k], p), p)
    # Koszul K_i = sum_j c_ij F_j, c antisymmetric degree d-4 ; Euler G z, deg G = d
    c = {}
    for i in range(6):
        for j in range(i + 1, 6):
            c[(i, j)] = {m: random.randrange(p) for m in monomials(d - 4)}
    K = [dict() for _ in range(6)]
    for (i, j), cij in c.items():
        K[i] = add(K[i], mul(cij, Fi[j], p), p); K[j] = sub(K[j], mul(cij, Fi[i], p), p)
    G = {m: random.randrange(p) for m in monomials(d)}
    Nfull = [add(add(N[i], K[i], p), mul(G, {tuple(1 if t == i else 0 for t in range(6)): 1}, p), p) for i in range(6)]
    # E1 check: sum N_i F_i = A D + 6 G F
    s = {}
    for i in range(6): s = add(s, mul(Nfull[i], Fi[i], p), p)
    e1 = sub(sub(s, mul(Ap, D, p), p), scale(mul(G, Fp, p), 6, p), p)
    # full E2 := D div N - N(D) - 6 G D
    divN = {}
    for i in range(6): divN = add(divN, pderiv(Nfull[i], i, p), p)
    E2 = mul(D, divN, p)
    for i in range(6): E2 = sub(E2, mul(Nfull[i], pderiv(D, i, p), p), p)
    E2 = sub(E2, scale(mul(G, D, p), 6, p), p)
    # relaxation value E(u) = D div N0 - N0(D)
    divN0 = {}
    for i in range(6): divN0 = add(divN0, pderiv(N[i], i, p), p)
    Eu = mul(D, divN0, p)
    for i in range(6): Eu = sub(Eu, mul(N[i], pderiv(D, i, p), p), p)
    diff = sub(E2, Eu, p)
    inJ = all(jacobian.in_jacobian(Fp, {e: v for e, v in diff.items() if core.yparity(e) == par}, p)[0]
              for par in (0, 1) if any(core.yparity(e) == par for e in diff))
    selftest = {'e1_identity_holds': not e1, 'E2full_minus_E_in_J': inJ}
    print('selftest', selftest, flush=True)

out = {'exp': 'hessian-relaxation', 'd': d, 'sector': args.sector, 'p': p, 'unknowns': n, 'with_d8': args.with_d8,
       'num_quadrics_dimR2d': nq, 'sym2_dim': sym2, 'quadric_span_rank': rank, 'sym2_minus_rank': sym2 - rank,
       'trivial_cone_dim_affine': cone, 'cone_members_satisfy': cone_ok, 'selftest': selftest,
       'seconds': round(time.time() - t0, 1)}
log(out)
# ---- write Groebner inputs (reduced to a basis of the quadric span) ------------------------------
tag = f'd{d}_{args.sector}{"_d8" if args.with_d8 else ""}{args.tag}'
Qr, rk2 = Qn.rref()
names = [f'u{k}' for k in range(n)]
polys = []
for i in range(rk2):
    terms = []
    for j, (k, l) in enumerate(pairs):
        v = int(Qr[i, j])
        if v:
            if v > p // 2: v -= p
            terms.append(f'{v}*{names[k]}*{names[l]}' if k != l else f'{v}*{names[k]}^2')
    polys.append('+'.join(terms).replace('+-', '-'))
with open(f'relax_{tag}.ms', 'w') as f:
    f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
with open(f'relax_{tag}.sing', 'w') as f:
    f.write(f'ring r = {p}, ({",".join(names)}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n')
    f.write('option(redSB);\nideal G = std(I);\nprint("DIM " + string(dim(G)));\nprint("MULT " + string(mult(G)));\n')
    f.write('print("HILB"); hilb(G);\nprint("NGENS " + string(size(G)));\nint k; int nlin = 0; for (k = 1; k <= size(G); k++) { if (deg(G[k]) == 1) { nlin++; } }\nprint("NLIN " + string(nlin));\nquit;\n')
json.dump({'basis': [list(map(str, b)) for b in basis], 'pairs': pairs, 'p': p, 'd': d, 'sector': args.sector},
          open(f'relax_{tag}_basis.json', 'w'))
print('wrote', f'relax_{tag}.ms', f'relax_{tag}.sing', 'rank', rank, 'of', sym2, round(time.time() - t0, 1), 's', flush=True)
