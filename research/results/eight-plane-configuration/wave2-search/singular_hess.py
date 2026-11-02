"""Exact Hessian-relaxation ideal at Jacobian-type degree d (parity parD): quadrics in w with
w^T Hess(F) w in Ker_{2d-6}^+ (= J_{2d-6}^+ for d <= 7; J + modular E1-kernel lift for d = 8, then mod p).
Writes Singular input wave2-search/hess_d{d}_par{parD}.sing computing sat(I, Euler-cone ideal)."""
import sys, time, json, pickle, numpy as np
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-search')
from w2ep import *
from flint import fmpq_mat, fmpq, nmod_mat
d, parD = int(sys.argv[1]), int(sys.argv[2])
p = int(sys.argv[3]) if len(sys.argv) > 3 else 0     # 0 = exact over Q
F, A, q = load()
fam = jacobian_family(F, A, d, parD); wdesc = fam['wdesc']; K = len(wdesc)
ST = space(2 * d - 6, 0)
Fi = [pderiv(F, i) for i in range(NV)]
Hij = {(i, j): pderiv(Fi[i], j) for i in range(NV) for j in range(NV)}
# exact quadratic map: for k <= l, coefficient vector of z^{m_k} z^{m_l} (H_{ik il} + H_{il ik} if k != l)
pairs = []
for k in range(K):
    for l in range(k, K):
        (i, mk), (j, ml) = wdesc[k], wdesc[l]
        mm = tuple(mk[t] + ml[t] for t in range(6))
        P = shift(Hij[(i, j)], mm, 1 if k == l else 2)
        pairs.append(((k, l), P))
# kernel subspace generators (exact over Q for J; mod p lift for extra kernel)
gens = []
for i in range(NV):
    for m in monomials(2 * d - 11, VARPAR[i]):
        gens.append(shift(Fi[i], m, 1))
extra = []
if 2 * d - 6 >= 10:
    assert p, 'extra kernel only available mod p'
    extra = pickle.load(open(f'{WD}/e1kernel_basis_d{2*d-6}_par0.pkl', 'rb'))
idx = {m: t for t, m in enumerate(ST.mons)}
n = ST.n
if p == 0:
    G = fmpq_mat(len(gens), n)
    for r_, P in enumerate(gens):
        for e, v in P.items(): G[r_, idx[e]] = fmpq(v.numerator, v.denominator)
    Gr, rk = G.rref()
    piv = []
    for r_ in range(rk):
        c = 0
        while Gr[r_, c] == 0: c += 1
        piv.append(c)
    nonpiv = [c for c in range(n) if c not in set(piv)]
    print(f'J_{2*d-6}^+ rank {rk}, quotient dim {len(nonpiv)}', flush=True)
    def normal_form(P):
        v = [fmpq(0)] * n
        for e, c in P.items(): v[idx[e]] = fmpq(c.numerator, c.denominator)
        for r_, c in enumerate(piv):
            if v[c] != 0:
                f = v[c] / Gr[r_, c]
                for t in range(n):
                    g = Gr[r_, t]
                    if g != 0: v[t] -= f * g
        return [v[t] for t in nonpiv]
else:
    rows = []
    for P in gens:
        row = [0] * n
        for e, v in core.reduce_mod(P, p).items(): row[idx[e]] = v
        rows.append(row)
    for P in extra:
        row = [0] * n
        for e, v in P.items(): row[idx[e]] = v % p
        rows.append(row)
    G = nmod_mat(len(rows), n, [x for row in rows for x in row], p)
    Gr, rk = G.rref()
    piv = []
    for r_ in range(rk):
        c = 0
        while int(Gr[r_, c]) == 0: c += 1
        piv.append(c)
    nonpiv = [c for c in range(n) if c not in set(piv)]
    print(f'Ker_{2*d-6}^+ rank mod p {rk}, quotient dim {len(nonpiv)}', flush=True)
    Grl = [[int(Gr[r_, t]) for t in range(n)] for r_ in range(rk)]
    def normal_form(P):
        v = [0] * n
        for e, c in core.reduce_mod(P, p).items(): v[idx[e]] = c
        for r_, c in enumerate(piv):
            if v[c]:
                f = v[c] * pow(Grl[r_][c], -1, p) % p
                row = Grl[r_]
                for t in range(n):
                    if row[t]: v[t] = (v[t] - f * row[t]) % p
        return [v[t] for t in nonpiv]
t0 = time.time()
quad = {}   # t -> list of (k, l, coeff)
for (k, l), P in pairs:
    nf = normal_form(P)
    for t, c in enumerate(nf):
        if c != 0:
            quad.setdefault(t, []).append((k, l, c))
print(f'{len(quad)} nonzero quadrics in {K} unknowns ({time.time()-t0:.1f}s)', flush=True)
# Euler cone: w = G z, G in S_{d-6}^{parD}: ideal of this linear subspace = all linear forms vanishing on it
Gm = monomials(d - 6, parD)
E = np.zeros((len(Gm), K))
for g, m in enumerate(Gm):
    for i in range(NV):
        mm = tuple(m[t] + (1 if t == i else 0) for t in range(6)); E[g, wdesc.index((i, mm))] = 1.0
# linear forms vanishing on the span: nullspace of E (rows) -> exact: E has 0/1 entries, use sympy-free integer nullspace via rref over Q
Ef = fmpq_mat(len(Gm), K, [int(x) for x in E.ravel()]) if len(Gm) else None
lin = []
if Ef is not None:
    # nullspace of Ef: vectors v with E v = 0
    Er, rk = Ef.rref()
    pivc = []
    for r_ in range(rk):
        c = 0
        while Er[r_, c] == 0: c += 1
        pivc.append(c)
    free = [c for c in range(K) if c not in set(pivc)]
    for fcol in free:
        v = [fmpq(0)] * K; v[fcol] = fmpq(1)
        for r_, c in enumerate(pivc): v[c] = -Er[r_, fcol] / Er[r_, c]
        lin.append(v)
def var(k): return f'w{k}'
with open(f'{SD_DIR}/hess_d{d}_par{parD}.sing', 'w') as f:
    f.write(f'ring R = {p}, ({",".join(var(k) for k in range(K))}), dp;\n')
    f.write('option(redSB);\nideal I;\n')
    for t, terms in quad.items():
        f.write('I[size(I)+1] = ' + '+'.join(f'({c})*{var(k)}*{var(l)}' for k, l, c in terms) + ';\n')
    f.write('ideal L;\n')
    for v in lin:
        f.write('L[size(L)+1] = ' + '+'.join(f'({c})*{var(k)}' for k, c in enumerate(v) if c != 0) + ';\n')
    f.write('int t0 = timer;\nideal SI = std(I);\n"dim(I) (affine, includes origin) =", dim(SI), " degree =", mult(SI), " time", timer - t0;\n')
    f.write('LIB "elim.lib";\nt0 = timer;\nideal S = sat(SI, L)[1];\n"saturation w.r.t. Euler cone: S == 1 ?", (S[1] == 1), " dim(S) =", dim(std(S)), " time", timer - t0;\n')
    f.write('"Euler cone dim (affine) =", dim(std(L));\nquit;\n')
print('written', f'{SD_DIR}/hess_d{d}_par{parD}.sing', 'with', len(quad), 'quadrics,', len(lin), 'linear forms (Euler cone codim)')
