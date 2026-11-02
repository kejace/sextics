"""Independent check of wave2-excl R2 (degree 6): D = sum g_i F_i, g in S_1^6 (36 unknowns w_{ik}), all characters jointly.
ob(g) = -[A sum_{ij} g_i g_j F_ij] in (S/J)_12.  Quadrics = coordinates of the class in a monomial basis of (S/J)_12
(all exponents <= 4 -- a basis iff good reduction, checked).  Euler line T = C.z is spurious (D = 6cF); ob is T-invariant,
so restrict to the hyperplane W = {w_00 = 0} (35 vars).  Then run Singular: std, dim, vdim, kbase (max degree) mod p.
If kbase has no monomial of degree >= 7 then I_p contains S_7(W) and, by rank_p <= rank_Q of Macaulay matrices, so does
I_Q: V(ob) = T over C, hence no degree-6 denominator (given E1 => D in J_6 in all characters)."""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000003
Fs, As, qv = load_rec(); F = modp(parse(Fs), p); A = modp(parse(As), p)
Fi = [deriv(F, i, p) for i in range(6)]
Fij = [[deriv(Fi[i], j, p) for j in range(6)] for i in range(6)]
t0 = time.time()
mons = monos(12); idx = {m: j for j, m in enumerate(mons)}
Jrows = [shift(Fi[i], m) for i in range(6) for m in monos(7)]
MJ = nmod_from_rows(Jrows, idx, len(mons), p)
Mr, rJ = MJ.rref()
std = [m for m in mons if max(m) <= 4]
print(f'rank J_12 = {rJ}, dim (S/J)_12 = {len(mons)-rJ}, Fermat-type monomials {len(std)}', flush=True)
# normal-form map: coordinates w.r.t. the basis 'std' of (S/J)_12.  Build the reduction: for a vector v in S_12, find
# coordinates c (len(std)) with v - sum c_m m in rowspace(J).  Do it by solving [J; std]^T x = v for all v at once:
# equivalently compute an rref of [Jrows; std-rows] and check pivots cover everything.
rows = Jrows + [{m: 1} for m in std]
Mfull = nmod_from_rows(rows, idx, len(mons), p)
assert Mfull.rank() == len(mons), 'std monomials do not complete J to a basis (bad reduction)'
# To get coordinates we solve the transposed system; use nullspace trick: coordinates of v = last len(std) entries of a
# solution of  [Jrows;std]^T x = v.  Precompute via rref of the augmented transposed matrix for all needed v.
pairs = [(i, k) for i in range(6) for k in range(6)]       # w_{ik}: g_i = sum_k w_{ik} z_k
zk = [tuple(1 if t == k else 0 for t in range(6)) for k in range(6)]
# class vectors v_{(ik),(jl)} = A z_k z_l F_ij  for (ik) <= (jl)
vecs = []; labels = []
for a in range(36):
    for b in range(a, 36):
        i, k = pairs[a]; j, l = pairs[b]
        if not Fij[i][j]: continue
        v = mul(A, shift(Fij[i][j], tuple(x + y for x, y in zip(zk[k], zk[l]))), p)
        vecs.append(v); labels.append((a, b))
print(f'{len(vecs)} nonzero class vectors to reduce ({time.time()-t0:.0f}s)', flush=True)
# Solve X * Mfull = V  (X: nvec x (nJ+nstd)), i.e. Mfull^T X^T = V^T.  Use flint solve via rref of [Mfull^T | V^T].
nJ = len(Jrows); ns = len(std); R = nJ + ns
MT = nmod_mat(len(mons), R + len(vecs), p)
for r, poly in enumerate(rows):
    for e, c in poly.items(): MT[idx[e], r] = c % p
for c, v in enumerate(vecs):
    for e, val in v.items(): MT[idx[e], R + c] = val
t1 = time.time()
MTr, rk = MT.rref()
assert rk == len(mons)
# since Mfull^T has full column rank (= len(mons) rows? no: it is len(mons) x R with R >= len(mons)); extract solution
# of the rows: after rref, rows correspond to pivots; the solution for free columns is non-unique only in the J part.
# Simpler: coordinates of v w.r.t. std are unique; recover them by reading the rref rows whose pivot is a std column.
coords = [[0] * ns for _ in vecs]
for r in range(rk):
    j = 0
    while MTr[r, j] == 0: j += 1
    if j >= nJ and j < R:
        for c in range(len(vecs)):
            coords[c][j - nJ] = int(MTr[r, R + c])
print(f'coordinates computed ({time.time()-t1:.0f}s)', flush=True)
# quadrics: for each basis element m of (S/J)_12, Q_m = sum_{(a<=b)} coords[c][m] * w_a w_b  ; restrict w_00 = 0 (a=0)
varnames = [f'w{i}{k}' for (i, k) in pairs]
quads = []
for mi in range(ns):
    terms = []
    for c, (a, b) in enumerate(labels):
        cf = coords[c][mi]
        if cf and a != 0 and b != 0:
            terms.append(f'{cf}*{varnames[a]}*{varnames[b]}')
    if terms: quads.append(' + '.join(terms))
print(f'{len(quads)} nonzero quadrics in 35 variables', flush=True)
with open(f'hess6_p{p}.sing', 'w') as f:
    f.write(f'ring r = {p}, ({",".join(varnames[1:])}), dp;\n')
    f.write('ideal I = ' + ',\n'.join(quads) + ';\n')
    f.write('int t0 = timer; option(redSB); ideal G = std(I); "GB size:", size(G), " time:", timer - t0;\n')
    f.write('"dim:", dim(G); "vdim:", vdim(G);\n')
    f.write('ideal kb = kbase(G); int i; int md = 0; for (i = 1; i <= size(kb); i++) { if (deg(kb[i]) > md) { md = deg(kb[i]); } }\n')
    f.write('"kbase size:", size(kb), " max degree of a standard monomial:", md;\n')
    f.write('intvec hf; int d; for (d = 0; d <= md + 1; d++) { hf[d+1] = size(kbase(G, d)); } "Hilbert function HF(0..md+1):", hf;\n')
    f.write('quit;\n')
print('Singular input written; running std ...', flush=True)
import subprocess
out = subprocess.run(['timeout', '3000', 'Singular', '-q', f'hess6_p{p}.sing'], capture_output=True, text=True)
print(out.stdout); print(out.stderr[-500:])
logjson(VLOG, {'check': 'hess6', 'p': p, 'n_quadrics': len(quads), 'singular': out.stdout.strip()[:600]})
