"""Exact (mod p) Hessian relaxation at a Jacobian-type degree as a system of quadrics, for Singular/msolve.

usage: python3 w2_relax_exact.py D PD P
Writes w2-relax-d{D}-p{PD}-{P}.sing (Singular driver), .ms (msolve input), .json (metadata).
Variables: coordinates a_k of D = sum_k a_k m_k F_{i_k} (basis jac_basis), with the 12 F-multiple directions
removed by setting a_{(0, X0*h)} = 0 (a complement of the F-multiple subspace).  Quadrics: coefficients of
NF_J(A * Hess_F(a,a)) on the standard monomials of (R/J)_{2D}^{odd}.  LS = linear forms cutting the spurious
component  D = g*h  (g = f(U) - f(V), h quadric of parity PD+1)."""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import w2lib as L
from hodge import core, jacobian
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, monomials
from fractions import Fraction as Fr
from flint import nmod_mat
d, pD, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
F, A, q = L.load()
Fi, Fij = L.Fparts(F)
PF = pD < 0   # parity-free (both sectors; mixed-parity D allowed)
basis = (L.jac_basis(d, 0) + L.jac_basis(d, 1)) if PF else L.jac_basis(d, pD); K = len(basis)
SB = L.space(2 * d, None if PF else 1)
std, red = jacobian.jacobian_ring_basis(F, 2 * d, p, parity=(None if PF else 1))
assert red.mons == SB.mons
nstd = len(std); stdidx = [SB.idx[m] for m in std]; piv = red.pivots
print(f'K={K} rows={SB.n} J-rank={red.rank} quotient={nstd} ({time.time()-t0:.1f}s)', flush=True)
# reducer rows restricted to standard columns: R_std (npiv x nstd)
stdpos = {SB.idx[m]: k for k, m in enumerate(std)}
Rflat = [0] * (len(piv) * nstd)
for r, row in enumerate(red.rows):
    base = r * nstd
    for e, v in row.items():
        j = SB.idx[e]
        k = stdpos.get(j)
        if k is not None:
            Rflat[base + k] = v
R_std = nmod_mat(len(piv), nstd, Rflat, p)
pivpos = {c: i for i, c in enumerate(piv)}
# pair polynomials A m_k m_l F_{ik il} mod p
pairs = [(k, l) for k in range(K) for l in range(k, K)]
Ap = core.reduce_mod(A, p)
Fijp = [[core.reduce_mod(Fij[i][j], p) for j in range(NV)] for i in range(NV)]
def pairvec(k, l):
    ik, mk, _ = basis[k]; il, ml, _ = basis[l]
    P = shift(shift(mul(Ap, Fijp[ik][il], p), mk, 1, p), ml, 1, p)
    return P
quad = np.zeros((len(pairs), nstd), dtype=np.int64)
CH = 1500
t1 = time.time()
for c0 in range(0, len(pairs), CH):
    chunk = pairs[c0:c0 + CH]
    Vp = [0] * (len(chunk) * len(piv)); Vs = [0] * (len(chunk) * nstd)
    for r, (k, l) in enumerate(chunk):
        P = pairvec(k, l)
        bp = r * len(piv); bs = r * nstd
        for e, v in P.items():
            j = SB.idx[e]
            if j in pivpos:
                Vp[bp + pivpos[j]] = v
            else:
                Vs[bs + stdpos[j]] = v
    Mp = nmod_mat(len(chunk), len(piv), Vp, p); Ms = nmod_mat(len(chunk), nstd, Vs, p)
    NF = Ms - Mp * R_std
    for r in range(len(chunk)):
        for k in range(nstd):
            quad[c0 + r, k] = int(NF[r, k])
    print(f'  pairs {c0+len(chunk)}/{len(pairs)} ({time.time()-t1:.0f}s)', flush=True)
# sanity: a spurious point must satisfy all quadrics: a for g*h
f = L.f_ternary(); g = sub(L.uv_to_xy(f), L.uv_to_xy(L.swapUV(f)))
S_D = L.space(d, None if PF else pD)
Mb = np.zeros((S_D.n, K), dtype=np.int64)
for k, (i0, m, Dk) in enumerate(basis):
    for e, c in core.reduce_mod(Dk, p).items():
        Mb[S_D.idx[e], k] = c
Mbn = nmod_mat(S_D.n, K, Mb.ravel().tolist(), p)
# drop coordinates (0, X0*h) for h in monomials(d-6, pD): complement of F-multiples
drop = []
pos = {(i0, m): k for k, (i0, m, _) in enumerate(basis)}
for h in monomials(d - 6, None if PF else pD):
    m = list(h); m[0] += 1; drop.append(pos[(0, tuple(m))])
keep = [k for k in range(K) if k not in set(drop)]
kpos = {k: i for i, k in enumerate(keep)}
def project_quotient(a):
    """a (mod p vector, length K) -> representative with dropped coords = 0 (a - sum 6 a_drop s_h)."""
    a = list(a)
    for h in monomials(d - 6, None if PF else pD):
        m = list(h); m[0] += 1; kd = pos[(0, tuple(m))]
        c = a[kd]
        if c:
            for i in range(NV):
                mm = list(h); mm[i] += 1
                a[pos[(i, tuple(mm))]] = (a[pos[(i, tuple(mm))]] - c) % p   # 6 * c * (1/6)
    assert all(a[k] == 0 for k in drop)
    return [a[k] for k in keep]
spur_pts = []
for h in monomials(d - 6, None if PF else (pD + 1) % 2):
    v = S_D.vec(mul(g, {h: Fr(1)}), dtype=np.int64, p=p)
    aug = nmod_mat(S_D.n, K + 1, np.concatenate([Mb, v.reshape(-1, 1)], axis=1).ravel().tolist(), p)
    Mr, rk = aug.rref()
    a = [0] * K
    for i in range(rk):
        j = 0
        while j <= K and int(Mr[i, j]) == 0: j += 1
        assert j < K, 'g*h not in span'
        a[j] = int(Mr[i, K])
    spur_pts.append(project_quotient(a))
# check quadrics vanish on spurious points
def eval_quads(apt):
    full = [0] * K
    for i, k in enumerate(keep): full[k] = apt[i]
    vals = np.zeros(nstd, dtype=object)
    for r, (k, l) in enumerate(pairs):
        c = (full[k] * full[l]) % p
        if c:
            vals += (c * quad[r].astype(object)) % p
    return [int(x) % p for x in vals]
chk = [max(eval_quads(a)) for a in spur_pts[:2]]
print('spurious points satisfy quadrics:', chk, flush=True)
# linear forms cutting the spurious component: nullspace of the (nspur x nkeep) matrix
nk = len(keep)
Sm = nmod_mat(len(spur_pts), nk, [x for a in spur_pts for x in a], p)
Nsp, nul = Sm.nullspace()      # columns = vectors v with Sm v = 0 : linear forms sum v_i x_i vanish on S
LS = []
for c in range(nul):
    LS.append([int(Nsp[i, c]) for i in range(nk)])
print(f'spurious component dim {len(spur_pts)} -> {len(LS)} linear forms', flush=True)
# write Singular
name = f'w2-relax-d{d}-p{pD}-{p}'
with open(name + '.sing', 'w') as fh:
    fh.write('LIB "elim.lib";\n')
    fh.write(f'ring r = {p}, (x(1..{nk})), dp;\n')
    fh.write('option(redSB);\n')
    fh.write('ideal I;\n')
    for s in range(nstd):
        terms = []
        for r, (k, l) in enumerate(pairs):
            c = int(quad[r, s])
            if c and k in kpos and l in kpos:
                if c > p // 2: c -= p
                terms.append(f'({c})*x({kpos[k]+1})*x({kpos[l]+1})')
        if terms:
            fh.write(f'I[{s+1}] = ' + '+'.join(terms) + ';\n')
    fh.write('ideal LS;\n')
    for i, v in enumerate(LS):
        terms = [f'({(c if c <= p//2 else c-p)})*x({j+1})' for j, c in enumerate(v) if c]
        fh.write(f'LS[{i+1}] = ' + '+'.join(terms) + ';\n')
    fh.write('''
int t = timer; system("--ticks-per-sec", 1);
"nvars =", nvars(r), " nquadrics =", ncols(I), " nlinforms =", ncols(LS);
ideal G = std(I);
"std(I): dim =", dim(G), " mult =", mult(G), " size =", size(G), " time =", timer - t;
intvec hv = hilb(G, 1); "hilbert numerator (first kind):"; hv;
t = timer;
ideal Isat = sat(I, LS);
ideal Gs = std(Isat);
"sat(I, LS): dim =", dim(Gs), " mult =", mult(Gs), " size =", size(Gs), " time =", timer - t;
"dim -1 means: the only solutions of the relaxation (mod p) are the spurious ones (empty projective residual)";
if (dim(Gs) >= 1) { "RESIDUAL COMPONENT PRESENT: dim", dim(Gs), "mult", mult(Gs); intvec hv2 = hilb(Gs, 1); hv2; }
quit;
''')
with open(name + '.ms', 'w') as fh:
    fh.write(','.join(f'x{i+1}' for i in range(nk)) + '\n' + str(p) + '\n')
    lines = []
    for s in range(nstd):
        terms = []
        for r, (k, l) in enumerate(pairs):
            c = int(quad[r, s])
            if c and k in kpos and l in kpos:
                terms.append(f'{c}*x{kpos[k]+1}*x{kpos[l]+1}')
        if terms: lines.append('+'.join(terms))
    fh.write(',\n'.join(lines) + '\n')
meta = {'d': d, 'pD': pD, 'p': p, 'K': K, 'nvars': nk, 'nquadrics': nstd, 'J_rank_2d': red.rank, 'spurious_dim': len(spur_pts),
        'n_linear_forms': len(LS), 'spurious_check_max_residue': chk, 'seconds': round(time.time() - t0, 1)}
json.dump(meta, open(name + '.json', 'w'), indent=1)
L.log_line(f'w2-relax-exact-d{d}.jsonl', {'experiment': 'hessian-relaxation-exact-setup', **meta})
print(json.dumps(meta), flush=True)
