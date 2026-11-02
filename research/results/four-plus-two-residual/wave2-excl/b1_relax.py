"""(b) E2-mod-J relaxation at degree d (Jacobian-type reduction with non-Jacobian E1-kernel part), mod p.
Unknowns v = (c, w): D(v) = sum_b c_b K_b + sum_i w_i F_i, with K_b the lifts of the E1 kernel modulo J(F)_d
(Thom-Sebastiani: K_b = (core kernel element of degree k) * (monomial lift of R(B)_{d-k})), w_i in R_{d-5} with the
Y-parity pattern of the chosen D-parity.  N(v) = sum_b c_b N_b + A w, N_b a particular solution of N_b(F) = A K_b.
E2 mod J:  Phi(v) = D(v) div N(v) - N(v)(D(v))  in (R/J)_{2d};  every certificate (after Euler normalization) satisfies
Phi(v) = 0 and the coprimality of D with F fails only on the trivial family L = {c = 0, w = l z}.
Output: msolve input file relax-d{d}-{mode}-p{p}.ms with the independent quadrics; metadata JSON.
Usage: python3 b1_relax.py d mode p   (mode in even|odd|mixed)"""
import sys, json, time, os; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
d = int(sys.argv[1]); mode = sys.argv[2]; p = int(sys.argv[3]) if len(sys.argv) > 3 else PRIMES[1]
tag = f'relax-d{d}-{mode}-p{p}'
T0 = time.time()
def tick(msg): print(f'[{round(time.time()-T0,1)}s] {msg}', flush=True)
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Fi = [pderiv(Fp, i, p) for i in range(6)]
Bp = reduce_mod(B, p); Bi = [pderiv(Bp, 2, p), pderiv(Bp, 5, p)]
parities = {'even': [0], 'odd': [1], 'mixed': [0, 1]}[mode]

# ---- exact core kernel bases (characteristic 0), reduced mod p
exact_dims = {}
try:
    for r in map(json.loads, open('receipts/excl-kernels.jsonl')):
        if r.get('exp') == 'core-kernel-exact': exact_dims[(r['d'], r['parity'])] = r['dim_modJ']
except FileNotFoundError: pass
if os.path.exists('corekernel-exact.json'):
    for key, r in json.load(open('corekernel-exact.json')).items():
        k_, par_ = map(int, key.split(',')); exact_dims[(k_, par_)] = r['kernel_mod_J']
def core_kernel_basis(k, par):
    r = core_kernel(k, par, p, want_basis=True)
    assert r['certified_upper_bound'], (k, par)
    bas = r.get('basis', [])
    if (k, par) in exact_dims:
        assert len(bas) == exact_dims[(k, par)], ('mod-p kernel dim != exact dim', k, par, len(bas), exact_dims[(k, par)])
        tick(f'core kernel k={k} par={par}: dim {len(bas)} mod p = exact dim (certified equal)')
    else:
        tick(f'core kernel k={k} par={par}: dim {len(bas)} mod p (upper bound; exact dim not available)')
    return bas
# ---- monomial lifts of a basis of R(B)_j (both parities): greedy independent modulo J(B)_j
def RB_lift(j, parB):
    mons = [e for e in monomials(j, parB) if all(e[i] == 0 for i in (0, 1, 3, 4))]
    if j < 5: return mons
    Jcols = []
    for Bd, i in zip(Bi, (2, 5)):
        for m in monomials(j - 5, (parB - (0 if i == 2 else 1)) % 2):
            if all(m[t] == 0 for t in (0, 1, 3, 4)): Jcols.append(shift(Bd, m, 1, p))
    cur = list(Jcols); currank = build(cur, p)[0].rank() if cur else 0
    out = []
    for m in mons:
        Mt, _ = build(cur + [{m: 1}], p)
        if Mt.rank() > currank: out.append(m); cur.append({m: 1}); currank += 1
    assert len(out) == dim_RB(j, parB), (j, parB, len(out), dim_RB(j, parB))
    return out
# ---- unknowns
Dlist = []; labels = []
for parD in parities:
    for k in range(6, d + 1):
        for eps in (0, 1):
            cb = core_kernel_basis(k, eps)
            if not cb: continue
            j = d - k; parB = (parD - eps) % 2
            for mB in RB_lift(j, parB):
                for bi, b in enumerate(cb):
                    Dlist.append(shift(b, mB, 1, p)); labels.append(('c', parD, k, eps, mB, bi))
n_c = len(Dlist)
for parD in parities:
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(d - 5, (parD - pi) % 2):
            Dlist.append(shift(Fi[i], m, 1, p)); labels.append(('w', parD, i, m))
n = len(Dlist); n_w = n - n_c
tick(f'd={d} mode={mode} p={p}: unknowns {n} = {n_c} kernel + {n_w} Jacobian')
# ---- N_alpha: particular solutions for kernel elements (multi-RHS), A*e_i*m for Jacobian part
Nlist = [None] * n
for parD in parities:
    idxs = [a for a in range(n_c) if labels[a][1] == parD]
    if not idxs: continue
    tpar = (1 + parD) % 2
    cols = []; info = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (tpar - pi) % 2):
            cols.append(shift(Fi[i], m, 1, p)); info.append((i, m))
    rhs = [mul(Ap, Dlist[a], p) for a in idxs]
    M, rows = build(cols + rhs, p)
    R, rk = M.rref()
    ncol = len(cols)
    pivcols = []
    jj = 0
    for i in range(rk):
        while int(R[i, jj]) == 0: jj += 1
        pivcols.append(jj); jj += 1
    assert all(pc < ncol for pc in pivcols), 'E1 inconsistent for a kernel element (should not happen)'
    assert rk == dim_JF(d + 6, tpar), (rk, dim_JF(d + 6, tpar))
    for t, a in enumerate(idxs):
        N = [dict() for _ in range(6)]
        for i, pc in enumerate(pivcols):
            v = int(R[i, ncol + t])
            if v:
                ii, m = info[pc]; N[ii][m] = v
        Nlist[a] = N
    tick(f'E1 particular solutions for parity {parD}: {len(idxs)} rhs, rank {rk}')
for a in range(n_c, n):
    _, parD, i, m = labels[a]
    N = [dict() for _ in range(6)]; N[i] = shift(Ap, m, 1, p); Nlist[a] = N
# sanity: N_alpha(F) = A D_alpha
for a in (0, n_c - 1, n_c, n - 1):
    if a < 0 or a >= n: continue
    lhs = vecfield_apply(Nlist[a], Fp, p); rhs = mul(Ap, Dlist[a], p)
    assert lhs == rhs, a
# ---- normal former in degree 2d
par2 = None if mode == 'mixed' else 1   # E2 of a pure-parity D is Y-odd
nf = NormalFormer(2 * d, p, par=par2)
assert nf.rk == (dim_JF(2 * d, 0) + dim_JF(2 * d, 1) if par2 is None else dim_JF(2 * d, 1)), 'J rank mod p != Hilbert value'
n_std = len(nf.std)
tick(f'NormalFormer deg {2*d}: {len(nf.mons)} monomials, rank J = {nf.rk} (certified), std = {n_std}')
# ---- bilinear products via multiplication matrices
mons_d = monomials(d); idx_d = {m: i for i, m in enumerate(mons_d)}
mons_dm1 = monomials(d - 1); idx_dm1 = {m: i for i, m in enumerate(mons_dm1)}
mons2d = nf.mons; idx2d = nf.idx
def colmat(polys, idx):
    R, C = len(idx), len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[idx[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p)
Dm = colmat(Dlist, idx_d)                       # mons_d x n
dD = [colmat([pderiv(Dl, i, p) for Dl in Dlist], idx_dm1) for i in range(6)]   # mons_{d-1} x n
def multmat(G, src_mons, src_idx):
    """matrix of multiplication by G: R_{deg src} -> R_{2d} in bases src_mons -> mons2d (dense)"""
    R, C = len(mons2d), len(src_mons)
    flat = [0] * (R * C)
    for j, m in enumerate(src_mons):
        for e, v in G.items():
            s = tuple(x + y for x, y in zip(e, m))
            r = idx2d.get(s)
            if r is not None: flat[r * C + j] = v
    return nmod_mat(R, C, flat, p)
# accumulate NF rows for ordered pairs: store as dict beta -> nmod_mat (n x n_std) with row alpha = NF(P_{alpha,beta})
NFrows = []
for b in range(n):
    Nb = Nlist[b]
    divb = divergence(Nb, p)
    E = multmat(divb, mons_d, idx_d) * Dm if divb else nmod_mat(len(mons2d), n, p)
    for i in range(6):
        if Nb[i]:
            E = E - multmat(Nb[i], mons_dm1, idx_dm1) * dD[i]
    # E: mons2d x n ; row alpha of NF = NF(column alpha)
    Et = E.transpose()
    NFrows.append(nf.nf_matrix(Et))
    if b % 10 == 0: tick(f'products beta={b}/{n}')
tick('all bilinear products reduced mod J')
# ---- quadrics: for each std monomial s, Q_s(v) = sum_{a<=b} q_{ab,s} v_a v_b
pairs = [(a, b) for a in range(n) for b in range(a, n)]
pidx = {ab: t for t, ab in enumerate(pairs)}
npairs = len(pairs)
# matrix QM: rows = std monomials, cols = pairs
flat = [0] * (n_std * npairs)
for b in range(n):
    M = NFrows[b]
    for a in range(n):
        t = pidx[(a, b)] if a <= b else pidx[(b, a)]
        for s in range(n_std):
            v = int(M[a, s])
            if v: flat[s * npairs + t] = (flat[s * npairs + t] + v) % p
QM = nmod_mat(n_std, npairs, flat, p)
Qr, qrk = QM.rref()
tick(f'quadric space: {n_std} std monomials, {npairs} pair-monomials, independent quadrics {qrk}')
# ---- write msolve input
vars_ = [f'v{a}' for a in range(n)]
with open(tag + '.ms', 'w') as f:
    f.write(','.join(vars_) + '\n'); f.write(f'{p}\n')
    lines = []
    for i in range(qrk):
        terms = []
        for t, (a, b) in enumerate(pairs):
            c = int(Qr[i, t])
            if c:
                terms.append(f'{c}*v{a}*v{b}' if a != b else f'{c}*v{a}^2')
        lines.append('+'.join(terms))
    f.write(',\n'.join(lines) + '\n')
# ---- Singular input for the same system (independent GB engine)
with open(tag + '.sing', 'w') as f:
    f.write(f'ring r = {p}, ({",".join(vars_)}), dp;\n')
    f.write('ideal I;\n')
    for i in range(qrk):
        terms = []
        for t, (a, b) in enumerate(pairs):
            c = int(Qr[i, t])
            if c: terms.append(f'{c}*v{a}*v{b}')
        f.write(f'I[{i+1}] = ' + '+'.join(terms) + ';\n')
    f.write('option(redSB); ideal G = std(I); G = std(G); print(size(G)); print(dim(G)); print(mult(G)); print(hilb(G)); write(":w ' + tag + '.gb.sing", G);\n$\n')
meta = dict(exp='b1-relax-build', d=d, mode=mode, p=p, n_unknowns=n, n_kernel=n_c, n_jacobian=n_w, labels=[[str(x) for x in l] for l in labels],
            n_std=n_std, n_pairs=npairs, n_quadrics=qrk, jacobian_rank_certified=True,
            trivial_family='c=0, w_i = l * z_i (Y-even part), l in R_{d-6}; projective dimension dim R_{d-6} - 1',
            seconds=round(time.time() - T0, 1))
json.dump(meta, open(tag + '.meta.json', 'w'), indent=1)
log({k: v for k, v in meta.items() if k != 'labels'})
