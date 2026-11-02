"""(b) Rigorous lower bound on deg D via the Jacobian-type reduction + Hessian relaxation.

For an E1-admissible denominator D = sum_i w_i F_i + t K (K = non-Jacobian E1-kernel element with a fixed
E1 solution N^K, A K = sum N^K_i F_i), every E1 solution is N_i = A w_i - H z_i/6 + Koszul + t N^K_i, and the
E2 condition  F | E/F  reduces modulo J(F) to the quadratic condition (all H / Koszul terms drop out)

   Phi(w,t) := t K (w.grad A + A div w) + t^2 (K div N^K - N^K(K)) - A Hess_F(w,w) - t (A w(K) + sum_ij N^K_i w_j F_ij)  in J(F)_{2d}

(w(P) = sum_i w_i d_i P).  The coordinates of Phi in R_{2d} = S/J are quadrics in (w,t); if they have no common
zero in P(W + <K>) over C, no denominator of degree d (of this parity) exists.  Mod-p ranks are lower bounds
for rational ranks, so (i) full rank of the quadric coefficient matrix mod p, or (ii) HF_p(k) = 0 for some k
(Groebner / Macaulay), are rigorous over Q and hence over C, given rank J_k mod p = CI value for k in {d, d+6, 2d}
(checked in a_verify.py, receipts a_rankJ.json).

Usage: python3 b_hessian.py d par [--fermat] [--prime P] [--kernel D7|D7Y|none]
"""
import sys, time, json, argparse, itertools, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual/wave2-excl')
from xlib import *
from hodge import jacobian, fixed_d
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('par', type=int)
ap.add_argument('--fermat', action='store_true'); ap.add_argument('--prime', type=int, default=P2)
ap.add_argument('--kernel', default='auto'); ap.add_argument('--groebner', default='msolve')
ap.add_argument('--threads', type=int, default=4); ap.add_argument('--selftest', action='store_true')
args = ap.parse_args()
d, par, p = args.d, args.par, args.prime
t00 = time.time()
def say(*a): print(f'[{time.time()-t00:7.1f}s]', *a, flush=True)

F, A, q = load()
if args.fermat:
    F = core.add(p_k(X, 6), p_k(Y, 6))
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
Fij = [[core.pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
dA = [core.pderiv(Ap, i, p) for i in range(6)]

# ---- Jacobian basis of J_d^{par}: u_a = m_a F_{i_a}
basis = []
for i in range(6):
    par_m = par if i < 3 else (par + 1) % 2
    for m in core.monomials(d - 5, par_m):
        basis.append((i, m))
n = len(basis)
# ---- kernel element
Kpoly = None; NK = None
if args.kernel == 'auto':
    args.kernel = {(7, 0): 'D7', (7, 1): 'D7Y'}.get((d, par), 'none')
if args.kernel != 'none' and not args.fermat:
    DX = core.mul(core.mul(core.sub(X[0], X[1]), core.sub(X[0], X[2])), core.sub(X[1], X[2]))
    DY = core.mul(core.mul(core.sub(Y[0], Y[1]), core.sub(Y[0], Y[2])), core.sub(Y[1], Y[2]))
    X2 = [core.power(v, 2) for v in X]; Y2 = [core.power(v, 2) for v in Y]
    Kpoly = core.reduce_mod(core.mul(DX, e_k(X2, 2)) if args.kernel == 'D7' else core.mul(DY, e_k(Y2, 2)), p)
    # E1 solution N^K with H = 0: A K = sum N^K_i F_i  (mod p), via rref of [J_{d+6} | A K]
    M, mons, nJ = jac_matrix(F, d + 6, (par + 1) % 2, p, extra_cols=[core.mul(Ap, Kpoly, p)])
    R, rk_ = M.rref()
    colinfo = []
    for i in range(6):
        pm = (par + 1) % 2 if i < 3 else par
        for m in core.monomials(d + 1, pm):
            colinfo.append((i, m))
    assert len(colinfo) == nJ
    NK = [{} for _ in range(6)]
    for r in range(rk_):
        j = 0
        while j < M.ncols() and int(R[r, j]) == 0: j += 1
        assert j < nJ, 'A*K not in J mod p: not an E1-kernel element'
        v = int(R[r, nJ])
        if v:
            i, m = colinfo[j]
            NK[i] = core.add(NK[i], {m: v}, p)
    chk = {}
    for i in range(6): chk = core.add(chk, core.mul(NK[i], Fi[i], p), p)
    assert core.sub(chk, core.mul(Ap, Kpoly, p), p) == {}, 'E1 solution check failed'
    say('kernel element', args.kernel, 'E1 solution found; nonzero N_i:', [len(x) for x in NK])
nv = n + (1 if Kpoly is not None else 0)
say(f'd={d} par={par} p={p} fermat={args.fermat}: n_J={n} n_vars={nv}')

# ---- reducer for R_{2d}^{odd}: left nullspace Y of J_{2d}^{odd}
deg2 = 2 * d
MJ, mons2, nJ2 = jac_matrix(F, deg2, 1, p)
ring_dim = len(mons2) - MJ.rank()
assert ring_dim == ci_ring_dim(deg2, 1), ('rank J_2d mod p not CI', ring_dim, ci_ring_dim(deg2, 1))
Xn, nullity = MJ.transpose().nullspace()
assert nullity == ring_dim
Yl = nmod_mat(ring_dim, len(mons2), p)
for k in range(ring_dim):
    for i in range(len(mons2)):
        v = int(Xn[i, k])
        if v: Yl[k, i] = v
say(f'R_{deg2}^odd dim {ring_dim}; left nullspace built')
codes2 = mono_codes(mons2); order2 = np.argsort(codes2); sc2 = codes2[order2]
def col_of(P):
    """row-index/value arrays of poly P in the S_{2d}^{odd} basis"""
    if not P: return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    E, C = poly_to_arrays(P, p)
    cc = (E * (32 ** np.arange(6, dtype=np.int64))).sum(axis=1)
    pos = np.searchsorted(sc2, cc)
    assert np.all(sc2[np.minimum(pos, len(sc2) - 1)] == cc)
    return order2[pos], C

# ---- quadric coefficient vectors: unordered pairs of variables 0..nv-1 (index nv-1 = t)
pairs = [(a, b) for a in range(nv) for b in range(a, nv)]
pidx = {pr: k for k, pr in enumerate(pairs)}
V = nmod_mat(len(mons2), len(pairs), p)
def addcol(k, P, scal=1):
    ri, vals = col_of(P)
    for r, v in zip(ri.tolist(), vals.tolist()):
        V[int(r), k] = (int(V[int(r), k]) + scal * int(v)) % p
t0 = time.time()
# -A Hess(w,w): ordered pairs (a,b) -> unordered coefficient: a==b: 1, a<b: 2
for a in range(n):
    ia, ma = basis[a]
    for b in range(a, n):
        ib, mb = basis[b]
        mm = tuple(x + y for x, y in zip(ma, mb))
        P = core.shift(Fij[ia][ib], mm, 1, p)
        P = core.mul(Ap, P, p)
        addcol(pidx[(a, b)], P, scal=(p - 1) * (1 if a == b else 2))
say('Hessian block built', round(time.time() - t0, 1), 's')
if Kpoly is not None:
    tI = nv - 1
    divNK = {}
    for i in range(6): divNK = core.add(divNK, core.pderiv(NK[i], i, p), p)
    NKK = {}
    for i in range(6): NKK = core.add(NKK, core.mul(NK[i], core.pderiv(Kpoly, i, p), p), p)
    ctt = core.sub(core.mul(Kpoly, divNK, p), NKK, p)
    addcol(pidx[(tI, tI)], ctt, 1)
    for a in range(n):
        ia, ma = basis[a]
        # t * w_a coefficient: K (m_a dA_ia + A d_ia m_a) - (A m_a d_ia K + m_a sum_i NK_i F_{i,ia})
        term = core.add(core.shift(dA[ia], ma, 1, p), core.mul(Ap, core.pderiv({ma: 1}, ia, p), p), p)
        c1 = core.mul(Kpoly, term, p)
        c2 = core.mul(Ap, core.shift(core.pderiv(Kpoly, ia, p), ma, 1, p), p)
        s = {}
        for i in range(6): s = core.add(s, core.mul(NK[i], Fij[i][ia], p), p)
        c2 = core.add(c2, core.shift(s, ma, 1, p), p)
        addcol(pidx[(a, tI)], core.sub(c1, c2, p), 1)
    say('kernel cross terms built')
Qm = Yl * V          # ring_dim x npairs
rQ = Qm.rank()
say(f'quadric coefficient matrix {ring_dim} x {len(pairs)}: rank {rQ} (full = {len(pairs)})')
rec = dict(experiment='b:hessian-relaxation', d=d, par=par, p=p, fermat=args.fermat, kernel=args.kernel, n_jacobian=n, n_vars=nv,
           ring_dim_2d=ring_dim, n_pairs=len(pairs), rank_quadrics=rQ, full_rank=(rQ == len(pairs)))

# ---- self test on Fermat: w = X_i/6 (D = p6X) must satisfy all quadrics
def eval_quadrics(wvec):
    val = nmod_mat(len(pairs), 1, p)
    for k, (a, b) in enumerate(pairs):
        val[k, 0] = (wvec[a] * wvec[b]) % p
    return Qm * val
if args.fermat and d == 6 and par == 0:
    inv6 = pow(6, -1, p)
    wv = [0] * nv
    for a, (i, m) in enumerate(basis):
        if i < 3 and m[i] == 1: wv[a] = inv6
    r = eval_quadrics(wv)
    z = all(int(r[k, 0]) == 0 for k in range(ring_dim))
    say('Fermat self-test: p6X direction satisfies all quadrics:', z)
    rec['selftest_p6X_in_V'] = z
    # random direction: full system vs relaxation
    rng = np.random.default_rng(0)
    for trial in range(2):
        wv = [int(x) for x in rng.integers(0, p, nv)]
        Dp = {}
        for a, (i, m) in enumerate(basis): Dp = core.add(Dp, core.shift(Fi[i], m, wv[a], p), p)
        r = eval_quadrics(wv); relax_ok = all(int(r[k, 0]) == 0 for k in range(ring_dim))
        Dq = {e: sym_int(c, p) for e, c in Dp.items()}
        res = fixed_d.solve_fixed_D_fast(F, A, q, Dq, p=p, want_exact=False)
        say(f'Fermat random w trial {trial}: relaxation holds={relax_ok}, full system consistent={res["consistent"]}')
        rec[f'selftest_random_{trial}'] = dict(relax=relax_ok, full=res['consistent'])
log(**rec)
receipt(f'b_hessian_d{d}_p{par}' + ('_fermat' if args.fermat else ''), rec)
if rQ == len(pairs):
    say('EXCLUDED: quadrics span Sym^2 => only w = 0; no denominator of degree', d, 'parity', par)
    sys.exit(0)

# ---- Groebner step: write quadrics, run msolve / Singular, test Krull dimension 0 (projective emptiness)
names = [f'w{a}' for a in range(n)] + (['t'] if Kpoly is not None else [])
# extract a basis of the row space of Qm (rank rQ) to reduce the number of quadrics
Qr, rk2 = Qm.rref()
polys = []
for r in range(rk2):
    terms = []
    for k, (a, b) in enumerate(pairs):
        c = int(Qr[r, k])
        if c:
            terms.append(f'{c}*{names[a]}*{names[b]}')
    polys.append('+'.join(terms))
tag = f'd{d}_p{par}' + ('_fermat' if args.fermat else '')
with open(f'{PK}/b_quadrics_{tag}.ms', 'w') as f:
    f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
with open(f'{PK}/b_quadrics_{tag}.sing', 'w') as f:
    f.write(f'ring r = {p}, ({",".join(names)}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n'
            'option(redSB);\nideal G = std(I);\nprint("KRULLDIM"); print(dim(G)); print("HILB"); print(hilb(G,1)); quit;\n')
say(f'wrote {len(polys)} independent quadrics in {nv} variables')
if args.groebner == 'singular':
    t0 = time.time()
    out = subprocess.run(['Singular', '-q', f'{PK}/b_quadrics_{tag}.sing'], capture_output=True, text=True, timeout=36000)
    say('Singular output:', out.stdout[-2000:], out.stderr[-500:])
    kd = None
    lines = out.stdout.strip().split('\n')
    if 'KRULLDIM' in lines: kd = int(lines[lines.index('KRULLDIM') + 1])
    rec2 = dict(experiment='b:groebner', d=d, par=par, p=p, kernel=args.kernel, tool='Singular', krull_dim=kd, seconds=round(time.time() - t0, 1))
    log(**rec2); receipt(f'b_groebner_{tag}', dict(rec2, stdout=out.stdout[-5000:]))
    say('EXCLUDED (Krull dim 0 => V_proj empty)' if kd == 0 else f'NOT excluded by relaxation: Krull dim {kd}')
else:
    t0 = time.time()
    outp = f'{PK}/b_msolve_{tag}.out'
    cmd = ['msolve', '-v', '1', '-t', str(args.threads), '-f', f'{PK}/b_quadrics_{tag}.ms', '-o', outp]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=36000)
    txt = open(outp).read() if os.path.exists(outp) else ''
    say('msolve stderr tail:', out.stderr[-1500:]); say('msolve output head:', txt[:500])
    rec2 = dict(experiment='b:groebner', d=d, par=par, p=p, kernel=args.kernel, tool='msolve', seconds=round(time.time() - t0, 1), output_head=txt[:300])
    log(**rec2); receipt(f'b_groebner_{tag}', dict(rec2, stderr=out.stderr[-5000:]))
