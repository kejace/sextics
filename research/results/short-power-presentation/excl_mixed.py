"""Wave-2 (excl lens): mixed-parity exclusion by the "annihilator functional" argument.

Setting: D = D^+ + D^- mixed, both parts nonzero (pure parts are handled by the pure relaxations).  With
x = (lam^+, w^+), y = (lam^-, w^-), the relaxed E2 condition splits as
    Xi^+(x) + Xi^-(y) in J_{2d}^-      (odd part, no cross terms)      and      Xi^x(x, y) in J_{2d}^+.
Let L^- = span of the coordinate polynomials of Xi^- modulo J_{2d}^-.  Every linear functional phi on S_{2d}^-
vanishing on J_{2d}^- and on L^- satisfies phi(Xi^+(x)) = 0 for every mixed solution.  This yields
dim L^+ - dim(L^+ ∩ L^-) independent quadrics in x alone; if their ideal is m-primary (mod a prime for which
rank_p [J ; L^- polys] = rank_Q [J ; L^- polys], checked exactly with fmpz_mat), then x = 0 for every mixed
solution over Q-bar, hence Xi^-(y) in J and y = 0 by the pure odd result.  Symmetrically with x <-> y.

Usage: python3 excl_mixed.py <d> [--pidx k] [--threads T] [--no-exact]
"""
import sys, time, json, os, subprocess, math, re
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
import excl_relax as R
from hodge import core
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, reduce_mod
from flint import nmod_mat, fmpz_mat
from fractions import Fraction

d = int(sys.argv[1])
pidx = int(sys.argv[sys.argv.index('--pidx') + 1]) if '--pidx' in sys.argv else 0
threads = sys.argv[sys.argv.index('--threads') + 1] if '--threads' in sys.argv else '4'
p = W.PRIMES30[pidx]
job = os.environ.get('W2JOB', f'excl_mixed_d{d}')
F, A, q, rec = W.load()
Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
t00 = time.time()
e = 2 * d
# ---- kernel elements (exact if available from excl_kernel_D7.json, else mod p only) -------------------------
kern = {0: [], 1: []}
exact_ok = True
if d >= 7:
    path = os.path.join(W.RES, 'excl_kernel_D7.json')
    if d == 7 and os.path.exists(path):
        ex = json.load(open(path))
        for par in [0, 1]:
            kern[par].append((core.parse(ex[str(par)]['D']), [core.parse(s) for s in ex[str(par)]['N']]))
    else:
        exact_ok = False
        for par in [0, 1]:
            info, basis = W.e1_kernel(F, A, d, par, p, want_basis=True)
            for Da in basis or []:
                Na, _ = W.membership_mod_p(Fp, mul(Ap, Da, p), p)
                kern[par].append((Da, Na))
print(f'kernel elements: even {len(kern[0])}, odd {len(kern[1])} (exact: {exact_ok})', flush=True)

# ---- pair polynomials of Xi^+ (x-pairs) and Xi^- (y-pairs), exact when possible ---------------------------
def pair_polys(par, modp):
    """list of (key, poly) for the pure-parity quadratic map Xi^{par}; keys ('ll',a,b), ('lw',a,(j,m)), ('ww',(j,m),(l,m2))."""
    kk = kern[par]
    if modp:
        Fx, Ax = Fp, Ap
        kk = [(reduce_mod(Da, p), [reduce_mod(n, p) for n in Na]) for Da, Na in kk]
        pp = p
    else:
        Fx, Ax, pp = F, A, None
    Fi = [pderiv(Fx, i, pp) for i in range(NV)]
    H = [[pderiv(Fi[j], l, pp) for l in range(NV)] for j in range(NV)]
    Ai = [pderiv(Ax, i, pp) for i in range(NV)]
    out = []
    divN = []
    for Da, Na in kk:
        dv = {}
        for i in range(NV):
            dv = add(dv, pderiv(Na[i], i, pp), pp)
        divN.append(dv)
    Dg = [[pderiv(Da, i, pp) for i in range(NV)] for Da, _ in kk]
    for a, (Da, Na) in enumerate(kk):
        for b, (Db, Nb) in enumerate(kk):
            B = mul(Da, divN[b], pp)
            for i in range(NV):
                B = sub(B, mul(Nb[i], Dg[a][i], pp), pp)
            out.append((('ll', min(a, b), max(a, b)), B))
    fb = [(j, m) for (j, m) in W.field_basis(d, par) if not (j == 0 and m[0] > 0)]   # Euler slice
    for a, (Da, Na) in enumerate(kk):
        for (j, m) in fb:
            HNj = {}
            for i in range(NV):
                HNj = add(HNj, mul(Na[i], H[i][j], pp), pp)
            P = shift(sub(sub(mul(Da, Ai[j], pp), mul(Ax, Dg[a][j], pp), pp), HNj, pp), m, 1, pp)
            out.append((('lw', a, (j, m)), P))
    AH = [[mul(Ax, H[j][l], pp) for l in range(NV)] for j in range(NV)]
    for u in range(len(fb)):
        j, m = fb[u]
        for v in range(u, len(fb)):
            l, m2 = fb[v]
            mm = tuple(a_ + b_ for a_, b_ in zip(m, m2))
            c = -1 if u == v else -2
            P = shift(AH[j][l], mm, (c % pp) if pp else c, pp)
            out.append((('ww', (j, m), (l, m2)), P))
    # merge duplicate keys (ll a,b and b,a)
    merged = {}
    for k, P in out:
        merged[k] = add(merged.get(k, {}), P, pp)
    return [(k, P) for k, P in merged.items() if P]

mons = monomials(e, 1); idx = {m: i for i, m in enumerate(mons)}
n = len(mons)
gens_p = [g[2] for g in W.jac_gens(Fp, e, 1, p)]
rJ_expected = W.dimJ(e, 1)
results = {}
for side, par_x, par_y in [('x', 0, 1), ('y', 1, 0)]:
    t0 = time.time()
    PX = pair_polys(par_x, True)
    PY = pair_polys(par_y, True)
    # only odd-parity parts matter (Xi^± live in S_{2d}^-); sanity: all pair polys here are odd
    for k, P in PX + PY:
        assert all(core.yparity(mo) == 1 for mo in P), 'parity bookkeeping error'
    MJ = W.rows_matrix(gens_p, idx, p)
    MY = W.rows_matrix([P for _, P in PY], idx, p)
    MJY = W.rows_matrix(gens_p + [P for _, P in PY], idx, p)
    rJ = MJ.rank(); rJY = MJY.rank()
    assert rJ == rJ_expected, (rJ, rJ_expected)
    # annihilator of [J ; Y-pairs] : right nullspace
    Phi, nphi = MJY.nullspace()        # n x nphi (columns = functionals), nphi = n - rJY
    Phi = nmod_mat(n, nphi, [int(Phi[i, k]) for i in range(n) for k in range(nphi)], p)
    MX = W.rows_matrix([P for _, P in PX], idx, p)
    Qm = MX * Phi                      # (#x-pairs) x nphi : column k = quadric phi_k(Xi^+)
    QT = Qm.transpose()
    Qr, rq = QT.rref()                 # rq independent quadrics in x
    dimLx = W.rows_matrix(gens_p + [P for _, P in PX], idx, p).rank() - rJ
    dimLy = rJY - rJ
    keysX = [k for k, _ in PX]
    # unknown indexing for the x-side
    unk = []
    for k in keysX:
        for part in (k[1:] if k[0] != 'ww' else k[1:]):
            pass
    names = {}
    def uname(t):
        if t not in names:
            names[t] = len(names)
        return names[t]
    eqs = []
    for i in range(rq):
        eq = {}
        for c, k in enumerate(keysX):
            v = int(Qr[i, c])
            if not v:
                continue
            if k[0] == 'll':
                a, b = uname(('lam', k[1])), uname(('lam', k[2]))
            elif k[0] == 'lw':
                a, b = uname(('lam', k[1])), uname(('w',) + k[2])
            else:
                a, b = uname(('w',) + k[1]), uname(('w',) + k[2])
            eq[(min(a, b), max(a, b))] = v
        eqs.append(eq)
    nu = len(names)
    allx = [('lam', a) for a in range(len(kern[par_x]))] + [('w',) + jm for jm in [(j, m) for (j, m) in W.field_basis(d, par_x) if not (j == 0 and m[0] > 0)]]
    invisible = [u for u in allx if u not in names]
    results.setdefault('_invisible', {})[side] = invisible
    print(f'[{side}] visible unknowns {nu} of {len(allx)}; invisible: {invisible}', flush=True)
    sysd = {'unknowns': [None] * nu, 'eqs': eqs}
    base = os.path.join(W.RES, f'excl_mixed_d{d}_{side}_p{pidx}')
    W.write_msolve(sysd, base + '.ms', p)
    W.write_singular(sysd, base + '.sing', p)
    print(f'[{side}] dim L_{par_x} = {dimLx}, dim L_{par_y} = {dimLy}, functionals {nphi}, independent quadrics in {side}: {rq} in {nu} unknowns [{time.time()-t0:.0f}s]', flush=True)
    ms = R.run_msolve(base + '.ms', base + '.ms.out', threads)
    mons_lead = R.parse_msolve_leading(base + '.ms.out')
    degs = [sum(k for _, k in m) for m in mons_lead]
    maxdeg = max(degs) if degs else 0
    tmax = 0
    for t in range(1, maxdeg + 4):
        cost = sum(math.comb(t - dm + nu - 1, nu - 1) for dm in degs if dm <= t)
        if cost > 3e7:
            break
        tmax = t
    h, _ = R.hilbert_from_leading(mons_lead, nu, tmax) if tmax >= 1 else (None, 0)
    sg = R.run_singular(base + '.sing', base + '.sing.out', nu)
    empty = (h is not None and 0 in h[1:]) or str(sg.get('dim', '')) == '0'
    results[side] = {'dimLx': dimLx, 'dimLy': dimLy, 'nphi': nphi, 'quadrics': rq, 'unknowns': nu, 'msolve_hilbert': h,
                     'singular_dim': sg.get('dim'), 'singular_vdim': sg.get('vdim'), 'singular_hilbert': sg.get('hilbert_function'),
                     'proj_empty': bool(empty), 'rank_JY_modp': rJY, 'seconds': round(time.time() - t0, 1)}
    print(f'[{side}] msolve h = {h}; singular dim = {sg.get("dim")} vdim = {sg.get("vdim")} -> proj_empty = {empty}', flush=True)
    W.wlog(exp='mixed-annihilator-relaxation', d=d, pars='m', side=side, prime=p, result='proj-empty' if empty else 'not-empty',
           **{k: v for k, v in results[side].items() if k not in ('proj_empty',)}, job=job,
           note=f'quadrics phi(Xi^{side}) for functionals phi killing J_{e}^- and the image span of the other parity; proj-empty => mixed solutions have {side}=0')
    if empty and '--no-exact' not in sys.argv and exact_ok:
        # exact certification: dim_Q L_{par_y} == dim_p L_{par_y}  (the functionals used vanish on J + L_{par_y})
        t1 = time.time()
        PYx = pair_polys(par_y, False)
        dep = W.exact_dependencies(F, [P for _, P in PYx], e, 1, primes=W.PRIMES31, max_primes=5)
        results[side]['exact_dependencies'] = dep
        W.wlog(exp='mixed-annihilator-exact-certificate', d=d, pars='m', side=side, prime=dep.get('primes'),
               result='exact' if dep.get('exact_equals_modp') else 'INCOMPLETE', **{k: v for k, v in dep.items() if k != 'primes'}, job=job,
               note='exact rational dependencies among the image polynomials of the other parity modulo J (verified with fmpq_mpoly); exact => the annihilator functionals are reductions of rational ones => mixed exclusion rigorous')
        print(f'[{side}] exact dependency certificate: {dep}', flush=True)

print(json.dumps(results, default=str)[:3000], flush=True)
print('done', round(time.time() - t00), 's')
