"""Non-Jacobian E1-kernel of the K3-sector target and fixed-D tests for it.

A_K3 (XY frame) = [(X0+Y0)^2(X1+Y1)^2(X2+Y2)^2 - (X0-Y0)^2(X1-Y1)^2(X2-Y2)^2]/64 = K3-sector part of A.
For the full A, [A] is injective on (R/J)_d for d<=9 (E1 kernel = J).  For A_K3 = eta(x)1 - 1(x)eta the kernel
is bigger (eta: R_{f,a} -> R_{f,a+6} has a kernel for a>=4), so non-Jacobian E1-admissible D exist from degree 8 on.
This script: (1) kernel dimensions of [A_K3] on (R/J)_d for d in DEGS, both parities (engine jacobian.mult_kernel);
(2) structured elements: ker(eta: R_4->R_10) =: K4 in the ternary ring; D_sym = sum_k g_k(U) g_k(V), D_alt = sum_k (-1)^k ...,
    and rank-1 g(U)g'(V) (reducible control); (3) random kernel elements; (4) fixed-D solve for target A_K3 (compressed mod p),
    plus the E1 check for the full A (expected to fail: not in J). Results appended to k3kernel-<tag>.jsonl.
usage: k3kernel.py <tag> <d> <nrandom> [p]
"""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from hodge import core, jacobian, fixed_d
from fractions import Fraction as Fr
tag = sys.argv[1]; d = int(sys.argv[2]); nrand = int(sys.argv[3]); p = int(sys.argv[4]) if len(sys.argv) > 4 else core.PRIMES[0]
OUT = f'k3kernel-{tag}.jsonl'
F, A, q, rec = core.load_instance('polarized-k3-remainder')
def var(i):
    e = [0]*6; e[i] = 1; return {tuple(e): Fr(1)}
X = [var(i) for i in range(3)]; Y = [var(3+i) for i in range(3)]
U = [core.add(X[i], Y[i]) for i in range(3)]; V = [core.sub(X[i], Y[i]) for i in range(3)]
def prodsq(L):
    out = {(0,)*6: Fr(1)}
    for l in L: out = core.mul(out, core.mul(l, l))
    return out
AK3 = core.scale(core.sub(prodsq(U), prodsq(V)), Fr(1, 64))
def uv_to_xy(P):
    imgs = U + V
    out = {}
    for e, c in P.items():
        term = {(0,)*6: Fr(c)}
        for i in range(6):
            if e[i]: term = core.mul(term, core.power(imgs[i], e[i]))
        out = core.add(out, term)
    return out
def emit(rec):
    rec['time'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['tag'] = tag
    with open(OUT, 'a') as fh: fh.write(json.dumps(rec) + '\n')
    print(json.dumps(rec)[:400], flush=True)

def test_D(label, D, family):
    t0 = time.time()
    par = core.parity_of(D)
    if par is None:
        emit({'exp': 'fixedD-AK3', 'label': label, 'family': family, 'd': core.degree(D), 'result': 'skipped (mixed parity)'}); return
    okK, infoK = fixed_d.e1_solvable(F, AK3, q, D, p=p)
    okA, infoA = fixed_d.e1_solvable(F, A, q, D, p=p)
    inJ, infoJ = jacobian.in_jacobian(F, D, p=p)
    base = {'exp': 'fixedD-AK3', 'label': label, 'family': family, 'd': core.degree(D), 'parity': par, 'prime': p,
            'D_in_J': bool(inJ), 'e1_AK3': bool(okK), 'e1_A': bool(okA), 'D_xy': core.to_str(D) if len(D) <= 60 else f'{len(D)} terms'}
    if not okK:
        base.update({'result': 'E1 fails for A_K3', 'seconds': round(time.time()-t0, 1)}); emit(base); return
    out = fixed_d.solve_fixed_D_fast(F, AK3, q, D, p=p, verbose=False, want_exact=True)
    base.update({'result': 'CONSISTENT (K3 target served!)' if out['consistent'] else 'inconsistent mod p',
                 'rank': out['rank'], 'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
                 'exact': out.get('exact'), 'seconds': round(time.time()-t0, 1)})
    if out['consistent'] and out.get('exact'):
        json.dump({'d': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q']),
                   'target': 'A_K3'}, open(f'WITNESS-AK3-{tag}-{label}.json', 'w'))
    emit(base)

# (1) kernel dimensions
kers = {}
for par in (0, 1):
    t0 = time.time()
    k = jacobian.mult_kernel(F, AK3, d, p=p, parity=par)
    kers[par] = k
    emit({'exp': 'e1kernel-AK3', 'd': d, 'parity': par, 'prime': p, 'ring_dim': k['ring_dim'], 'jacobian_dim': k['jacobian_dim'],
          'kernel_mod_J': k['dim'], 'kernel_total': k['kernel_dim'], 'seconds': round(time.time()-t0, 1)})
# (2) structured: K4 = ker(eta: R_4 -> R_10) in the ternary Jacobian ring (U,V frame tensor tools)
if d == 8:
    from uvring import TernaryJacobianRing, load_f, tern_monomials, tmul
    from flint import nmod_mat
    f = load_f(p); R = TernaryJacobianRing(f, p=p, maxdeg=12)
    eta = {(2, 2, 2): 1}
    m4 = tern_monomials(4); vecs = [R.nf(tmul({m: 1}, eta, p), 10) for m in m4]
    M = nmod_mat(len(m4), R.dim(10), [x for v in vecs for x in v], p)
    K, nul = M.transpose().nullspace()   # vectors in the monomial basis of S_4 (15 monomials) with eta*g in J(f)
    # K is (15 x nul): columns = kernel vectors g (coefficients on m4)
    gs = []
    for k in range(nul):
        g = {m4[i]: int(K[i, k]) for i in range(len(m4)) if int(K[i, k])}
        gs.append(g)
    emit({'exp': 'ker-eta-R4', 'prime': p, 'dim_ker_eta_S4_to_R10': nul, 'note': 'quartics g with eta*g in J(f)_10 (includes J(f)_4 = 0, so = ker on R_4)'})
    # lift to small integers via symmetric residues (these are mod-p vectors; we use them as mod-p candidates after
    # rational reconstruction when possible)
    def lift(g):
        out = {}
        for m, v in g.items():
            r = core.ratrecon(v, p)
            if r is None: return None
            out[m] = r
        return out
    gl = [lift(g) for g in gs]
    okl = [g for g in gl if g is not None]
    emit({'exp': 'ker-eta-R4-lift', 'lifted': len(okl), 'of': len(gl)})
    def toU(g): return {(m[0], m[1], m[2], 0, 0, 0): c for m, c in g.items()}
    def toV(g): return {(0, 0, 0, m[0], m[1], m[2]): c for m, c in g.items()}
    if okl:
        Dsym = {}; Dalt = {}
        for k, g in enumerate(okl):
            t = core.mul(toU(g), toV(g))
            Dsym = core.add(Dsym, t); Dalt = core.add(Dalt, core.scale(t, (-1)**k))
        test_D('K4-sym-pairing', uv_to_xy(Dsym), 'ker(eta)_4 (x) ker(eta)_4 sum g_k(U)g_k(V)')
        test_D('K4-alt-pairing', uv_to_xy(Dalt), 'ker(eta)_4 (x) ker(eta)_4 alternating')
        if len(okl) >= 2:
            g0, g1 = okl[0], okl[1]
            D01 = core.add(core.mul(toU(g0), toV(g1)), core.mul(toU(g1), toV(g0)))
            test_D('K4-g0g1+g1g0', uv_to_xy(D01), 'symmetrized rank-2 tensor in ker(eta)_4 (x) ker(eta)_4')
            test_D('K4-g0g0', uv_to_xy(core.mul(toU(g0), toV(g0))), 'rank-1 g0(U)g0(V) (reducible control)')
# (3) random kernel elements (mod p representatives in standard monomials; lifted to small integers when possible)
random.seed(12345)
for par in (0, 1):
    basis = kers[par]['basis']
    if not basis: continue
    for t in range(nrand):
        coeffs = [random.randrange(-3, 4) for _ in basis]
        D = {}
        for c, b in zip(coeffs, basis):
            if c: D = core.add(D, core.scale({e: (v if v <= p//2 else v - p) for e, v in b.items()}, c))
        D = {e: Fr(v) for e, v in D.items() if v}
        test_D(f'random-ker-par{par}-{t}', D, f'random small-integer combination of {len(basis)} kernel reps (mod p lifted symmetric)')
emit({'exp': 'done', 'd': d})
