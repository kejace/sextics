"""Self-test of w2lib on accepted witnesses of the shape F = f(x) + c Y2^6, A free of Y2, D free of Y2.

1. reduce the accepted 6-variable witness to a 5-variable solution (necessity direction) and check (*) exactly;
2. check the Operator5 assembly (exact residual ~ 0 at the reduced D), the regularised inner solve, the gradient;
3. check that D lies in the float E1 kernel, then run a few varpro restarts and try to recover a solution exactly.
"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import numpy as np
from fractions import Fraction
import w2lib as W
from w2lib import core, mul, add, sub, scale, pderiv, shift

def reduce_witness(prob, D6, N6):
    """u0 = average over mu6(Y2) of N6/D6, then Y2^6 -> -f/c.  Requires D6 free of Y2."""
    assert all(e[5] == 0 for e in D6)
    N5 = []
    for i in range(5):
        acc = {}
        for e, v in N6[i].items():
            if e[5] % 6:
                continue
            j = e[5] // 6
            mon = {e[:5] + (0,): Fraction(v)}
            acc = add(acc, mul(mon, core.power(scale(prob.f, -1 / prob.c), j)))
        N5.append(acc)
    return N5

def identity5(prob, D, N, T):
    divN = {}; ND = {}
    for i in range(5):
        divN = add(divN, pderiv(N[i], i)); ND = add(ND, mul(N[i], pderiv(D, i)))
    BD = sub(mul(D, divN), ND)
    Nf = {}
    for i in range(5):
        Nf = add(Nf, mul(N[i], prob.fi[i]))
    return add(sub(mul(T, mul(D, D)), mul(D, Nf)), scale(mul(prob.f, BD), Fraction(6, 5)))

def run(iid, restarts=12, maxiter=300):
    print(f'=== {iid} ===', flush=True)
    prob = W.Problem.from_instance(iid)
    d6, D6, N6, Q6 = core.load_witness(iid)
    D = {e: Fraction(v) for e, v in D6.items()}
    N5 = reduce_witness(prob, D, N6)
    res = identity5(prob, D, N5, prob.A)
    print(f'  reduction: deg D = {core.degree(D)}, 5-var identity residual terms = {len(res)} (0 expected)')
    pD, sD = W.parity_of(D, W.rpar), (W.parity_of(D, W.spar) if prob.use_sigma else None)
    print(f'  sector of D: rho={pD} sigma={sD}; use_sigma={prob.use_sigma}; A parities rho={prob.rA} sigma={prob.sA}')
    d = core.degree(D)
    op = W.Operator5(prob, d, pD, sD, prob.A, verbose=True)
    cstar = W.vec_from_poly(op.Dm, D)
    vp = W.VarPro5(op, np.eye(op.nD))
    rel, x, r = vp.exact_residual(cstar)
    print(f'  exact residual at reduced D: rel={rel:.3e} (expect < 1e-10)')
    # inner solve check and gradient check at a random point
    rng = np.random.default_rng(1)
    c = rng.standard_normal(op.nD); c /= np.linalg.norm(c)
    M, b = op.assemble(c)
    x1, r1, lam = vp.inner(M, b)
    x2, *_ = W.sla.lstsq(M, b, cond=1e-10, lapack_driver='gelsd')
    r2 = M @ x2 - b
    print(f'  inner solve: |r_reg|={np.linalg.norm(r1):.6e} |r_svd|={np.linalg.norm(r2):.6e} lam={lam:.2e} (M {M.shape})')
    u = rng.standard_normal(op.nD)
    f0, g0 = vp.value_and_grad(u)
    h = 1e-6; errs = []
    for k in rng.choice(op.nD, size=min(6, op.nD), replace=False):
        du = np.zeros(op.nD); du[k] = h
        fp, _ = vp.value_and_grad(u + du); fm, _ = vp.value_and_grad(u - du)
        fd = (fp - fm) / (2 * h)
        errs.append(abs(fd - g0[k]) / max(abs(g0[k]), 1e-12))
    print(f'  gradient check: max rel err over 6 coords = {max(errs):.2e}')
    # E1 kernel and family
    Kb, info = W.e1_kernel(prob, d, pD, sD, verbose=False)
    print(f'  E1 kernel: {info}')
    proj = Kb @ (Kb.T @ cstar)
    print(f'  D in float kernel? ||c* - proj|| / ||c*|| = {np.linalg.norm(cstar - proj) / np.linalg.norm(cstar):.2e}')
    Bf = W.fmult_family(prob, d, pD, sD)
    B = W.complement_within(Kb, Bf) if Bf.shape[1] else Kb
    print(f'  family dim K = {B.shape[1]} (kernel {Kb.shape[1]} minus f-multiples {Bf.shape[1]})')
    vp = W.VarPro5(op, B)
    t0 = time.time()
    best = None
    hits = 0
    for s in range(restarts):
        out = vp.restart(1000 + s, maxiter=maxiter, init='gauss' if s % 2 == 0 else 'sparse')
        if out['rel'] < 1e-8:
            hits += 1
        if best is None or out['rel'] < best['rel']:
            best = out
        print(f"    restart {s}: rel={out['rel']:.3e} iters={out['iters']} nfev={out['nfev']} {out['secs']}s {out['msg']}")
    print(f'  {hits}/{restarts} restarts reached rel < 1e-8 in {time.time()-t0:.1f}s; best rel = {best["rel"]:.3e}')
    # try to confirm the best candidate exactly
    cands = W.rationalize(best['c'])
    print(f'  rational candidates: {len(cands)}')
    for fr in cands[:4]:
        Dc = W.poly_from_vec(op.Dm, fr)
        st, inf = W.exact_confirm(prob, Dc, pD, sD, prob.A)
        print(f'    confirm {core.to_str(Dc)[:70]} -> {st} {dict((k, v) for k, v in inf.items() if k in ("rank", "nullity", "verified6", "primes_used", "secs"))}')
        if st == 'solved':
            break
    # also the known solution through the exact confirmer
    st, inf = W.exact_confirm(prob, D, pD, sD, prob.A)
    print(f'  exact_confirm on the reduced accepted D -> {st}, verified6={inf.get("verified6")}')
    sa = W.served_analysis(prob, D, pD, sD)
    print(f'  served_analysis on accepted D: {sa}')

if __name__ == '__main__':
    ids = sys.argv[1:] or ['family-hesse-web']
    for iid in ids:
        run(iid)
