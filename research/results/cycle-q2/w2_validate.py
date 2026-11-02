"""Validation of w2lib.KoszulVarPro on accepted q=2 witnesses (block-q2, fermat-q2) and on cycle-q2 d=6 kernels."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
from fractions import Fraction
import w2lib as W
from hodge import core

def known(iid):
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    return F, A, q, D

def run_known(iid, n_random=3, restarts=12, seed=1):
    print(f'=== {iid} ===', flush=True)
    F, A, q, D = known(iid)
    S = W.Sectors(F, A)
    kappa = S.key(D); d = core.degree(D)
    print(f'group order {S.order}, gens {S.ng}, key(D)={kappa}, d={d}, parity(D)={core.parity_of(D)}', flush=True)
    t0 = time.time()
    Dm, Kb, info = W.e1_kernel_numeric(S, d, kappa, q)
    print('E1 kernel numeric:', info, f'({time.time()-t0:.1f}s)', flush=True)
    kmod = W.e1_kernel_modp(S, d, kappa, q)
    print('E1 kernel mod p: dimS', kmod['dimS'], 'dimK', kmod['dimK'], 'dimJ', kmod['dimJ'], flush=True)
    # is D in the kernel?
    dv = W.vecs_from_dicts(Dm, [D])[:, 0]
    proj = Kb @ (Kb.T @ dv)
    print('D in E1 kernel (numeric rel err):', np.linalg.norm(proj - dv) / np.linalg.norm(dv), flush=True)
    # family: D plus random kernel elements
    rng = np.random.default_rng(seed)
    Bcols = [dv / np.linalg.norm(dv)]
    for _ in range(min(n_random, Kb.shape[1] - 1)):
        v = Kb @ rng.standard_normal(Kb.shape[1]); Bcols.append(v / np.linalg.norm(v))
    Bmat = np.array(Bcols).T
    t0 = time.time()
    kvp = W.KoszulVarPro(S, d, kappa, Dm, Bmat, q=q, verbose=True)
    a0 = np.zeros(kvp.K, dtype=complex); a0[0] = 1
    f, x, r, b = kvp.value(a0)
    print(f'f at known D: {f:.3e}, |b|={np.linalg.norm(b):.3e}, rel={f/max(np.linalg.norm(b)**2,1e-300):.3e} ({time.time()-t0:.1f}s)', flush=True)
    # gradient check at a random point (complex)
    z0 = rng.standard_normal(2 * kvp.K); z0 /= np.linalg.norm(z0)
    f0, g0 = kvp.value_and_grad_real(z0)
    eps = 1e-6; num = np.zeros_like(z0)
    for i in range(len(z0)):
        zp = z0.copy(); zp[i] += eps; zm = z0.copy(); zm[i] -= eps
        num[i] = (kvp.value_and_grad_real(zp)[0] - kvp.value_and_grad_real(zm)[0]) / (2 * eps)
    print(f'gradient check: |g-num|/|g| = {np.linalg.norm(g0-num)/max(np.linalg.norm(g0),1e-300):.2e}', flush=True)
    t0 = time.time()
    cands = kvp.search(n_restarts=restarts, seed=seed, maxiter=300, verbose=False)
    best = cands[0]
    print(f'search: best rel={best["rel"]:.3e} f={best["f"]:.3e} restart {best["restart"]} iters {best["iters"]} ({time.time()-t0:.1f}s); '
          f'all rels: {[f"{c["rel"]:.1e}" for c in sorted(cands, key=lambda c: c["restart"])]}', flush=True)
    dvv = kvp.D_coeffs(best['a']); v, imrel = W.normalize_phase(dvv)
    print('imag rel of best D:', imrel, flush=True)
    for fr in W.rational_candidates(v.real):
        Dr = W.exact_D_from_fracs(Dm, fr)
        print('candidate D:', core.to_str(Dr)[:120], flush=True)
        res = W.exact_confirm(F, A, q, Dr)
        print('exact confirm:', {k: v for k, v in res.items() if k not in ('N', 'Q')}, flush=True)
        break

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'block'
    if which in ('block', 'all'):
        run_known('block-q2')
    if which in ('fermat', 'all'):
        run_known('fermat-q2')
    if which in ('cycle', 'all'):
        F, A, q, rec = core.load_instance('cycle-q2')
        S = W.Sectors(F, A)
        print('cycle-q2 group order', S.order, 'gens', S.ng, 'kA', S.kA, flush=True)
        for kappa in (21, 96, 597):
            t0 = time.time()
            Dm, Kb, info = W.e1_kernel_numeric(S, 6, kappa, q)
            km = W.e1_kernel_modp(S, 6, kappa, q)
            print(f'key {kappa}: numeric {info}; modp dimS {km["dimS"]} dimK {km["dimK"]} dimJ {km["dimJ"]} ({time.time()-t0:.1f}s)', flush=True)
            if Kb.shape[1] == 0: continue
            t0 = time.time()
            kvp = W.KoszulVarPro(S, 6, kappa, Dm, Kb, q=q, verbose=True)
            cands = kvp.search(n_restarts=6, seed=3, maxiter=200)
            print(f'  search best rel {cands[0]["rel"]:.3e} ({time.time()-t0:.1f}s) rels {[f"{c["rel"]:.1e}" for c in cands]}', flush=True)
