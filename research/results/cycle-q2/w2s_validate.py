"""Validation of w2s_lib on accepted q=2 witnesses (block-q2: family mode; fermat-q2: wspace/hess modes)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
import w2s_lib as L
from fractions import Fraction
from hodge import core

def check_instance(iid, Dstr, modes):
    I = L.Inst(iid)
    print(iid, 'group order', I.n, 'elements', len(I.elements), 'gens', I.gens, 'separate', I.gens_separate)
    D = core.parse(Dstr)
    d = core.degree(D)
    keyD = I.key(next(iter(D)))
    assert all(I.key(e) == keyD for e in D), 'D not in one sector'
    print(' D =', Dstr, 'd =', d, 'key', keyD, 'parity', core.parity_of(D))
    t = time.time()
    kd = L.sector_kernel(I, d, keyD, parity=core.parity_of(D), verbose=True)
    print(' kernel dims', {k: v for k, v in kd.items() if k in ('nmons', 'kernel_dim', 'jac_dim', 'nonjac_dim', 'ring_dim', 'primes', 'shape', 'seconds')})
    # is D in the kernel span?  (exact: solve in rref basis: coefficients = free coordinates)
    piv = kd['pivots']; mons = kd['mons']
    # express D via free coordinates
    coeff = [D.get(mons[p], Fraction(0)) for p in piv]
    Dre = {}
    for c, B in zip(coeff, kd['basis']):
        if c: Dre = core.add(Dre, core.scale(B, c))
    print(' D reconstructed from kernel basis equals D:', Dre == {e: v for e, v in D.items() if v != 0})
    for mode in modes:
        if mode == 'family':
            svp = L.SectorVarPro(I, keyD, mode='family', basis=kd['basis'], parity=core.parity_of(D), verbose=True)
            a_true = np.array([float(c) for c in coeff])
        else:
            wb = L.wspace_basis(I, d, keyD, parity=core.parity_of(D))
            print(' wbasis size', len(wb))
            svp = L.SectorVarPro(I, keyD, mode=mode, wbasis=wb, parity=core.parity_of(D), verbose=True)
            # find w with sum m F_i = D : least squares on Dk
            mons2 = {}
            for B in svp.Dk + [D]:
                for e in B: mons2.setdefault(e, len(mons2))
            M = np.zeros((len(mons2), svp.K)); b = np.zeros(len(mons2))
            for k, B in enumerate(svp.Dk):
                for e, v in B.items(): M[mons2[e], k] = float(v)
            for e, v in D.items(): b[mons2[e]] = float(v)
            a_true, *_ = np.linalg.lstsq(M, b, rcond=None)
            print(' w-fit residual', np.linalg.norm(M @ a_true - b))
        f, x, r, bb = svp.value(a_true)
        print(f' mode={mode}: f at true D = {f:.3e} rel = {f/bb:.3e}  (bb={bb:.3e})')
        t = time.time()
        rng = np.random.default_rng(1)
        best = None
        for it in range(6):
            c = svp.one_restart(rng, maxiter=300)
            print(f'   restart {it}: rel={c["rel"]:.3e} iters={c["iters"]} spurious={c["spurious"]}')
            if not c['spurious'] and (best is None or c['rel'] < best['rel']): best = c
        print(f'   best rel {best["rel"]:.3e} in {time.time()-t:.1f}s')
        if best['rel'] < 1e-10:
            cands = L.reconstruct_candidates(best['a'])
            for fr, err in cands[:6]:
                Dx = svp.D_exact(fr)
                st, data = L.exact_confirm(I, Dx)
                print('   reconstruct err', err, '->', st, core.to_str(Dx)[:100])
                if st == 'solved': break

check_instance('block-q2', 'X0**6 + 2*X0**4*Y0**2 - X0**2*Y0**4 + Y0**6', ['family'])
check_instance('fermat-q2', '6*X1**6 + 6*X2**6', ['hess', 'wspace', 'family'])
