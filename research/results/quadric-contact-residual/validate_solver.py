"""Validate SectorSystem (trivial character) against known witnesses of solved instances."""
import sys, time
from common import *
triv = lambda e: (0, 0, 0, 0)
for iid in ['collision-h3', 'hesse-triangle', 'ternary-collision-cb-d2']:
    try:
        Fi, Ai, qi, _ = core.load_instance(iid)
        d, Dw, Nw, Qw = core.load_witness(iid)
    except Exception as ex:
        print(iid, 'skip', ex); continue
    if qi != 1: print(iid, 'q!=1 skip'); continue
    t0 = time.time()
    S = SectorSystem(Dw, F_=Fi, A_=Ai, chi_fn=triv, verbose=True)
    res = S.solve(verbose=True)
    N, Q = S.unpack(res['x'])
    Fp = core.reduce_mod(Fi, S.p); Ap = core.reduce_mod(Ai, S.p); Dp = core.reduce_mod(Dw, S.p)
    E = core.defect(Fp, Ap, Dp, 1, N, S.p)
    Rr = mul(power(Fp, 2, S.p), Q, S.p)
    print(iid, 'd', d, 'consistent', res['consistent'], 'defect ok mod p:', core.sub(E, Rr, S.p) == {}, f'({time.time()-t0:.1f}s)')
    # perturbed D must be inconsistent
    D2 = dict(Dw); e = next(iter(D2)); D2[e] = D2[e] + 1
    S2 = SectorSystem(D2, F_=Fi, A_=Ai, chi_fn=triv); r2 = S2.solve()
    print(iid, 'perturbed consistent =', r2['consistent'])
