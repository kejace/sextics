"""Calibration of the served-space machinery on solved instances with target X0^3Y0^3."""
import sys, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
for iid in ['family-hesse-web', 'family-block-cube', 'f3-block-x0-quartic-strict']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D)
    # reduce mod J_6^- : rank of J rows + basis
    from hodge.core import monomials, shift, pderiv, NV
    P = core.PRIMES[0]
    Fp = core.reduce_mod(F, P); Fi = [pderiv(Fp, i, P) for i in range(NV)]
    mons = monomials(6, 1); idx = {m: i for i, m in enumerate(mons)}
    Jrows = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(1, (1 - pi) % 2):
            Jrows.append(shift(Fi[i], m, 1, P))
    from flint import nmod_mat
    def rank_of(rows):
        R, C = len(rows), len(mons); flat = [0]*(R*C)
        for r, poly in enumerate(rows):
            for e, v in poly.items(): flat[r*C + idx[e]] = v
        return nmod_mat(R, C, flat, P).rank()
    rJ = rank_of(Jrows); rT = rank_of(Jrows + basis)
    # Y2-even/odd sector split of served (projection valid if D even in Y2)
    ev = [ {e: v for e, v in b.items() if e[5] % 2 == 0} for b in basis]
    od = [ {e: v for e, v in b.items() if e[5] % 2 == 1} for b in basis]
    rE = rank_of(Jrows + ev) - rJ; rO = rank_of(Jrows + od) - rJ
    print(f'{iid}: d={d} served raw={dim} modJ={rT-rJ} (Y2-even sector {rE}, Y2-odd {rO}) ({time.time()-t0:.0f}s)', flush=True)
