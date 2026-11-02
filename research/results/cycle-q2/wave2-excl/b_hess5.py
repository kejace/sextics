"""E4: degree-5 exclusion via the Hessian relaxation.
If deg D = 5 then (Ann_5 = 0) D = sum c_i F_i, N = (1/2) c A + Koszul, and E2 mod J forces
sum_{i,j} c_i c_j [A F_{ij}] = 0 in R_16.  If the 21 classes [A F_{ij}] (i<=j) are linearly independent
mod p (hence over Q), then c_i c_j = 0 for all i,j, so c = 0: no certificate of degree 5 at all
(any parity, any coefficients).  We compute the rank sector-wise (classes are semi-invariant)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
for p in (P2, P1):
    t0 = time.time()
    JS = JacobianSectors(F2, p)
    Fi = [pderiv(F2, i) for i in range(6)]
    classes = {}
    for i in range(6):
        for j in range(i, 6):
            H = pmul(A, pderiv(Fi[i], j))
            if not H:
                classes[(i, j)] = None; continue
            key = key_of_poly(H, JS.gens)
            classes.setdefault(key, []).append(((i, j), H))
    total_rank = 0; detail = {}
    nclasses = 0
    for key, lst in classes.items():
        if key is None or lst is None:
            continue
        polys = [H for _, H in lst]
        C = JS.normal_form_coords(polys, 16, key)
        rk = C.rank()
        total_rank += rk; nclasses += len(polys)
        detail[','.join(map(str, key))] = dict(pairs=[list(ij) for ij, _ in lst], rank=rk, ring_dim=JS.ring_dim(16, key))
        print(key, [ij for ij, _ in lst], 'rank', rk, '/', len(polys), flush=True)
    zero = [str(k) for k, v in classes.items() if v is None]
    rec = dict(exp='E4-deg5-hessian-classes', engine=ENGINE_TAG, p=p, n_classes=nclasses, zero_classes=zero, total_rank=total_rank,
               independent=(total_rank == 21), detail=detail, seconds=round(time.time() - t0, 1),
               conclusion=('RIGOROUS: no certificate with deg D = 5 (any parity/coefficients)' if total_rank == 21 else 'classes dependent: residual quadric system in c remains'))
    print(json.dumps(rec)[:600])
    xlib.log(rec)
