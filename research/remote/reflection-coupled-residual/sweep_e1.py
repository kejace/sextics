"""Degree sweep: kernel of multiplication by [A] on R_d^{even} (Y-even, degree d) -> R_{d+6}^{odd}."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
from rcr import *
out = {}
for d in range(1, 9):
    t0 = time.time()
    Jd = JacobianDegree(d, 0)          # even part in degree d
    Jd6 = JacobianDegree(d + 6, 1)     # odd part in degree d+6
    # multiplication by A on classes: build matrix from R_d^even basis (nonpivot monomials) to R_{d+6}^odd coords
    cols = []
    for j in Jd.nonpiv:
        m = Jd.mons[j]
        cols.append(Jd6.normal_form(shift(Ap, m, 1, P)))
    n_in, n_out = len(cols), Jd6.quot_dim
    flat = [0] * (n_out * n_in)
    for c, v in enumerate(cols):
        for r in range(n_out):
            flat[r * n_in + c] = v[r]
    M = nmod_mat(n_out, n_in, flat, P)
    X, nullity = M.nullspace()
    rank = n_in - nullity
    # also: kernel restricted to Y2-even (r2-invariant) classes
    ker_polys = []
    for k in range(nullity):
        poly = {}
        for c, j in enumerate(Jd.nonpiv):
            v = int(X[c, k])
            if v: poly[Jd.mons[j]] = v
        ker_polys.append(poly)
    print(f'd={d}: dim R_d^even={n_in} dim R_(d+6)^odd={n_out} rank={rank} KERNEL dim={nullity}  ({time.time()-t0:.1f}s)', flush=True)
    for kp in ker_polys[:12]:
        s = core.to_str({e: (v if v < P//2 else v - P) for e, v in kp.items()})
        print('   ker:', s[:200])
    out[d] = {'dim_in': n_in, 'dim_out': n_out, 'rank': rank, 'nullity': nullity,
              'kernel': [core.to_str({e: (v if v < P//2 else v - P) for e, v in kp.items()}) for kp in ker_polys]}
    json.dump(out, open('sweep_e1.json', 'w'), indent=1)
