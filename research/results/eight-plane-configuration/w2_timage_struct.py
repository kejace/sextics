"""target_image for structured Jacobian-type denominators (exact rational D)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from hodge import jacobian
Q0 = core.parse('X0**2+Y0**2'); Q1 = core.parse('X1**2+Y1**2'); Q2 = core.parse('X2**2+Y2**2')
cands = {'Q0*Q1*dF/dX0 (d=8, vanishes on all eight planes)': mul(mul(Q0, Q1), Fi[0])}
for name, D in cands.items():
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D, p=PRIME)
    inJ = sum(1 for b in basis if jacobian.in_jacobian(F, b, p=PRIME)[0])
    rec = {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'target_image of structured Jacobian-type D (mod p)', 'D': name, 'd': core.degree(D),
           'target_image_dim': int(dim), 'basis_elements_in_J6': inJ, 'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True)
    jlog(f'{WD}/wave2-search.log.jsonl', rec)
