"""target_image for D in J (pass E1): simplicity predicts dim 18 (only J_6^-) or 226 (all)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
F, A, q = load()
Fi = [pderiv(F, i) for i in range(NV)]
X0 = core.parse('X0'); Y0 = core.parse('Y0')
def EX(P_): return add(add(mul(core.parse('X0'), pderiv(P_, 0)), mul(core.parse('X1'), pderiv(P_, 1))), mul(core.parse('X2'), pderiv(P_, 2)))
cands = [('F_X0', Fi[0]), ('E_X(F)', EX(F)), ('X0*F_X0 + Y0*F_Y0', add(mul(X0, Fi[0]), mul(Y0, Fi[3]))),
         ('X0^2*F_X0', mul(core.parse('X0**2'), Fi[0]))]
for nm, D in cands:
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D)
    inJ6 = 0
    print(f'{nm}: d={core.degree(D)} target_image dim = {dim}  rows={S.nrows} cols={S.ncols} ({time.time()-t0:.1f}s)', flush=True)
    with open('/tmp/hodge/results/dense-rational/timage.jsonl', 'a') as f:
        f.write(json.dumps({'D': nm, 'd': core.degree(D), 'target_image_dim': dim, 'seconds': round(time.time()-t0,1), 'engine': 2}) + '\n')
