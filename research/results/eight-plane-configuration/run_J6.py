import sys, time
sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration')
from lib import *
import lib
_log = lib.log
def log2(e): e['engine'] = 2; _log(e)
lib.log = log2
Fi = [core.pderiv(F, i) for i in range(6)]
def vf(coeffs):  # coeffs: list of (poly, i) meaning poly * dF/dZ_i
    D = {}
    for c, i in coeffs:
        D = core.add(D, core.mul(c, Fi[i]))
    return D
cands = [
 ('(X0 dX0 + Y0 dY0)F', 'block-radial', vf([(X0,0),(Y0,3)])),
 ('(X1 dX1 + Y1 dY1)F', 'block-radial', vf([(X1,1),(Y1,4)])),
 ('(X2 dX2 + Y2 dY2)F', 'block-radial', vf([(X2,2),(Y2,5)])),
 ('eps_X(F)', 'block-euler', vf([(X0,0),(X1,1),(X2,2)])),
 ('X0*dF/dX0', 'J6-monomial', vf([(X0,0)])),
 ('Y0*dF/dY0', 'J6-monomial', vf([(Y0,3)])),
 ('X1*dF/dX0 - X0*dF/dX1', 'J6-rotation', core.sub(core.mul(X1,Fi[0]), core.mul(X0,Fi[1]))),
 ('Y0*dF/dX0 - X0*dF/dY0', 'J6-rotation-odd', core.sub(core.mul(Y0,Fi[0]), core.mul(X0,Fi[3]))),
 ('X0*dF/dY0', 'J6-monomial-odd', vf([(X0,3)])),
 ('Y0*dF/dX0', 'J6-monomial-odd', vf([(Y0,0)])),
]
for name, fam, D in cands:
    e = test_D(D, name, fam, note='D in J(F)_6 (E1 trivial); tangency locus of a linear vector field', fast=True, want_exact=True)
    print(name, {k: e[k] for k in ('d','e1','consistent','rank','nullity','rows','cols','seconds') if k in e}, flush=True)
