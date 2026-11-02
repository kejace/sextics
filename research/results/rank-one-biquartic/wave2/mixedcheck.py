"""Exact (mod p) identification of the degree-9 Y-even admissible space: W_9 = J_9 + (Y-base-type part) + K_7 . S_2(Y)?"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
from flint import nmod_mat
F, A, q = load(); p = core.PRIMES[0]
S = Space(9, 0)
def rows_of(polys):
    flat = [0] * (len(polys) * S.n)
    for r_, P_ in enumerate(polys):
        for e, v in core.reduce_mod(P_, p).items():
            flat[r_ * S.n + S.idx[e]] = v
    return flat, len(polys)
Jg = [g[2] for g in jacobian.jacobian_generators(F, 9, p, 0)]
k9 = json.load(open(f'{OUT}/e1kernel_d9_p0.json')); K9 = [core.parse(b, p) for b in k9['basis_modp']]
k7 = json.load(open(f'{OUT}/e1kernel_d7_p0.json')); K7 = [core.parse(b, p) for b in k7['basis_modp']]
Hs = [core.parse(h) for h in ('Y0**2', 'Y0*Y1', 'Y0*Y2', 'Y1**2', 'Y1*Y2', 'Y2**2')]
prod = [mul(Ki, core.reduce_mod(H, p), p) for Ki in K7 for H in Hs]           # 18 products (3 kernel x 6 quadrics)
yb = ybase_polys(F, 9, 0)
def rank(polys):
    flat, n = rows_of(polys); return nmod_mat(n, S.n, flat, p).rank()
rJ = rank(Jg); rW = rank(Jg + K9); rJP = rank(Jg + prod); rJY = rank(Jg + yb); rJYP = rank(Jg + yb + prod)
# Y-base part of W: dim (span(yb) cap W) = rank(yb) + rank(W) - rank(yb + W)
rY = rank(yb); rYW = rank(yb + Jg + K9)
rJYK = rank(Jg + yb + K9)
print(json.dumps({'batch': 'mixed-identification', 'dim_J+Ybase': rJY, 'new_from_K7S2Y_beyond_J+Ybase': rJYP - rJY, 'dim_J+Ybase+K9': rJYK, 'K7S2Y_inside_J+Ybase+W': rank(Jg + yb + K9 + prod) == rJYK, 'd': 9, 'p': p, 'dim_J9': rJ, 'dim_W9': rW, 'dim_J+K7S2Y': rJP, 'new_from_K7S2Y': rJP - rJ,
                  'dim_Ybase_cap_W': rY + rW - rYW, 'dim_J+YbaseW+K7S2Y': rJYP, 'W_spanned_by_J_Ybase_K7S2Y': rJYP == rW and rank(Jg + yb + prod + K9) == rW}))
