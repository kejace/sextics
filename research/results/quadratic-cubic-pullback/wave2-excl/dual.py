"""Dual E1 count: dim ker([A]: R_d^{parD} -> R_{d+6}^{1-parD}) = dim R_d^{parD} - rank_Q([A]) and rank_Q >= rank_p, so
   ker_Q <= dim R_d - rank_p =: ker_p ; by Gorenstein duality of R (socle degree 24) rank([A]: R_d -> R_{d+6}) =
   rank([A]: R_{18-d} -> R_{24-d}) (parities swapped), which is cheaper for d >= 10.  Lower bound ker >= dim R_d - dim R_{d+6}.
Usage: python3 dual.py d parD p"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *
d = int(sys.argv[1]); parD = int(sys.argv[2]); p = int(sys.argv[3])
F, A, q = load(); S, R = hilbert()
t0 = time.time()
e = 18 - d; parE = (parD + 1) % 2          # source R_e^{parE}, target R_{e+6}^{parD}
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
parT = (parE + 1) % 2
Emons = monomials(e, parE)
cols = [shift(Ap, m, 1, p) for m in Emons]
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(e + 1, (parT - pi) % 2):
        cols.append(shift(Fi[i], m, 1, p))
rowsl = monomials(e + 6, parT); M = mat_from_polys(cols, rowsl, p).transpose()
rk = M.rank()
hilbJ = S[e + 6][parT] - R[e + 6][parT]
rankA = rk - hilbJ
ker_upper = R[e][parE] - rankA            # = ker of [A] on R_e^{parE} (mod p)  -- but we want ker on R_d^{parD}
# rank([A]: R_d^{parD} -> R_{d+6}) = rank([A]: R_e^{parE} -> R_{e+6}^{parD}) by Gorenstein duality
ker_d = R[d][parD] - rankA
gor = max(0, R[d][parD] - R[d + 6][(parD + 1) % 2])
rec = {'exp': 'e1kernel_dual', 'd': d, 'parD': parD, 'p': p, 'via_e': e, 'rows': M.nrows(), 'cols': M.ncols(),
       'rank_stacked': rk, 'hilbert_J_target': hilbJ, 'rank_A_mod_p': rankA, 'dimR_d': R[d][parD],
       'kernel_mod_J_upper_p': ker_d, 'gorenstein_lower_bound': gor,
       'rigorous': ('dim(K_Q/J_d) = %d exactly' % ker_d) if ker_d == gor else ('dim(K_Q/J_d) <= %d' % ker_d),
       'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True); log(rec); receipt(f'dual_d{d}_par{parD}_p{p}.json', rec)
