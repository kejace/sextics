import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
iid = 'coupled-pair-rational-d4-strict'
F, A, q, rec = core.load_instance(iid)
d, D, N, Q = core.load_witness(iid)
N2 = normalize_lift(F, A, D, N)
prob = E2Float(F, A, d, parity_of(D), verbose=False)
Dv = prob.SD.vec(D); Nv = prob.Nvec(N2)
fam = E1Family(prob)
co = fam.Dbasis @ Dv
N0 = fam.Nbasis.T @ co
diff = N0 - Nv
# E1 of diff
ST = fam.ST
Jm = sp.hstack([mult_op(prob.Fi[i], prob.SN[i], ST) for i in range(6)]).toarray()
print('|Jm diff| =', np.linalg.norm(Jm @ diff), ' |Jm N0 - A D| =', np.linalg.norm(Jm @ N0 - mult_op(A, prob.SD, ST).toarray() @ Dv), ' |Jm Nv - A D|', np.linalg.norm(Jm @ Nv - mult_op(A, prob.SD, ST).toarray() @ Dv))
K = prob.Kmap.toarray()
c, res, rk, sv = np.linalg.lstsq(K, diff, rcond=None)
print('Koszul fit of N0 - N2: residual', np.linalg.norm(K @ c - diff), ' |diff|', np.linalg.norm(diff), ' rank K', rk)
# all 15 constant koszul fields (ignoring parity)
cols = []
for i in range(6):
    for j in range(i+1, 6):
        Nk = [{} for _ in range(6)]
        Nk[i] = prob.Fi[j]; Nk[j] = scale(prob.Fi[i], -1)
        try:
            cols.append(prob.Nvec(Nk))
        except KeyError as ex:
            print('pair', i, j, 'parity-illegal (KeyError)')
K15 = np.array(cols).T
c, res, rk, sv = np.linalg.lstsq(K15, diff, rcond=None)
print('fit with all legal-space constant Koszul fields:', np.linalg.norm(K15 @ c - diff), 'ncols', K15.shape[1])
print('prob.pc', prob.pc)
print('per-pair Sc dims', {ij: prob.Sc[ij].n for ij in prob.pairs})
