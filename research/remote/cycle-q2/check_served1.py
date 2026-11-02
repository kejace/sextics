import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
from flint import nmod_mat
F, A, q, rec = core.load_instance('cycle-q2')
p = 2147483647
Gw = eq.symmetry_group(F); gens = eq.generators(Gw)
kA = eq.key_of_poly(A, gens)
Fp = core.reduce_mod(F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
# full D=1 system, no symmetry: columns N_i m (deg 7), image 2 N.gradF - F div N
cols = []
for i in range(6):
    for m in core.monomials(7):
        Nm = {m: 1}
        v = core.scale(core.mul(Nm, Fi[i], p), 2, p)
        dN = core.pderiv(Nm, i, p)
        v = core.sub(v, core.mul(Fp, dN, p), p)
        cols.append(v)
nN = len(cols)
allT = eq.monomial_exps(12); kT = eq.keys_of(allT, gens)
Tsel = [tuple(int(x) for x in e) for e in allT[kT == kA]]
for m in Tsel: cols.append({m: 1})
rows = {}
for c in cols:
    for e in c: rows.setdefault(e, len(rows))
M = np.zeros((len(rows), len(cols)), dtype=np.int64)
for j, c in enumerate(cols):
    for e, v in c.items(): M[rows[e], j] = v % p
t0 = time.time()
r_all = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p).rank()
r_N = nmod_mat(M.shape[0], nN, M[:, :nN].ravel().tolist(), p).rank()
print('no-symmetry D=1: N cols', nN, 'targets', len(Tsel), 'rank N', r_N, 'rank all', r_all, 'served =', r_all - r_N, f'{time.time()-t0:.0f}s')
# sector-only J_12 and the reduction map rank
