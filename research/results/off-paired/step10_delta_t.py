"""Step 10: full fixed-D solve for Delta_t = 4Qx^3 + 27 f6^2 g6 (degree 18, G-invariant) with a chunked, memory-light
row compression (random combinations of the 44k sparse rows, computed in blocks) and a generator-built nmod_mat."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import oplib
from hodge import core, fixed_d
from hodge.core import mul, add, scale, power, sub
import numpy as np, scipy.sparse as sp
from flint import nmod_mat
F, A, q = oplib.load(); p = oplib.P0
f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6'); Qx = sub(sub(F, f6), g6)
Dt = add(scale(power(Qx, 3), 4), scale(mul(mul(f6, f6), g6), 27))
t0 = time.time()
S = oplib.ReducedSystem(F, A, q, Dt, p, verbose=True)
M = fixed_d._sparse_matrix(S)            # nrows x (C+1), int64 residues
C = S.ncols; k = min(S.nrows, C + 1 + 40)
MT = M.T.tocsr(); del M
rng = np.random.default_rng(1)
blk = 1024
out = np.empty((k, C + 1), dtype=np.int64)
for s in range(0, k, blk):
    e = min(k, s + blk)
    R = rng.integers(0, p, size=(e - s, S.nrows), dtype=np.int64)
    lo = (MT @ (R & 0xFFFF).T) % p
    hi = (MT @ (R >> 16).T) % p
    out[s:e] = ((lo + (hi * (1 << 16)) % p) % p).T
    del R, lo, hi
t1 = time.time(); print(f'compressed k={k} C={C} in {t1-t0:.0f}s', file=sys.stderr, flush=True)
del MT
def gen():
    for i in range(k):
        for v in out[i].tolist(): yield v
Mn = nmod_mat(k, C + 1, gen(), p); del out
t2 = time.time(); print(f'nmod_mat built {t2-t1:.0f}s', file=sys.stderr, flush=True)
Mr, rank = Mn.rref(); del Mn
t3 = time.time()
pivots = []
for i in range(rank):
    j = 0
    while j <= C and int(Mr[i, j]) == 0: j += 1
    pivots.append(j)
consistent = not (pivots and pivots[-1] == C)
rank0 = rank - (0 if consistent else 1)
print(f'RESULT consistent={consistent} rank={rank0} nullity={C-rank0} rref {t3-t2:.0f}s total {t3-t0:.0f}s', file=sys.stderr, flush=True)
oplib.log('4*Qx**3 + 27*f6**2*g6 (cubic-cover discriminant)', 18, 'cubic-cover-discriminant', True, consistent, rank0, C - rank0, t3 - t0,
          f'full G-reduced compressed solve rows={S.nrows} cols={S.ncols} k={k}')
if consistent:
    print('!!!!! CONSISTENT — attempting exact lift', file=sys.stderr, flush=True)
    st, data = oplib.lift_solution(F, A, q, Dt, verbose=True)
    print(st, file=sys.stderr, flush=True)
    if st == 'solved':
        json.dump({'D': core.to_str(Dt), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open('witness-candidate.json', 'w'))
