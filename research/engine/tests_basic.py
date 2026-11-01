import sys, time
sys.path.insert(0, '.')
from hodge import core, fixed_d
from flint import nmod_mat

# API sanity
M = nmod_mat(2, 3, [1,2,3,4,5,6], 101)
print('rref ->', M.rref())
print('nullspace ->', M.nullspace())
print('rank ->', M.rank())

for iid in ['collision-h3', 'cone-ladder-linear', 'hesse-triangle', 'block-q2', 'family-sylvester-norm']:
    t0 = time.time()
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    print(f'{iid}: q={q} d={d} verify={ok} {info} ({time.time()-t0:.1f}s)')
