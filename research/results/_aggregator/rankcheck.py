"""Confirm the engine rank-reporting convention for inconsistent fixed-D systems:
compare fixed_d's reported rank/nullity with an independently computed rank of the coefficient matrix (no RHS)."""
import sys, json, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
import inspect
p = 2147483563
F, A, q, rec = core.load_instance('cyclic-triple-residual', p)
D = core.pderiv(F, 0, p)   # dF/dX0, degree 5, in J -> E1 holds, known inconsistent
t0=time.time()
print('ENGINE_VERSION', fixed_d.ENGINE_VERSION)
names=[n for n,_ in inspect.getmembers(fixed_d, inspect.isfunction)]
print('fixed_d functions:', names)
res = fixed_d.solve_fixed_D(F, A, q, D, primes=[p]) if 'primes' in inspect.signature(fixed_d.solve_fixed_D).parameters else fixed_d.solve_fixed_D(F, A, q, D)
print('solve_fixed_D ->', str(res)[:600])
# independent: build the system and compute rank(M) and rank([M|b]) directly
S = fixed_d.FixedDSystem(F, A, q, D, p)
attrs=[a for a in dir(S) if not a.startswith('_')]
print('FixedDSystem attrs:', attrs)
try:
    M, b = S.matrix(), S.rhs()
except Exception as e:
    print('matrix/rhs API differs:', e); M=None
if M is not None:
    from flint import nmod_mat
    r = M.rank()
    # augmented
    R,C = M.nrows(), M.ncols()
    aug = nmod_mat(R, C+1, [0]*(R*(C+1)), p)
    for i in range(R):
        for j in range(C):
            v=int(M[i,j])
            if v: aug[i,j]=v
        bv=int(b[i,0]) if hasattr(b,'nrows') else int(b[i])
        if bv: aug[i,C]=bv
    ra = aug.rank()
    print(json.dumps({'rank_M': r, 'rank_aug': ra, 'ncols': C, 'nullity_M': C-r, 'seconds': round(time.time()-t0,1)}))
