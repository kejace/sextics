import lib, time
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
# Hilbert function of S/J(F) in degrees 0..24 (mod p) + A in J?
def jac_dim(deg):
    cols = []
    for i in range(6):
        for m in core.monomials(deg - 5):
            cols.append(core.shift(Fi[i], m, 1, p))
    mons = core.monomials(deg)
    idx = {m: k for k, m in enumerate(mons)}
    if not cols: return len(mons), 0
    flat = [0] * (len(mons) * len(cols))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * len(cols) + j] = v
    Mx = nmod_mat(len(mons), len(cols), flat, p)
    r = Mx.rank()
    return len(mons) - r, r
t0 = time.time()
hil = []
for deg in range(0, 13):
    dimq, r = jac_dim(deg)
    hil.append(dimq)
print('Hilbert function of Jacobian ring (deg 0..12):', hil, f'({time.time()-t0:.1f}s)')
ok, info = fixed_d.e1_solvable(F, A, q, core.parse('1'))
print('A in J(F)?', ok, info)
res = run_D(core.parse('1'), '1', 'sanity', note='D=1: is the class zero?')
show(res)
