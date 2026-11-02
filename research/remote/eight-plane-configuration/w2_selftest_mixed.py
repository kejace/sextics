import sys, time, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from w2lib import _col
random.seed(2)
fam = JacobianFamily(6, None, want_C=False)
print('mixed build: K', fam.K, 'dimZ', fam.dimZ, 'rows', fam.rows, 'nc', fam.nc, 'nQ', fam.nQ)
w = [Fraction(random.randint(-3, 3)) for _ in range(fam.K)]
cvec = [Fraction(random.randint(-3, 3)) for _ in range(fam.nc)]
Qvec = [Fraction(random.randint(-3, 3)) for _ in range(fam.nQ)]
Dex = fam.w_to_D(w)
N = [{} for _ in range(NV)]
for i in range(NV):
    for k, m in enumerate(fam.Ws[i].mons):
        c = w[fam.wbase[i] + k]
        if c: N[i] = add(N[i], scale(shift(A, m), c))
for (i, j), (cb, Cs) in fam.Cs.items():
    for k, m in enumerate(Cs.mons):
        c = cvec[cb + k]
        if c:
            N[i] = add(N[i], scale(shift(Fi[j], m), c)); N[j] = add(N[j], scale(shift(Fi[i], m), -c))
Qex = {m: Qvec[k] for k, m in enumerate(fam.Qs.mons) if Qvec[k]}
NF = {}
for i in range(NV): NF = add(NF, mul(N[i], Fi[i]))
assert not sub(NF, mul(A, Dex))
rex = sub(beta_exact(Dex, N), mul(F, Qex)); rvec = fam.Rs.vec(rex)
a = np.array([float(x) for x in w]); c = np.array([float(x) for x in cvec]); Qv = np.array([float(x) for x in Qvec])
D = _col(fam.Dfull @ a); Om = fam.pat.omega(D)
rnum = Om @ (_col(fam.Nfull @ a) + fam.Kmat @ c) - fam.TF @ Qv
print('mixed residual assembly max abs diff', np.max(np.abs(rnum - rvec)), 'scale', np.max(np.abs(rvec)))
t0 = time.time(); info, a1, o = fam.restart(seed=1, max_nfev=80); print('LM', {k: info[k] for k in ('rel', 'nfev', 'status', 'seconds')})
