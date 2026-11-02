import sys, time, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from w2lib import _col, _dense
random.seed(1)
for pD in (0, 1):
    t0 = time.time()
    fam = JacobianFamily(6, pD, want_C=True)
    print('build', time.time()-t0, 'dimZ', fam.dimZ, 'C dim', fam.C.shape[1])
    # random rational w, c, Q
    w = [Fraction(random.randint(-3, 3)) for _ in range(fam.K)]
    cvec = [Fraction(random.randint(-3, 3)) for _ in range(fam.nc)]
    Qvec = [Fraction(random.randint(-3, 3)) for _ in range(fam.nQ)]
    Dex = fam.w_to_D(w)
    # N exact = A w + K_c
    N = [{} for _ in range(NV)]
    for i in range(NV):
        for k, m in enumerate(fam.Ws[i].mons):
            c = w[fam.wbase[i] + k]
            if c: N[i] = add(N[i], scale(shift(A, m), c))
    for (i, j), (cb, Cs) in fam.Cs.items():
        for k, m in enumerate(Cs.mons):
            c = cvec[cb + k]
            if c:
                N[i] = add(N[i], scale(shift(Fi[j], m), c))
                N[j] = add(N[j], scale(shift(Fi[i], m), -c))
    Qex = {m: Qvec[k] for k, m in enumerate(fam.Qs.mons) if Qvec[k]}
    # E1 check: N(F) == A D exactly
    NF = {}
    for i in range(NV): NF = add(NF, mul(N[i], Fi[i]))
    assert not sub(NF, mul(A, Dex)), 'E1 normalization broken'
    rex = sub(beta_exact(Dex, N), mul(F, Qex))
    rvec = fam.Rs.vec(rex)
    a = np.array([float(x) for x in w]); c = np.array([float(x) for x in cvec]); Qv = np.array([float(x) for x in Qvec])
    D = _col(fam.Dfull @ a); Om = fam.pat.omega(D)
    rnum = Om @ (_col(fam.Nfull @ a) + fam.Kmat @ c) - fam.TF @ Qv
    print('pD', pD, 'residual assembly max abs diff', np.max(np.abs(rnum - rvec)), 'scale', np.max(np.abs(rvec)))
    # Jacobian finite-difference check
    rng = np.random.default_rng(0)
    u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
    J = fam.jac_u(u)
    r0 = fam.residual_u(u)
    h = 1e-6
    errs = []
    for k in range(3):
        du = np.zeros(fam.dimP); du[k] = h
        fd = (fam.residual_u(u + du) - fam.residual_u(u - du)) / (2*h)
        errs.append(np.linalg.norm(fd - J[:, k]) / max(np.linalg.norm(fd), 1e-300))
    print('pD', pD, 'jac fd rel errs', errs)
    # R-part check: C^T applied to J-elements should vanish
    print('pD', pD, '|C^T (Om Kmat c)|', np.linalg.norm(fam.C.T @ (Om @ (fam.Kmat @ c))), '|C^T TF Q|', np.linalg.norm(fam.C.T @ (fam.TF @ Qv)))
    t0 = time.time()
    info, a1, o = fam.restart(seed=1, max_nfev=60)
    print('pD', pD, 'LM restart', info, 'time', round(time.time()-t0, 2))
