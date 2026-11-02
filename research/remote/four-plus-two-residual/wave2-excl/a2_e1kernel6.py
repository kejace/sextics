"""(a) Re-verify the 6-variable E1 kernel {D : A D in J(F)} / J(F)_d at a NEW prime, both parities, with J-rank certification."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
p = int(sys.argv[1]) if len(sys.argv) > 1 else PRIMES[3]
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 7
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Fi = [pderiv(Fp, i, p) for i in range(6)]
for d in range(1, dmax + 1):
    for parD in (0, 1):
        t0 = time.time()
        monsD = monomials(d, parD)
        tpar = (1 + parD) % 2
        Jcols = []
        for i in range(6):
            for m2 in monomials(d + 1, (tpar - (0 if i < 3 else 1)) % 2):
                Jcols.append(shift(Fi[i], m2, 1, p))
        Acols = [shift(Ap, m, 1, p) for m in monsD]
        MJ, _ = build(Jcols, p); rJ = MJ.rank()
        MAJ, _ = build(Acols + Jcols, p); rAJ = MAJ.rank()
        dimKer = len(monsD) - (rAJ - rJ)
        Jd = []
        if d >= 5:
            for i in range(6):
                for m3 in monomials(d - 5, (parD - (0 if i < 3 else 1)) % 2):
                    Jd.append(shift(Fi[i], m3, 1, p))
        dimJd = build(Jd, p)[0].rank() if Jd else 0
        out = dict(exp='a2-e1kernel6-verify', d=d, parD=parD, p=p, n_mons=len(monsD), kernel_dim_p=dimKer,
                   rank_J_d6_p=rJ, dim_J_d6_exact=dim_JF(d + 6, tpar), rank_J_d_p=dimJd, dim_J_d_exact=dim_JF(d, parD),
                   kernel_mod_J_p=dimKer - dimJd, rows=MAJ.nrows(), cols=MAJ.ncols(), seconds=round(time.time() - t0, 1))
        out['certified_upper_bound'] = (rJ == out['dim_J_d6_exact']) and (dimJd == out['dim_J_d_exact'])
        log(out)
