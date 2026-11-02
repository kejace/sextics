"""Referee replay of E1 by Klein character at an unused prime.
Usage: python3 ref_e1.py <dmin> <dmax> <p>
E1-admissible space K_d^chi = {D in S_d^chi : A*D in J}; kernel_mod_J = dim K - dim J_d^chi.
Good reduction: dim (S/J)_n^chi mod p == Fermat count (monomials with all exponents <= 4) in degrees d and d+6."""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
from flint import nmod_mat
dmin, dmax, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
S = R.Setup(p)
OUT = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee/ref_e1.jsonl'
for d in range(dmin, dmax + 1):
    for c in R.CHARS:
        t0 = time.time()
        cols = R.monomials(d, c)
        Rd = R.Reducer(S, d, c)
        tc = R.cadd(c, R.CHI_A)
        Rt = R.Reducer(S, d + 6, tc)
        Arows = [R.shift(S.Ap, m, p) for m in cols]
        rank_A_mod_J = Rt.rank_mod_J(Arows) if Arows else 0
        ker = len(cols) - rank_A_mod_J
        R.receipt(OUT, dict(exp='ref-e1-kernel', d=d, chi=list(c), p=p, dim_Sd=len(cols), rankJ_d=Rd.rank, dim_RJ_d=Rd.dimRJ,
                            fermat_d=Rd.fermat, good_d=Rd.good, target_chi=list(tc), dim_S_target=Rt.dimS, rankJ_target=Rt.rank,
                            dim_RJ_target=Rt.dimRJ, fermat_target=Rt.fermat, good_target=Rt.good, rank_A_mod_J=rank_A_mod_J,
                            e1_kernel_dim=ker, kernel_mod_J=ker - Rd.rank, seconds=round(time.time() - t0, 1)))
print('REF-E1 DONE', flush=True)
