import time, sys
from w2lib import *
d = int(sys.argv[1])
t0 = time.time()
colon = dimJ = dimS = R6 = 0
for c in sorted({chi(e) for e in monomials(d)}):
    t1 = time.time()
    mons = [e for e in monomials(d) if chi(e) == c]
    S6 = sector(d + 6, chi_add(c, CHI_A))
    NF = S6.nf_matrix([S6.vec(shift(Ap, e, 1, p)) for e in mons])
    nul = len(mons) - NF.rank()
    Sd = sector(d, c)
    colon += nul; dimJ += Sd.rank; dimS += len(mons); R6 += S6.dimR
    print(c, len(mons), Sd.rank, nul, S6.dimR, f'{time.time()-t1:.0f}s', flush=True)
print(f'd={d}: dim S_d={dimS}, dim J_d={dimJ}, dim (J:A)_d={colon}, non-Jacobian={colon-dimJ}, dim R_(d+6)={R6}, dim S_d - dim R_(d+6)={dimS - R6} ({time.time()-t0:.0f}s)', flush=True)
log('E1-colon', d=d, dim_colon=colon, dim_J=dimJ, n_extra=colon - dimJ, dim_S=dimS, dim_R_d6=R6, seconds=round(time.time() - t0, 1),
    note='dims only; dossier claims (J:A)_d = S_d for d >= 13 "since R_19 = 0" but R_19 = 246; prediction by Gorenstein duality: dim (J:A)_13 = 8568 - 246 = 8322')
