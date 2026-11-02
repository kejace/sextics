"""Certify that the Jacobian ideal of F mod p has the characteristic-zero Hilbert function in every degree m <= mmax
(sum over G_aff character sectors of rank J_m[c]).  Equality in all degrees <= 25 shows J mod p is a complete intersection
(F mod p smooth), and in a given degree m it makes modular kernel/colon dimensions upper bounds for the rational ones and
modular non-membership in J_m a proof of rational non-membership.  Usage: python3 jachilb.py [mmax]"""
import sys, time, json
from w2lib import *
HF = [1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, 1686, 1506, 1246, 951, 666, 426, 246, 126, 56, 21, 6, 1, 0, 0]
mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
for m in range(5, mmax + 1):
    t0 = time.time()
    chars = sorted({chi(e) for e in monomials(m)})
    rk = 0; per = {}
    for c in chars:
        S = Sector(m, c)
        rk += S.rank; per[str(c)] = (S.n, S.rank)
    exp = len(monomials(m)) - HF[m]
    log('jacobian-hilbert', m=m, rank_p=rk, expected_rank=exp, ok=(rk == exp), dimS=len(monomials(m)), n_sectors=len(chars), seconds=round(time.time() - t0, 1),
        note='rank of J_m mod p summed over G_aff sectors vs dim S_m - HF_R(m) (regular-sequence Hilbert function)')
print('DONE', file=sys.stderr)
