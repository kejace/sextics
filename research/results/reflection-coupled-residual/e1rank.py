"""E1-kernel dimension at degree d for arbitrary (F, A): dim{[D] in R_d^{Y-even} : A D in J(F)} (one prime).
Usage: python3 e1rank.py <d> <A-expr> [F-expr | 'instance'] [tag]"""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv
from flint import nmod_mat
P = core.PRIMES[0]
d = int(sys.argv[1]); A = core.parse(sys.argv[2])
Fsrc = sys.argv[3] if len(sys.argv) > 3 else 'instance'
tag = sys.argv[4] if len(sys.argv) > 4 else 'x'
if Fsrc == 'instance':
    F, _, q, rec = core.load_instance('reflection-coupled-residual')
else:
    F = core.parse(Fsrc)
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]
parA = core.parity_of(A); parD = 0; parT = (parD + parA) % 2
def rank_of(rows, mons):
    idx = {m: i for i, m in enumerate(mons)}
    R, C = len(rows), len(mons); flat = [0] * (R * C)
    for r, poly in enumerate(rows):
        base = r * C
        for e, v in poly.items(): flat[base + idx[e]] = v
    return nmod_mat(R, C, flat, P).rank()
t0 = time.time()
monsT = monomials(d + 6, parT); monsD = monomials(d, parD)
Jrows = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d + 1, (parT - pi) % 2): Jrows.append(shift(Fi[i], m, 1, P))
Jd = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    if d >= 5:
        for m in monomials(d - 5, (parD - pi) % 2): Jd.append(shift(Fi[i], m, 1, P))
rJ = rank_of(Jrows, monsT); rJd = rank_of(Jd, monsD) if Jd else 0
dimRd = len(monsD) - rJd; dimRT = len(monsT) - rJ
rTot = rank_of(Jrows + [shift(Ap, m, 1, P) for m in monsD], monsT)
ker = dimRd - (rTot - rJ)
# smoothness check of F: Hilbert function of R at degree d+6 must equal the Fermat count
import itertools
fermat = sum(1 for e in itertools.product(range(5), repeat=6) if sum(e) == d + 6 and (sum(e[3:]) % 2) == parT)
print(f'[{tag}] d={d} A={sys.argv[2]} dim R_d^even={dimRd} dim R_(d+6)={dimRT} (Fermat count {fermat}, smooth={dimRT==fermat}) KER={ker} ({time.time()-t0:.0f}s)', flush=True)
