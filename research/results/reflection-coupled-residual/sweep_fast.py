"""Fast E1-kernel dimension sweep: dim Ann(A)_d in the Jacobian ring (Y-even classes D of degree d with A*D in J).
Rank-based (flint nmod_mat.rank), no explicit normal forms.  Usage: python3 sweep_fast.py <target-expr> <dmin> <dmax> <tag>"""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv
from flint import nmod_mat
P = core.PRIMES[0]
F, _, q, rec = core.load_instance('reflection-coupled-residual')
A = core.parse(sys.argv[1]); dmin, dmax = int(sys.argv[2]), int(sys.argv[3]); tag = sys.argv[4]
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]
parA = core.parity_of(A)

def rank_of(rows, mons):
    idx = {m: i for i, m in enumerate(mons)}
    R, C = len(rows), len(mons)
    flat = [0] * (R * C)
    for r, poly in enumerate(rows):
        base = r * C
        for e, v in poly.items():
            flat[base + idx[e]] = v
    return nmod_mat(R, C, flat, P).rank()

out = {}
for d in range(dmin, dmax + 1):
    t0 = time.time()
    parD = 0  # Y-even D
    parT = (parD + parA) % 2
    monsT = monomials(d + 6, parT)
    Jrows = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (parT - pi) % 2):
            Jrows.append(shift(Fi[i], m, 1, P))
    rJ = rank_of(Jrows, monsT)
    Jd_rows = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        if d - 5 >= 0:
            for m in monomials(d - 5, (parD - pi) % 2):
                Jd_rows.append(shift(Fi[i], m, 1, P))
    monsD = monomials(d, parD)
    rJd = rank_of(Jd_rows, monsD) if Jd_rows else 0
    dimRd = len(monsD) - rJd
    dimRT = len(monsT) - rJ
    Arows = [shift(Ap, m, 1, P) for m in monsD]
    rTot = rank_of(Jrows + Arows, monsT)
    rank_mult = rTot - rJ
    ker = dimRd - rank_mult
    generic = max(0, dimRd - dimRT)
    print(f'[{tag}] d={d}: dim R_d^even={dimRd} dim R_(d+6)^{"odd" if parT else "even"}={dimRT} rank(A.)={rank_mult} '
          f'KER dim Ann(A)_d={ker} (generic prediction {generic}) |S_T|={len(monsT)} Jrows={len(Jrows)} ({time.time()-t0:.0f}s)', flush=True)
    out[d] = dict(dimRd=dimRd, dimRT=dimRT, rank=rank_mult, ker=ker, generic=generic, seconds=round(time.time()-t0))
    json.dump(out, open(f'sweep_fast_{tag}.json', 'w'), indent=1)
