"""Is the E1 kernel generated in degree 7?  Compare span(E1_7 * S_1) + J_8 with E1_8 (both parities), mod p, from saved bases."""
import json, time
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
def rank_of(polys, mons):
    idx = {m: k for k, m in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for j, c in enumerate(polys):
        for e, v in c.items(): flat[j * len(mons) + idx[e]] = v
    return nmod_mat(len(polys), len(mons), flat, p).rank()
E7 = [core.parse(s, mod=p) for s in json.load(open(f'{WORK}/e1kernel_d7_par0.json'))['basis_modp']]
print('E1_7 even basis size', len(E7))
for par in (0, 1):
    E8 = [core.parse(s, mod=p) for s in json.load(open(f'{WORK}/e1kernel_d8_par{par}.json'))['basis_modp']]
    mons = core.monomials(8, par)
    J8 = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(3, (par - pi) % 2):
            J8.append(core.shift(Fi[i], m, 1, p))
    lin = [3, 4, 5] if par == 1 else [0, 1, 2]
    prods = [core.shift(e, tuple(1 if k == i else 0 for k in range(6)), 1, p) for e in E7 for i in lin]
    rJ = rank_of(J8, mons); rP = rank_of(J8 + prods, mons); rE = rank_of(E8, mons); rAll = rank_of(J8 + prods + E8, mons)
    print(f'parity {par}: dim J_8={rJ}, dim(J_8 + E1_7*S_1)={rP}, dim E1_8={rE}, dim(sum)={rAll}  -> generated in degree 7: {rP == rE == rAll}')
