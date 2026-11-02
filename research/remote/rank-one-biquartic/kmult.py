"""Which multipliers m kill the degree-7 kernel K_7 (mod J)?  count = 3 - (rank[J_{7+deg m} | m*E1_7] - rank J_{7+deg m})."""
import json
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
E7 = [core.parse(s, mod=p) for s in json.load(open(f'{WORK}/e1kernel_d7_par0.json'))['basis_modp']]
def rank_of(polys, mons):
    idx = {m: k for k, m in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for j, c in enumerate(polys):
        for e, v in c.items(): flat[j * len(mons) + idx[e]] = v
    return nmod_mat(len(polys), len(mons), flat, p).rank()
mults = [('X0', P('X0')), ('X1', P('X1')), ('X0^2', P('X0**2')), ('X0*X1', P('X0*X1')), ('L', L), ('X0^3', P('X0**3')), ('X0^2*X1', P('X0**2*X1')),
         ('X0*X1*X2', P('X0*X1*X2')), ('X0*L', core.mul(P('X0'), L)), ('Y0', P('Y0')), ('Y0^2', P('Y0**2')), ('M', M), ('X0*Y0', P('X0*Y0')),
         ('X0^2*Y0', P('X0**2*Y0')), ('X0*Y0^2', P('X0*Y0**2')), ('X0*M', core.mul(P('X0'), M)), ('Y0^3', P('Y0**3')), ('X0^2*Y0^2', P('X0**2*Y0**2')),
         ('L*M', core.mul(L, M)), ('X0^3*Y0', P('X0**3*Y0')), ('X0^3*Y0^2', P('X0**3*Y0**2')), ('X0^2*Y0^3', P('X0**2*Y0**3')), ('X0*Y0^3', P('X0*Y0**3')),
         ('f', f), ('g', g), ('L^2', core.mul(L, L)), ('M^2', core.mul(M, M)), ('L^2*M', Pc), ('F_3 (dF/dY0)', core.pderiv(F, 3)), ('F_0 (dF/dX0)', core.pderiv(F, 0))]
for name, m in mults:
    mp = core.reduce_mod(m, p); dm = core.degree(m); par = (0 + core.parity_of(m)) % 2
    mons = core.monomials(7 + dm, par)
    J = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for mm in core.monomials(7 + dm - 5, (par - pi) % 2):
            J.append(core.shift(Fi[i], mm, 1, p))
    prods = [core.mul(e, mp, p) for e in E7]
    rJ = rank_of(J, mons); rJE = rank_of(J + prods, mons)
    print(f'{name:14s}: kernel directions killed = {3 - (rJE - rJ)} of 3', flush=True)
