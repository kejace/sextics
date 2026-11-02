"""How much of the saved E1 kernel (degrees 7, 8) is base-type: Y-shape P(X)M^a g^b, X-shape Q(Y)L^a f^b, or both, modulo J."""
import json, time, sys
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
def xmons(k, ys=False):
    if ys: return [m for m in core.monomials(k) if m[0] == m[1] == m[2] == 0]
    return [m for m in core.monomials(k) if m[3] == m[4] == m[5] == 0]
def base(d, side, par):
    Q2, Q6 = (M, g) if side == 'Y' else (L, f)
    out = []
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0: continue
            W = core.power(Q2, a) if a else {(0,)*6: Fraction(1)}
            if b: W = core.mul(W, Q6)
            Wp = core.reduce_mod(W, p)
            for m in xmons(k, ys=(side == 'X')):
                P_ = core.shift(Wp, m, 1, p)
                if core.parity_of(P_) == par: out.append(P_)
    return out
def rank_of(polys, mons):
    idx = {m: k for k, m in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for j, c in enumerate(polys):
        for e, v in c.items(): flat[j * len(mons) + idx[e]] = v
    return nmod_mat(len(polys), len(mons), flat, p).rank() if polys else 0
for d, par in [(7, 0), (8, 0), (8, 1)]:
    E = [core.parse(s, mod=p) for s in json.load(open(f'{WORK}/e1kernel_d{d}_par{par}.json'))['basis_modp']]
    mons = core.monomials(d, par)
    J = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(d - 5, (par - pi) % 2):
            J.append(core.shift(Fi[i], m, 1, p))
    rJ = rank_of(J, mons); rE = rank_of(E + J, mons)
    BY = base(d, 'Y', par); BX = base(d, 'X', par)
    for name, B in [('Y-shape P(X)M^a g^b', BY), ('X-shape Q(Y)L^a f^b', BX), ('both', BY + BX)]:
        rJB = rank_of(J + B, mons); rEJB = rank_of(E + J + B, mons)
        print(f'd={d} par={par} {name:22s}: |B|={len(B)}, dim J={rJ}, dim E1={rE}, dim(J+B)={rJB}, dim(E1+J+B)={rEJB}; '
              f'kernel dirs NOT in J+B: {rEJB - rJB}; base-type kernel dirs: {(rE - rJ) - (rEJB - rJB)}', flush=True)
