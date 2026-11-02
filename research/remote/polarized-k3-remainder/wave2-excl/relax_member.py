"""Exact membership test  B_D(T a) in L_{2d}  for explicit a (D = sum P_i a_i), per character block, mod p and over Q.
Used to certify explicit points of the L-relaxation over Q (e.g. the family D = f(U) h) and to re-verify the thread's
DEGREE_EIGHT.md claims (even-V-degree quartic multiples of f_U (E_U T - 3T) lie in L_16; some odd ones do not).
usage: relax_member.py <label> <d> "<a-spec>" [--target=full|k3] [--prime=P] [--out=FILE]
   a-spec as in fixedD_koszul.py (python expression -> list of 6 dict polys; helpers U, V, z, mon, f_U, f_V, mul, add, scale, Fr).
"""
import sys, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fractions import Fraction as Fr
from flint import nmod_mat, fmpz_mat
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, shift, power
from relax import G, P, target, vdeg, z, mon, poly7f, L_generators
NV = 6
U = [z(i) for i in range(3)]; V = [z(3+i) for i in range(3)]
f_U = scale(poly7f(0), Fr(1, 7)); f_V = scale(poly7f(3), Fr(1, 7))

def B_of(a, T):
    D = {}
    for i in range(NV):
        if a[i]: D = add(D, mul(P[i], a[i]))
    Ta = [mul(T, a[i]) if a[i] else {} for i in range(NV)]
    divTa = {}
    for i in range(NV):
        if Ta[i]: divTa = add(divTa, pderiv(Ta[i], i))
    B = mul(D, divTa)
    for i in range(NV):
        if Ta[i]: B = sub(B, mul(Ta[i], pderiv(D, i)))
    return D, B

def member(label, d, a, kind='full', p=65521, exact=True, verbose=True, gens_cache={}):
    t0 = time.time()
    T = target(kind)
    D, B = B_of(a, T)
    deg2 = 2 * d
    if deg2 not in gens_cache: gens_cache[deg2] = L_generators(deg2)
    gens = gens_cache[deg2]
    chs = sorted({vdeg(e) % 6 for e in B})
    res = {'label': label, 'd': d, 'target': kind, 'prime': p, 'B_terms': len(B), 'B_chars': chs, 'blocks': []}
    inL_p = True; inL_Q = True
    for r in chs:
        mons = [m for m in monomials(deg2) if vdeg(m) % 6 == r]
        idx = {m: i for i, m in enumerate(mons)}
        gr = [g for (c, g) in gens if c == r]
        Br = {e: v for e, v in B.items() if vdeg(e) % 6 == r}
        den = 1
        for v in Br.values():
            fv = Fr(v); den = den * fv.denominator // math.gcd(den, fv.denominator)
        Bi = {e: int(Fr(v) * den) for e, v in Br.items()}
        mr, nr = len(mons), len(gr)
        flat = [0] * (mr * (nr + 1))
        for j, g in enumerate(gr):
            for e, v in g.items(): flat[idx[e] * (nr + 1) + j] = int(v)
        for e, v in Bi.items(): flat[idx[e] * (nr + 1) + nr] = v
        Mp = nmod_mat(mr, nr + 1, [v % p for v in flat], p); ra_p = Mp.rank()
        M0p = nmod_mat(mr, nr, [flat[i * (nr + 1) + j] % p for i in range(mr) for j in range(nr)], p); r_p = M0p.rank()
        blk = {'char': r, 'rows': mr, 'L_cols': nr, 'L_rank_p': r_p, 'aug_rank_p': ra_p, 'in_L_p': r_p == ra_p}
        inL_p &= blk['in_L_p']
        if exact:
            Mz = fmpz_mat(mr, nr + 1, flat); ra_Q = Mz.rank()
            M0z = fmpz_mat(mr, nr, [flat[i * (nr + 1) + j] for i in range(mr) for j in range(nr)]); r_Q = M0z.rank()
            blk.update({'L_rank_Q': r_Q, 'aug_rank_Q': ra_Q, 'in_L_Q': r_Q == ra_Q}); inL_Q &= blk['in_L_Q']
        res['blocks'].append(blk)
        if verbose: print('[member-block]', json.dumps(blk), flush=True)
    res['in_L_p'] = inL_p
    if exact: res['in_L_Q'] = inL_Q
    res['seconds'] = round(time.time() - t0, 1)
    return res

if __name__ == '__main__':
    label = sys.argv[1]; d = int(sys.argv[2]); aspec = sys.argv[3]
    kind = 'full'; p = 65521; out = None; exact = '--modular-only' not in sys.argv
    for s in sys.argv[4:]:
        if s.startswith('--target='): kind = s.split('=')[1]
        if s.startswith('--prime='): p = int(s.split('=')[1])
        if s.startswith('--out='): out = s.split('=')[1]
    a = eval(aspec)
    res = member(label, d, a, kind, p, exact)
    print('[member]', json.dumps({k: v for k, v in res.items() if k != 'blocks'}), flush=True)
    if out:
        with open(out, 'a') as fh: fh.write(json.dumps(res) + '\n')
