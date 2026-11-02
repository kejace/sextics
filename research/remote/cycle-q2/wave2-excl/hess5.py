"""Hessian relaxation at deg D = 5 (lens excl, lane cycle-q2).

Rigorous chain: E1 forces A D in J; Ann_5([A]) = 0 (ann.sing) forces D = d_w F = sum_i w_i F_i, w in C^6.
With N0 = (A/2) w the E1 equation holds exactly; E2 reads  H(w) + (Koszul terms in J) = 2 F^2 Q,
H(w) = (d_w F)(d_w A) - A (d_w^2 F)  (degree 16).  Hence H(w) in J_16 is necessary.
H(w) = sum_{i<=j} w_i w_j H_ij.  If the 21 classes [H_ij] are linearly independent in R_16 = S_16/J_16
then w_i w_j = 0 for all i,j, i.e. w = 0: no certificate with deg D = 5 exists.
Independence is tested sector-wise (G-characters) mod two primes; dim J_16 = 19398 is forced by the regular
sequence (smoothness), so rank_p(J gens + H) = dim J_16 + 21 transfers to characteristic 0.
"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials
from flint import nmod_mat

def key_of(m, gens):
    return tuple(sum(a*b for a, b in zip(m, g)) % 6 for g in gens)

def poly_key(P, gens):
    ks = {key_of(m, gens) for m in P}
    assert len(ks) == 1, ks
    return ks.pop()

def main():
    t0 = time.time()
    F, A, q, G = w2lib.load()
    gens = G.generators
    F2 = scale(F, 2)
    Fi = [pderiv(F2, i) for i in range(6)]
    Ai = [pderiv(A, i) for i in range(6)]
    H = {}
    for i in range(6):
        for j in range(i, 6):
            if i == j:
                H[(i, i)] = add(mul(Fi[i], Ai[i]), scale(mul(A, pderiv(Fi[i], i)), -1))
            else:
                H[(i, j)] = add(add(mul(Fi[i], Ai[j]), mul(Fi[j], Ai[i])), scale(mul(A, pderiv(Fi[i], j)), -2))
    for k, P in H.items():
        assert P and core.degree(P) == 16, k
    keys = {k: poly_key(P, gens) for k, P in H.items()}
    S16 = monomials(16); S11 = monomials(11)
    bykey16 = {}
    for m in S16:
        bykey16.setdefault(key_of(m, gens), []).append(m)
    gens_by_key = {}
    for m in S11:
        for k in range(6):
            P = shift(Fi[k], m, 1)
            kk = key_of(next(iter(P)), gens)
            gens_by_key.setdefault(kk, []).append(P)
    out = {'primes': [w2lib.P1, w2lib.P2], 'sectors': {}, 'dimS16': len(S16), 'H_keys': {str(k): list(v) for k, v in keys.items()}}
    total_J = {p: 0 for p in (w2lib.P1, w2lib.P2)}
    ok_all = True
    for kk, mons in sorted(bykey16.items()):
        idx = {m: i for i, m in enumerate(mons)}
        Hs = [k for k in H if keys[k] == kk]
        gl = gens_by_key.get(kk, [])
        rec = {'dimS': len(mons), 'nJgens': len(gl), 'H': [list(k) for k in Hs]}
        for p in (w2lib.P1, w2lib.P2):
            R = len(gl) + len(Hs); C = len(mons)
            flat = [0] * (R * C)
            for r, P in enumerate(gl + [H[k] for k in Hs]):
                base = r * C
                for e, v in P.items():
                    flat[base + idx[e]] = int(v) % p
            M = nmod_mat(R, C, flat, p)
            rkJ = nmod_mat(len(gl), C, flat[:len(gl) * C], p).rank() if gl else 0
            rk = M.rank()
            rec[f'rankJ_p{p}'] = rkJ; rec[f'rankJH_p{p}'] = rk
            total_J[p] += rkJ
            if Hs:
                rec[f'independent_p{p}'] = (rk == rkJ + len(Hs))
                ok_all &= (rk == rkJ + len(Hs))
        out['sectors'][str(kk)] = rec
    out['total_dimJ16'] = total_J; out['expected_dimJ16_regular_sequence'] = len(S16) - 951
    out['all_H_independent_mod_J'] = ok_all
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open(f'{w2lib.PACKET}/hess5.receipt.json', 'w'), indent=1)
    w2lib.log(exp='hess5', family='all deg-5 denominators D = d_w F (forced by Ann_5 = 0)', d=5, prime=[w2lib.P1, w2lib.P2],
              result=('EXCLUDED: 21 classes [H_ij] independent in R_16 => no certificate with deg D = 5' if ok_all else 'NOT decided: some H_ij dependent mod J_16'),
              rank={'dimJ16': total_J, 'expected': len(S16) - 951}, seconds=out['seconds'], job='excl-cycle-q2-hess5',
              note='sector-wise ranks in wave2-excl/hess5.receipt.json; transfers to char 0 since dim_Q J_16 = dim_p J_16 (regular sequence) and rank_p <= rank_Q')
    print(json.dumps({k: v for k, v in out.items() if k != 'sectors'}, indent=1))
    for kk, rec in out['sectors'].items():
        if rec['H']:
            print(kk, rec)

if __name__ == '__main__':
    main()
