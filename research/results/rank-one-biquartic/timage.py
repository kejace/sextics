"""target_image for a list of denominators; saves served-space bases mod p and logs dims."""
import sys, time, json, pickle
import lib
from lib import *
from flint import nmod_mat
p = core.PRIMES[0]
Ap = core.reduce_mod(A, p)
odd6 = core.monomials(6, 1)
idx6 = {m: k for k, m in enumerate(odd6)}

def in_span(basis, target):
    """Is target (poly mod p) in the span of basis (list of polys mod p)?"""
    if not basis: return False
    R = len(odd6); C = len(basis)
    flat = [0] * (R * (C + 1))
    for j, b in enumerate(basis):
        for e, v in b.items(): flat[idx6[e] * (C + 1) + j] = v
    for e, v in target.items(): flat[idx6[e] * (C + 1) + C] = v
    Mx = nmod_mat(R, C + 1, flat, p)
    M0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p)
    return Mx.rank() == M0.rank()

def jac6_basis():
    Fp = core.reduce_mod(F, p)
    out = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(1, (1 - pi) % 2):
            out.append(core.shift(core.pderiv(Fp, i, p), m, 1, p))
    return out

def run_timage(D, expr, family, note=''):
    t0 = time.time()
    try:
        dim, basis, S = fixed_d.target_image(F, q, D)
        hasA = in_span(basis, Ap)
        res = {'D': expr, 'd': core.degree(D), 'family': family, 'e1': None, 'consistent': hasA,
               'rank': dim, 'nullity': dim - 18, 'seconds': round(time.time() - t0, 1),
               'note': f'target_image: served dim={dim} (geometric {dim-18}); A in span={hasA}. ' + note}
        with open(f'{WORK}/timage_{family}_{abs(hash(expr))%10**8}.pkl', 'wb') as fh:
            pickle.dump({'D': expr, 'dim': dim, 'basis': basis}, fh)
    except Exception as ex:
        res = {'D': expr, 'd': core.degree(D), 'family': family, 'e1': None, 'consistent': None, 'rank': -1,
               'nullity': -1, 'seconds': round(time.time() - t0, 1), 'note': f'target_image error {ex!r}'}
    log(res); show(res)
    return res

if __name__ == '__main__':
    name = sys.argv[1]
    Fi = [core.pderiv(F, i) for i in range(6)]
    B = {}
    B['T1'] = [(P('Y0'), 'Y0', 'timage', 'calibration: hyperplane should serve only J_6'),
               (P('X0'), 'X0', 'timage', ''),
               (L, 'L', 'timage', ''), (M, 'M', 'timage', ''),
               (f, 'f(X)', 'timage', 'geometric prediction: 3 (p_g of double octic plane)'),
               (g, 'g(Y)', 'timage', 'geometric prediction: 0'),
               (RX, 'eps_X F', 'timage', 'ramification divisor R'),
               (core.add(f, g), 'f+g', 'timage', ''),
               (Pc, 'L^2 M', 'timage', '')]
    B['T2'] = [(Fi[0], 'F_0', 'timage', ''), (Fi[1], 'F_1', 'timage', ''), (Fi[2], 'F_2', 'timage', ''),
               (Fi[3], 'F_3', 'timage', 'Y-odd D'), (Fi[4], 'F_4', 'timage', 'Y-odd D'), (Fi[5], 'F_5', 'timage', 'Y-odd D'),
               (core.mul(L, M), 'L M', 'timage', ''),
               (core.mul(f, core.scale(Pc, 1)), 'f * L^2M', 'timage', 'd=12 too big? skip if slow')]
    for D, expr, fam, note in B[name]:
        if core.degree(D) > 8: continue
        run_timage(D, expr, fam, note)
