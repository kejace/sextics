"""Sealed exclusion of structured candidate denominators (goal c): D8 products, polar products, structural pieces of F,
coordinate multiples, up to degree 12 (sector solver + modular-rank argument with exact null-vector sealing).
usage: python3 families.py threads [maxdeg]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from multiprocessing import Pool
threads = int(sys.argv[1]); maxdeg = int(sys.argv[2]) if len(sys.argv) > 2 else 12
S = Setup(); F, A = S.F, S.A
P = core.parse; M = core.mul
t = Fraction(5, 7)
X0, X1, X2, Y0, Y1, Y2 = [P(v) for v in core.VARS]
C = P('X0**3 + X0*Y1*Y2 + X1*Y0*Y1')
H = P('X0**6 + X1**6 + Y0**6 + Y1**6 + Y2**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0**2*Y0*Y1*Y2**2/7 + 5*X0*X1*Y1**3*Y2/7')
Z = P('X2**3'); Zp = core.add(core.scale(Z, 2), core.scale(C, t)); Hp = core.sub(H, core.scale(M(C, C), t * t / 4))
Fi = [core.pderiv(F, i) for i in range(6)]
D8 = P(json.load(open(WD + '/kernel_d8_0_4_0.json'))['D'])
cands = []
def add(D, name, fam):
    if core.degree(D) <= maxdeg and core.parity_of(D) is not None: cands.append((D, name, fam))
# D8 times coordinates / squares / structural pieces (degrees 9..12)
for i, v in enumerate(core.VARS):
    add(M(P(v), D8), f'{v}*D8', 'D8-products'); add(M(M(P(v), P(v)), D8), f'{v}^2*D8', 'D8-products')
for i in range(6):
    for j in range(i + 1, 6):
        add(M(M(P(core.VARS[i]), P(core.VARS[j])), D8), f'{core.VARS[i]}*{core.VARS[j]}*D8', 'D8-products')
add(M(C, D8), 'C*D8', 'D8-products'); add(M(Zp, D8), "Z'*D8", 'D8-products'); add(M(Z, D8), 'X2^3*D8', 'D8-products')
add(M(X2, M(X2, M(X2, M(X2, D8)))), 'X2^4*D8', 'D8-products')
for i in range(6): add(M(Fi[i], D8), f'F_{core.VARS[i]}*D8', 'D8-products')   # degree 13 > 12 skipped unless maxdeg >= 13
# polar products and structural pieces of F (degrees 8..12)
for i, j in [(0, 3), (2, 2), (0, 0), (3, 3), (1, 4), (2, 5), (0, 2), (2, 3), (1, 1), (4, 4), (5, 5), (0, 1), (3, 4), (3, 5), (4, 5)]:
    add(M(Fi[i], Fi[j]), f'F_{core.VARS[i]}*F_{core.VARS[j]}', 'polar-products')
for D, name, fam in [(M(M(X2, X2), H), 'X2^2*H', 'pieces-8'), (M(M(X2, X2), Hp), "X2^2*H'", 'pieces-8'), (M(M(X2, X2), M(C, Zp)), "X2^2*C*Z'", 'pieces-8'),
                     (M(C, M(C, M(X2, X2))), 'X2^2*C^2', 'pieces-8'), (M(Zp, M(Zp, M(X2, X2))), "X2^2*Z'^2", 'pieces-8'),
                     (M(Z, H), 'X2^3*H', 'pieces-9'), (M(Z, Hp), "X2^3*H'", 'pieces-9'), (M(C, Hp), "C*H'", 'pieces-9'), (M(Zp, Hp), "Z'*H'", 'pieces-9'),
                     (M(C, M(C, C)), 'C^3', 'pieces-9'), (M(Zp, M(Zp, Zp)), "Z'^3", 'pieces-9'), (M(C, M(C, Zp)), "C^2*Z'", 'pieces-9'), (M(C, M(Zp, Zp)), "C*Z'^2", 'pieces-9'),
                     (M(H, M(C, Zp)), "H*C*Z'", 'pieces-12'), (M(H, Hp), "H*H'", 'pieces-12'), (M(H, M(C, C)), 'H*C^2', 'pieces-12'), (M(Hp, M(C, Zp)), "H'*C*Z'", 'pieces-12'),
                     (M(H, M(Zp, Zp)), "H*Z'^2", 'pieces-12'), (M(H, H), 'H^2', 'pieces-12'), (M(Hp, Hp), "H'^2", 'pieces-12'), (M(C, M(C, M(C, C))), 'C^4', 'pieces-12'),
                     (M(M(X2, X2), M(X2, M(X2, M(X2, X2)))), 'X2^6', 'coordinate-powers'), (M(M(X0, X0), M(X0, M(X0, M(X0, X0)))), 'X0^6', 'coordinate-powers'),
                     (M(M(Y0, Y0), M(Y0, M(Y0, M(Y0, Y0)))), 'Y0^6', 'coordinate-powers'), (M(M(Y1, Y1), M(Y1, M(Y1, M(Y1, Y1)))), 'Y1^6', 'coordinate-powers'),
                     (M(M(Y2, Y2), M(Y2, M(Y2, M(Y2, Y2)))), 'Y2^6', 'coordinate-powers'), (M(M(X1, X1), M(X1, M(X1, M(X1, X1)))), 'X1^6', 'coordinate-powers')]:
    add(D, name, fam)
# Hessian determinant is degree 24 (too big); gradient-related: sum_i F_i^2 (degree 10), sum X_i^2 F_i^2? skip
add(core.add(core.add(M(Fi[0], Fi[0]), M(Fi[1], Fi[1])), M(Fi[2], Fi[2])), 'sum_X F_i^2', 'polar-products')
add(core.add(core.add(M(Fi[3], Fi[3]), M(Fi[4], Fi[4])), M(Fi[5], Fi[5])), 'sum_Y F_i^2', 'polar-products')
print(f'{len(cands)} candidates (degrees {sorted(set(core.degree(D) for D,_,_ in cands))})', flush=True)
def run(k):
    D, name, fam = cands[k]
    try:
        out = sealed_exclusion(F, A, 1, D)
    except Exception as ex:
        return k, {'error': str(ex)[:300]}
    return k, out
with Pool(threads) as pool:
    for k, out in pool.imap_unordered(run, range(len(cands))):
        D, name, fam = cands[k]
        if 'error' in out:
            log({'exp': 'structured-family-sealed-test', 'D': name, 'd': core.degree(D), 'family': fam, 'prime': None, 'result': 'ERROR ' + out['error'], 'rank': None, 'seconds': 0, 'job': 'excl-ctr-families'})
            print(name, 'ERROR', out['error'], flush=True); continue
        res = ('SEALED INCONSISTENT over Q' if out['sealed_inconsistent'] else ('CONSISTENT mod p ***' if out['consistent_mod_p'] else 'inconsistent mod p, not sealed'))
        log({'exp': 'structured-family-sealed-test', 'D': name, 'd': core.degree(D), 'family': fam, 'prime': ','.join(map(str, out['primes'])), 'result': res,
             'rank': out['rank'][0], 'seconds': out['seconds'], 'job': 'excl-ctr-families',
             'note': f"sector system rows={out['rows']} cols={out['cols']} |G|={out['G_order']} nullity={out['nullity']} rank_aug={out['rank_aug']} verified_null={out['verified_null']}"})
        print(f"{name:22s} d={core.degree(D):2d} {res:32s} rank={out['rank']} aug={out['rank_aug']} null={out['nullity']} sealed={out['verified_null']} {out['seconds']}s", flush=True)
