"""Transported 'contact' denominators (Chapter 2a packets) for four-plus-two-residual, mod p.
For the binary sextics f(X0,X1), g(Y0,Y1) of the core H = f + g + K find all F_p-rational contact presentations
f = C_f^3 + r_f^2 (C quadratic, r cubic; 40 over C), likewise g, by a vectorized scan over (c0,c1,c2) in F_p^3,
then test the 2a-type denominators (r_f^2 + r_g^2, r_f^2 + g, f + r_g^2, C_f^3 + C_g^3, products, +-B, +-H, ...)
for E1 (A D in J(F) mod p) and, when E1 passes, for the full sector-restricted fixed-D system mod p.
usage: python3 fpt_contact.py [p]   (p = 3 mod 4 prime, default 1019)"""
import os, sys, json, time, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '2'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
import numpy as np
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1019
assert p % 4 == 3
OUT = f'{LANE}/wave2-search/contact-p{p}.jsonl'
inst = Instance()
F, A, q = inst.F, inst.A, inst.q
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
inv7 = pow(7, -1, p); inv2 = pow(2, -1, p)
sq = np.full(p, -1, dtype=np.int64)
ar = np.arange(p, dtype=np.int64)
cand = np.array([pow(int(a), (p + 1) // 4, p) for a in range(p)], dtype=np.int64)
ok = (cand * cand) % p == ar
sq[ok] = cand[ok]

def presentations(s):
    """all (c0,c1,c2,r0,r1,r2,r3) mod p with sum s_k X^(6-k) Y^k = C^3 + r^2."""
    s = [x % p for x in s]
    sols = []
    C1, C2 = np.meshgrid(ar, ar, indexing='ij'); C1 = C1.ravel(); C2 = C2.ravel()
    for c0 in range(p):
        t = (s[0] - pow(c0, 3, p)) % p
        if sq[t] < 0:
            continue
        for r0 in ({0} if t == 0 else {int(sq[t]), (-int(sq[t])) % p}):
            if r0:
                i2r0 = pow(2 * r0, -1, p)
                r1 = ((s[1] - 3 * c0 * c0 % p * C1) % p) * i2r0 % p
                r2 = ((s[2] - 3 * c0 * c0 % p * C2 - 3 * c0 * (C1 * C1 % p) - r1 * r1) % p) * i2r0 % p
                r3 = ((s[3] - C1 * C1 % p * C1 - 6 * c0 * (C1 * C2 % p) - 2 * (r1 * r2 % p)) % p) * i2r0 % p
                e4 = (3 * c0 * (C2 * C2 % p) + 3 * (C1 * C1 % p) * C2 + r2 * r2 + 2 * (r1 * r3 % p) - s[4]) % p
                e5 = (3 * C1 * (C2 * C2 % p) + 2 * (r2 * r3 % p) - s[5]) % p
                e6 = (C2 * C2 % p * C2 + r3 * r3 - s[6]) % p
                hit = np.nonzero((e4 == 0) & (e5 == 0) & (e6 == 0))[0]
                for h in hit:
                    sols.append((c0, int(C1[h]), int(C2[h]), r0, int(r1[h]), int(r2[h]), int(r3[h])))
            else:
                # r0 = 0: e1: s1 = 3 c0^2 c1
                if c0 == 0:
                    if s[1] % p: continue
                    c1s = range(p)
                else:
                    c1s = [s[1] * pow(3 * c0 * c0, -1, p) % p]
                for c1 in c1s:
                    tt = (s[2] - 3 * c0 * c0 % p * ar - 3 * c0 * c1 * c1) % p   # r1^2 as function of c2
                    r1 = sq[tt]
                    for sign in (1, -1):
                        R1 = (sign * r1) % p
                        valid = r1 >= 0
                        # e3: 2 r1 r2 = s3 - c1^3 - 6 c0 c1 c2
                        rhs3 = (s[3] - pow(c1, 3, p) - 6 * c0 * c1 * ar) % p
                        with np.errstate(all='ignore'):
                            inv = np.array([pow(int(2 * x % p), -1, p) if x % p else 0 for x in R1], dtype=np.int64)
                        r2 = rhs3 * inv % p
                        zero_r1 = (R1 == 0)
                        # e4, e5, e6 with r3 from e5: 2 r2 r3 = s5 - 3 c1 c2^2  (needs r2 != 0)
                        rhs5 = (s[5] - 3 * c1 * (ar * ar % p)) % p
                        inv2r2 = np.array([pow(int(2 * x % p), -1, p) if x % p else 0 for x in r2], dtype=np.int64)
                        r3 = rhs5 * inv2r2 % p
                        e4 = (3 * c0 * (ar * ar % p) + 3 * c1 * c1 % p * ar + r2 * r2 + 2 * (R1 * r3 % p) - s[4]) % p
                        e6 = (ar * ar % p * ar + r3 * r3 - s[6]) % p
                        good = valid & (~zero_r1) & (r2 != 0) & (e4 == 0) & (e6 == 0) & (((2 * r2 * r3 + 3 * c1 * (ar * ar % p) - s[5]) % p) == 0)
                        for c2 in np.nonzero(good)[0]:
                            sols.append((c0, c1, int(c2), 0, int(R1[c2]), int(r2[c2]), int(r3[c2])))
    # dedupe
    return sorted(set(sols))

def binpoly(coeffs, vi, vj, deg):
    """sum coeffs[k] z_i^(deg-k) z_j^k as a 6-var poly mod p"""
    P = {}
    for k, c in enumerate(coeffs):
        if c % p:
            e = [0] * 6; e[vi] = deg - k; e[vj] = k
            P[tuple(e)] = c % p
    return P

fco = [1, 5 * inv7, 0, 0, 0, 10 * inv7, 1]
gco = [1, 5 * inv7, 0, 0, 0, 15 * inv7, 1]
t0 = time.time()
Pf = presentations(fco); print(f'f: {len(Pf)} F_p-rational presentations ({time.time()-t0:.0f}s)', flush=True)
t0 = time.time()
Pg = presentations(gco); print(f'g: {len(Pg)} F_p-rational presentations ({time.time()-t0:.0f}s)', flush=True)
log_line(OUT, {'kind': 'contact-presentations', 'p': p, 'n_f': len(Pf), 'n_g': len(Pg), 'f': Pf[:60], 'g': Pg[:60]})
# verify a presentation
def check(pres, co, vi, vj):
    c0, c1, c2, r0, r1, r2, r3 = pres
    C = binpoly([c0, c1, c2], vi, vj, 2); r = binpoly([r0, r1, r2, r3], vi, vj, 3)
    S = add(power(C, 3, p), mul(r, r, p), p)
    return S == binpoly(co, vi, vj, 6), C, r
nbad = 0
for pres in Pf:
    okk, _, _ = check(pres, fco, 0, 1); nbad += (not okk)
for pres in Pg:
    okk, _, _ = check(pres, gco, 3, 4); nbad += (not okk)
print('verification failures:', nbad, flush=True)
Hp = {e: c for e, c in Fp.items() if e[2] == 0 and e[5] == 0}
Bp = {e: c for e, c in Fp.items() if e[2] or e[5]}
fP = binpoly(fco, 0, 1, 6); gP = binpoly(gco, 3, 4, 6)

def e1_modp(D):
    return jacobian.in_jacobian(Fp, mul(Ap, D, p), p)[0]

def full_test(D, name, family):
    t0 = time.time()
    e1 = e1_modp(D)
    rec = {'kind': 'contact-test', 'p': p, 'D': name, 'family': family, 'd': core.degree(D), 'e1': bool(e1), 'engine': ENGINE}
    if e1:
        try:
            S, subg = symmetry.restricted_system(Fp, Ap, q, D, p, inst.G)
            r = S.solve()
            rec.update({'consistent': bool(r['consistent']), 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols, 'subgroup': len(subg)})
            if r['consistent']:
                print('!!!!! CONSISTENT mod p:', name, flush=True)
        except Exception as ex:
            rec['error'] = f'{ex.__class__.__name__}: {ex}'
    rec['seconds'] = round(time.time() - t0, 1)
    log_line(OUT, rec)
    print({k: rec[k] for k in rec if k not in ('kind',)}, flush=True)
    return rec

tests = []
for i, pf in enumerate(Pf[:12]):
    _, Cf, rf = check(pf, fco, 0, 1)
    rf2 = mul(rf, rf, p); Cf3 = power(Cf, 3, p)
    tests.append((rf2, f'r_f^2 (f-sheet {i})', 'contact-single'))
    tests.append((add(rf2, gP, p), f'r_f^2 + g (f-sheet {i})', 'contact-mixed'))
    tests.append((add(Cf3, gP, p), f'C_f^3 + g (f-sheet {i})', 'contact-mixed'))
    for jg, pg in enumerate(Pg[:12]):
        _, Cg, rg = check(pg, gco, 3, 4)
        rg2 = mul(rg, rg, p); Cg3 = power(Cg, 3, p)
        if i == 0:
            tests.append((rg2, f'r_g^2 (g-sheet {jg})', 'contact-single'))
            tests.append((add(fP, rg2, p), f'f + r_g^2 (g-sheet {jg})', 'contact-mixed'))
        for lam, lname in ((1, '+'), (p - 1, '-')):
            tests.append((add(rf2, scale(rg2, lam, p), p), f'r_f^2 {lname} r_g^2 (sheets {i},{jg})', 'contact-2a'))
        tests.append((add(Cf3, Cg3, p), f'C_f^3 + C_g^3 (sheets {i},{jg})', 'contact-2a'))
        tests.append((mul(rf, rg, p), f'r_f r_g (sheets {i},{jg})', 'contact-product'))
        tests.append((add(add(rf2, rg2, p), Bp, p), f'r_f^2 + r_g^2 + B (sheets {i},{jg})', 'contact-2a+B'))
        tests.append((add(add(rf2, rg2, p), Hp, p), f'r_f^2 + r_g^2 + H (sheets {i},{jg})', 'contact-2a+H'))
print(f'{len(tests)} candidate denominators; E1 screen first', flush=True)
n_e1 = 0
for D, name, fam in tests:
    if not D or core.parity_of(D) is None:
        log_line(OUT, {'kind': 'contact-test', 'p': p, 'D': name, 'family': fam, 'note': 'zero or mixed parity', 'engine': ENGINE}); continue
    rec = full_test(D, name, fam)
    n_e1 += bool(rec.get('e1'))
print(f'done: {n_e1} of {len(tests)} pass E1', flush=True)
