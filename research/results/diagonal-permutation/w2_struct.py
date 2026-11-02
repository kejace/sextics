"""Structured candidate denominators (degrees 6..12, Y-even) suggested by the lane's geometry, with the E1 test done
by exact membership mod p in the E1-admissible space  Adm_d = J(F)_d^+ (+ D8 at d=8) (+ kernel basis at d>=10 from
w2-kernel-d<d>-par0.json when available).  Admissible candidates get the full fixed-D solve (engine, compressed, exact).
usage: python3 w2_struct.py [dmin dmax]       appends to w2-struct.jsonl
Families: cones over plane curves g(X), g(Y) (vertex plane), cones over surfaces in 4 linear forms (vertex line),
discriminants/Hessians of restrictions, block radial derivatives, transposition-adapted products, companion-type
denominators transported from accepted witnesses (Fermat cofactors X^4 - X^2 Y^2 + Y^4 etc.), S3-orbit sums."""
import os, sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hodge import core, fixed_d, jacobian
from hodge.core import monomials, NV, pderiv, shift, mul, add, scale, parse, to_str, power
from flint import nmod_mat
from fractions import Fraction
LANE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(LANE, 'w2-struct.jsonl')
dmin = int(sys.argv[1]) if len(sys.argv) > 1 else 6
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
F, A, q, rec = core.load_instance('diagonal-permutation')
p = core.PRIMES[0]
PAR = [0, 0, 0, 1, 1, 1]
X = [parse(v) for v in ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']]
def log(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(OUT, 'a') as fh: fh.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:400], flush=True)

# ------------------------------------------------------------------ candidate generation (strings -> polys) ----
def S3_orbit_sum(P):
    from w2sym import sig_poly
    perms = [(0,1,2),(1,0,2),(0,2,1),(2,1,0),(1,2,0),(2,0,1)]
    out = {}
    for s in perms:
        Q = {}
        for e, c in P.items():
            e2 = (e[s[0]], e[s[1]], e[s[2]], e[3+s[0]], e[3+s[1]], e[3+s[2]])
            Q[e2] = Q.get(e2, 0) + c
        out = add(out, Q)
    return {e: v for e, v in out.items() if v}
cands = []   # (name, family, poly)
def addc(name, fam, P):
    if not P or not core.is_homogeneous(P) or core.parity_of(P) != 0: return
    d = core.degree(P)
    if dmin <= d <= dmax: cands.append((name, fam, P))
Fi = [pderiv(F, i) for i in range(NV)]
# ternary forms in X and in Y
def tern(vars3, name):
    a, b, c = vars3
    forms = {
        'fermat6': add(add(power(a,6), power(b,6)), power(c,6)),
        'fermat4': add(add(power(a,4), power(b,4)), power(c,4)),
        'fermat3': add(add(power(a,3), power(b,3)), power(c,3)),
        'e3': mul(mul(a,b),c), 'e2': add(add(mul(a,b), mul(b,c)), mul(a,c)), 'p2': add(add(mul(a,a), mul(b,b)), mul(c,c)),
        'hesse3': add(add(add(power(a,3), power(b,3)), power(c,3)), scale(mul(mul(a,b),c), Fraction(-3))),
        'vandermonde': mul(mul(add(a, scale(b,-1)), add(b, scale(c,-1))), add(a, scale(c,-1))),
    }
    out = {}
    for k, P in forms.items(): out[f'{name}:{k}'] = P
    out[f'{name}:fermat3^2'] = mul(forms['fermat3'], forms['fermat3'])
    out[f'{name}:hesse3^2'] = mul(forms['hesse3'], forms['hesse3'])
    out[f'{name}:vandermonde^2'] = mul(forms['vandermonde'], forms['vandermonde'])
    out[f'{name}:e3^2'] = mul(forms['e3'], forms['e3']); out[f'{name}:e3^3'] = power(forms['e3'], 3); out[f'{name}:e3^4'] = power(forms['e3'], 4)
    out[f'{name}:e2^3'] = power(forms['e2'], 3); out[f'{name}:p2^3'] = power(forms['p2'], 3); out[f'{name}:p2^4'] = power(forms['p2'], 4); out[f'{name}:p2^5'] = power(forms['p2'], 5)
    out[f'{name}:fermat4*p2'] = mul(forms['fermat4'], forms['p2']); out[f'{name}:fermat6*p2'] = mul(forms['fermat6'], forms['p2'])
    out[f'{name}:fermat6*p2^2'] = mul(forms['fermat6'], power(forms['p2'], 2)); out[f'{name}:fermat6^2'] = mul(forms['fermat6'], forms['fermat6'])
    out[f'{name}:fermat6*fermat4'] = mul(forms['fermat6'], forms['fermat4']); out[f'{name}:fermat4^2'] = power(forms['fermat4'], 2); out[f'{name}:fermat4^3'] = power(forms['fermat4'], 3)
    out[f'{name}:fermat3^4'] = power(forms['fermat3'], 4); out[f'{name}:fermat3^3'] = power(forms['fermat3'], 3)
    # discriminant-type: product over the 3 differences squared times p2
    out[f'{name}:vandermonde^2*p2'] = mul(out[f'{name}:vandermonde^2'], forms['p2']); out[f'{name}:vandermonde^2*p2^2'] = mul(out[f'{name}:vandermonde^2'], power(forms['p2'],2))
    out[f'{name}:vandermonde^4'] = power(forms['vandermonde'], 4)
    return out
for k, P in tern(X[:3], 'coneX').items(): addc(k, 'cone-vertex-plane-Y', P)
for k, P in tern(X[3:], 'coneY').items(): addc(k, 'cone-vertex-plane-X', P)
# mixed: g(X) * h(Y) products (Y-even requires even Y-degree)
addc('coneX:fermat6+coneY:fermat6', 'sum-cones', add(tern(X[:3],'a')['a:fermat6'], tern(X[3:],'b')['b:fermat6']))
addc('coneX:fermat6-coneY:fermat6', 'sum-cones', add(tern(X[:3],'a')['a:fermat6'], scale(tern(X[3:],'b')['b:fermat6'], -1)))
addc('coneX:fermat3^2 + coneY:fermat3^2', 'sum-cones', add(tern(X[:3],'a')['a:fermat3^2'], tern(X[3:],'b')['b:fermat3^2']))
addc('e3(X)^2 + e3(Y)^2', 'sum-cones', add(tern(X[:3],'a')['a:e3^2'], tern(X[3:],'b')['b:e3^2']))
addc('e3(X)^2 - e3(Y)^2', 'sum-cones', add(tern(X[:3],'a')['a:e3^2'], scale(tern(X[3:],'b')['b:e3^2'], -1)))
addc('e3(X)*e3(Y)*(e3X*e3Y)', 'sum-cones', power(mul(tern(X[:3],'a')['a:e3'], tern(X[3:],'b')['b:e3']), 2))
# cones with vertex a line: forms in 4 linear forms (X0,X1,Y0,Y1), (X0,X2,Y0,Y2), (X0,Y0,X1,Y1) with Fermat / products
def quat(vs, name):
    out = {}
    out[f'{name}:fermat4'] = {}
    for v in vs: out[f'{name}:fermat4'] = add(out[f'{name}:fermat4'], power(v, 4))
    out[f'{name}:fermat6'] = {}
    for v in vs: out[f'{name}:fermat6'] = add(out[f'{name}:fermat6'], power(v, 6))
    out[f'{name}:prod2'] = power(mul(mul(vs[0], vs[1]), mul(vs[2], vs[3])), 2)
    out[f'{name}:fermat4^2'] = power(out[f'{name}:fermat4'], 2)
    out[f'{name}:fermat6*fermat4'] = mul(out[f'{name}:fermat6'], out[f'{name}:fermat4'])
    out[f'{name}:fermat4^3'] = power(out[f'{name}:fermat4'], 3)
    return out
for vs, nm in [([X[0],X[1],X[3],X[4]], 'L01'), ([X[0],X[2],X[3],X[5]], 'L02'), ([X[1],X[2],X[4],X[5]], 'L12')]:
    for k, P in quat(vs, nm).items(): addc(k, 'cone-vertex-line', P)
# pair blocks T_i and A-related
T = [mul(power(X[i],3), power(X[3+i],3)) for i in range(3)]
B2 = add(add(T[0], T[1]), scale(T[2], -2))
addc('A^2', 'target-related', mul(A, A)); addc('B2^2', 'target-related', mul(B2, B2)); addc('A*B2', 'target-related', mul(A, B2))
addc('sum T_i^2', 'target-related', add(add(mul(T[0],T[0]), mul(T[1],T[1])), mul(T[2],T[2])))
addc('(sum T_i)^2', 'target-related', power(add(add(T[0],T[1]),T[2]), 2))
# Fermat cofactors (companion denominators from fermat-pair-rational-d4-strict) and their orbit sums / products
cof = [add(add(power(X[i],4), scale(mul(power(X[i],2), power(X[3+i],2)), -1)), power(X[3+i],4)) for i in range(3)]
addc('cof0*cof1', 'fermat-cofactor', mul(cof[0], cof[1])); addc('cof0*cof1*cof2', 'fermat-cofactor', mul(mul(cof[0], cof[1]), cof[2]))
addc('cof0^2*cof1^2... (cof0 cof1)^2', 'fermat-cofactor', power(mul(cof[0], cof[1]), 2))
addc('cof0*cof1*p2(X)', 'fermat-cofactor', mul(mul(cof[0], cof[1]), tern(X[:3],'a')['a:p2']))
addc('sum cof_i cof_j', 'fermat-cofactor', add(add(mul(cof[0],cof[1]), mul(cof[1],cof[2])), mul(cof[0],cof[2])))
for i in range(3):
    addc(f'cof{i}^2', 'fermat-cofactor', mul(cof[i], cof[i])); addc(f'cof{i}^3', 'fermat-cofactor', power(cof[i], 3))
sq = [add(power(X[i],2), power(X[3+i],2)) for i in range(3)]       # X_i^2 + Y_i^2
dq = [add(power(X[i],2), scale(power(X[3+i],2), -1)) for i in range(3)]
addc('prod (X_i^2+Y_i^2)', 'pair-quadrics', mul(mul(sq[0], sq[1]), sq[2])); addc('prod (X_i^2-Y_i^2)', 'pair-quadrics', mul(mul(dq[0], dq[1]), dq[2]))
addc('prod (X_i^2+Y_i^2)^2', 'pair-quadrics', power(mul(mul(sq[0], sq[1]), sq[2]), 2)); addc('prod (X_i^2-Y_i^2)^2', 'pair-quadrics', power(mul(mul(dq[0], dq[1]), dq[2]), 2))
addc('(X0^2+Y0^2)^3 (X1^2+Y1^2)^3', 'pair-quadrics', mul(power(sq[0],3), power(sq[1],3)))
# block radial derivatives and their products (Jacobian-type, d=6,12) and polar products
EX = {}; EY = {}
for i in range(3): EX = add(EX, mul(X[i], Fi[i])); EY = add(EY, mul(X[3+i], Fi[3+i]))
addc('E_X(F)^2', 'block-radial', mul(EX, EX)); addc('E_Y(F)^2', 'block-radial', mul(EY, EY)); addc('E_X(F) E_Y(F)', 'block-radial', mul(EX, EY))
addc('E_X(F)*p2(X)^2', 'block-radial', mul(EX, power(tern(X[:3],'a')['a:p2'], 2))); addc('E_X(F)*fermat4(X)', 'block-radial', mul(EX, tern(X[:3],'a')['a:fermat4']))
addc('E_X(F)*e3(X)^2', 'block-radial', mul(EX, tern(X[:3],'a')['a:e3^2'])); addc('E_Y(F)*e3(Y)^2', 'block-radial', mul(EY, tern(X[3:],'a')['a:e3^2']))
addc('E_X(F)*fermat6(X)', 'block-radial', mul(EX, tern(X[:3],'a')['a:fermat6'])); addc('E_Y(F)*fermat6(Y)', 'block-radial', mul(EY, tern(X[3:],'a')['a:fermat6']))
addc('F_X0*F_X1 (polar product)', 'polar-products', mul(Fi[0], Fi[1])); addc('F_Y0*F_Y1', 'polar-products', mul(Fi[3], Fi[4])); addc('F_X0*F_Y0', 'polar-products', mul(Fi[0], Fi[3]))
addc('F_X0^2', 'polar-products', mul(Fi[0], Fi[0])); addc('F_X2*F_Y2', 'polar-products', mul(Fi[2], Fi[5])); addc('(F_X0-F_X1)^2', 'polar-products', power(add(Fi[0], scale(Fi[1],-1)), 2))
addc('(F_X0-F_X1)(F_Y0-F_Y1)', 'polar-products', mul(add(Fi[0], scale(Fi[1],-1)), add(Fi[3], scale(Fi[4],-1))))
addc('sum_i F_Xi^2', 'polar-products', add(add(mul(Fi[0],Fi[0]), mul(Fi[1],Fi[1])), mul(Fi[2],Fi[2])))
addc('sum_i F_Yi^2', 'polar-products', add(add(mul(Fi[3],Fi[3]), mul(Fi[4],Fi[4])), mul(Fi[5],Fi[5])))
addc('F_X0 F_X1 F_X2 ... (F_X0 F_X1)(F_X2^2)?', 'polar-products', mul(mul(Fi[0], Fi[1]), mul(Fi[2], Fi[2])))
# transposition-adapted: products over the three transposition quadrics
r = [add(X[0], scale(X[1],-1)), add(X[1], scale(X[2],-1)), add(X[0], scale(X[2],-1))]
s = [add(X[3], scale(X[4],-1)), add(X[4], scale(X[5],-1)), add(X[3], scale(X[5],-1))]
qp = {}; qm = {}
for k in range(3):
    qp[k] = add(mul(r[k], r[k]), mul(s[k], s[k])); qm[k] = add(mul(r[k], r[k]), scale(mul(s[k], s[k]), -1))
addc('prod (r^2+s^2)', 'transposition', mul(mul(qp[0], qp[1]), qp[2])); addc('prod (r^2-s^2)', 'transposition', mul(mul(qm[0], qm[1]), qm[2]))
addc('prod (r^2+s^2)^2', 'transposition', power(mul(mul(qp[0], qp[1]), qp[2]), 2)); addc('prod(r^2+s^2) prod(r^2-s^2)', 'transposition', mul(mul(mul(qp[0], qp[1]), qp[2]), mul(mul(qm[0], qm[1]), qm[2])))
addc('Delta_X^2 Delta_Y^2', 'transposition', mul(power(mul(mul(r[0],r[1]),r[2]),2), power(mul(mul(s[0],s[1]),s[2]),2)))
addc('(r0^2+s0^2)^3 * (r0 s0)^... (r0^2+s0^2)^4', 'transposition', power(qp[0], 4)); addc('(r0^2+s0^2)^3', 'transposition', power(qp[0], 3))
addc('(r0 s0)^2 (r1 s1)^2 (r2 s2)^2', 'transposition', power(mul(mul(mul(r[0],s[0]), mul(r[1],s[1])), mul(r[2],s[2])), 2))
# Hessian determinants of restrictions: Hess of F restricted to the (X0,X1,Y0,Y1) block etc. (degree 16: too big) -> skip; use det of 2x2 Hessian blocks
H2 = []
for (i, j) in [(0, 3), (0, 1), (3, 4), (0, 4)]:
    Fij = pderiv(Fi[i], j); Fii = pderiv(Fi[i], i); Fjj = pderiv(Fi[j], j)
    H2.append(((i, j), add(mul(Fii, Fjj), scale(mul(Fij, Fij), -1))))
for (i, j), P in H2: addc(f'2x2 Hessian det ({i},{j})', 'hessian-blocks', P)
# S3-orbit sums of a few non-symmetric ones
addc('S3-sum cof0^2', 'S3-orbit', S3_orbit_sum(mul(cof[0], cof[0]))); addc('S3-sum F_X0^2', 'S3-orbit', S3_orbit_sum(mul(Fi[0], Fi[0])))
addc('S3-sum (X0^2+Y0^2)^3(X1^2+Y1^2)^3', 'S3-orbit', S3_orbit_sum(mul(power(sq[0],3), power(sq[1],3))))
addc('S3-sum L01:fermat4^2', 'S3-orbit', S3_orbit_sum(power(quat([X[0],X[1],X[3],X[4]],'a')['a:fermat4'], 2)))
# the coupling orbits O1,O2,O3 and products
O1 = {}; O2 = {}; O3 = {}
for (i, j, k) in [(0,1,2),(1,0,2),(0,2,1),(2,0,1),(1,2,0),(2,1,0)]:
    O1 = add(O1, mul(mul(power(X[i],3), X[j]), mul(X[3+i], X[3+j]))); O2 = add(O2, mul(mul(power(X[i],3), X[j]), power(X[3+k],2))); O3 = add(O3, mul(mul(power(X[3+i],3), X[3+j]), power(X[k],2)))
for nm, P in [('O1', O1), ('O2', O2), ('O3', O3)]:
    addc(nm, 'coupling-orbits', P); addc(nm + '^2', 'coupling-orbits', mul(P, P))
addc('O1*O2', 'coupling-orbits', mul(O1, O2)); addc('O1*O3', 'coupling-orbits', mul(O1, O3)); addc('O2*O3', 'coupling-orbits', mul(O2, O3))
addc('F - Fermat part (coupling)', 'coupling-orbits', add(add(scale(O1, Fraction(5,7)), scale(O2, Fraction(5,11))), scale(O3, Fraction(5,13))))
print(f'{len(cands)} candidates in degrees {dmin}..{dmax}', flush=True)

# ------------------------------------------------------------------ admissibility test mod p -------------------
Fp = core.reduce_mod(F, p); Fip = [pderiv(Fp, i, p) for i in range(NV)]
adm_cache = {}
def admissible_space(d):
    """rows (mod p) spanning Adm_d inside S_d^+ : J_d^+ generators + D8 (d=8) + kernel basis file (if present)."""
    if d in adm_cache: return adm_cache[d]
    mons = monomials(d, 0); idx = {m: j for j, m in enumerate(mons)}
    rows = []
    for i in range(NV):
        for m in monomials(d - 5, PAR[i]) if d >= 5 else []:
            v = [0] * len(mons)
            for e, c in shift(Fip[i], m, 1, p).items(): v[idx[e]] = c
            rows.append(v)
    extra_note = ''
    if d == 8:
        D8 = core.reduce_mod(core.parse(json.load(open(os.path.join(LANE, 'D8_e1kernel.json')))['D']), p)
        v = [0] * len(mons)
        for e, c in D8.items(): v[idx[e]] = c
        rows.append(v); extra_note = '+D8'
    kf = os.path.join(LANE, f'w2-kernel-d{d}-par0.json')
    if os.path.exists(kf):
        kz = json.load(open(kf))
        if kz['p'] == p:
            for b in kz['basis_modp']:
                v = [0] * len(mons)
                for ks, c in b.items(): v[idx[tuple(int(x) for x in ks.split(','))]] = c
                rows.append(v)
            extra_note += f'+kernel({len(kz["basis_modp"])})'
        else:
            extra_note += '+kernel-file-prime-mismatch'
    elif d >= 10:
        extra_note += '+KERNEL-MISSING(J-only test: a non-admissible verdict is NOT rigorous)'
    M = nmod_mat(len(rows), len(mons), sum(rows, []), p) if rows else None
    rk = M.rank() if M is not None else 0
    adm_cache[d] = (mons, idx, rows, rk, extra_note)
    return adm_cache[d]
def is_admissible(D):
    d = core.degree(D)
    mons, idx, rows, rk, note = admissible_space(d)
    Dp = core.reduce_mod(D, p)
    v = [0] * len(mons)
    for e, c in Dp.items(): v[idx[e]] = c
    rk2 = nmod_mat(len(rows) + 1, len(mons), sum(rows, []) + v, p).rank()
    return rk2 == rk, note
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(json.loads(line)['D'])
        except Exception: pass
for name, fam, P in cands:
    if name in done: continue
    t0 = time.time()
    d = core.degree(P)
    adm, note = is_admissible(P)
    e = {'D': name, 'd': d, 'family': fam, 'terms': len(P), 'e1': adm, 'adm_note': note, 'consistent': None, 'rank': -1, 'nullity': -1}
    divF = False
    if d >= 6:
        divF, _ = jacobian.divisible_by(F, P)
        e['F_divides_D'] = divF
    if adm and not divF:
        try:
            out = fixed_d.solve_fixed_D_fast(F, A, q, P, want_exact=True)
            e.update({'consistent': bool(out['consistent']), 'rank': int(out['rank']), 'nullity': int(out['nullity']), 'rows': out['rows'], 'cols': out['cols']})
            if out['consistent'] and out.get('exact'):
                w = {'d': d, 'D': to_str(P), 'N': [to_str(n) for n in out['N']], 'Q': to_str(out['Q'])}
                json.dump(w, open(os.path.join(LANE, 'witness-candidate.json'), 'w'), indent=1)
                e['note'] = 'EXACT WITNESS'
        except Exception as ex:
            e['note'] = f'solver error {ex!r}'
    e['seconds'] = round(time.time() - t0, 1)
    log(e)
print('STRUCT DONE', flush=True)
