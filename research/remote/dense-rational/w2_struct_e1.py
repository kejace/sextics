"""Structured linear families vs the E1-admissible space at degree d (mod p).
For a family V (Y-even forms of degree d spanned by a monomial set or explicit polynomials) compute
   dim (V cap K_full)   where K_full = {D : A D in J(F)}   (exact mod p; = dim over Q for all but finitely many p)
   dim (V cap J_d)      (trivial Jacobian part)
   dim (V cap K_full) - dim (V cap J_d)  = non-Jacobian E1-admissible directions inside V.
Reducers: normal forms modulo J_{d+6}^- and J_d^+ (hodge.jacobian._Reducer, rref mod p).
Families (geometric meaning in the report): quaternary cones (vertex a line), ternary cones (vertex a plane),
binary-block forms, forms in (X-block) x (Y-block) products, block products of F's bidegree parts, |grad F|^2, ...
Usage: python3 w2_struct_e1.py d [p]"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2lib import *
from flint import nmod_mat
d = int(sys.argv[1]); p = int(sys.argv[2]) if len(sys.argv) > 2 else core.PRIMES[0]
F, A, q = load()
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
t0 = time.time()
red_hi = jacobian._Reducer(F, d + 6, p, 1)      # normal forms in S_{d+6}^- modulo J
print(f'[reducer d+6={d+6} odd] J dim {red_hi.jacobian_dim} ring dim {red_hi.ring_dim} ({time.time()-t0:.1f}s)', flush=True)
t1 = time.time()
red_lo = jacobian._Reducer(F, d, p, 0)          # normal forms in S_d^+ modulo J
print(f'[reducer d={d} even] J dim {red_lo.jacobian_dim} ring dim {red_lo.ring_dim} ({time.time()-t1:.1f}s)', flush=True)
hi_idx = {m: i for i, m in enumerate(red_hi.standard)}
lo_idx = {m: i for i, m in enumerate(red_lo.standard)}


def nullity_of(polys, reducer, idx, mulA):
    """dim of {c : sum c_k P_k reduces to 0}, i.e. dim(span P cap J) (mulA: reduce A*P instead)."""
    n = len(polys)
    if n == 0:
        return 0, 0
    rows = []
    for P in polys:
        Q = mul(Ap, P, p) if mulA else P
        nf = reducer.normal_form(Q)
        rows.append([0] * len(idx) if not nf else None)
        if nf:
            r = [0] * len(idx)
            for e, v in nf.items():
                r[idx[e]] = v
            rows[-1] = r
    M = nmod_mat(n, len(idx), [v for r in rows for v in r], p)
    rk = M.rank()
    return n - rk, rk


def span_dim(polys, mons):
    idx = {m: i for i, m in enumerate(mons)}
    M = nmod_mat(len(polys), len(mons), [0] * (len(polys) * len(mons)), p)
    flat = [0] * (len(polys) * len(mons))
    for r, P in enumerate(polys):
        for e, v in P.items():
            flat[r * len(mons) + idx[e]] = v
    return nmod_mat(len(polys), len(mons), flat, p).rank()


mons_d = monomials(d, 0)
fams = {}
# monomial-closed families: forms in a subset of variables (Y-even)
subsets = {'quaternary X0X1X2Y0 (vertex line Y1=Y2=0... i.e. {X=0,Y0=0})': (0, 1, 2, 3),
           'quaternary X0X1Y0Y1': (0, 1, 3, 4), 'quaternary X0Y0Y1Y2 (vertex line {X1=X2=0,Y=0})': (0, 3, 4, 5),
           'quaternary X1X2Y1Y2': (1, 2, 4, 5),
           'ternary X (vertex plane P_Y)': (0, 1, 2), 'ternary Y (vertex plane P_X)': (3, 4, 5),
           'ternary X0X1Y0': (0, 1, 3), 'ternary X0Y0Y1': (0, 3, 4),
           'binary X0Y0': (0, 3), 'binary X0X1': (0, 1), 'binary Y0Y1': (3, 4),
           'quinary X0X1X2Y0Y1 (cone over a 3-fold, vertex a point)': (0, 1, 2, 3, 4),
           'quinary X0X1Y0Y1Y2': (0, 1, 3, 4, 5)}
for name, vs in subsets.items():
    fams[name] = [{m: 1} for m in mons_d if all(m[i] == 0 for i in range(NV) if i not in vs)]
# products of block forms: (ternary in X of degree a) * (ternary in Y of degree d-a), a even
for a in range(0, d + 1, 2):
    fam = []
    for m in mons_d:
        if sum(m[:3]) == a:
            fam.append({m: 1})
    fams[f'bidegree ({a},{d-a}) forms'] = fam
# explicit structured polynomials
Fi = [pderiv(F, i) for i in range(NV)]
blocks = {}
for m, c in F.items():
    k = sum(m[3:]); blocks.setdefault(k, {})[m] = c
B60, B42, B24, B06 = blocks[0], blocks[2], blocks[4], blocks[6]
named = {}
if d == 10:
    named['|grad F|^2 = sum F_i^2'] = sum_ = None
    s = {}
    for i in range(NV):
        s = add(s, mul(Fi[i], Fi[i]))
    named['|grad F|^2 = sum F_i^2'] = s
    named['sum_{i<3} F_i^2 (X-part)'] = add(add(mul(Fi[0], Fi[0]), mul(Fi[1], Fi[1])), mul(Fi[2], Fi[2]))
    named['sum_{i>=3} F_i^2 (Y-part)'] = add(add(mul(Fi[3], Fi[3]), mul(Fi[4], Fi[4])), mul(Fi[5], Fi[5]))
    named['F_X0^2'] = mul(Fi[0], Fi[0]); named['F_Y0^2'] = mul(Fi[3], Fi[3]); named['F_X0 F_X1'] = mul(Fi[0], Fi[1])
    named['F_Y0 F_Y1'] = mul(Fi[3], Fi[4])
    named['X0^4 B60'] = mul(core.parse('X0**4'), B60); named['Y0^4 B06'] = mul(core.parse('Y0**4'), B06)
    named['X0^2 Y0^2 B42'] = mul(core.parse('X0**2*Y0**2'), B42)
    named['(X0^2+X1^2+X2^2)^2 B60'] = mul(power(core.parse('X0**2+X1**2+X2**2'), 2), B60)
    named['(Y0^2+Y1^2+Y2^2)^2 B06'] = mul(power(core.parse('Y0**2+Y1**2+Y2**2'), 2), B06)
    named['X0^4 F_X0^... no: X0^3 F_X0 X0 = X0^4 F_X0 (deg 9)'] = None
    named = {k: v for k, v in named.items() if v}
if d == 12:
    prods = {'B60*B06': (B60, B06), 'B60*B24': (B60, B24), 'B42*B06': (B42, B06), 'B42*B24': (B42, B24),
             'B60*B42': (B60, B42), 'B24*B06': (B24, B06), 'B60^2': (B60, B60), 'B42^2': (B42, B42),
             'B24^2': (B24, B24), 'B06^2': (B06, B06)}
    for k, (P, Q) in prods.items():
        named[k] = mul(P, Q)
    named['(B60+B06)^2'] = mul(add(B60, B06), add(B60, B06))
    named['(B42+B24)^2'] = mul(add(B42, B24), add(B42, B24))
    named['B60*B06 - B42*B24'] = sub(mul(B60, B06), mul(B42, B24))
    named['Hess 2x2 minor (X0X0,Y0Y0)^... deg 8: skip'] = None
    named = {k: v for k, v in named.items() if v}
    named['span of the 10 block products'] = [mul(P, Q) for (P, Q) in prods.values()]
if d == 11:
    named['X0 * |grad F|^2 ... odd: skip'] = None
    named = {k: v for k, v in named.items() if v}
for k, v in named.items():
    fams[k] = v if isinstance(v, list) else [v]

out = []
for name, polys in fams.items():
    polys = [core.reduce_mod(P, p) for P in polys]
    polys = [P for P in polys if P]
    t2 = time.time()
    dimV = span_dim(polys, mons_d) if polys else 0
    nullK, _ = nullity_of(polys, red_hi, hi_idx, True)      # dim (V cap K_full) (V assumed independent? use span dim)
    nullJ, _ = nullity_of(polys, red_lo, lo_idx, False)     # dim (V cap J_d)
    # correct for dependent generators: nullity counts relations; subtract (len - dimV)
    dep = len(polys) - dimV
    rec = {'batch': 'structured-family-E1', 'family': name, 'd': d, 'sector': 'even', 'p': p, 'dim_family': dimV,
           'dim_cap_E1': nullK - dep, 'dim_cap_J': nullJ - dep, 'nonjacobian_E1_dirs': (nullK - dep) - (nullJ - dep),
           'seconds': round(time.time() - t2, 1),
           'note': 'dim(V cap {D: AD in J}) and dim(V cap J_d) mod p; nonjacobian = E1-admissible directions in V not of Jacobian type'}
    print(json.dumps(rec), flush=True)
    log_batch(rec); out.append(rec)
with open(f'{LANE}/w2_struct_e1_d{d}.json', 'w') as f:
    json.dump(out, f, indent=1)
