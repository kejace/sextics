#!/usr/bin/env python3
"""w2_struct.py -- wave 2, item 4 (linear part): intersections of structured subspaces W of S_d^even with the exact
E1-admissible space K_d (mod p, basis from w2_kernel_<d>_0.json), and the non-Jacobian part of the intersection.

E1 for a denominator D is exactly D in K_d, so dim(W cap K_d) = 0 means: no element of the family W passes E1 (mod p).
For each W: dim W, dim(W cap K_d), dim((W cap K_d) + J_d)/J_d, and explicit non-J elements (saved to w2_struct_d<d>_cands.jsonl
for the full fixed-D solve).  Light: all matrices are <= ~2000 x |S_d^even|.
"""
import sys, os, json, time, itertools, argparse, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from flint import nmod_mat
from hodge import core, jacobian
from hodge.core import NV, monomials, mul, add, scale, pderiv, power, reduce_mod, shift

LANE = 'quadratic-cubic-pullback'
ap = argparse.ArgumentParser(); ap.add_argument('--d', type=int, default=10); ap.add_argument('--kernel', default=None)
ap.add_argument('--mirror', default=f'/tmp/hodge/remote/{LANE}'); ap.add_argument('--out', default=None)
a = ap.parse_args()
d = a.d
kpath = a.kernel or f'{a.mirror}/w2_kernel_{d}_0.json'
K = json.load(open(kpath)); p = K['p']
Dm = [tuple(m) for m in K['monomials']]; idx = {m: i for i, m in enumerate(Dm)}; nD = len(Dm)
kern = [{tuple(e): c for e, c in v} for v in K['kernel_basis']]
print(f'[struct d={d}] kernel basis {len(kern)} vectors over {nD} monomials, K/J = {len(K["quotient_basis"])}')
F, A, q, _ = core.load_instance(LANE)
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
red = jacobian._Reducer(F, d, p, 0)

def rows_of(polys):
    flat = []
    for P in polys:
        row = [0] * nD
        for e, c in P.items():
            row[idx[e]] = c % p
        flat.extend(row)
    return nmod_mat(len(polys), nD, flat, p)

def rank(polys):
    return rows_of(polys).rank() if polys else 0

Kmat = rows_of(kern); rK = Kmat.rank()
assert rK == len(kern)

def intersect(name, polys, out_f):
    """polys: exact or mod-p polys spanning W (reduced mod p here)."""
    t0 = time.time()
    Wp = [reduce_mod(P, p) if any(isinstance(c, Fraction) or isinstance(c, int) and c < 0 for c in P.values()) else dict(P) for P in polys]
    Wp = [{e: c % p for e, c in P.items() if c % p} for P in Wp]
    Wp = [P for P in Wp if P]
    if not Wp:
        print(f'  {name}: empty'); return
    rW = rank(Wp)
    stacked = rows_of(Wp + kern); r = stacked.rank()
    dim_int = rW + rK - r
    rec = {'lane': LANE, 'item': 'struct-intersection', 'd': d, 'family': name, 'dimW': rW, 'n_gens': len(Wp),
           'dim_W_cap_K': dim_int, 'p': p}
    nonJ = []
    if dim_int > 0:
        # explicit intersection: nullspace of [W^T | -K^T] (nD x (nW + nK)) -> combinations
        Wm = rows_of(Wp); 
        M = nmod_mat(nD, len(Wp) + len(kern), [0] * (nD * (len(Wp) + len(kern))), p)
        flat = [0] * (nD * (len(Wp) + len(kern)))
        C = len(Wp) + len(kern)
        for j, P in enumerate(Wp):
            for e, c in P.items():
                flat[idx[e] * C + j] = c
        for j, P in enumerate(kern):
            for e, c in P.items():
                flat[idx[e] * C + len(Wp) + j] = (p - c) % p
        M = nmod_mat(nD, C, flat, p)
        X, nul = M.nullspace()
        elems = []
        for k in range(nul):
            coeffs = [int(X[j, k]) for j in range(len(Wp))]
            Dk = {}
            for cj, P in zip(coeffs, Wp):
                if cj:
                    Dk = add(Dk, scale(P, cj, p), p)
            if Dk:
                elems.append(Dk)
        # independent ones
        if elems:
            Em = rows_of(elems); Er, er = Em.rref()
            elems2 = []
            for i in range(er):
                elems2.append({Dm[j]: int(Er[i, j]) for j in range(nD) if int(Er[i, j])})
            elems = elems2
        nf = [red.normal_form(E) for E in elems]
        nf_nz = [v for v in nf if v]
        rec['dim_int_mod_J'] = rank(nf_nz) if nf_nz else 0
        for E, v in zip(elems, nf):
            if v:
                nonJ.append(E)
        rec['n_nonJ_elements'] = len(nonJ)
        for E in nonJ[:50]:
            out_f.write(json.dumps({'family': name, 'd': d, 'p': p,
                                    'D_modp': [[list(e), c] for e, c in jacobian.symmetric(E, p).items()]}) + '\n')
    rec['seconds'] = round(time.time() - t0, 2)
    print('  ' + json.dumps(rec)); return rec

# ---------------------------------------------------------------------------------------------------------------------
# structured families
Qs = [core.parse(s) for s in ['X0**2 + 5*(X1*X2 + Y0*Y1)/7', 'X1**2 + 5*(X0*X2 + Y1*Y2)/11', 'X2**2 + 5*(X0*X1 + Y0*Y2)/13',
                              'Y0**2 + 5*(X0*X1 + Y1*Y2)/17', 'Y1**2 + 5*(X1*X2 + Y0*Y2)/19', 'Y2**2 + 5*(X0*X2 + Y0*Y1)/23']]
chk = {}
for Qi in Qs: chk = add(chk, power(Qi, 3))
assert add(F, scale(chk, -1)) == {}, 'Q presentation mismatch'
delta_src = open('/tmp/hodge/results/quadratic-cubic-pullback/delta.txt').read().strip()
delta = core.parse(delta_src)
assert core.degree(delta) == 6 and core.parity_of(delta) == 1
Fi = [pderiv(F, i) for i in range(NV)]

def Qmonos(deg, S=range(6)):
    out = []
    for alpha in itertools.combinations_with_replacement(list(S), deg):
        P = {(0,) * 6: Fraction(1)}
        for i in alpha:
            P = mul(P, Qs[i])
        out.append(P)
    return out

fams = []
if d % 2 == 0:
    fams.append((f'C[Q]_{d//2} (pullbacks of degree-{d//2} hypersurface sections of the cubic)', Qmonos(d // 2)))
    for S in itertools.combinations(range(6), 3):
        fams.append((f'C[Q_S]_{d//2} S={S}', Qmonos(d // 2, S)))
    fams.append((f'C[Q_0,Q_1,Q_2]_a * C[Q_3,Q_4,Q_5]_b bihomogeneous products (a+b={d//2})',
                 [mul(P1, P2) for k in range(1, d // 2) for P1 in Qmonos(k, [0, 1, 2]) for P2 in Qmonos(d // 2 - k, [3, 4, 5])]))
# delta * S_{d-6}^odd
if d >= 7:
    fams.append((f'delta * S_{d-6}^odd', [mul(delta, {m: Fraction(1)}) for m in monomials(d - 6, 1)]))
if d >= 12:
    fams.append(('delta^2 * S_{d-12}^even', [mul(mul(delta, delta), {m: Fraction(1)}) for m in monomials(d - 12, 0)]))
# cones with vertex a plane: polynomials in X only / in Y only (any Y-even degree)
fams.append(('C[X0,X1,X2]_d (cones, vertex the Y-plane)', [{m: Fraction(1)} for m in monomials(d, 0) if m[3] == m[4] == m[5] == 0]))
if d % 2 == 0:
    fams.append(('C[Y0,Y1,Y2]_d (cones, vertex the X-plane)', [{m: Fraction(1)} for m in monomials(d, 0) if m[0] == m[1] == m[2] == 0]))
# cones with vertex a line: polynomials in 4 coordinates (all 15 choices, Y-even monomials)
for V in itertools.combinations(range(6), 4):
    comp = [i for i in range(6) if i not in V]
    fams.append((f'C[z_{V}]_d (cone, vertex line z_{comp[0]}=z_{comp[1]}=0)',
                 [{m: Fraction(1)} for m in monomials(d, 0) if all(m[i] == 0 for i in comp)]))
# bigraded pieces: X-degree a, Y-degree b (Y-even) -- cones over curves in the join sense (products of X-forms and Y-forms)
for b in range(0, d + 1, 2):
    fams.append((f'bidegree ({d-b},{b}) monomials', [{m: Fraction(1)} for m in monomials(d, 0) if m[3] + m[4] + m[5] == b]))
# E_I(F) * S_{d-6}  (block radial derivatives times anything: reducible, but tests the polar + multiplier shape; in J)
# polars times quadrics etc. are in J -> skip.  Products of two Jacobian elements (Y-even): F_i F_j * S_{d-10}
if d >= 10:
    fams.append(('F_i F_j * S_{d-10} (products of two polars)',
                 [mul(mul(Fi[i], Fi[j]), {m: Fraction(1)}) for i in range(6) for j in range(i, 6)
                  for m in monomials(d - 10, (0 if (i < 3) == (j < 3) else 1))]))
# A-multiples and A-related: A * S_{d-6}^odd (D divisible by the target itself)
fams.append((f'A * S_{d-6}^odd', [mul(A, {m: Fraction(1)}) for m in monomials(d - 6, 1)]))
# Hessian of the cubic-in-Q structure: det of Hess_Z (sum Q_i^3) is degree 24; skip.  Discriminant-type: skip.
# --- union families and shape-symmetry invariant subspaces (dossier E3)
def union(names):
    out = []
    for nm, polys in fams:
        if any(nm.startswith(s) for s in names):
            out.extend(polys)
    return out
if d % 2 == 0:
    fams.append(('UNION C[Q]_k + delta*S + A*S + F_iF_j*S', union([f'C[Q]_{d//2}', 'delta * S', 'A * S', 'F_i F_j'])))
    fams.append(('UNION C[Q]_k + bihomogeneous + all C[Q_S]_k', union([f'C[Q]_{d//2}', 'C[Q_0,Q_1,Q_2]', 'C[Q_S]'])))
fams.append(('UNION all cones with vertex a line', union(['C[z_'])))
fams.append(('UNION cones vertex a plane + extreme bidegrees', union(['C[X0,X1,X2]', 'C[Y0,Y1,Y2]', f'bidegree ({d-2},2)', f'bidegree (2,{d-2})'])))
def orbit_sums(perms):
    seen = set(); out = []
    for m in monomials(d, 0):
        if m in seen: continue
        orb = {m}; frontier = [m]
        while frontier:
            x = frontier.pop()
            for pm in perms:
                y = tuple(x[pm[i]] for i in range(6))
                if y not in orb: orb.add(y); frontier.append(y)
        seen |= orb
        P = {}
        for y in orb: P[y] = Fraction(1)
        out.append(P)
    return out
rho = (1, 2, 0, 4, 5, 3)     # X0->X1->X2->X0, Y0->Y1->Y2->Y0  (as index permutation of exponents)
sigma = (3, 4, 5, 0, 1, 2)   # X_i <-> Y_i
fams.append(('rho-invariant (cyclic shift of both blocks) Y-even polynomials', orbit_sums([rho])))
fams.append(('sigma-invariant (X_i <-> Y_i) Y-even polynomials', orbit_sums([sigma])))
fams.append(('(rho,sigma)-invariant Y-even polynomials', orbit_sums([rho, sigma])))
tau_x = (1, 0, 2, 4, 3, 5)   # X0<->X1, Y0<->Y1 (a transposition of the shape symmetry S3)
fams.append(('S3 x (X<->Y) shape-symmetric Y-even polynomials', orbit_sums([rho, sigma, tau_x])))
out_path = a.out or f'/tmp/hodge/results/{LANE}/w2_struct_d{d}_cands.jsonl'
recs = []
with open(out_path, 'a') as out_f:
    for name, polys in fams:
        r = intersect(name, polys, out_f)
        if r: recs.append(r)
with open(f'/tmp/hodge/results/{LANE}/w2_struct.jsonl', 'a') as f:
    for r in recs:
        f.write(json.dumps(r) + '\n')
print('done', len(recs), 'families')
