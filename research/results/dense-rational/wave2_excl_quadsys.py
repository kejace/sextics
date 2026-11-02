#!/usr/bin/env python3
"""Hessian-relaxation (necessary) quadratic system for Jacobian-type denominators D = w(F) = sum_i w_i F_i of
degree d in {6,7} (e = d-5 = deg w_i), written as an msolve input file (homogeneous, no slice).

Sectors: even  (D Y-even: w_i has Y-parity p_i, p_i = 0 for X_i, 1 for Y_i)
         odd   (D Y-odd : w_i has Y-parity p_i + 1)
         mixed (no parity restriction: all 6*dim S_e coefficients; D = D_even + D_odd both nonzero allowed)
Condition: w^T Hess(F) w in J_{2d-6}  (both Y-parities of S_{2d-6} in the mixed sector; in a pure sector w^T Hess w is Y-even).
Why necessary (q = 1, A = X0^3 Y0^3): E1 forces D in J_d for d <= 9 (E1-kernel sweep, both parities), so D = w(F);
every E1 numerator is N = A w + K grad F (K skew: Koszul syzygies); E2 gives B_D(N) = -A w^T Hess(F) w + (element of J) = F C,
so A w^T Hess w in J_{2d}, and injectivity of [A] on R_{2d-6}^{+/-} (2d-6 <= 8) gives w^T Hess w in J_{2d-6}.
Gauge: w -> w + l*Euler (l in S_{d-6}, parity of D) preserves the condition (Euler: Hess(F) z = 5 grad F); removed by
zeroing the pivot coordinates of the gauge subspace. A legal D is never a gauge element (it would be l*F, gcd(D,F) != 1),
and w(F) = 0 forces w = 0 (no syzygies of grad F below degree 5).
Rigor: the equations are p-integral reductions of the rational ones (the pivot minor of J_{2d-6} is a p-unit because the
modular rank attains the exact Koszul count 6*dim S_{2d-11}); a projective variety over Q with a point over C has a point
over F_p-bar; hence "no nonzero solution mod p" (leading ideal contains a pure power of every variable) proves
"no nonzero solution over C", i.e. no certificate with a degree-d denominator in this sector.
Usage: python3 wave2_excl_quadsys.py d sector prime
Writes wave2-excl/quadsys_d{d}_{sector}_p{prime}.ms and ..._receipt.json
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
t0 = time.time()
d = int(sys.argv[1]); sector = sys.argv[2]; p = int(sys.argv[3])
assert d in (6, 7) and sector in ('even', 'odd', 'mixed')
e = d - 5
F, A, q = load()
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
basis = []
for i in range(NV):
    pars = {'even': [par_i(i)], 'odd': [(par_i(i) + 1) % 2], 'mixed': [0, 1]}[sector]
    for par in pars:
        for m in monomials(e, par):
            basis.append((i, m))
nvar = len(basis); bidx = {b: k for k, b in enumerate(basis)}
gpars = {'even': [0], 'odd': [1], 'mixed': [0, 1]}[sector]
gauge = []
for gp in gpars:
    for l in monomials(d - 6, gp):
        v = [0] * nvar
        for i in range(NV):
            ei = [0] * NV; ei[i] = 1
            m = tuple(a + b for a, b in zip(l, ei))
            v[bidx[(i, m)]] = 1
        gauge.append(v)
ng = len(gauge)
deg = 2 * d - 6
tpars = [0] if sector != 'mixed' else [0, 1]
mons = []
for tp in tpars: mons += monomials(deg, tp)
idx = {m: j for j, m in enumerate(mons)}; nm = len(mons)
gens = []
for i in range(NV):
    for tp in tpars:
        for m in monomials(deg - 5, (tp - par_i(i)) % 2):
            gens.append(shift(Fi[i], m, 1, p))
G = nmod_mat(len(gens), nm, p)
for r, g in enumerate(gens):
    for mm, v in g.items(): G[r, idx[mm]] = v
Gr, rk = G.rref()
expected = len(gens)          # no syzygies of grad F below degree 10 (deg - 5 <= 3): dim J_deg = 6 dim S_{deg-5} (per parity)
assert rk == expected, (rk, expected)
pivots = []
for r in range(rk):
    j = 0
    while int(Gr[r, j]) == 0: j += 1
    pivots.append(j)
pivset = set(pivots); nonpiv = [j for j in range(nm) if j not in pivset]
print(f'd={d} sector={sector} p={p}: nvar={nvar} gauge={ng} dim J_{deg} = {rk} (expected {expected}) dim quotient = {len(nonpiv)}', flush=True)
pairs = []; rowsV = []
for k in range(nvar):
    ik, mk = basis[k]
    for l in range(k, nvar):
        il, ml = basis[l]
        Pm = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1, p)
        if Pm: pairs.append((k, l)); rowsV.append(Pm)
npairs = len(pairs)
V = nmod_mat(npairs, nm, p)
for r, Pm in enumerate(rowsV):
    for mm, c in Pm.items(): V[r, idx[mm]] = c % p
VP = nmod_mat(npairs, rk, [int(V[r, j]) for r in range(npairs) for j in pivots], p)
NF = V - VP * Gr
Qcoef = {}
for r, (k, l) in enumerate(pairs):
    nf = [int(NF[r, j]) for j in nonpiv]
    if any(nf): Qcoef[(k, l)] = nf
Gm = nmod_mat(ng, nvar, [x for v in gauge for x in v], p) if ng else None
gpiv = []
if ng:
    Gmr, grk = Gm.rref(); assert grk == ng
    for r in range(grk):
        j = 0
        while int(Gmr[r, j]) == 0: j += 1
        gpiv.append(j)
keep = [k for k in range(nvar) if k not in set(gpiv)]
kpos = {k: t for t, k in enumerate(keep)}
varnames = [f'c{k}' for k in keep]
eqs = []
for r in range(len(nonpiv)):
    terms = []
    for (k, l), nf in Qcoef.items():
        if k not in kpos or l not in kpos: continue
        c = nf[r]
        if k != l: c = (2 * c) % p
        if c: terms.append(f'{c}*c{k}*c{l}' if k != l else f'{c}*c{k}^2')
    if terms: eqs.append('+'.join(terms))
fn = f'{OUT}/quadsys_d{d}_{sector}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n'); f.write(f'{p}\n'); f.write(',\n'.join(eqs) + '\n')
receipt = {'exp': 'hessian-relaxation-system', 'd': d, 'sector': sector, 'prime': p, 'nvar_w': nvar, 'gauge_dim': ng,
           'unknowns': len(varnames), 'target_degree': deg, 'target_parities': tpars, 'dim_J_target': rk,
           'dim_J_expected': expected, 'dim_quotient': len(nonpiv), 'equations': len(eqs), 'file': fn,
           'basis': [[i, list(m)] for i, m in basis], 'kept_vars': keep, 'gauge_pivots': gpiv,
           'seconds': round(time.time() - t0, 1)}
with open(fn.replace('.ms', '_receipt.json'), 'w') as f: json.dump(receipt, f)
log2(**{k: v for k, v in receipt.items() if k not in ('basis',)}, note='msolve input written; run msolve -g 1 and wave2_excl_ltcheck.py')
print(f'wrote {fn}: {len(varnames)} vars, {len(eqs)} eqs ({time.time()-t0:.1f}s)', flush=True)
