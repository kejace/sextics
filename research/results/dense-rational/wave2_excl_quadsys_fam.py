#!/usr/bin/env python3
"""Hessian-relaxation quadratic system restricted to a structured linear FAMILY of polar fields w (D = w(F)), for
d in {7, 8}, as an msolve input.  Necessary condition for a certificate with D = w(F) (E1 forces D in J_d for d <= 9):
   d = 7:  w^T Hess(F) w in J_8          ([A] injective on R_8^{+/-})
   d = 8:  w^T Hess(F) w in (J:A)_10^+   (quotient dim 762 - 294 = 468; (J:A)_10^+ = K_full(10) from the kernel sweep)
Families (W_fam a linear subspace of the sector's w-space W; a legal D = w(F) with w in W_fam is excluded when the
restricted homogeneous system has only the trivial solution mod p — rigorous over C by the projective reduction argument,
after removing W_fam ∩ gauge):
  full       all coefficients (reference; for d = 8 only for LM)
  blockdiag  w_Xi in C[X0,X1,X2]_e,  w_Yj in C[Y0,Y1,Y2]_e           (pure even at d=8; mixed at d=7)
  cross      w_Xi in C[Y]_e,         w_Yj in C[X]_e                   (pure odd at d=8; mixed at d=7)
  xonly      w_Y = 0  /  yonly  w_X = 0
  weuler     w_i = z_i g_i  (weighted Euler fields; gauge g_i = common l)
  gradient   w = grad g, g of degree e+1
  pair-i-j / comp-i-j-k-...   only the listed components nonzero (all monomials)
Usage: python3 wave2_excl_quadsys_fam.py d sector family prime   (sector even|odd|mixed)
Writes wave2-excl/fam_d{d}_{sector}_{family}_p{prime}.ms + _receipt.json
"""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
t0 = time.time()
d = int(sys.argv[1]); sector = sys.argv[2]; fam = sys.argv[3]; p = int(sys.argv[4])
assert d in (7, 8) and sector in ('even', 'odd', 'mixed')
if d == 8: assert sector != 'mixed', 'mixed d=8 would need (J:A)_10^- as well'
e = d - 5
F, A, q = load()
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
pd = {'even': 0, 'odd': 1, 'mixed': None}[sector]
basis = []
for i in range(NV):
    pars = [0, 1] if pd is None else [(par_i(i) + pd) % 2]
    for par in pars:
        for m in monomials(e, par): basis.append((i, m))
nvar = len(basis); bidx = {b: k for k, b in enumerate(basis)}
# gauge vectors l*Euler in full coordinates
gpars = [0, 1] if pd is None else [pd]
gauge = []
for gp in gpars:
    for l in monomials(d - 6, gp):
        v = [0] * nvar
        for i in range(NV):
            ei = [0] * NV; ei[i] = 1
            v[bidx[(i, tuple(a + b for a, b in zip(l, ei)))]] = 1
        gauge.append(v)
# family basis vectors (lists of length nvar)
def unit(i, m):
    v = [0] * nvar; v[bidx[(i, m)]] = 1; return v
Xsupp = lambda m: m[3] + m[4] + m[5] == 0
Ysupp = lambda m: m[0] + m[1] + m[2] == 0
famvecs = []
if fam == 'full': famvecs = [unit(i, m) for (i, m) in basis]
elif fam == 'blockdiag': famvecs = [unit(i, m) for (i, m) in basis if (Xsupp(m) if i < 3 else Ysupp(m))]
elif fam == 'cross': famvecs = [unit(i, m) for (i, m) in basis if (Ysupp(m) if i < 3 else Xsupp(m))]
elif fam == 'xonly': famvecs = [unit(i, m) for (i, m) in basis if i < 3]
elif fam == 'yonly': famvecs = [unit(i, m) for (i, m) in basis if i >= 3]
elif fam == 'weuler':
    for i in range(NV):
        ei = [0] * NV; ei[i] = 1
        for gp in ([0, 1] if pd is None else [pd]):
            for g in monomials(e - 1, gp):
                m = tuple(a + b for a, b in zip(g, ei))
                if (i, m) in bidx: famvecs.append(unit(i, m))
elif fam == 'gradient':
    for gp in ([0, 1] if pd is None else [pd]):
        for g in monomials(e + 1, gp):
            v = [0] * nvar
            for i in range(NV):
                if g[i] == 0: continue
                m = list(g); m[i] -= 1; m = tuple(m)
                v[bidx[(i, m)]] = (v[bidx[(i, m)]] + g[i]) % p
            famvecs.append(v)
elif fam.startswith('pair-') or fam.startswith('comp-'):
    comps = set(int(x) for x in fam.split('-')[1:])
    famvecs = [unit(i, m) for (i, m) in basis if i in comps]
else: raise SystemExit('unknown family')
nfam = len(famvecs)
Vf = nmod_mat(nvar, nfam, [famvecs[c][r] for r in range(nvar) for c in range(nfam)], p)
assert Vf.rank() == nfam
# gauge inside the family: u with Vf u = g  for gauge vectors g (solve exactly mod p)
gin = []
for g in gauge:
    M = nmod_mat(nvar, nfam + 1, [famvecs[c][r] if c < nfam else g[r] for r in range(nvar) for c in range(nfam + 1)], p)
    Mr, rk = M.rref()
    if rk == nfam:        # consistent (no pivot in last column) iff rank unchanged
        u = [0] * nfam
        for r in range(rk):
            j = 0
            while int(Mr[r, j]) == 0: j += 1
            u[j] = int(Mr[r, nfam])
        gin.append(u)
    elif rk == nfam + 1: pass
    else: raise RuntimeError
# some gauge vectors may be inside the family only as combinations (e.g. gradient family); handle via subspace intersection
if gauge:
    # intersection of span(gauge) with span(fam): kernel of [Vf | -G]
    G = nmod_mat(nvar, len(gauge), [gauge[c][r] for r in range(nvar) for c in range(len(gauge))], p)
    Mk = nmod_mat(nvar, nfam + len(gauge), [int(Vf[r, c]) if c < nfam else (p - int(G[r, c - nfam])) % p for r in range(nvar) for c in range(nfam + len(gauge))], p)
    K, nk = Mk.nullspace()
    gin = []
    for k in range(nk):
        gin.append([int(K[j, k]) for j in range(nfam)])
Gin = nmod_mat(max(len(gin), 1), nfam, [x for u in gin for x in u] if gin else [0] * nfam, p)
Ginr, grk = Gin.rref()
gpiv = []
for r in range(grk):
    j = 0
    while int(Ginr[r, j]) == 0: j += 1
    gpiv.append(j)
keep = [c for c in range(nfam) if c not in set(gpiv)]
# quotient-defining space Z in the target degree (J_{2d-6} for d=7; (J:A)_10^+ for d=8)
deg = 2 * d - 6
tpars = [0] if pd is not None else [0, 1]
mons = []
for tp in tpars: mons += monomials(deg, tp)
idx = {m: j for j, m in enumerate(mons)}; nm = len(mons)
if d == 7:
    zrows = []
    for i in range(NV):
        for tp in tpars:
            for m in monomials(deg - 5, (tp - par_i(i)) % 2): zrows.append(shift(Fi[i], m, 1, p))
    expected = len(zrows)
else:
    bf = f'{OUT}/e1kernel_d10_par0_p{p}_basis.json'
    data = json.load(open(bf)); assert data['p'] == p and data['d'] == 10 and data['parD'] == 0
    zrows = [{tuple(int(x) for x in k.split()): v for k, v in b.items()} for b in data['basis']]
    expected = 1044
Z = nmod_mat(len(zrows), nm, p)
for r, g in enumerate(zrows):
    for mm, v in g.items(): Z[r, idx[mm]] = v
Zr, rk = Z.rref(); assert rk == expected, (rk, expected)
pivots = []
for r in range(rk):
    j = 0
    while int(Zr[r, j]) == 0: j += 1
    pivots.append(j)
pivset = set(pivots); nonpiv = [j for j in range(nm) if j not in pivset]
neq = len(nonpiv)
print(f'd={d} {sector} fam={fam}: nvar={nvar} nfam={nfam} gauge_in_fam={grk} unknowns={len(keep)} quotient dim={neq} ({time.time()-t0:.1f}s)', flush=True)
# normal forms of all pair products restricted later; build V (pairs x mons) only for pairs touching the family support
supp = sorted({r for v in famvecs for r in range(nvar) if v[r]})
pairs = []; rowsV = []
for a in range(len(supp)):
    k = supp[a]; ik, mk = basis[k]
    for b in range(a, len(supp)):
        l = supp[b]; il, ml = basis[l]
        Pm = shift(H[ik][il], tuple(x + y for x, y in zip(mk, ml)), 1, p)
        if Pm: pairs.append((k, l)); rowsV.append(Pm)
V = nmod_mat(len(pairs), nm, p)
for r, Pm in enumerate(rowsV):
    for mm, c in Pm.items(): V[r, idx[mm]] = c % p
VP = nmod_mat(len(pairs), rk, [int(V[r, j]) for r in range(len(pairs)) for j in pivots], p)
NF = V - VP * Zr
# for each equation r: full quadratic form Q_r on supp coordinates, restrict with Vf (rows supp), gauge-fix
sidx = {k: t for t, k in enumerate(supp)}
Vs = nmod_mat(len(supp), nfam, [famvecs[c][k] for k in supp for c in range(nfam)], p)
VsT = Vs.transpose()
eqs = []
for rr, j in enumerate(nonpiv):
    col = [int(NF[r, j]) for r in range(len(pairs))]
    if not any(col): continue
    Qr = nmod_mat(len(supp), len(supp), p)
    for r, (k, l) in enumerate(pairs):
        c = col[r]
        if c:
            Qr[sidx[k], sidx[l]] = (int(Qr[sidx[k], sidx[l]]) + c) % p
            if k != l: Qr[sidx[l], sidx[k]] = (int(Qr[sidx[l], sidx[k]]) + c) % p
    Rr = VsT * Qr * Vs
    terms = []
    for a in range(nfam):
        if a not in set(keep): continue
        for b in range(a, nfam):
            if b not in set(keep): continue
            c = int(Rr[a, b]) if a == b else (int(Rr[a, b]) + int(Rr[b, a])) % p
            if c: terms.append(f'{c}*c{a}*c{b}' if a != b else f'{c}*c{a}^2')
    if terms: eqs.append('+'.join(terms))
varnames = [f'c{a}' for a in keep]
fn = f'{OUT}/fam_d{d}_{sector}_{fam}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n'); f.write(f'{p}\n'); f.write(',\n'.join(eqs) + '\n')
receipt = {'exp': 'hessian-relaxation-family-system', 'd': d, 'sector': sector, 'family': fam, 'prime': p, 'nvar_full': nvar,
           'family_dim': nfam, 'gauge_in_family': grk, 'unknowns': len(varnames), 'quotient': ('J_8' if d == 7 else '(J:A)_10^+'),
           'quotient_dim': neq, 'equations': len(eqs), 'file': fn, 'seconds': round(time.time() - t0, 1)}
with open(fn.replace('.ms', '_receipt.json'), 'w') as f: json.dump(receipt, f)
log2(**receipt, note='msolve input written')
print(json.dumps(receipt), flush=True)
