#!/usr/bin/env python3
"""Y-degree strata of the degree-d Hessian relaxation (VM excl agent).

Full system (sector even/odd, d = 7): w in (S_2)^6 with parities, condition Q(w) := w^T Hess(F) w in J_8^+ (483 quadrics, 60 unknowns).
J is not bigraded, but two *projections* of the condition are valid necessary conditions:
  restriction Y0 (set Y = 0):  Q(w)|_{Y=0} in J_8|_{Y=0} = J(a)_8,  a = F|_{Y=0} (plane sextic), only Y-free w-coefficients survive;
  restriction X0 (set X = 0):  Q(w)|_{X=0} in J(dY)_8,  dY = F|_{X=0}.
(Because (m F_i)|_{Y=0} = m|_{Y=0} a_i and F_{Y_j}|_{Y=0} = 0.)  Each gives 15 quadrics (dim R(a)_8 = 15 for a smooth plane sextic).
Even sector: Y0 -> unknowns w_X^{(0)} (18, gauge l*X with l in <X>: 3) ; X0 -> unknowns w_X^{(2)} (18, no gauge).
Odd sector:  Y0 -> unknowns w_Y^{(0)} (18, no gauge)                ; X0 -> unknowns w_Y^{(2)} (18, gauge l*Y: 3).
Also 'none' = full system (same as wave2_excl_quadsys.py, re-derived here independently) for cross-checks.
Usage: python3 w2x_strata.py d sector {Y0|X0|none} prime [--run]   -> wave2-excl/strata_d{d}_{sector}_{restr}_p{prime}.{ms,sing}
"""
import sys, time, json, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2x_log import *
t0 = time.time()
d = int(sys.argv[1]); sector = sys.argv[2]; restr = sys.argv[3]; p = int(sys.argv[4]); run = '--run' in sys.argv
e = d - 5
F, A, q = load()
Fp = core.reduce_mod(F, p)
def restrict(P):
    if restr == 'Y0': return {m: v for m, v in P.items() if m[3] + m[4] + m[5] == 0}
    if restr == 'X0': return {m: v for m, v in P.items() if m[0] + m[1] + m[2] == 0}
    return dict(P)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[restrict(pderiv(Fi[i], j, p)) for j in range(NV)] for i in range(NV)]
Fir = [restrict(Fi[i]) for i in range(NV)]
pd = 0 if sector == 'even' else 1
def keep_mon(m):
    if restr == 'Y0': return m[3] + m[4] + m[5] == 0
    if restr == 'X0': return m[0] + m[1] + m[2] == 0
    return True
basis = [(i, m) for i in range(NV) for m in monomials(e, (pd + par_i(i)) % 2) if keep_mon(m)]
nvar = len(basis); bidx = {b: k for k, b in enumerate(basis)}
# gauge vectors (restricted)
gauge = []
for l in monomials(d - 6, pd):
    v = [0] * nvar; nz = False
    for i in range(NV):
        m = tuple(l[t] + (1 if t == i else 0) for t in range(NV))
        if (i, m) in bidx: v[bidx[(i, m)]] = 1; nz = True
    if nz: gauge.append(v)
deg = 2 * d - 6
mons = [m for m in monomials(deg, 0) if keep_mon(m)]
idx = {m: j for j, m in enumerate(mons)}; nm = len(mons)
gens = []
for i in range(NV):
    if not Fir[i]: continue
    for m in monomials(deg - 5, par_i(i)):
        if keep_mon(m): gens.append(shift(Fir[i], m, 1, p))
G = nmod_mat(len(gens), nm, p)
for r, g in enumerate(gens):
    for mm, v in g.items(): G[r, idx[mm]] = v
Gr, rk = G.rref()
pivots = []
for r in range(rk):
    j = 0
    while int(Gr[r, j]) == 0: j += 1
    pivots.append(j)
pivset = set(pivots); nonpiv = [j for j in range(nm) if j not in pivset]
print(f'd={d} sector={sector} restr={restr} p={p}: nvar={nvar} gauge={len(gauge)} target monomials {nm}, J-gens {len(gens)} rank {rk}, quotient {len(nonpiv)}', flush=True)
pairs = []; rowsV = []
for k in range(nvar):
    ik, mk = basis[k]
    for l in range(k, nvar):
        il, ml = basis[l]
        Pm = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1, p)
        Pm = {m: v for m, v in Pm.items() if m in idx}
        if Pm: pairs.append((k, l)); rowsV.append(Pm)
V = nmod_mat(len(pairs), nm, p)
for r, Pm in enumerate(rowsV):
    for mm, c in Pm.items(): V[r, idx[mm]] = c % p
VP = nmod_mat(len(pairs), rk, [int(V[r, j]) for r in range(len(pairs)) for j in pivots], p)
NF = V - VP * Gr
Qcoef = {}
for r, (k, l) in enumerate(pairs):
    nf = [int(NF[r, j]) for j in nonpiv]
    if any(nf): Qcoef[(k, l)] = nf
# gauge check: each equation must be invariant: Q(w + g) - Q(w) = 2 B(w,g) + Q(g) must vanish identically mod J => check B(e_k, g) = 0 for all k and Q(g)=0
viol = 0
for g in gauge:
    for r in range(len(nonpiv)):
        # linear form in w: sum_{k,l} Qcoef[(k,l)][r] * (e_k g_l + e_l g_k)
        lin = [0] * nvar
        for (k, l), nf in Qcoef.items():
            c = nf[r]
            if not c: continue
            if k == l: lin[k] = (lin[k] + 2 * c * g[k]) % p
            else:
                lin[k] = (lin[k] + c * g[l]) % p; lin[l] = (lin[l] + c * g[k]) % p
        viol += sum(1 for x in lin if x)
gpiv = []
if gauge:
    Gm = nmod_mat(len(gauge), nvar, [x for v in gauge for x in v], p)
    Gmr, grk = Gm.rref(); assert grk == len(gauge)
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
tag = f'strata_d{d}_{sector}_{restr}_p{p}'
fn = f'{OUT}/{tag}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n'); f.write(f'{p}\n'); f.write(',\n'.join(eqs) + '\n')
sing = f'{OUT}/{tag}.sing'
with open(sing, 'w') as f:
    f.write(f'ring R = {p}, ({",".join(varnames)}), dp;\noption(redSB);\nideal I;\n')
    for eq in eqs: f.write('I[size(I)+1] = ' + eq + ';\n')
    f.write('int t = timer;\nideal G = std(I);\nprint("dim = " + string(dim(G)));\nprint("mult = " + string(mult(G)));\n')
    f.write('print("time = " + string(timer - t));\nintvec hv = hilb(G, 1);\nprint("hilb1 = " + string(hv));\nintvec hv2 = hilb(G, 2);\nprint("hilb2 = " + string(hv2));\n')
    f.write('if (dim(G) > 1) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); int ii; for (ii = 1; ii <= size(L); ii++) { print("component " + string(ii) + ": dim " + string(dim(std(L[ii]))) + ", mult " + string(mult(std(L[ii]))) + ", ngens " + string(size(L[ii]))); } L; }\nquit;\n')
receipt = {'exp': 'hessian-relaxation-strata', 'd': d, 'sector': sector, 'restriction': restr, 'prime': p, 'nvar_w': nvar,
           'gauge_dim': len(gauge), 'gauge_invariance_violations': viol, 'unknowns': len(varnames), 'target_degree': deg,
           'target_monomials': nm, 'J_gens': len(gens), 'rank_J': rk, 'dim_quotient': len(nonpiv), 'equations': len(eqs),
           'file': fn, 'basis': [[i, list(m)] for i, m in basis], 'kept_vars': keep, 'gauge_pivots': gpiv, 'seconds': round(time.time() - t0, 1)}
with open(fn.replace('.ms', '_receipt.json'), 'w') as f: json.dump(receipt, f)
print(f'wrote {fn}: {len(varnames)} vars, {len(eqs)} eqs, gauge violations {viol} ({time.time()-t0:.1f}s)', flush=True)
rec = {k: v for k, v in receipt.items() if k != 'basis'}
if run:
    t1 = time.time()
    r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=36000)
    out = r.stdout + r.stderr
    with open(f'{OUT}/{tag}.singular.out', 'w') as f: f.write(out)
    print(out[-4000:], flush=True)
    for key in ('dim', 'mult', 'hilb1', 'hilb2', 'ncomp'):
        mt = re.search(rf'^{key} = (.*)$', out, re.M)
        if mt: rec['singular_' + key] = mt.group(1).strip()
    rec['singular_seconds'] = round(time.time() - t1, 1)
    rec['note'] = 'Singular dim is affine dim of the cone: dim 0 <=> projectively empty (then rigorous over C by the Hilbert-function semicontinuity argument); dim 1 <=> finitely many projective points (mult = degree)'
log2(**rec)
