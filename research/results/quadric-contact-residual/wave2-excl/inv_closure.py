"""Gamma-invariant family (invariant under G_aff and the exchanges s0: X0<->Y0, r1: X1<->Y1) at even degree d:
basis = orbit sums of chi=0 monomials; E1 kernel; gauge (s*F, s invariant of degree d-6) first; Jacobian part (symmetrized
m*F_i); extras (lifted).  Hessian relaxation + msolve -g 1 verdict.  Usage: inv_closure.py d [threads] [timeout]"""
import sys, time, json, os
from w2lib import *
from family import e1_profile_sector
from msverdict import verdict
d = int(sys.argv[1]); threads = int(sys.argv[2]) if len(sys.argv) > 2 else 4; tmo = int(sys.argv[3]) if len(sys.argv) > 3 else 1800
OUT = PK + '/inv'
os.makedirs(OUT, exist_ok=True)
ID = (0, 1, 2, 3, 4, 5); S0P = (3, 1, 2, 0, 4, 5); R1P = (0, 4, 2, 3, 1, 5); S0R1 = (3, 4, 2, 0, 1, 5)
GRP = [ID, S0P, R1P, S0R1]
INV4 = pow(4, p - 2, p); INV6 = pow(6, p - 2, p)
def perm_poly(P, g):
    return {tuple(e[g[j]] for j in range(NV)): c for e, c in P.items()}
def sym_poly(P):
    out = {}
    for g in GRP: out = add(out, perm_poly(P, g), p)
    return scale(out, INV4, p)
def perm_field(N, g):
    return [perm_poly(N[g[i]], g) for i in range(NV)]
def sym_field(N):
    out = [{} for _ in range(NV)]
    for g in GRP:
        gN = perm_field(N, g); out = [add(out[i], gN[i], p) for i in range(NV)]
    return [scale(c, INV4, p) for c in out]
assert sym_poly(Fp) == Fp
def inv_basis(deg):
    seen, basis = set(), []
    for e in monomials(deg):
        if chi(e) != CHI0 or e in seen: continue
        orb = sorted({tuple(e[g[j]] for j in range(NV)) for g in GRP}); seen.update(orb)
        basis.append({m: 1 for m in orb})
    return basis
t0 = time.time()
basis = inv_basis(d); nb = len(basis)
bidx = {}
for k, B in enumerate(basis):
    for e in B: bidx[e] = k
def coords(P):
    v = [0] * nb
    for e, c in P.items(): v[bidx[e]] = (v[bidx[e]] + c) % p
    return v
S6 = sector(d + 6, CHI_A)
NF = S6.nf_matrix([S6.vec(mul(Ap, B, p)) for B in basis])
X, nul = NF.transpose().nullspace()
Ds, Ns, vecs = [], [], []
def try_add(D, N):
    v = coords(D)
    test = nmod_mat(len(vecs) + 1, nb, [x for u in vecs + [v] for x in u], p)
    if test.rank() == len(vecs) + 1:
        vecs.append(v); Ds.append(D); Ns.append(N); return True
    return False
ng = 0
if d >= 6:
    for s in inv_basis(d - 6):
        N = [mul(Ap, scale(mul(s, {tuple(1 if j == i else 0 for j in range(NV)): 1}, p), INV6, p), p) for i in range(NV)]
        if try_add(mul(s, Fp, p), N): ng += 1
nj = 0
for mm in monomials(d - 5):
    for i in range(NV):
        if chi_add(chi(mm), chi_neg(CHI_Z[i])) != CHI0: continue
        C = [{} for _ in range(NV)]; C[i] = {mm: 1}
        Cs = sym_field(C)
        D = {}
        for j in range(NV): D = add(D, mul(Cs[j], Fi[j], p), p)
        if D and try_add(D, [mul(Ap, c, p) for c in Cs]): nj += 1
ne = 0
# extras: kernel vectors not in the span so far
for j in range(nul):
    v = [int(X[i, j]) for i in range(nb)]
    P = {}
    for k, cv in enumerate(v):
        if cv: P = add(P, scale(basis[k], cv, p), p)
    test = nmod_mat(len(vecs) + 1, nb, [x for u in vecs + [v] for x in u], p)
    if test.rank() == len(vecs) + 1:
        N = S6.lift_field(mul(Ap, P, p)); N = sym_field(N)
        assert sub(vfield(N, Fp), mul(Ap, P, p), p) == {}
        vecs.append(v); Ds.append(P); Ns.append(N); ne += 1
n = len(Ds)
info = {'d': d, 'dim_inv_Sd': nb, 'dim_E1_kernel': nul, 'n_gauge': ng, 'n_jac': nj, 'n_extra': ne, 'n_params': n}
print(info, file=sys.stderr, flush=True)
quads, qi = quadrics_from_family(Ds, Ns, d, verbose=False)
keep = list(range(ng, n)); m = n - ng
quads2 = []
for cf in quads:
    cf2 = {(k - ng, l - ng): v for (k, l), v in cf.items() if k >= ng and l >= ng}
    if cf2: quads2.append(cf2)
rank, forced, pairs = quad_rank_and_squares(quads2, m)
prefix = OUT + f'/inv_d{d}'
neq, nv = write_systems(quads2, m, prefix)
json.dump({'p': p, 'chi': [[0, 0, 0, 0]] * m, 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads2], 'D': [core.to_str(D) for D in Ds[ng:]]}, open(prefix + '.quads.json', 'w'))
info.update({'n_unknowns': m, 'n_quadrics': len(quads2), 'quadric_rank': rank, 'n_pair_monomials': len(pairs), 'forced_zero_by_squares': forced})
if neq:
    v = verdict(prefix + '.ms', threads=threads, timeout=tmo, verbose=1, logfile=prefix + '.msolve.log')
    info['msolve'] = {k: val for k, val in v.items() if k not in ('log_tail', 'pure_powers')}
    info['result'] = 'EXCLUDED (only c=0)' if v.get('verdict', '').startswith('ONLY-ORIGIN') else v.get('verdict')
else:
    info['result'] = 'NO EQUATIONS'
info['seconds'] = round(time.time() - t0, 1)
log('closure-Gamma-invariant', **info, note='Gamma = G_aff x| <s0, r1>; gauge-fixed Hessian relaxation of the Gamma-invariant E1 kernel; wave1 (p=2^31-1): d=8 20 params, 3 gauge, 52 quadrics rank 35')
print(json.dumps(info))
