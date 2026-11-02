"""Hessian relaxation of the polar family D = w . grad F, deg w = e = d-5, for quadratic-cubic-pullback.
E2 => A (w^T Hess(F) w) in J_{2d};  if [A] is injective on R_{2d-6} (both parities) this gives
      phi(w) := w^T Hess(F) w in J_{2d-6}   (a system of quadratic forms in the coefficients of w).
Gauge: w -> w + l*Z (l in S_{e-1}) leaves phi mod J unchanged (Euler); gauge-fixed complement W.
Certificate of emptiness of the projective zero set in P(W): Macaulay matrix in degree delta has full
column rank mod p  =>  (w)^delta subset of I over Q  =>  only w = 0 (rigorous over Q for homogeneous systems).
Usage: python3 hess.py d mode prime [delta]      mode in {mixed, even, odd}  (parity of D)
Writes receipt hess_d{d}_{mode}_p{p}.json and (for msolve) hess_d{d}_{mode}_p{p}.ms (homogeneous input)."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *
from math import comb
d = int(sys.argv[1]); mode = sys.argv[2]; p = int(sys.argv[3]); delta_ = int(sys.argv[4]) if len(sys.argv) > 4 else None
e = d - 5
F, A, q = load(); S, R = hilbert()
t0 = time.time()
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
# basis of w: component i gets monomials m of degree e with parity depending on mode
def comp_pars(mode):
    if mode == 'mixed': return [None] * NV
    pD = 0 if mode == 'even' else 1
    # w_i F_i must have parity pD; F_i has parity 0 (i<3) / 1 (i>=3)
    return [(pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
pars = comp_pars(mode)
basis = []
for i in range(NV):
    for m in monomials(e, pars[i]):
        basis.append((i, m))
nvar = len(basis); bidx = {b: k for k, b in enumerate(basis)}
# gauge vectors: l in S_{e-1} with l*Z_i of the right parity for every i (mixed: all l; pure: parity of l fixed)
gauge = []
if e >= 1:
    lpar = None if mode == 'mixed' else (pars[0] - 0) % 2   # l*X0 must have parity pars[0] -> l parity = pars[0]
    for l in monomials(e - 1, lpar):
        v = [0] * nvar
        for i in range(NV):
            ei = [0] * NV; ei[i] = 1
            lm = tuple(a + b for a, b in zip(l, ei))
            if (i, lm) in bidx: v[bidx[(i, lm)]] = 1
        gauge.append(v)
ng = len(gauge)
# gauge fix: rref of gauge matrix, drop pivot variables
if ng:
    Gm = nmod_mat(ng, nvar, [x for v in gauge for x in v], p)
    Gmr, grk = Gm.rref()
    gpiv = []
    for r in range(grk):
        j = 0
        while int(Gmr[r, j]) == 0: j += 1
        gpiv.append(j)
    assert grk == ng
else:
    gpiv = []
keep = [k for k in range(nvar) if k not in set(gpiv)]
nk = len(keep); kpos = {k: t for t, k in enumerate(keep)}
# target space R_{2d-6} with parity split; normal forms via mod-p reducer per parity
deg = 2 * d - 6
reds = {par: jacobian._Reducer(F, deg, p, par) for par in (0, 1)}
hilbJ = {par: S[deg][par] - R[deg][par] for par in (0, 1)}
assert all(reds[par].rank == hilbJ[par] for par in (0, 1)), 'J-block rank mod p must equal Hilbert value'
std = {par: reds[par].standard for par in (0, 1)}
sidx = {par: {m: j for j, m in enumerate(std[par])} for par in (0, 1)}
# quadratic forms: for each pair (k<=l) of kept variables, normal form of H[ik][il] * m_k m_l, split by parity
t1 = time.time()
# vectorized normal form per parity: build V matrix (pairs x monomials), NF = V - V[:,piv] * Gr
pairs = [(k, l) for a, k in enumerate(keep) for l in keep[a:]]
rowsV = {0: [], 1: []}; pairpar = []
for (k, l) in pairs:
    ik, mk = basis[k]; il, ml = basis[l]
    P_ = shift(H[ik][il], tuple(a + b for a, b in zip(mk, ml)), 1, p)
    par = core.parity_of(P_) if P_ else None
    pairpar.append(par)
    if P_:
        rowsV[par].append((len(pairpar) - 1, P_))
# in mixed mode a pair poly has pure parity anyway (monomials times Hessian entries of fixed parity)
eqs = {}   # (par, row index in R_deg^par) -> dict {pair index: coeff}
for par in (0, 1):
    if not rowsV[par]: continue
    red = reds[par]; mons = red.mons; idx = red.idx
    npr = len(rowsV[par]); nm = len(mons)
    flat = [0] * (npr * nm)
    for r, (pi_, P_) in enumerate(rowsV[par]):
        base = r * nm
        for mm, c in P_.items(): flat[base + idx[mm]] = c
    V = nmod_mat(npr, nm, flat, p); del flat
    piv = red.pivots
    # Gr: rank x nm rref rows
    Grflat = [0] * (red.rank * nm)
    for r, row in enumerate(red.rows):
        for mm, c in row.items(): Grflat[r * nm + idx[mm]] = c
    Gr = nmod_mat(red.rank, nm, Grflat, p); del Grflat
    VP = nmod_mat(npr, red.rank, [int(V[r, j]) for r in range(npr) for j in piv], p)
    NF = V - VP * Gr
    nonpiv = [j for j in range(nm) if j not in set(piv)]
    for r, (pi_, P_) in enumerate(rowsV[par]):
        for jj, j in enumerate(nonpiv):
            c = int(NF[r, j])
            if c: eqs.setdefault((par, jj), {})[pi_] = c
t2 = time.time()
neq = len(eqs)
print(f'[hess d={d} {mode} p={p}] nvar={nvar} gauge={ng} kept={nk}; target R_{deg} dims {R[deg]}; '
      f'{neq} nonzero quadratic equations ({t2-t0:.1f}s)', flush=True)
eqlist = list(eqs.values())
# generic-Hilbert prediction of the certificate degree
def generic_delta(n, m):
    # smallest delta with coefficient of t^delta in (1-t^2)^m/(1-t)^n  <= 0
    for dl in range(2, 12):
        c = sum((-1) ** k * comb(m, k) * comb(n - 1 + dl - 2 * k, n - 1) for k in range(0, dl // 2 + 1))
        if c <= 0: return dl
    return None
gd = generic_delta(nk, neq)
rec = {'exp': 'hessian_relaxation', 'd': d, 'mode': mode, 'p': p, 'nvar': nvar, 'gauge_dim': ng, 'kept': nk,
       'target_deg': deg, 'R_target_dims': R[deg], 'J_rank_attained': True, 'n_equations': neq,
       'generic_certificate_degree': gd, 'seconds_build': round(t2 - t0, 1)}
# write msolve homogeneous input (variables c_k for kept k)
varnames = [f'c{k}' for k in keep]
fn = f'{PK}/hess_d{d}_{mode}_p{p}.ms'
with open(fn, 'w') as f:
    f.write(','.join(varnames) + '\n'); f.write(f'{p}\n')
    lines = []
    for eq in eqlist:
        terms = []
        for pi_, c in eq.items():
            k, l = pairs[pi_]
            terms.append(f'{c}*c{k}*c{l}' if k != l else f'{c}*c{k}^2')
        lines.append('+'.join(terms))
    f.write(',\n'.join(lines) + '\n')
rec['msolve_file'] = fn
# Macaulay certificate at degree delta (dense, only when feasible)
if delta_ is None: delta_ = gd
ncols = comb(nk + delta_ - 1, delta_)
rec['macaulay_delta'] = delta_; rec['macaulay_cols'] = ncols; rec['macaulay_rows'] = neq * comb(nk + delta_ - 3, delta_ - 2)
if ncols * rec['macaulay_rows'] <= 3e8:
    t3 = time.time()
    # monomials in kept variables of degree delta and delta-2
    def kmons(dd):
        return list(itertools.combinations_with_replacement(range(nk), dd))
    colm = kmons(delta_); cidx = {m: j for j, m in enumerate(colm)}
    mulm = kmons(delta_ - 2)
    Rr = len(mulm) * neq
    flat = [0] * (Rr * ncols)
    r = 0
    for mm in mulm:
        for eq in eqlist:
            base = r * ncols
            for pi_, c in eq.items():
                k, l = pairs[pi_]
                mon = tuple(sorted(mm + (kpos[k], kpos[l])))
                j = cidx[mon]
                flat[base + j] = (flat[base + j] + c) % p
            r += 1
    Mm = nmod_mat(Rr, ncols, flat, p); del flat
    rk = Mm.rank()
    rec['macaulay_rank'] = rk; rec['macaulay_full_column_rank'] = (rk == ncols)
    rec['seconds_macaulay'] = round(time.time() - t3, 1)
    rec['conclusion'] = ('projective zero set of the Hessian relaxation is EMPTY (rigorous over Q): no degree-%d polar '
                         'denominator in mode %s passes E2' % (d, mode)) if rk == ncols else \
                        'Macaulay rank deficient at this degree (no conclusion at delta=%d)' % delta_
else:
    rec['conclusion'] = 'Macaulay matrix too large for dense rank; use msolve on the .ms file'
rec['seconds'] = round(time.time() - t0, 1)
print(json.dumps(rec), flush=True)
log(rec); receipt(f'hess_d{d}_{mode}_p{p}.json', rec)
