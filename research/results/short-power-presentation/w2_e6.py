import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_families import form_polys
from collections import Counter
F, A, q, rec = load(); FP = form_polys()
def residual_space(lname, d, par, p):
    lp = reduce_mod(FP[lname], p); pare = (par + core.parity_of(lp)) % 2
    gens = jacobian.jacobian_generators(F, d + 1, p, pare); Dm = monomials(d, par)
    cols = [g[2] for g in gens] + [shift(lp, m, 1, p) for m in Dm]
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols); flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p); X, nul = M.nullspace(); nJ = len(gens)
    Ds = [{Dm[j]: int(X[nJ + j, c]) for j in range(len(Dm)) if int(X[nJ + j, c])} for c in range(nul)]
    Ds = [D for D in Ds if D]
    jref = JRef(F, d, par, p, verbose=False)
    nfs = [nf for nf in jref.nf_polys(Ds) if nf]
    sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}
    Q = nmod_mat(len(nfs), len(sm), [nf.get(jref.mons[sm[t]], 0) for nf in nfs for t in range(len(sm))], p)
    Qr, rk = Q.rref()
    return [{jref.mons[sm[t]]: int(Qr[i, t]) for t in range(len(sm)) if int(Qr[i, t])} for i in range(rk)], jref
out = {}
for lname in ('L', 'M'):
    lists = {}
    for p in PGB[:2] + core.PRIMES[:2]:
        reps, _ = residual_space(lname, 6, 0, p); lists[p] = reps
        rat, info = reconstruct_polys(lists) if len(lists) >= 2 else (None, {})
        if rat: break
    print(f'Ann({lname})_6^+ : dim {len(lists[p])}, reconstruction {"ok" if rat else "failed"} {info}', flush=True)
    if rat:
        E = rat[0]
        print(f'  terms {len(E)}, Y-degrees {sorted(Counter(sum(e[3:]) for e in E).items())}, X-degrees {sorted(Counter(sum(e[:3]) for e in E).items())}', flush=True)
        dY = [core.pderiv(E, i) for i in range(3, 6)]; dX = [core.pderiv(E, i) for i in range(3)]
        print('  depends on Y only via M:', core.sub(dY[1], core.scale(dY[0], 3)) == {} and core.sub(dY[2], core.scale(dY[0], 5)) == {},
              '; on X only via L:', core.sub(dX[1], core.scale(dX[0], 2)) == {} and core.sub(dX[2], core.scale(dX[0], 3)) == {}, flush=True)
        okE1, _ = fixed_d.e1_solvable(F, A, 1, E)
        print('  E1 for A:', okE1, flush=True)
        out[lname] = core.to_str(E)
json.dump(out, open('E6.json', 'w'))
# is D7 ≡ (a0 X0 + a1 X1 + a2 X2) * E6 mod J_7 ?  and D7' ≡ (b Y) * E6' ?
p = PGB[0]
D7 = reduce_mod(core.parse(json.load(open('kernel_d7_p0.json'))['basis_rational'][0]), p)
if 'L' in out:
    E6 = reduce_mod(core.parse(out['L']), p)
    jref7 = JRef(F, 7, 0, p, verbose=False); sm = jref7.std; sidx = {jref7.mons[j]: t for t, j in enumerate(sm)}
    cands = [shift(E6, tuple(1 if t == i else 0 for t in range(NV)), 1, p) for i in range(3)]
    if 'M' in out:
        E6M = reduce_mod(core.parse(out['M']), p)
        cands += [shift(E6M, tuple(1 if t == i else 0 for t in range(NV)), 1, p) for i in range(3)]
    nfs = jref7.nf_polys(cands + [D7])
    def vec(nf): return [nf.get(jref7.mons[sm[t]], 0) for t in range(len(sm))]
    Mx = nmod_mat(len(nfs), len(sm), [x for nf in nfs for x in vec(nf)], p)
    r_c = nmod_mat(len(nfs) - 1, len(sm), [x for nf in nfs[:-1] for x in vec(nf)], p).rank()
    print(f'rank(X_i E6, X_i E6M) = {r_c}; with D7: {Mx.rank()} -> D7 in span(X_i*E6, X_i*E6M) mod J_7: {Mx.rank() == r_c}', flush=True)
