"""Polynomial (degree-1 in the family parameters) Farkas certificates for a linear family of denominators
D(a) = sum_k a_k D_k lying in ONE character sector (lens excl, lane cycle-q2).

The sector system is M(a) x = b(a) with M(a) = sum_k a_k M_k (N-columns) | M_Q (Q-columns, constant),
b(a) = -A D(a)^2.  We look for y(a) = sum_k a_k y_k with y(a)^T M(a) == 0 identically in a:
    y_k^T M_Q = 0 (all k),  y_k^T M_k = 0,  y_k^T M_l + y_l^T M_k = 0 (k < l).
Each solution gives the cubic form phi_y(a) = y(a)^T b(a); every D(a) with phi_y(a) != 0 admits no certificate.
Output: cubic forms (mod p) and a Singular script intersecting V(phi_y) with the relaxation variety of the family.
usage: python3 wave2-excl/polyfarkas.py --family trivial6 [--prime P]
"""
import sys, json, time, argparse, os, subprocess, re, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hodge import core, symmetry
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod
from flint import nmod_mat

def family_polys(name, F):
    X = ['X0','X1','X2','Y0','Y1','Y2']
    if name == 'trivial6':
        return [mul(core.parse(X[i]), pderiv(F, i)) for i in range(6)], 'D(a) = sum_i a_i z_i F_i (trivial sector, deg 6)'
    if name == 'test-FX0':
        return [pderiv(F, 0)], 'test: single D = F_X0 (k = 1)'
    if name == 'test-FX0-FX1pair':
        return [pderiv(F, 0), mul(core.parse('X1'), pderiv(F, 1))], 'test: mixed-degree (invalid) family'
    raise ValueError(name)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--family', default='trivial6'); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--tag', default=None)
    a = ap.parse_args(); p = a.prime
    t0 = time.time()
    F, A, q, G = w2lib.load()
    Ds, desc = family_polys(a.family, F)
    k = len(Ds)
    tag = a.tag or f'{a.family}_p{p}'
    allowed, sub = symmetry.sector_filter(F, A, q, Ds[0], G)
    systems = [symmetry.RestrictedFixedDSystem(F, A, q, D, p, allowed) for D in Ds]
    cols = systems[0].cols
    for S in systems[1:]:
        assert S.cols == cols, 'column sets differ within the family'
    rows = {}
    for S in systems:
        for col in S.colpolys:
            for e in col: rows.setdefault(e, len(rows))
        for e in S.rhs: rows.setdefault(e, len(rows))
    R = len(rows); n = len(cols)
    Nidx = [j for j, c in enumerate(cols) if c[0] == 'N']; Qidx = [j for j, c in enumerate(cols) if c[0] == 'Q']
    print(f'family {a.family}: k={k}, rows R={R}, cols n={n} (N {len(Nidx)}, Q {len(Qidx)}), stabilizer {len(sub)}', flush=True)
    # column vectors as dict row->val for each system
    colvec = [[{rows[e]: v for e, v in col.items()} for col in S.colpolys] for S in systems]
    # equations: unknown index (kk, r) -> kk*R + r
    eqs = []   # list of dict unknown->coef
    for j in Qidx:
        for kk in range(k):
            eqs.append({kk * R + r: v for r, v in colvec[kk][j].items()})   # y_kk . colQ_j(kk)  (Q columns identical across the family)
    for j in Nidx:
        for kk in range(k):
            for ll in range(kk, k):
                eq = {}
                if kk == ll:
                    for r, v in colvec[kk][j].items(): eq[kk * R + r] = v
                else:
                    for r, v in colvec[ll][j].items(): eq[kk * R + r] = (eq.get(kk * R + r, 0) + v) % p
                    for r, v in colvec[kk][j].items(): eq[ll * R + r] = (eq.get(ll * R + r, 0) + v) % p
                eq = {u: v for u, v in eq.items() if v}
                if eq: eqs.append(eq)
    E, U = len(eqs), k * R
    print(f'equations {E} x unknowns {U}  [{time.time()-t0:.1f}s]', flush=True)
    flat = [0] * (E * U)
    for i, eq in enumerate(eqs):
        base = i * U
        for u, v in eq.items(): flat[base + u] = v
    M = nmod_mat(E, U, flat, p)
    del flat
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    print(f'degree-1 polynomial Farkas space: nullity {nullity}  [nullspace {t2-t1:.1f}s]', flush=True)
    # rhs vectors b_{jl} = -A D_j D_l  (rows)
    Ap = reduce_mod(A, p)
    Dp = [reduce_mod(D, p) for D in Ds]
    bvec = {}
    for j in range(k):
        for l in range(j, k):
            P = scale(mul(Ap, mul(Dp[j], Dp[l], p), p), p - 1, p)
            bvec[(j, l)] = {rows[e]: v for e, v in P.items() if e in rows}
            assert all(e in rows for e in P), 'rhs monomial outside row set?!'
    # cubic forms phi_y(a) = sum_kk a_kk y_kk . b(a),  b(a) = sum_{j<=l} (1 or 2) a_j a_l b_{jl}
    cubics = []
    vars_ = [f'a{i}' for i in range(k)]
    for c in range(nullity):
        y = [[int(X[kk * R + r, c]) for r in range(R)] for kk in range(k)]
        coef = {}
        for kk in range(k):
            for (j, l), bv in bvec.items():
                s = sum(y[kk][r] * v for r, v in bv.items()) % p
                if not s: continue
                mult = 1 if j == l else 2
                key = tuple(sorted((kk, j, l)))
                coef[key] = (coef.get(key, 0) + mult * s) % p
        coef = {kk: v for kk, v in coef.items() if v}
        if coef:
            cubics.append('+'.join(f'{v}*' + '*'.join(vars_[i] for i in kk) for kk, v in sorted(coef.items())))
    print(f'nonzero cubic forms: {len(cubics)}  [{time.time()-t0:.1f}s]', flush=True)
    os.makedirs(f'{w2lib.PACKET}/polyfarkas', exist_ok=True)
    with open(f'{w2lib.PACKET}/polyfarkas/{tag}.cubics.txt', 'w') as f:
        f.write('\n'.join(cubics) + '\n')
    # Singular: relaxation ideal of the family (recomputed here from the explicit quadrics for trivial6) + cubics
    sing = f'{w2lib.PACKET}/polyfarkas/{tag}.sing'
    with open(sing, 'w') as f:
        f.write(f'ring r = {p}, ({",".join(vars_)}), dp;\n')
        if a.family == 'trivial6':
            f.write('ideal V0 = (3*a0+a1+a3+a4)^2 - 6*(3*a0^2+a1^2+a3^2+a4^2), (3*a1+a2+a4+a5)^2 - 6*(3*a1^2+a2^2+a4^2+a5^2), (3*a2+a0+a5+a3)^2 - 6*(3*a2^2+a0^2+a5^2+a3^2);\n')
        else:
            f.write('ideal V0 = 0;\n')
        f.write('ideal C = ' + (',\n'.join(cubics) if cubics else '0') + ';\n')
        f.write('ideal W = std(V0 + C); "DIMW", dim(W); if (dim(W) <= 1) { "VDIMW", vdim(W); }\n')
        f.write('ideal Eu = a0-a5,a1-a5,a2-a5,a3-a5,a4-a5; "EULER_IN_W (0 = yes)", size(reduce(W, std(Eu)));\n')
        f.write('LIB "elim.lib"; def WS = sat(W, Eu); ideal Wsat; if (typeof(WS) == "list") { Wsat = WS[1]; } else { Wsat = WS; }\n')
        f.write('Wsat = std(Wsat); "DIM_W_minus_Euler", dim(Wsat); if (dim(Wsat) <= 1) { "VDIM_W_minus_Euler", vdim(Wsat); }\nquit;\n')
    r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=1800)
    print(r.stdout[-1500:], flush=True)
    dimW = re.search(r'DIMW (\d+)', r.stdout); vd = re.search(r'VDIMW (-?\d+)', r.stdout); dimWE = re.search(r'DIM_W_minus_Euler (-?\d+)', r.stdout); vdE = re.search(r'VDIM_W_minus_Euler (-?\d+)', r.stdout)
    res = {'exp': 'polyfarkas-deg1', 'family': desc, 'tag': tag, 'd': core.degree(Ds[0]), 'prime': p, 'rows': R, 'cols': n, 'equations': E, 'unknowns': U,
           'nullity': nullity, 'ncubics': len(cubics), 'dimW_affine': int(dimW.group(1)) if dimW else None, 'vdimW': int(vd.group(1)) if vd else None,
           'dimW_minus_Euler_affine': int(dimWE.group(1)) if dimWE else None, 'vdimW_minus_Euler': int(vdE.group(1)) if vdE else None,
           'seconds': round(time.time() - t0, 1), 'job': os.environ.get('W2JOB', 'polyfarkas')}
    res['result'] = (f"W = V0 cap V(cubics): affine dim {res['dimW_affine']}; after removing the Euler point (D = F): affine dim {res['dimW_minus_Euler_affine']}"
                     + (' => family EXCLUDED mod p except D = F' if res['dimW_minus_Euler_affine'] == 0 else ''))
    json.dump(res, open(f'{w2lib.PACKET}/polyfarkas/{tag}.receipt.json', 'w'), indent=1)
    w2lib.log(**res)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main()
