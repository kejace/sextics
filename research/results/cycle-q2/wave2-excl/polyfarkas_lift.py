"""Exact (characteristic-0) degree-1 polynomial Farkas certificates for a one-sector linear family D(a) = sum_k a_k D_k.

Builds the EXACT integer system E (equations: coefficients of a_k a_l in y(a)^T M(a), y(a) = sum_k a_k y_k) from the
exact integer columns (F2 = 2F, integral D_k), finds its pivot structure mod p, then lifts several nullspace vectors
exactly (one Dixon solve per vector: square pivot system, free coordinates = unit vectors), verifies y(a)^T M(a) == 0
over Z, and writes the cubic forms phi_y(a) = y(a)^T b(a) over Z plus a Singular script for W = V0 cap V(phi_y) over Q.
usage: python3 wave2-excl/polyfarkas_lift.py --family trivial6 [--nvec 6] [--prime P]
"""
import sys, json, time, argparse, os, subprocess, re, hashlib
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, polyfarkas
from hodge import core, symmetry
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod
from flint import nmod_mat, fmpz_mat, fmpq_mat
from fractions import Fraction

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--family', default='trivial6'); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--nvec', type=int, default=6); ap.add_argument('--tag', default=None)
    a = ap.parse_args(); p = a.prime
    t0 = time.time()
    F, A, q, G = w2lib.load()
    Ds, desc = polyfarkas.family_polys(a.family, F)
    k = len(Ds); tag = a.tag or f'{a.family}_lift'
    allowed, sub = symmetry.sector_filter(F, A, q, Ds[0], G)
    S0 = symmetry.RestrictedFixedDSystem(F, A, q, Ds[0], p, allowed)
    cols = S0.cols
    # exact integer columns for each D_k (same column list)
    exact = [w2lib.exact_columns(F, A, q, D, cols) for D in Ds]   # (colpolys, rhs, Dz)
    rows = {}
    for colpolys, rhs, Dz in exact:
        for col in colpolys:
            for e in col: rows.setdefault(e, len(rows))
        for e in rhs: rows.setdefault(e, len(rows))
    R = len(rows); n = len(cols)
    Nidx = [j for j, c in enumerate(cols) if c[0] == 'N']; Qidx = [j for j, c in enumerate(cols) if c[0] == 'Q']
    colvec = [[{rows[e]: v for e, v in col.items()} for col in colpolys] for colpolys, _, _ in exact]
    print(f'{desc}: k={k} R={R} n={n}', flush=True)
    eqs = []
    for j in Qidx:
        for kk in range(k):
            eqs.append({kk * R + r: v for r, v in colvec[kk][j].items()})
    for j in Nidx:
        for kk in range(k):
            for ll in range(kk, k):
                eq = {}
                if kk == ll:
                    for r, v in colvec[kk][j].items(): eq[kk * R + r] = v
                else:
                    for r, v in colvec[ll][j].items(): eq[kk * R + r] = eq.get(kk * R + r, 0) + v
                    for r, v in colvec[kk][j].items(): eq[ll * R + r] = eq.get(ll * R + r, 0) + v
                eq = {u: v for u, v in eq.items() if v}
                if eq: eqs.append(eq)
    E, U = len(eqs), k * R
    print(f'exact system {E} x {U}  [{time.time()-t0:.1f}s]', flush=True)
    # mod p: rank, pivot columns (rref of E) and pivot rows (rref of E^T)
    flat = [0] * (E * U)
    for i, eq in enumerate(eqs):
        base = i * U
        for u, v in eq.items(): flat[base + u] = v % p
    Mp = nmod_mat(E, U, flat, p)
    Mr, rank = Mp.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j < U and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    free = [j for j in range(U) if j not in set(piv)]
    # mod-p nullspace vectors (one per free column, free coordinate = 1) and their cubic forms: keep only those with phi != 0
    X, nullity = Mp.nullspace()
    assert nullity == U - rank
    bmod = {}
    Dzs_ = [ex[2] for ex in exact]
    for j in range(k):
        for l in range(j, k):
            P = scale(mul(A, mul(Dzs_[j], Dzs_[l])), -1)
            bmod[(j, l)] = {rows[e]: int(v) % p for e, v in P.items()}
    good = []
    for c in range(nullity):
        yc = [int(X[u, c]) for u in range(U)]
        nz = False
        for kk in range(k):
            for (j, l), bv in bmod.items():
                if sum(yc[kk * R + r] * v for r, v in bv.items()) % p:
                    nz = True; break
            if nz: break
        if nz: good.append(c)
    print(f'nullspace vectors with nonzero cubic (mod p): {len(good)} of {nullity}', flush=True)
    # map nullspace column c -> its free coordinate (the unique free coordinate with entry 1 and others 0)
    goodfree = []
    for c in good:
        fc = [u for u in free if int(X[u, c]) != 0]
        goodfree.append(fc[0] if len(fc) == 1 else None)
    del X
    del Mr
    print(f'rank {rank}, nullity {U - rank}  [{time.time()-t0:.1f}s]', flush=True)
    flatT = [0] * (U * E)
    for i, eq in enumerate(eqs):
        for u, v in eq.items(): flatT[u * E + i] = v % p
    MT = nmod_mat(U, E, flatT, p)
    MTr, rankT = MT.rref()
    assert rankT == rank
    prow = []
    for i in range(rankT):
        j = 0
        while j < E and int(MTr[i, j]) == 0: j += 1
        prow.append(j)
    del MTr, MT, Mp
    print(f'pivot rows found  [{time.time()-t0:.1f}s]', flush=True)
    # square exact system: E[prow, piv] y_piv = - E[prow, free] f
    pivpos = {u: i for i, u in enumerate(piv)}
    Esq = fmpz_mat(rank, rank)
    for i, ri in enumerate(prow):
        for u, v in eqs[ri].items():
            if u in pivpos: Esq[i, pivpos[u]] = v
    sel = [fc for fc in goodfree if fc is not None]
    # spread the selection over the list
    if len(sel) > a.nvec:
        step = len(sel) // a.nvec
        sel = [sel[i * step] for i in range(a.nvec)]
    nvec = len(sel)
    os.makedirs(f'{w2lib.PACKET}/polyfarkas', exist_ok=True)
    Ap = A; Dzs = [ex[2] for ex in exact]
    bvec = {}
    for j in range(k):
        for l in range(j, k):
            P = scale(mul(Ap, mul(Dzs[j], Dzs[l])), -1)
            bvec[(j, l)] = {rows[e]: int(v) for e, v in P.items()}
    cubics = []; vectors = []
    vars_ = [f'a{i}' for i in range(k)]
    for t in range(nvec):
        fcoord = sel[t]
        rhs = fmpz_mat(rank, 1)
        for i, ri in enumerate(prow):
            v = eqs[ri].get(fcoord, 0)
            if v: rhs[i, 0] = -v
        t1 = time.time()
        ysol = Esq.solve(rhs)
        y = [Fraction(0)] * U
        y[fcoord] = Fraction(1)
        for i, u in enumerate(piv):
            y[u] = Fraction(int(ysol[i, 0].p), int(ysol[i, 0].q))
        # exact verification of ALL equations
        bad = 0
        for eq in eqs:
            s = Fraction(0)
            for u, v in eq.items():
                if y[u]: s += y[u] * v
            if s != 0: bad += 1
        maxbits = max(max(v.numerator.bit_length(), v.denominator.bit_length()) for v in y)
        print(f'vector {t}: solve {time.time()-t1:.1f}s, bad equations {bad}, max bits {maxbits}', flush=True)
        if bad: continue
        coef = {}
        for kk in range(k):
            for (j, l), bv in bvec.items():
                s = Fraction(0)
                for r, v in bv.items():
                    if y[kk * R + r]: s += y[kk * R + r] * v
                if s == 0: continue
                mult = 1 if j == l else 2
                key = tuple(sorted((kk, j, l)))
                coef[key] = coef.get(key, Fraction(0)) + mult * s
        coef = {kk: v for kk, v in coef.items() if v}
        if coef:
            den = 1
            for v in coef.values():
                den = den * v.denominator // __import__('math').gcd(den, v.denominator)
            cubics.append('+'.join(f'({int(v*den)})*' + '*'.join(vars_[i] for i in kk) for kk, v in sorted(coef.items())))
            vectors.append({'free_coord': fcoord, 'sha256': hashlib.sha256('\n'.join(str(v) for v in y).encode()).hexdigest(), 'maxbits': maxbits})
    with open(f'{w2lib.PACKET}/polyfarkas/{tag}.cubics_Q.txt', 'w') as f:
        f.write('\n'.join(cubics) + '\n')
    sing = f'{w2lib.PACKET}/polyfarkas/{tag}_Q.sing'
    with open(sing, 'w') as f:
        f.write(f'LIB "elim.lib";\nring r = 0, ({",".join(vars_)}), dp;\n')
        if a.family == 'trivial6':
            f.write('ideal V0 = (3*a0+a1+a3+a4)^2 - 6*(3*a0^2+a1^2+a3^2+a4^2), (3*a1+a2+a4+a5)^2 - 6*(3*a1^2+a2^2+a4^2+a5^2), (3*a2+a0+a5+a3)^2 - 6*(3*a2^2+a0^2+a5^2+a3^2);\n')
        else:
            f.write('ideal V0 = 0;\n')
        f.write('ideal C = ' + (',\n'.join(cubics) if cubics else '0') + ';\n')
        f.write('ideal W = std(V0 + C); "DIMW", dim(W); if (dim(W) <= 1) { "VDIMW", vdim(W); }\n')
        if a.family == 'trivial6':
            f.write('ideal Eu = a0-a5,a1-a5,a2-a5,a3-a5,a4-a5; "EULER_IN_W (0 = yes)", size(reduce(W, std(Eu)));\n')
            f.write('def WS = sat(W, Eu); ideal Wsat; if (typeof(WS) == "list") { Wsat = WS[1]; } else { Wsat = WS; }\nWsat = std(Wsat); "DIM_W_minus_Euler", dim(Wsat); if (dim(Wsat) <= 1) { "VDIM_W_minus_Euler", vdim(Wsat); }\n')
        f.write('quit;\n')
    r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=3600)
    print(r.stdout[-1200:], flush=True)
    dimW = re.search(r'DIMW (\d+)', r.stdout); dimWE = re.search(r'DIM_W_minus_Euler (-?\d+)', r.stdout); vdE = re.search(r'VDIM_W_minus_Euler (-?\d+)', r.stdout)
    res = {'exp': 'polyfarkas-deg1-exact', 'family': desc, 'tag': tag, 'd': core.degree(Ds[0]), 'prime': p, 'rows': R, 'cols': n, 'equations': E, 'unknowns': U,
           'rank_mod_p': rank, 'nullity_mod_p': U - rank, 'vectors_lifted': len(vectors), 'vectors': vectors, 'ncubics_Q': len(cubics),
           'dimW_affine_Q': int(dimW.group(1)) if dimW else None, 'dimW_minus_Euler_affine_Q': int(dimWE.group(1)) if dimWE else None,
           'vdimW_minus_Euler_Q': int(vdE.group(1)) if vdE else None, 'seconds': round(time.time() - t0, 1), 'job': os.environ.get('W2JOB', 'polyfarkas_lift')}
    res['result'] = (f"exact degree-1 polynomial Farkas vectors over Q: {len(vectors)}; W = V0 cap V(phi): affine dim {res['dimW_affine_Q']}; minus Euler point: affine dim {res['dimW_minus_Euler_affine_Q']}"
                     + (' => family EXCLUDED over C except D = F (rigorous)' if res['dimW_minus_Euler_affine_Q'] == 0 else ''))
    json.dump(res, open(f'{w2lib.PACKET}/polyfarkas/{tag}.receipt.json', 'w'), indent=1)
    w2lib.log(**res)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main()
