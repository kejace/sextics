"""Reduced (Koszul-parametrized) certificate system for a one-sector linear family D(a) = sum_k a_k D_k, and
polynomial Farkas certificates of degree 1, 2, 3 in a on it (lens excl, lane cycle-q2).

Every E1 solution is N = N0(D) + (R' z)/12 + K/2 with K = alpha.grad F (alpha antisymmetric, deg d+2) and the Euler
part drops out of E identically, so a certificate exists iff
        H(D) + D div K - K.grad D = 2 F^2 Q,      H(D) = D div N0 - N0.grad D        (degree 2d+6),
with N0 a fixed E1 solution (for Jacobian-type D_k = z_k F_k: N0_k = (A/2) z_k e_k; otherwise from an E1 solve mod p).
Unknowns x = (alpha_{ij,m}, Q_m'); M(a) = sum_k a_k M_k is linear in a, b(a) = -H(D(a)) quadratic.
Polynomial Farkas of degree t: y(a) = sum_{|e|=t} a^e y_e with y(a)^T M(a) == 0 identically; phi_y(a) = y(a)^T b(a)
(degree t+2).  Every a with phi_y(a) != 0 for some y admits no certificate (mod p statement; exact lift separately).
usage: python3 wave2-excl/redsys.py --family trivial6 --degree 2 [--prime P] [--check-degree1]
"""
import sys, json, time, argparse, os, subprocess, re, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, polyfarkas, hess6
from hess6 import key_of, poly_key
from hodge import core, symmetry
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod, power
from flint import nmod_mat

def build(F, A, q, Ds, p, G, N0s=None):
    gens = G.generators
    F2 = scale(F, 2); Fp = reduce_mod(F2, p); Fi = [pderiv(Fp, i, p) for i in range(6)]
    Ap = reduce_mod(A, p)
    k = len(Ds); d = core.degree(Ds[0])
    Dp = [reduce_mod(D, p) for D in Ds]
    kA = poly_key(A, gens); kF = poly_key(F, gens); kD = poly_key(Ds[0], gens)
    wz = [key_of(tuple(1 if t == i else 0 for t in range(6)), gens) for i in range(6)]
    if N0s is None:
        raise ValueError('N0s required')
    # H bilinear parts: H(D(a)) = sum_k a_k^2 H_kk + sum_{k<l} 2 a_k a_l H_kl   (H_kl symmetric part)
    def Hpart(Da, Nb, Db, Na):
        out = {}; divNb = {}; divNa = {}
        for t in range(6):
            divNb = add(divNb, pderiv(Nb[t], t, p), p); divNa = add(divNa, pderiv(Na[t], t, p), p)
        out = add(mul(Da, divNb, p), mul(Db, divNa, p), p)
        for t in range(6):
            out = core.sub(out, mul(Nb[t], pderiv(Da, t, p), p), p); out = core.sub(out, mul(Na[t], pderiv(Db, t, p), p), p)
        return scale(out, pow(2, -1, p), p)
    Hb = {}
    for a_ in range(k):
        for b_ in range(a_, k):
            Hb[(a_, b_)] = Hpart(Dp[a_], N0s[b_], Dp[b_], N0s[a_])
    degE = 2 * d + 6
    kH = tuple((x + 2 * y - z) % 6 for x, y, z in zip(kA, kD, kF))
    # columns: Koszul pairs (i<j, m in S_{d+2}) with key(m) = kA + kD - 2kF + w_i + w_j ; Q: m' in S_{2d-6} with key = kA + 2kD - 3kF
    cols = []  # (kind, i, j, m)
    colvec = [[] for _ in range(k)]   # per family member: dict row-> val  (filled after row indexing)
    kQ = tuple((x + 2 * y - 3 * z) % 6 for x, y, z in zip(kA, kD, kF))
    Fq = scale(mul(Fp, Fp, p), p - 2, p)   # -2 F^2
    polys_k = [[] for _ in range(k)]; polys_Q = []
    for i in range(6):
        for j in range(i + 1, 6):
            need = tuple((x + y - 2 * z + u + v) % 6 for x, y, z, u, v in zip(kA, kD, kF, wz[i], wz[j]))
            for m in monomials(d + 2):
                if key_of(m, gens) != need: continue
                mFj = shift(Fi[j], m, 1, p); mFi = shift(Fi[i], m, 1, p)
                dmFj = pderiv(mFj, i, p); dmFi = pderiv(mFi, j, p)
                cols.append(('K', i, j, m))
                for kk in range(k):
                    D = Dp[kk]
                    col = core.sub(mul(D, core.sub(dmFj, dmFi, p), p), core.sub(mul(mFj, pderiv(D, i, p), p), mul(mFi, pderiv(D, j, p), p), p), p)
                    polys_k[kk].append(col)
    nK = len(cols)
    for m in monomials(2 * d - 6):
        if key_of(m, gens) != kQ: continue
        cols.append(('Q', -1, -1, m)); polys_Q.append(shift(Fq, m, 1, p))
    rows = {}
    for kk in range(k):
        for P in polys_k[kk]:
            for e in P: rows.setdefault(e, len(rows))
    for P in polys_Q:
        for e in P: rows.setdefault(e, len(rows))
    for P in Hb.values():
        for e in P: rows.setdefault(e, len(rows))
    R = len(rows)
    for kk in range(k):
        colvec[kk] = [{rows[e]: v for e, v in P.items()} for P in polys_k[kk]]
    Qvec = [{rows[e]: v for e, v in P.items()} for P in polys_Q]
    bvec = {key: {rows[e]: (p - v) % p for e, v in P.items()} for key, P in Hb.items()}   # b(a) = -H(a)
    return {'rows': rows, 'R': R, 'cols': cols, 'nK': nK, 'nQ': len(polys_Q), 'colvec': colvec, 'Qvec': Qvec, 'bvec': bvec, 'k': k, 'p': p, 'd': d}

def monos(k, t):
    return list(itertools.combinations_with_replacement(range(k), t))

def polyfarkas_degree(sysd, t):
    """degree-t polynomial Farkas space; returns (nullity, list of forms phi_y as dict monomial->coef, timings)."""
    p = sysd['p']; k = sysd['k']; R = sysd['R']; colvec = sysd['colvec']; Qvec = sysd['Qvec']; bvec = sysd['bvec']
    ms = monos(k, t); midx = {m: i for i, m in enumerate(ms)}
    U = len(ms) * R
    eqs = []
    # Q columns: y_e . Qcol = 0 for all e
    for jq, qv in enumerate(Qvec):
        for e in ms:
            eqs.append({midx[e] * R + r: v for r, v in qv.items()})
    # K columns: coefficient of a^f (|f| = t+1) in sum_e a^e y_e . (sum_l a_l col_l) = sum over (e, l) with e + l = f
    fs = monos(k, t + 1)
    for jc in range(sysd['nK']):
        for f in fs:
            eq = {}
            for pos in range(len(f)):
                l = f[pos]; e = f[:pos] + f[pos + 1:]
                if pos > 0 and f[pos] == f[pos - 1]: continue   # avoid duplicate splits
                base = midx[e] * R
                for r, v in colvec[l][jc].items():
                    eq[base + r] = (eq.get(base + r, 0) + v) % p
            eq = {u: v for u, v in eq.items() if v}
            if eq: eqs.append(eq)
    E = len(eqs)
    t0 = time.time()
    flat = [0] * (E * U)
    for i, eq in enumerate(eqs):
        base = i * U
        for u, v in eq.items(): flat[base + u] = v
    M = nmod_mat(E, U, flat, p); del flat
    X, nullity = M.nullspace(); del M
    t1 = time.time()
    forms = []
    for c in range(nullity):
        y = [int(X[u, c]) for u in range(U)]
        coef = {}
        for e in ms:
            base = midx[e] * R
            if not any(y[base:base + R]): continue
            for (j, l), bv in bvec.items():
                s = 0
                for r, v in bv.items():
                    if y[base + r]: s += y[base + r] * v
                s %= p
                if not s: continue
                mult = 1 if j == l else 2
                key = tuple(sorted(e + (j, l)))
                coef[key] = (coef.get(key, 0) + mult * s) % p
        coef = {kk: v for kk, v in coef.items() if v}
        if coef: forms.append(coef)
    return nullity, forms, {'E': E, 'U': U, 'nullspace_s': round(t1 - t0, 1), 'forms_s': round(time.time() - t1, 1)}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--family', default='trivial6'); ap.add_argument('--degree', type=int, default=2)
    ap.add_argument('--prime', type=int, default=w2lib.P1); ap.add_argument('--check-degree1', action='store_true'); ap.add_argument('--tag', default=None)
    a = ap.parse_args(); p = a.prime
    t0 = time.time()
    F, A, q, G = w2lib.load()
    Ds, desc = polyfarkas.family_polys(a.family, F)
    k = len(Ds)
    if a.family == 'trivial6':
        Z = [core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
        Ap = reduce_mod(A, p)
        N0s = [[reduce_mod(scale(mul(Ap, Z[i], p), pow(2, -1, p), p), p) if t == i else {} for t in range(6)] for i in range(6)]   # N0_i = (A/2) z_i e_i  (for D_i = z_i F_i)
    else:
        raise SystemExit('only trivial6 implemented here')
    sysd = build(F, A, q, Ds, p, G, N0s)
    print(f"reduced system: rows {sysd['R']}, K-cols {sysd['nK']}, Q-cols {sysd['nQ']}  [{time.time()-t0:.1f}s]", flush=True)
    tag = a.tag or f'red_{a.family}_deg{a.degree}_p{p}'
    os.makedirs(f'{w2lib.PACKET}/polyfarkas', exist_ok=True)
    degs = ([1] if a.check_degree1 else []) + [a.degree]
    vars_ = [f'a{i}' for i in range(k)]
    results = {}
    for t in degs:
        nullity, forms, tm = polyfarkas_degree(sysd, t)
        polys = ['+'.join(f'{v}*' + '*'.join(vars_[i] for i in kk) for kk, v in sorted(f.items())) for f in forms]
        sing = f'{w2lib.PACKET}/polyfarkas/{tag}_t{t}.sing'
        with open(sing, 'w') as f:
            f.write(f'LIB "elim.lib";\nring r = {p}, ({",".join(vars_)}), dp;\n')
            f.write('ideal V0 = (3*a0+a1+a3+a4)^2 - 6*(3*a0^2+a1^2+a3^2+a4^2), (3*a1+a2+a4+a5)^2 - 6*(3*a1^2+a2^2+a4^2+a5^2), (3*a2+a0+a5+a3)^2 - 6*(3*a2^2+a0^2+a5^2+a3^2);\n')
            f.write('ideal C = ' + (',\n'.join(polys) if polys else '0') + ';\n')
            f.write('ideal W = std(V0 + C); "DIMW", dim(W); if (dim(W) <= 1) { "VDIMW", vdim(W); }\n')
            f.write('ideal Eu = a0-a5,a1-a5,a2-a5,a3-a5,a4-a5; "EULER_IN_W (0 = yes)", size(reduce(W, std(Eu)));\n')
            f.write('def WS = sat(W, Eu); ideal Wsat; if (typeof(WS) == "list") { Wsat = WS[1]; } else { Wsat = WS; }\nWsat = std(Wsat); "DIM_W_minus_Euler", dim(Wsat); if (dim(Wsat) <= 1) { "VDIM_W_minus_Euler", vdim(Wsat); }\nquit;\n')
        r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=3600)
        dimW = re.search(r'DIMW (\d+)', r.stdout); dimWE = re.search(r'DIM_W_minus_Euler (-?\d+)', r.stdout); vdE = re.search(r'VDIM_W_minus_Euler (-?\d+)', r.stdout)
        res = {'degree': t, 'nullity': nullity, 'nforms': len(forms), 'dimW_affine': int(dimW.group(1)) if dimW else None,
               'dimW_minus_Euler_affine': int(dimWE.group(1)) if dimWE else None, 'vdimW_minus_Euler': int(vdE.group(1)) if vdE else None, **tm,
               'singular_tail': r.stdout[-300:]}
        results[t] = res
        print(json.dumps(res), flush=True)
        w2lib.log(exp=f'redsys-polyfarkas-deg{t}', family=desc, tag=tag, d=sysd['d'], prime=p, rows=sysd['R'], cols=sysd['nK'] + sysd['nQ'],
                  result=f"degree-{t} polynomial Farkas on the reduced (Koszul) system: nullity {nullity}, {len(forms)} nonzero forms of degree {t+2}; W = V0 cap V(forms): affine dim {res['dimW_affine']}; minus Euler point: affine dim {res['dimW_minus_Euler_affine']}" + (' => trivial-sector family EXCLUDED mod p except D = F' if res['dimW_minus_Euler_affine'] == 0 else ''),
                  rank={'E': tm['E'], 'U': tm['U']}, seconds=round(time.time() - t0, 1), job=os.environ.get('W2JOB', 'redsys'))
    json.dump({'family': desc, 'prime': p, 'rows': sysd['R'], 'nK': sysd['nK'], 'nQ': sysd['nQ'], 'results': results, 'seconds': round(time.time() - t0, 1)},
              open(f'{w2lib.PACKET}/polyfarkas/{tag}.receipt.json', 'w'), indent=1)
    print('done', round(time.time() - t0, 1))

if __name__ == '__main__':
    main()
