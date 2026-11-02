"""EXACT (characteristic-0) degree-1 polynomial Farkas certificates on the reduced (Koszul-parametrized) system of the
trivial-sector family D(a) = sum_k a_k z_k F_k (deg 6), lane cycle-q2, lens excl.

Rigorous chain for this family: (1) [exact algebra] a certificate with denominator D(a) forces H(D(a)) in J_18, i.e.
a in V0 = {(a.e(m))^2 = 6 sum_i e_i(m) a_i^2, m the three couplings}; (2) [exact linear algebra] every solution N of E1
is A w + K/2 + Euler gauge, so a certificate exists iff the reduced system M(a) x = b(a) (x = Koszul coefficients alpha
and Q) is consistent; (3) exact rational vectors y(a) = sum_k a_k y_k with y(a)^T M(a) == 0 identically give forms
phi_y(a) = y(a)^T b(a) of degree 3; any a with phi_y(a) != 0 admits no certificate; (4) Singular over Q computes
W = V0 cap V(phi_y).  If W = {Euler point} the family is excluded over C except D = 6F (illegal, gcd(D,F) = F).
usage: python3 wave2-excl/redsys_lift.py [--nvec 40] [--prime P]
"""
import sys, json, time, argparse, os, subprocess, re, itertools, hashlib, math
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hess6 import key_of, poly_key
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials
from fractions import Fraction
from flint import nmod_mat, fmpz_mat

def build_exact(F, A, G):
    """exact reduced system for D(a) = sum a_k z_k F2_k (F2 = 2F), N0 = A w: rows (monomials of degree 18 in the H-sector),
    K-columns (i<j, m) per family member, Q-columns, rhs blocks b_{kl}; all entries integers after row scaling."""
    gens = G.generators
    F2 = scale(F, 2); Fi = [pderiv(F2, i) for i in range(6)]
    Z = [core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
    Ds = [mul(Z[i], Fi[i]) for i in range(6)]                      # D_k = z_k F2_k (integral)
    N0s = [[mul(A, Z[i]) if t == i else {} for t in range(6)] for i in range(6)]   # N0_k = A z_k e_k :  N0(F2) = A D  (E1 for F2)
    k = 6; d = 6
    kA = poly_key(A, gens); kF = poly_key(F, gens); kD = poly_key(Ds[0], gens)
    wz = [key_of(tuple(1 if t == i else 0 for t in range(6)), gens) for i in range(6)]
    def Hpart(Da, Nb, Db, Na):
        divNb = {}; divNa = {}
        for t in range(6):
            divNb = add(divNb, pderiv(Nb[t], t)); divNa = add(divNa, pderiv(Na[t], t))
        out = add(mul(Da, divNb), mul(Db, divNa))
        for t in range(6):
            out = core.sub(out, mul(Nb[t], pderiv(Da, t))); out = core.sub(out, mul(Na[t], pderiv(Db, t)))
        return scale(out, Fraction(1, 2))
    # reduced identity with F2: E = A D^2 - 2 D N(F) + F(D div N - N.grad D) with N = N0 + K/2 (+Euler), 2N(F) = N(F2) = A D for N0
    # => E = F (D div N - N.grad D - D R') ... in terms of F2: E = (F2/2) [ H(D,N0) + (1/2)(D div K - K.grad D) ] = F^3 Q = (F2^3/8) Q
    # multiply by 8/F2:  4 H + 2 (D div K - K.grad D) = F2^2 Q.   Columns: K_{ij,m} -> 2 (D div K - K.grad D); Q_m -> -F2^2 m; rhs: -4 H.
    Hb = {}
    for a_ in range(k):
        for b_ in range(a_, k):
            Hb[(a_, b_)] = scale(Hpart(Ds[a_], N0s[b_], Ds[b_], N0s[a_]), -4)
    cols = []; polys_k = [[] for _ in range(k)]; polys_Q = []
    for i in range(6):
        for j in range(i + 1, 6):
            need = tuple((x + y - 2 * z + u + v) % 6 for x, y, z, u, v in zip(kA, kD, kF, wz[i], wz[j]))
            for m in monomials(d + 2):
                if key_of(m, gens) != need: continue
                mFj = shift(Fi[j], m, 1); mFi = shift(Fi[i], m, 1)
                dmFj = pderiv(mFj, i); dmFi = pderiv(mFi, j)
                cols.append(('K', i, j, m))
                for kk in range(k):
                    D = Ds[kk]
                    col = core.sub(mul(D, core.sub(dmFj, dmFi)), core.sub(mul(mFj, pderiv(D, i)), mul(mFi, pderiv(D, j))))
                    polys_k[kk].append(scale(col, 2))
    nK = len(cols)
    kQ = tuple((x + 2 * y - 3 * z) % 6 for x, y, z in zip(kA, kD, kF))
    F2sq = scale(mul(F2, F2), -1)
    for m in monomials(2 * d - 6):
        if key_of(m, gens) != kQ: continue
        cols.append(('Q', -1, -1, m)); polys_Q.append(shift(F2sq, m, 1))
    rows = {}
    for kk in range(k):
        for P in polys_k[kk]:
            for e in P: rows.setdefault(e, len(rows))
    for P in polys_Q:
        for e in P: rows.setdefault(e, len(rows))
    for P in Hb.values():
        for e in P: rows.setdefault(e, len(rows))
    R = len(rows)
    # row scaling to integers
    den = [1] * R
    def upd(P):
        for e, v in P.items():
            fr = Fraction(v); den[rows[e]] = den[rows[e]] * fr.denominator // math.gcd(den[rows[e]], fr.denominator)
    for kk in range(k):
        for P in polys_k[kk]: upd(P)
    for P in polys_Q: upd(P)
    for P in Hb.values(): upd(P)
    def vec(P):
        return {rows[e]: int(Fraction(v) * den[rows[e]]) for e, v in P.items() if Fraction(v) != 0}
    colvec = [[vec(P) for P in polys_k[kk]] for kk in range(k)]
    Qvec = [vec(P) for P in polys_Q]
    bvec = {key: vec(P) for key, P in Hb.items()}
    return {'rows': rows, 'R': R, 'cols': cols, 'nK': nK, 'nQ': len(polys_Q), 'colvec': colvec, 'Qvec': Qvec, 'bvec': bvec, 'k': k, 'd': d, 'Ds': Ds}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--nvec', type=int, default=40); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--tag', default='red_trivial6_exact')
    a = ap.parse_args(); p = a.prime
    t0 = time.time()
    F, A, q, G = w2lib.load()
    S = build_exact(F, A, G)
    k, R = S['k'], S['R']; colvec, Qvec, bvec = S['colvec'], S['Qvec'], S['bvec']
    print(f"exact reduced system: rows {R}, K-cols {S['nK']}, Q-cols {S['nQ']}  [{time.time()-t0:.1f}s]", flush=True)
    # sanity: at a random integer point, the reduced system must agree with the engine's full sector system (consistency)
    # degree-1 Farkas equations: y(a) = sum_k a_k y_k ; unknown index kk*R + r
    eqs = []
    for qv in Qvec:
        for kk in range(k):
            eqs.append({kk * R + r: v for r, v in qv.items()})
    for jc in range(S['nK']):
        for kk in range(k):
            for ll in range(kk, k):
                eq = {}
                if kk == ll:
                    for r, v in colvec[kk][jc].items(): eq[kk * R + r] = v
                else:
                    for r, v in colvec[ll][jc].items(): eq[kk * R + r] = eq.get(kk * R + r, 0) + v
                    for r, v in colvec[kk][jc].items(): eq[ll * R + r] = eq.get(ll * R + r, 0) + v
                eq = {u: v for u, v in eq.items() if v}
                if eq: eqs.append(eq)
    E, U = len(eqs), k * R
    print(f'exact Farkas equations {E} x {U}', flush=True)
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
    X, nullity = Mp.nullspace()
    print(f'rank {rank}, nullity {nullity}  [{time.time()-t0:.1f}s]', flush=True)
    # which nullspace vectors give a nonzero cubic mod p
    bmod = {key: {r: v % p for r, v in bv.items()} for key, bv in bvec.items()}
    good = []
    for c in range(nullity):
        yc = [int(X[u, c]) for u in range(U)]
        nz = any(sum(yc[kk * R + r] * v for r, v in bv.items()) % p for kk in range(k) for bv in bmod.values())
        if nz:
            fc = [u for u in free if int(X[u, c]) != 0]
            if len(fc) == 1: good.append(fc[0])
    print(f'vectors with nonzero cubic mod p: {len(good)} of {nullity}', flush=True)
    # pivot rows via rref of transpose
    flatT = [0] * (U * E)
    for i, eq in enumerate(eqs):
        for u, v in eq.items(): flatT[u * E + i] = v % p
    MTr, rankT = nmod_mat(U, E, flatT, p).rref(); assert rankT == rank
    prow = []
    for i in range(rankT):
        j = 0
        while j < E and int(MTr[i, j]) == 0: j += 1
        prow.append(j)
    pivpos = {u: i for i, u in enumerate(piv)}
    Esq = fmpz_mat(rank, rank)
    for i, ri in enumerate(prow):
        for u, v in eqs[ri].items():
            if u in pivpos: Esq[i, pivpos[u]] = v
    sel = good if len(good) <= a.nvec else [good[(i * len(good)) // a.nvec] for i in range(a.nvec)]
    vars_ = [f'a{i}' for i in range(k)]
    cubics = []; vinfo = []
    os.makedirs(f'{w2lib.PACKET}/polyfarkas', exist_ok=True)
    ckpt = f'{w2lib.PACKET}/polyfarkas/{a.tag}.cubics_Q.jsonl'
    done_fc = set()
    if os.path.exists(ckpt):
        for line in open(ckpt):
            try:
                e = json.loads(line); cubics.append(e['cubic']); vinfo.append(e['info']); done_fc.add(e['info']['free_coord'])
            except Exception: pass
        print(f'checkpoint: {len(cubics)} lifted cubics loaded', flush=True)
    def singular_W(cubics):
        sing = f'{w2lib.PACKET}/polyfarkas/{a.tag}_Q.sing'
        with open(sing, 'w') as f:
            f.write('LIB "elim.lib";\nring r = 0, (a0,a1,a2,a3,a4,a5), dp;\n')
            f.write('ideal V0 = (3*a0+a1+a3+a4)^2 - 6*(3*a0^2+a1^2+a3^2+a4^2), (3*a1+a2+a4+a5)^2 - 6*(3*a1^2+a2^2+a4^2+a5^2), (3*a2+a0+a5+a3)^2 - 6*(3*a2^2+a0^2+a5^2+a3^2);\n')
            f.write('ideal C = ' + (',\n'.join(cubics) if cubics else '0') + ';\n')
            f.write('ideal Conly = std(C); "DIM_C_only", dim(Conly); if (dim(Conly) <= 1) { "VDIM_C_only", vdim(Conly); }\n')
            f.write('ideal W = std(V0 + C); "DIMW", dim(W); if (dim(W) <= 1) { "VDIMW", vdim(W); }\n')
            f.write('ideal Eu = a0-a5,a1-a5,a2-a5,a3-a5,a4-a5; "EULER_IN_W (0 = yes)", size(reduce(W, std(Eu)));\n')
            f.write('def WS = sat(W, Eu); ideal Wsat; if (typeof(WS) == "list") { Wsat = WS[1]; } else { Wsat = WS; }\nWsat = std(Wsat); "DIM_W_minus_Euler", dim(Wsat); if (dim(Wsat) <= 1) { "VDIM_W_minus_Euler", vdim(Wsat); }\n')
            f.write('def CS = sat(Conly, Eu); ideal Csat; if (typeof(CS) == "list") { Csat = CS[1]; } else { Csat = CS; }\nCsat = std(Csat); "DIM_C_minus_Euler", dim(Csat);\nquit;\n')
        r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=3600)
        g = lambda pat: (lambda m: int(m.group(1)) if m else None)(re.search(pat, r.stdout))
        return r.stdout, g
    if cubics:
        out, g = singular_W(cubics)
        print('checkpoint Singular:', out.replace(chr(10), ' | ')[-300:], flush=True)
        if g(r'DIM_W_minus_Euler (-?\d+)') in (-1, 0):
            sel = []   # already decided
    for t, fcoord in enumerate(sel):
        if fcoord in done_fc: continue
        rhs = fmpz_mat(rank, 1)
        for i, ri in enumerate(prow):
            v = eqs[ri].get(fcoord, 0)
            if v: rhs[i, 0] = -v
        ysol = Esq.solve(rhs)
        y = [Fraction(0)] * U; y[fcoord] = Fraction(1)
        for i, u in enumerate(piv): y[u] = Fraction(int(ysol[i, 0].p), int(ysol[i, 0].q))
        bad = 0
        for eq in eqs:
            s_ = Fraction(0)
            for u, v in eq.items():
                if y[u]: s_ += y[u] * v
            if s_ != 0: bad += 1
        maxbits = max(max(v.numerator.bit_length(), v.denominator.bit_length()) for v in y)
        if bad:
            print(f'vector {t}: {bad} bad equations -- skipped', flush=True); continue
        coef = {}
        for kk in range(k):
            for (j, l), bv in bvec.items():
                s_ = Fraction(0)
                for r, v in bv.items():
                    if y[kk * R + r]: s_ += y[kk * R + r] * v
                if s_ == 0: continue
                mult = 1 if j == l else 2
                key = tuple(sorted((kk, j, l)))
                coef[key] = coef.get(key, Fraction(0)) + mult * s_
        coef = {kk: v for kk, v in coef.items() if v}
        if not coef: continue
        den = 1
        for v in coef.values(): den = den * v.denominator // math.gcd(den, v.denominator)
        cub = '+'.join(f'({int(v * den)})*' + '*'.join(vars_[i] for i in kk) for kk, v in sorted(coef.items()))
        info = {'free_coord': fcoord, 'maxbits': maxbits, 'sha256_y': hashlib.sha256('\n'.join(str(v) for v in y).encode()).hexdigest()}
        cubics.append(cub); vinfo.append(info)
        with open(ckpt, 'a') as f: f.write(json.dumps({'cubic': cub, 'info': info}) + '\n')
        print(f'vector {t}: verified exactly (max bits {maxbits}), {len(cubics)} cubics  [{time.time()-t0:.1f}s]', flush=True)
        if len(cubics) % 5 == 0:
            out, g = singular_W(cubics)
            print('interim Singular:', out.replace(chr(10), ' | ')[-200:], flush=True)
            if g(r'DIM_W_minus_Euler (-?\d+)') in (-1, 0):
                print('EARLY STOP: W minus Euler already empty over Q', flush=True); break
    with open(f'{w2lib.PACKET}/polyfarkas/{a.tag}.cubics_Q.txt', 'w') as f: f.write('\n'.join(cubics) + '\n')
    out, g = singular_W(cubics)
    print(out[-1500:], flush=True)
    res = {'exp': 'redsys-polyfarkas-deg1-EXACT', 'family': 'D(a) = sum_k a_k z_k F_k (trivial sector, deg 6)', 'tag': a.tag, 'd': 6, 'prime': p, 'rows': R,
           'cols': S['nK'] + S['nQ'], 'equations': E, 'unknowns': U, 'rank_mod_p': rank, 'nullity_mod_p': nullity, 'vectors_lifted_verified': len(cubics),
           'maxbits': max((v['maxbits'] for v in vinfo), default=None), 'sha256_vectors': [v['sha256_y'] for v in vinfo],
           'dimW_affine_Q': g(r'DIMW (\d+)'), 'vdimW_Q': g(r'VDIMW (-?\d+)'), 'dimW_minus_Euler_affine_Q': g(r'DIM_W_minus_Euler (-?\d+)'),
           'vdimW_minus_Euler_Q': g(r'VDIM_W_minus_Euler (-?\d+)'), 'dim_cubics_only_Q': g(r'DIM_C_only (\d+)'), 'dim_cubics_only_minus_Euler_Q': g(r'DIM_C_minus_Euler (-?\d+)'),
           'seconds': round(time.time() - t0, 1), 'job': os.environ.get('W2JOB', 'redsys_lift')}
    res['result'] = (f"EXACT degree-1 polynomial Farkas on the reduced system: {len(cubics)} verified rational vectors; over Q: W = V0 cap V(phi) affine dim {res['dimW_affine_Q']}, "
                     f"minus Euler point {res['dimW_minus_Euler_affine_Q']}; V(phi) alone: dim {res['dim_cubics_only_Q']}, minus Euler {res['dim_cubics_only_minus_Euler_Q']}"
                     + (' => trivial-sector family EXCLUDED over C except D = 6F (illegal): RIGOROUS' if res['dimW_minus_Euler_affine_Q'] in (-1, 0) and res['dimW_minus_Euler_affine_Q'] is not None else ''))
    json.dump(res, open(f'{w2lib.PACKET}/polyfarkas/{a.tag}.receipt.json', 'w'), indent=1)
    w2lib.log(**res)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main()
