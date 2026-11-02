"""EXACT (characteristic-0) relaxation of one degree-6 character sector: E1 kernel over Q, exact N0, exact H(D_a,D_b),
exact normal forms modulo J_18 (fmpq_mat rref of the sector Jacobian generators), quadrics over Q, Singular dim over Q.
A result 'affine dim 0' here excludes the whole sector family over Q-bar WITHOUT the reduction argument T2.
usage: python3 wave2-excl/sector_exact.py --sector 0,0,1
"""
import sys, json, time, argparse, subprocess, re, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hess6 import key_of, poly_key
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials
from fractions import Fraction
from flint import fmpq_mat, fmpq

def rref_rows(polys, mons):
    """exact rref of the polynomials (rows) over the monomial basis; returns (pivots, rows as dicts over mons, rank)."""
    idx = {m: i for i, m in enumerate(mons)}
    R, C = len(polys), len(mons)
    M = fmpq_mat(R, C)
    for r, P in enumerate(polys):
        for e, v in P.items():
            fr = Fraction(v); M[r, idx[e]] = fmpq(fr.numerator, fr.denominator)
    Mr, rank = M.rref()
    piv = []; rows = []
    for i in range(rank):
        j = 0
        while j < C and Mr[i, j] == 0: j += 1
        piv.append(j)
        rows.append({mons[k]: Fraction(int(Mr[i, k].p), int(Mr[i, k].q)) for k in range(j, C) if Mr[i, k] != 0})
    return piv, rows, rank

class ExactReducer:
    def __init__(self, Fi, deg, key, gens):
        self.mons = [m for m in monomials(deg) if key_of(m, gens) == key]
        gl = []
        for k in range(6):
            kk = key_of(next(iter(Fi[k])), gens)
            need = tuple((a - b) % 6 for a, b in zip(key, kk))
            for m in monomials(deg - 5):
                if key_of(m, gens) == need: gl.append(shift(Fi[k], m, 1))
        self.piv, self.rows, self.rank = rref_rows(gl, self.mons) if gl else ([], [], 0)
        ps = set(self.piv); self.standard = [m for j, m in enumerate(self.mons) if j not in ps]
    def normal_form(self, P):
        P = {e: Fraction(v) for e, v in P.items() if Fraction(v) != 0}
        for j, row in zip(self.piv, self.rows):
            m = self.mons[j]; c = P.get(m, 0)
            if c:
                for e, v in row.items():
                    nv = P.get(e, 0) - c * v
                    if nv: P[e] = nv
                    else: P.pop(e, None)
        return P

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sector', required=True)
    a = ap.parse_args(); kD = tuple(int(x) for x in a.sector.split(','))
    t0 = time.time()
    F, A, q, G = w2lib.load(); gens = G.generators
    F2 = scale(F, 2); Fi = [pderiv(F2, i) for i in range(6)]
    kA = poly_key(A, gens); d = 6
    Dmons = [m for m in monomials(d) if key_of(m, gens) == kD]
    kAD = tuple((x + y) % 6 for x, y in zip(kA, kD))
    R18 = ExactReducer(Fi, 18, kAD, gens)
    print(f'sector {kD}: {len(Dmons)} monomials; J_18 sector: rank {R18.rank}, standard {len(R18.standard)}  [{time.time()-t0:.1f}s]', flush=True)
    # E1 kernel over Q: D with [A D] = 0 in R_18
    imgs = [R18.normal_form(shift(A, m, 1)) for m in Dmons]
    std = R18.standard; sidx = {m: i for i, m in enumerate(std)}
    M = fmpq_mat(len(std), len(Dmons))
    for j, nf in enumerate(imgs):
        for e, v in nf.items(): M[sidx[e], j] = fmpq(v.numerator, v.denominator)
    # nullspace over Q via rref of M
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < len(Dmons) and Mr[i, j] == 0: j += 1
        piv.append(j)
    free = [j for j in range(len(Dmons)) if j not in set(piv)]
    kervecs = []
    for fj in free:
        v = {Dmons[fj]: Fraction(1)}
        for i, pj in enumerate(piv):
            c = Mr[i, fj]
            if c != 0: v[Dmons[pj]] = -Fraction(int(c.p), int(c.q))
        kervecs.append(v)
    k = len(kervecs)
    print(f'E1 kernel over Q: dim {k}  [{time.time()-t0:.1f}s]', flush=True)
    if k < 2:
        print('nothing to relax (dim < 2)'); return
    # exact N0: solve N(F2) = A D over Q (sector-restricted columns m F2_t with key = kAD)
    def solve_E1(AD):
        cols = []; colp = []
        for t in range(6):
            kk = key_of(next(iter(Fi[t])), gens)
            need = tuple((x - y) % 6 for x, y in zip(kAD, kk))
            for m in monomials(13):
                if key_of(m, gens) == need: cols.append((t, m)); colp.append(shift(Fi[t], m, 1))
        rows = {}
        for P in colp:
            for e in P: rows.setdefault(e, len(rows))
        for e in AD: rows.setdefault(e, len(rows))
        Rn, Cc = len(rows), len(cols)
        Mq = fmpq_mat(Rn, Cc + 1)
        for j, P in enumerate(colp):
            for e, v in P.items(): fr = Fraction(v); Mq[rows[e], j] = fmpq(fr.numerator, fr.denominator)
        for e, v in AD.items(): fr = Fraction(v); Mq[rows[e], Cc] = fmpq(fr.numerator, fr.denominator)
        Mqr, rk = Mq.rref()
        pv = []
        for i in range(rk):
            j = 0
            while j <= Cc and Mqr[i, j] == 0: j += 1
            pv.append(j)
        assert pv[-1] != Cc, 'E1 inconsistent over Q?!'
        N = [{} for _ in range(6)]
        for i, j in enumerate(pv):
            val = Mqr[i, Cc]
            if val != 0:
                t, m = cols[j]; N[t][m] = Fraction(int(val.p), int(val.q))
        return N
    Ns = [solve_E1(mul(A, v)) for v in kervecs]
    print(f'exact E1 solutions found  [{time.time()-t0:.1f}s]', flush=True)
    def Hpart(Da, Nb, Db, Na):
        divNb = {}; divNa = {}
        for t in range(6):
            divNb = add(divNb, pderiv(Nb[t], t)); divNa = add(divNa, pderiv(Na[t], t))
        out = add(mul(Da, divNb), mul(Db, divNa))
        for t in range(6):
            out = core.sub(out, mul(Nb[t], pderiv(Da, t))); out = core.sub(out, mul(Na[t], pderiv(Db, t)))
        return scale(out, Fraction(1, 2))
    redH = {}; quad = {}
    for a_ in range(k):
        for b_ in range(a_, k):
            Hab = Hpart(kervecs[a_], Ns[b_], kervecs[b_], Ns[a_])
            Hab = {e: Fraction(v) for e, v in Hab.items() if Fraction(v) != 0}
            if not Hab: continue
            kH = poly_key(Hab, gens)
            if kH not in redH: redH[kH] = ExactReducer(Fi, 18, kH, gens)
            nf = redH[kH].normal_form(Hab)
            for e, v in nf.items():
                quad.setdefault((kH, e), {})[(a_, b_)] = v if a_ == b_ else 2 * v
    vars_ = [f'c{i+1}' for i in range(k)]
    polys = []
    for terms in quad.values():
        den = 1
        import math
        for v in terms.values(): den = den * v.denominator // math.gcd(den, v.denominator)
        polys.append('+'.join(f'({int(v*den)})*{vars_[x]}*{vars_[y]}' if x != y else f'({int(v*den)})*{vars_[x]}^2' for (x, y), v in sorted(terms.items())))
    sing = f'{w2lib.PACKET}/sectors/exact_d6_k{"".join(map(str,kD))}.sing'
    with open(sing, 'w') as f:
        f.write(f'ring r = 0, ({",".join(vars_)}), dp;\nideal I = ' + (',\n'.join(polys) if polys else '0') + ';\nideal G = std(I); "DIM", dim(G); if (dim(G) <= 1) { "VDIM", vdim(G); }\nquit;\n')
    r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=600)
    dm = re.search(r'DIM (\d+)', r.stdout)
    res = {'exp': 'sector-d6-exact', 'd': 6, 'prime': None, 'sector': list(kD), 'dimK_Q': k, 'nquadrics_Q': len(polys), 'relax_dim_affine_Q': int(dm.group(1)) if dm else None,
           'seconds': round(time.time() - t0, 1), 'job': 'local',
           'result': ('sector family EXCLUDED over Q-bar by EXACT relaxation (affine dim 0, no reduction argument)' if dm and dm.group(1) == '0' else f'exact relaxation affine dim {dm.group(1) if dm else None}')}
    w2lib.log(**res); print(json.dumps(res))

if __name__ == '__main__':
    main()
