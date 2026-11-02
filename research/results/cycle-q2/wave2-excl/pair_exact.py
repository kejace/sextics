"""EXACT relaxation (over Q) for the union of the trivial sector with one nontrivial sector at degree 6 (spot check of subsets6).
usage: python3 wave2-excl/pair_exact.py --sector 0,0,1"""
import sys, json, time, argparse, subprocess, re, math
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, sector_exact
from sector_exact import ExactReducer, rref_rows
from hess6 import key_of, poly_key
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials
from fractions import Fraction
from flint import fmpq_mat, fmpq

ap = argparse.ArgumentParser(); ap.add_argument('--sector', required=True); a = ap.parse_args()
kD2 = tuple(int(x) for x in a.sector.split(','))
t0 = time.time()
F, A, q, G = w2lib.load(); gens = G.generators
F2 = scale(F, 2); Fi = [pderiv(F2, i) for i in range(6)]; kA = poly_key(A, gens)
Z = [core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
# trivial sector basis: D_i = z_i F2_i with N0_i = A z_i e_i (exact E1 solutions for F2)
fam = [(mul(Z[i], Fi[i]), [mul(A, Z[i]) if t == i else {} for t in range(6)]) for i in range(6)]
# nontrivial sector: exact E1 kernel and E1 solutions (reuse sector_exact machinery)
Dmons = [m for m in monomials(6) if key_of(m, gens) == kD2]
kAD = tuple((x + y) % 6 for x, y in zip(kA, kD2))
R18 = ExactReducer(Fi, 18, kAD, gens)
imgs = [R18.normal_form(shift(A, m, 1)) for m in Dmons]
std = R18.standard; sidx = {m: i for i, m in enumerate(std)}
M = fmpq_mat(len(std), len(Dmons))
for j, nf in enumerate(imgs):
    for e, v in nf.items(): M[sidx[e], j] = fmpq(v.numerator, v.denominator)
Mr, rk = M.rref(); piv = []
for i in range(rk):
    j = 0
    while j < len(Dmons) and Mr[i, j] == 0: j += 1
    piv.append(j)
free = [j for j in range(len(Dmons)) if j not in set(piv)]
for fj in free:
    v = {Dmons[fj]: Fraction(1)}
    for i, pj in enumerate(piv):
        c = Mr[i, fj]
        if c != 0: v[Dmons[pj]] = -Fraction(int(c.p), int(c.q))
    # exact E1 solve
    AD = mul(A, v); cols = []; colp = []
    for t in range(6):
        kk = key_of(next(iter(Fi[t])), gens); need = tuple((x - y) % 6 for x, y in zip(kAD, kk))
        for m in monomials(13):
            if key_of(m, gens) == need: cols.append((t, m)); colp.append(shift(Fi[t], m, 1))
    rows = {}
    for P in colp:
        for e in P: rows.setdefault(e, len(rows))
    for e in AD: rows.setdefault(e, len(rows))
    Rn, Cc = len(rows), len(cols); Mq = fmpq_mat(Rn, Cc + 1)
    for j, P in enumerate(colp):
        for e, val in P.items(): fr = Fraction(val); Mq[rows[e], j] = fmpq(fr.numerator, fr.denominator)
    for e, val in AD.items(): fr = Fraction(val); Mq[rows[e], Cc] = fmpq(fr.numerator, fr.denominator)
    Mqr, rk2 = Mq.rref(); pv = []
    for i in range(rk2):
        j = 0
        while j <= Cc and Mqr[i, j] == 0: j += 1
        pv.append(j)
    assert pv[-1] != Cc
    N = [{} for _ in range(6)]
    for i, j in enumerate(pv):
        val = Mqr[i, Cc]
        if val != 0: t, m = cols[j]; N[t][m] = Fraction(int(val.p), int(val.q))
    fam.append((v, N))
k = len(fam)
print(f'union family: {k} = 6 (trivial) + {k-6} (sector {kD2})  [{time.time()-t0:.1f}s]', flush=True)
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
        Hab = Hpart(fam[a_][0], fam[b_][1], fam[b_][0], fam[a_][1])
        Hab = {e: Fraction(v) for e, v in Hab.items() if Fraction(v) != 0}
        if not Hab: continue
        kH = poly_key(Hab, gens)
        if kH not in redH: redH[kH] = ExactReducer(Fi, 18, kH, gens)
        nf = redH[kH].normal_form(Hab)
        for e, v in nf.items(): quad.setdefault((kH, e), {})[(a_, b_)] = v if a_ == b_ else 2 * v
vars_ = [f'c{i+1}' for i in range(k)]
polys = []
for terms in quad.values():
    den = 1
    for v in terms.values(): den = den * v.denominator // math.gcd(den, v.denominator)
    polys.append('+'.join(f'({int(v*den)})*{vars_[x]}*{vars_[y]}' if x != y else f'({int(v*den)})*{vars_[x]}^2' for (x, y), v in sorted(terms.items())))
sing = f'{w2lib.PACKET}/sectors/exact_pair_k{"".join(map(str,kD2))}.sing'
nontriv = '*'.join(vars_[6:])
open(sing, 'w').write(f'LIB "elim.lib";\nring r = 0, ({",".join(vars_)}), dp;\nideal I = ' + ',\n'.join(polys) + f';\nideal G = std(I); "DIM", dim(G);\npoly f = {nontriv};\ndef S = sat(I, f); ideal Sat; if (typeof(S) == "list") {{ Sat = S[1]; }} else {{ Sat = S; }}\n"SATDIM", dim(std(Sat));\nquit;\n')
r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=600)
dm = re.search(r'DIM (\d+)', r.stdout); sd = re.search(r'SATDIM (-?\d+)', r.stdout)
res = {'exp': 'pair-d6-exact', 'd': 6, 'prime': None, 'sectors': [[0,0,0], list(kD2)], 'dim_union': k, 'nquadrics_Q': len(polys), 'dim_affine_Q': int(dm.group(1)) if dm else None,
       'dim_off_V0_Q': int(sd.group(1)) if sd else None, 'seconds': round(time.time()-t0, 1), 'job': 'local',
       'result': f'exact over Q: union of trivial sector and {kD2}: relaxation variety affine dim {dm.group(1) if dm else None}; points with nontrivial coordinates nonzero: affine dim {sd.group(1) if sd else None}' + (' (none => contained in V0, EXACT)' if sd and sd.group(1) == '-1' else '')}
w2lib.log(**res); print(json.dumps(res))
