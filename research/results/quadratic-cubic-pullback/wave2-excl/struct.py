"""Structured subspaces W of S_d^{par} intersected with the E1 kernel K_d (mod p), compared with J_d:
   reports dim W, dim(W ∩ K), dim(W ∩ J) and the number of 'new' E1-admissible directions dim(W∩K) - dim(W∩J).
   Families (d = 10, Y-even unless noted): quintics in Q (252); delta * S_4^odd (60); P(Q_S) deg-5 for 20 triples;
   products H(Q_S)*H'(Q_S') (deg 3 x deg 2 in complementary triples); delta*Q_i*(odd quadric); (even quartic)*(Q-cubic);
   Fermat-type blocks X0^10+..., sums of 10th powers; F_part products; tau-invariant 'shape symmetric' sums.
Usage: python3 struct.py basisfile.json"""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *
data = json.load(open(sys.argv[1])); p = data['p']; d = data['d']; parD = data['parD']
basis = [{tuple(int(x) for x in k.split()): v for k, v in b.items()} for b in data['basis']]
F, A, q = load(); Fp = core.reduce_mod(F, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
Q = [core.reduce_mod(q_, p) for q_ in quadrics()]
dl = delta(p)
Dmons = monomials(d, parD); idx = {m: j for j, m in enumerate(Dmons)}; nD = len(Dmons)
def mat(polys):
    return mat_from_polys(polys, Dmons, p)
def rank(polys): return mat(polys).rank() if polys else 0
Jgens = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d - 5, (parD - pi) % 2): Jgens.append(shift(Fi[i], m, 1, p))
rK = rank(basis); rJ = rank(Jgens)
print(f'd={d} par={parD} p={p}: dim K={rK}, dim J_d={rJ}, kernel mod J={rK-rJ}', flush=True)
def polyQ(exps):
    out = {(0,)*NV: 1}
    for i, e in enumerate(exps):
        for _ in range(e): out = mul(out, Q[i], p)
    return out
fams = {}
# quintics in Q
fams['C[Q]_5'] = [polyQ(e) for e in monomials(5)]  # exps over 6 'variables' = Q indices
# delta * S_4^odd
fams['delta*S4odd'] = [mul(dl, {m: 1}, p) for m in monomials(4, 1)]
# cubics in Q times quadrics in Q is inside C[Q]_5; triples:
for Sset in itertools.combinations(range(6), 3):
    exps5 = [e for e in monomials(5) if all(e[i] == 0 for i in range(6) if i not in Sset)]
    fams[f'C[Q_{"".join(map(str,Sset))}]_5'] = [polyQ(e) for e in exps5]
# delta * Q_i * odd quadric  (subset of delta*S4odd) skip; (even quartic in Z) * (cubic in Q): big (66*56); sample: quartic monomials X0^4.. times Q cubics
# sums of 10th powers / Fermat-like blocks
fams['Z^10 monomials even'] = [{m: 1} for m in monomials(10, 0) if max(m) == 10 or sorted(m, reverse=True)[1] == 0]
# F-part products: F_{60}*F_{40}? degree 10 = products of bigraded pieces? (Fpart(6)=F_60 deg 6) * quartics: skip
# products of two polars F_i F_j * quadric -> in J; skip
# delta^2 is degree 12; skip. delta * (Q-quadric) = delta*S4 subset covered.
# E_X(F)-type: Euler images: E_I(F) * S_4 for I = {0,1,2}: subset of J? E_I(F) = sum_{i in I} Z_i F_i in J. skip.
res = {'d': d, 'parD': parD, 'p': p, 'dimK': rK, 'dimJ': rJ, 'families': {}}
for nm, polys in fams.items():
    polys = [P_ for P_ in polys if P_ and core.parity_of(P_) == parD]
    if not polys: continue
    rW = rank(polys); rKW = rank(basis + polys); rJW = rank(Jgens + polys)
    inter_K = rK + rW - rKW; inter_J = rJ + rW - rJW
    r = {'dimW': rW, 'dim W∩K': inter_K, 'dim W∩J': inter_J, 'new': inter_K - inter_J}
    res['families'][nm] = r
    print(json.dumps({'family': nm, **r}), flush=True)
    log({'exp': 'struct_intersection', 'd': d, 'parD': parD, 'p': p, 'family': nm, **r})
receipt(f'struct_d{d}_par{parD}_p{p}.json', res)
