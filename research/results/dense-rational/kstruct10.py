"""Structure of K_full(d) = {D : A D in J} (mod p): dimension of K_full ∩ g*S_{d-deg g} for small forms g,
   i.e. kernel elements divisible by g; compared with the trivial part J_d ∩ g S.  Usage: python3 kstruct10.py basisfile"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
data = json.load(open(sys.argv[1])); p = data['p']; d = data['d']; parD = data['parD']
basis = [{tuple(int(x) for x in k.split()): v for k, v in b.items()} for b in data['basis']]
F, A, q = load(); Fp = core.reduce_mod(F, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
Dmons = monomials(d, parD); idx = {m: j for j, m in enumerate(Dmons)}; nD = len(Dmons)
def mat(polys):
    M = nmod_mat(len(polys), nD, p)
    for r, P_ in enumerate(polys):
        for mm, v in P_.items(): M[r, idx[mm]] = v
    return M
def rank(polys): return mat(polys).rank() if polys else 0
Jgens = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(d - 5, (parD - pi) % 2): Jgens.append(shift(Fi[i], m, 1, p))
rK = rank(basis); rJ = rank(Jgens)
print(f'd={d} par={parD}: dim K_full={rK}, dim J_d={rJ}, kernel mod J={rK-rJ}', flush=True)
tests = {'X0': core.parse('X0'), 'Y0': core.parse('Y0'), 'X1': core.parse('X1'), 'Y1': core.parse('Y1'),
         'X0^2': core.parse('X0**2'), 'X0*Y0': core.parse('X0*Y0'), 'Y0^2': core.parse('Y0**2'),
         'X0^2+X1^2+X2^2': core.parse('X0**2+X1**2+X2**2'), 'X0+X1': core.parse('X0+X1'), 'X0+Y0': core.parse('X0+Y0')}
for nm, g in tests.items():
    gp = core.reduce_mod(g, p); e = core.degree(g); pg = core.parity_of(g)
    gS = [mul(gp, {m: 1}, p) for m in monomials(d - e, (parD - pg) % 2)]
    rg = rank(gS); rKg = rank(basis + gS); rJg = rank(Jgens + gS)
    inter_K = rK + rg - rKg; inter_J = rJ + rg - rJg
    print(json.dumps({'g': nm, 'dim g*S': rg, 'dim K∩gS': inter_K, 'dim J∩gS': inter_J, 'new (mod J) divisible by g': inter_K - inter_J}), flush=True)
