#!/usr/bin/env python3
"""Structure of the E1 kernel K_full(d) = {D in S_d^{parD} : A D in J} (mod p) from a saved basis:
 dim of K_full ∩ g*S_{d-deg g} versus J_d ∩ g*S for test forms g (kernel elements divisible by g), including g = A
 (is A*u in (J:A), i.e. A^2 u in J?), and the dimension of K_full ∩ (D's with at most t monomials)? (not meaningful)
Usage: python3 wave2_excl_kstruct.py basis.json"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
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
    for m in monomials(d - 5, (parD - par_i(i)) % 2): Jgens.append(shift(Fi[i], m, 1, p))
rK = rank(basis); rJ = rank(Jgens)
out = {'exp': 'e1-kernel-structure', 'd': d, 'parD': parD, 'prime': p, 'dim_Kfull': rK, 'dim_J': rJ, 'kernel_mod_J': rK - rJ, 'tests': []}
print(json.dumps({k: v for k, v in out.items() if k != 'tests'}), flush=True)
tests = {'X0': 'X0', 'Y0': 'Y0', 'X1': 'X1', 'Y1': 'Y1', 'X2': 'X2', 'Y2': 'Y2', 'X0+X1': 'X0+X1', 'X0+Y0': 'X0+Y0', 'X0+Y1': 'X0+Y1',
         'X0^2': 'X0**2', 'X0*Y0': 'X0*Y0', 'Y0^2': 'Y0**2', 'X0^2+X1^2+X2^2': 'X0**2+X1**2+X2**2', 'Y0^2+Y1^2+Y2^2': 'Y0**2+Y1**2+Y2**2',
         'X0^2+Y0^2': 'X0**2+Y0**2', 'X0^3': 'X0**3', 'Y0^3': 'Y0**3', 'X0^3*Y0^3=A': 'X0**3*Y0**3', 'X0^2*Y0^2': 'X0**2*Y0**2',
         'X0^3+X1^3+X2^3': 'X0**3+X1**3+X2**3', 'X0^4': 'X0**4', 'X0^6': 'X0**6', 'Y0^6': 'Y0**6', 'X0^5': 'X0**5', 'X0^4*Y0^2': 'X0**4*Y0**2'}
for nm, ex in tests.items():
    g = core.parse(ex); gp = core.reduce_mod(g, p); e = core.degree(g); pg = core.parity_of(g)
    if e > d: continue
    gS = [mul(gp, {m: 1}, p) for m in monomials(d - e, (parD - pg) % 2)]
    if not gS: continue
    rg = rank(gS); rKg = rank(basis + gS); rJg = rank(Jgens + gS)
    rec = {'g': nm, 'deg_g': e, 'dim_gS': rg, 'dim_K_cap_gS': rK + rg - rKg, 'dim_J_cap_gS': rJ + rg - rJg, 'new_mod_J_divisible_by_g': (rK + rg - rKg) - (rJ + rg - rJg)}
    out['tests'].append(rec); print(json.dumps(rec), flush=True)
log2(**{k: v for k, v in out.items() if k != 'tests'}, tests=out['tests'])
with open(sys.argv[1].replace('_basis.json', '_kstruct.json'), 'w') as f: json.dump(out, f, indent=1)
