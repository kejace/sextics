"""E1: independent re-verification of the symmetry/character facts of dossier §1.5 (A.1)."""
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
Gw = diagonal_group(F2)                       # includes scalars (k,..,k)
n_all = len(Gw)
# projective group: quotient by scalars (6 of them)
proj = {tuple(int(x) for x in (w - w[5]) % 6) for w in Gw}     # normalize w(Y2)=0
nG = len(proj)
# residue character chi_A(w) = sum (e_i+1) w_i mod 6, with e the exponent of A (the volume weight sum w_i included)
eA = list(A.keys())[0]
def chi(w): return int(sum((eA[i] + 1) * int(w[i]) for i in range(6)) % 6)
vals = sorted({int(chi(w)) for w in proj})
K = [w for w in proj if chi(w) == 0]
tau = (0, 0, 0, 3, 3, 3); g = (2, 0, 0, 2, 4, 2); h = (2, 1, 0, 3, 2, 1); gam = (2, 2, 0, 4, 0, 0)
s0 = (3, 0, 0, 3, 0, 0); s1 = (0, 3, 0, 0, 3, 0)
def inG(w): return tuple((w[i] - w[5]) % 6 for i in range(6)) in proj
def order(w):
    k = 1; x = w
    while any(((x[i] - x[5]) % 6) for i in range(6)):
        x = tuple((a + b) % 6 for a, b in zip(x, w)); k += 1
    return k
orders = {}
for w in proj:
    orders[order(w)] = orders.get(order(w), 0) + 1
ordersK = {}
for w in K:
    ordersK[order(w)] = ordersK.get(order(w), 0) + 1
# permutations preserving the support of F up to diagonal rescaling: support must map to support
supp = set(F2.keys())
perms = []
for pi in itertools.permutations(range(6)):
    img = {tuple(e[pi[i]] for i in range(6)) for e in supp}
    if img == supp:
        perms.append(pi)
# is R: (X_i,Y_i)->(X_{i+1},Y_{i+1}) exact on F2 (coefficients)?
R = (1, 2, 0, 4, 5, 3)   # new index of variable i: X0->X1 etc.
FR = {tuple(e[R.index(i)] for i in range(6)): c for e, c in F2.items()}
# R(A) characters
def permpoly(P, pi):
    return {tuple(e[pi.index(i)] for i in range(6)): c for e, c in P.items()}
RA = permpoly(A, R); R2A = permpoly(RA, R)
def chiP(P, w):
    e = list(P.keys())[0]; return int(sum((e[i] + 1) * int(w[i]) for i in range(6)) % 6)
diff1 = sorted({(chiP(RA, w) - chi(w)) % 6 for w in proj}); diff2 = sorted({(chiP(R2A, w) - chi(w)) % 6 for w in proj})
rec = dict(exp='E1-group-characters', engine=ENGINE_TAG, group_with_scalars=int(n_all), projective_order=nG, chiA_values=vals, kernel_order=len(K),
           tau_in_K=(chi(tau) == 0), chi_h=chi(h), chi_g=chi(g), gamma_in_K=(chi(gam) == 0), s0_in_K=chi(s0) == 0, s1_in_K=chi(s1) == 0,
           tau_in_G=inG(tau), g_in_G=inG(g), h_in_G=inG(h), gamma_in_G=inG(gam), element_orders_G=orders, element_orders_K=ordersK,
           support_perms=len(perms), R_preserves_F_exactly=(FR == F2), chi_RA_minus_chiA_values=diff1, chi_R2A_minus_chiA_values=diff2,
           deg12_monomials_with_chiA=None)
# count degree-12 monomials with the same G-character as A (full group with scalars: same key)
gens = generators(Gw); kA = key_of_poly(A, gens)
rec['deg12_monomials_with_chiA'] = len(sector_monos(12, gens, kA))
rec['deg12_Yodd_monomials_with_chiA'] = len([m for m in sector_monos(12, gens, kA) if (m[3] + m[4] + m[5]) % 2 == 1])
print(json.dumps(rec, indent=1))
xlib.log(rec)
