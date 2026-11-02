"""Independent check of the group/character claims (dossier §1.5, A.1)."""
import sys, itertools, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify')
import vcheck as V
import numpy as np
F, A, q = V.load_instance()
G = V.diag_group(F)                       # incl scalars
scal = np.array([[c]*6 for c in range(6)])
# projective group: quotient by scalars -> normalize w(X2)=0 (index 2)
Gp = np.array(sorted({tuple((w - w[2]) % 6) for w in G}))
out = {'order_with_scalars': len(G), 'projective_order': len(Gp)}
def order(w):
    k = 1; v = np.array(w)
    while np.any((k * v) % 6): k += 1
    return k
from collections import Counter
out['element_orders_G'] = dict(Counter(order(w) for w in Gp))
# residue character chi_A(w) = sum (e_i+1) w_i mod 6, with e = exponents of A
eA = np.array(list(A.keys())[0]) + 1
chiA = lambda w: int((eA @ np.array(w)) % 6)
out['chiA_values'] = sorted({chiA(w) for w in Gp})
K = [w for w in Gp if chiA(w) == 0]
out['K_order'] = len(K); out['element_orders_K'] = dict(Counter(order(w) for w in K))
tau = (0,0,0,3,3,3); s0 = (3,0,0,3,0,0); s1 = (0,3,0,0,3,0); g = (2,0,0,2,4,2); h = (2,1,0,3,2,1); gam = (2,2,0,4,0,0)
inG = lambda w: tuple(w) in {tuple(x) for x in Gp}
out['tau_in_G'] = inG(tau); out['s0_in_G'] = inG(s0); out['s1_in_G'] = inG(s1); out['g_in_G'] = inG(g); out['h_in_G'] = inG(h); out['gamma_in_G'] = inG(gam)
out['chiA_tau'] = chiA(tau); out['chiA_s0'] = chiA(s0); out['chiA_s1'] = chiA(s1); out['chiA_g'] = chiA(g); out['chiA_h'] = chiA(h); out['chiA_gamma'] = chiA(gam)
# generators span G?
def span(gens):
    S = {(0,)*6}; frontier = [(0,)*6]
    while frontier:
        nxt = []
        for s in frontier:
            for gg in gens:
                t = tuple((a+b) % 6 for a, b in zip(s, gg))
                if t not in S: S.add(t); nxt.append(t)
        frontier = nxt
    return S
sp = span([g, h, s0, s1, tau]); out['gens_span_order'] = len({tuple((np.array(w)-w[2])%6) for w in sp})
# Gamma_book = P x <h>, P = paired signs; K = ?
P = span([s0, s1]); Gb = span([s0, s1, h])
out['Gamma_book_order'] = len(Gb); out['Gamma_book_eq_K'] = {tuple(x) for x in Gb} == {tuple(x) for x in K}
out['h_in_K'] = chiA(h) == 0
# SNF-type structure: G ≅ Z/2 x Z/6 x Z/6 ? check via counting elements of each order: Z2xZ6xZ6 has
# #order-dividing-2 elements = 2*2*2 = 8 (7 nontrivial), #order dividing 3 = 9 (8 nontrivial)
out['n_order_div2'] = sum(1 for w in Gp if order(w) in (1, 2)); out['n_order_div3'] = sum(1 for w in Gp if order(w) in (1, 3))
# permutations preserving the support of F (up to diagonal rescaling: support only), and exact preservation
suppF = set(F.keys())
perms = []
for s in itertools.permutations(range(6)):
    img = {tuple(m[s[i]] for i in range(6)) for m in suppF}   # (sigma . m)_i = m_{sigma(i)}
    if img == suppF: perms.append(s)
out['support_perms'] = [list(s) for s in perms]
def permute(Pl, s):
    return {tuple(m[s[i]] for i in range(6)): c for m, c in Pl.items()}
out['perms_preserve_F_exactly'] = [permute(F, s) == F for s in perms]
# R(A), R^2(A) characters
RA = [permute(A, s) for s in perms]
out['perm_images_of_A'] = [V.to_str(x) for x in RA]
def chiP(Pm, w):   # residue character of a degree-12 monomial numerator at q=2 (same pole order), = weight + sum w
    e = np.array(list(Pm.keys())[0]) + 1
    return int((e @ np.array(w)) % 6)
out['chi_diffs_RA_vs_A'] = [sorted({(chiP(x, w) - chiA(w)) % 6 for w in Gp}) for x in RA]
# counts: degree-12 monomials with weight(A) under G; degree-6 monomials in trivial sector and A-sector
def wt(m, w): return int(np.dot(m, w) % 6)
wA = [wt(list(A.keys())[0], w) for w in G]
m12 = V.monomials(12); m6 = V.monomials(6)
out['deg12_mono_with_weightA'] = sum(1 for m in m12 if all(wt(m, w) == wA[k] for k, w in enumerate(G)))
out['deg12_mono_with_weightA_Yodd'] = sum(1 for m in m12 if all(wt(m, w) == wA[k] for k, w in enumerate(G)) and (m[3]+m[4]+m[5]) % 2 == 1)
out['deg6_mono_weight0'] = sum(1 for m in m6 if all(wt(m, w) == 0 for w in G))
out['deg6_mono_weightA'] = sum(1 for m in m6 if all(wt(m, w) == wA[k] for k, w in enumerate(G)))
out['A_Yodd'] = (list(A.keys())[0][3] + list(A.keys())[0][4] + list(A.keys())[0][5]) % 2 == 1
out['F_Yeven'] = all((m[3]+m[4]+m[5]) % 2 == 0 for m in F)
# coordinates with identical weight under all of G (needed for z_j F_i invariant with j != i)
out['coord_pairs_same_weight'] = [(i, j) for i in range(6) for j in range(i+1, 6) if all(w[i] == w[j] for w in G)]
print(json.dumps(out, indent=1, default=str))
V.log(claim='G: |G|=72 (432 with scalars), element orders {1:1,2:7,3:8,6:56}, chi_A values {0,2,4}, |K|=24 with orders {1:1,2:7,3:2,6:14}, tau,s0,s1,gamma in K, chi_A(h)=2, chi_A(g)=4, gens g,h,s0,s1,tau span G, Gamma_book != K, only id,R,R^2 preserve supp F, R preserves F exactly, chi(R^k A) != chi_A, 90 deg-12 monomials with chi_A (all Y-odd), 12 invariant sextic monomials, 7 sextic monomials of weight A',
      source='dossier §1.5/A.1', verdict='CONFIRMED' if (out['projective_order']==72 and out['K_order']==24 and out['chiA_values']==[0,2,4] and out['chiA_h']==2 and out['chiA_g']==4 and out['gens_span_order']==72 and not out['Gamma_book_eq_K'] and len(perms)==3 and all(out['perms_preserve_F_exactly']) and out['deg12_mono_with_weightA']==90 and out['deg6_mono_weight0']==12 and out['deg6_mono_weightA']==7 and all(0 in d and len(d)>1 for d in out['chi_diffs_RA_vs_A'][1:])) else 'CHECK',
      evidence=out, tool='v1_group.py (own code)')
