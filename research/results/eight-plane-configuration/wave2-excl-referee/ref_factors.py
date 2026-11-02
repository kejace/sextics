"""Sum-principle bookkeeping replay at a NEW prime with the no-drop check recorded:
Ann_R(D)_6^odd = {c in S_6^- : c*D in J_{6+e}} / J_6^-.  dim_p = nC - (rank_p[D m | J_{6+e}] - rank_p J_{6+e}) - 18.
Rigorous over Q when rank_p J_{6+e}^{par} = exact Koszul value (same projected-kernel argument as E1)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl-referee')
import ref_lib as L
from flint import nmod_mat
from fractions import Fraction
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000009
F, A, q = L.load(); Fp = L.red(F, p)
Fi = [L.pderiv(Fp, i, p) for i in range(L.NV)]
KD = json.load(open('/tmp/hodge/results/eight-plane-configuration/K_decomp.json'))
K = [L.parse_poly(s) for s in KD['K']]; H = [L.parse_poly(s) for s in KD['H']]
def P(s): return L.parse_poly(s)
Q = [P(f'X{i}**2+Y{i}**2') for i in range(3)]
# exact check of the decomposition F = sum Q_i K_i and K_i = diag_i + (5/p_i) H_i
S = {}
for i in range(3): S = L.add(S, L.mul(Q[i], K[i]))
okF = (L.add(S, F, s=-1) == {})
pr = [7, 11, 13]; okK = all(L.add(K[i], L.add(P(f'X{i}**4-X{i}**2*Y{i}**2+Y{i}**4'), {e: Fraction(5, pr[i]) * c for e, c in H[i].items()}), s=-1) == {} for i in range(3))
L.log({'check': 'exact: F = sum Q_i K_i and K_i = diag_i + (5/p_i) H_i from K_decomp.json', 'F_eq_sumQK': okF, 'K_eq_diag_plus_5H_over_p': okK, 'H_terms': [len(h) for h in H]})
facs = {}
for i in range(3):
    facs[f'Q{i}'] = Q[i]; facs[f'K{i}'] = K[i]; facs[f'H{i}'] = H[i]; facs[f'diag{i}'] = P(f'X{i}**4-X{i}**2*Y{i}**2+Y{i}**4'); facs[f'Q{i}^2'] = L.mul(Q[i], Q[i])
    facs[f'X{i}Y{i}'] = P(f'X{i}*Y{i}'); facs[f'X{i}^2-Y{i}^2'] = P(f'X{i}**2-Y{i}**2')
for i in range(3):
    for j in range(i + 1, 3):
        facs[f'Q{i}Q{j}'] = L.mul(Q[i], Q[j]); facs[f'Q{i}+Q{j}'] = L.add(Q[i], Q[j]); facs[f'Q{i}-Q{j}'] = L.add(Q[i], Q[j], s=-1)
facs['Q0+Q1+Q2'] = L.add(L.add(Q[0], Q[1]), Q[2]); facs['Q0Q1Q2'] = L.mul(L.mul(Q[0], Q[1]), Q[2])
facs['X0^2'] = P('X0**2'); facs['X0X1'] = P('X0*X1'); facs['X0Y1'] = P('X0*Y1')
facs['Fermat'] = P('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
mons6 = L.monomials(6, 1)
summary = []
for name, D in facs.items():
    t0 = time.time()
    Dp = L.red(D, p); e = L.deg(D); parD = L.parity(D); parp = (1 + parD) % 2
    mons = L.monomials(6 + e, parp); idx = {m: i for i, m in enumerate(mons)}
    gens = L.J_generators(Fi, 6 + e, parp, p)
    M = nmod_mat(len(gens) + len(mons6), len(mons), p)
    for r, G in enumerate(gens):
        for ee, v in G.items(): M[r, idx[ee]] = v
    for r, m in enumerate(mons6):
        for ee, v in L.shift(Dp, m, 1, p).items(): M[len(gens) + r, idx[ee]] = v
    rfull = M.rank()
    for r, m in enumerate(mons6):
        for ee, v in L.shift(Dp, m, 1, p).items(): M[len(gens) + r, idx[ee]] = 0
    rJ = M.rank(); del M
    kerdim = len(mons6) - (rfull - rJ); ann = kerdim - 18
    exact = L.J_dim(6 + e, parp)
    rec = {'check': f'factor annihilator {name}', 'deg': e, 'parity': parD, 'p': p, 'rank_J_p': rJ, 'exact_J': exact, 'nodrop': rJ == exact, 'kerdim_p': kerdim, 'Ann_R(D)_6odd_dim': ann,
           'rigorous_over_Q': (rJ == exact and ann == 0), 'seconds': round(time.time() - t0, 1)}
    L.log(rec); summary.append((name, ann, rJ == exact))
L.log({'check': 'factor annihilator summary', 'n_factors': len(summary), 'all_zero_and_nodrop': all(a == 0 and nd for _, a, nd in summary), 'nonzero': [n for n, a, nd in summary if a or not nd]})
