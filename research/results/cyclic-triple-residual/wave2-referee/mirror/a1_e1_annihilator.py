"""(a)/(b) Exact E1 annihilator table over Q by G-sectors: dim {D in S_d : A D in J(F)} / J(F)_d.
Also recomputed mod a second prime.  Rigorous statement produced: for every d in the range,
   Ann([A])_d = (J(F):A)_d / J(F)_d  has the stated dimension over Q (exact fmpz ranks per sector).
Usage: python3 a1_e1_annihilator.py dmin dmax [exact|modp]
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from lib import *

dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else 'exact'
p = None if mode == 'exact' else int(mode)
F, A, q, Fz, Az = load()
G = group_elements(Fz)
G = semi_invariant_subgroup(G, [Az])
gens = generators(G); chi = make_chi(gens); chiA = chi_of(Az, chi)
dF = [ideriv(Fz, i) for i in range(6)]
print(f'|G|={len(G)} gens={gens} chiA={chiA} mode={mode}', flush=True)
out = {'mode': mode, 'prime': p, 'group_order': len(G), 'generators': gens, 'chiA': chiA, 'degrees': {}}
for d in range(dmin, dmax + 1):
    t0 = time.time()
    bych = {}
    for m in mons(d): bych.setdefault(chi(m), []).append(m)
    tot_ker = 0; tot_J = 0; tot_JD = 0; sectors = {}
    for c, Dm in sorted(bych.items()):
        target = cadd(chiA, c)
        colsA = [{tuple(a+b for a, b in zip(m, e)): v for e, v in Az.items()} for m in Dm]
        colsJ = []
        for i in range(6):
            ci = cadd(target, chi(E6[i]))
            for m in mons(d + 1):
                if chi(m) != ci: continue
                col = {}
                for e, v in dF[i].items():
                    k = tuple(a+b for a, b in zip(m, e)); col[k] = col.get(k, 0) + v
                colsJ.append(col)
        rM = rank_of_columns(colsA + colsJ, p)
        rJ = rank_of_columns(colsJ, p)
        # dim {D in sector : A D in J_{d+6}} = dim S_d^chi - rank M + rank Jcols
        kerA = len(Dm) - rM + rJ
        colsJd = []
        for i in range(6):
            ci = cadd(c, chi(E6[i]))
            for m in mons(d - 5):
                if chi(m) != ci: continue
                col = {}
                for e, v in dF[i].items():
                    k = tuple(a+b for a, b in zip(m, e)); col[k] = col.get(k, 0) + v
                colsJd.append(col)
        rJd = rank_of_columns(colsJd, p) if colsJd else 0
        sectors[str(c)] = {'dimS': len(Dm), 'kerA_in_S': kerA, 'dimJ_d': rJd, 'ann_mod_J': kerA - rJd, 'dimJ_dplus6_sector': rJ, 'rankM': rM}
        tot_ker += kerA; tot_J += rJd; tot_JD += rJ
    rec = {'d': d, 'sum_kerA_in_S': tot_ker, 'sum_dimJ_d': tot_J, 'ci_dimJ_d': ci_dim_J(d), 'sum_dimJ_dplus6': tot_JD, 'ci_dimJ_dplus6': ci_dim_J(d+6),
           'ann_mod_J_total': tot_ker - tot_J, 'seconds': round(time.time() - t0, 1), 'sectors': sectors}
    out['degrees'][d] = rec
    print(f'd={d}: Ann([A])_d mod J = {tot_ker - tot_J}  (kerA_in_S={tot_ker}, dimJ_d={tot_J} ci={ci_dim_J(d)}, dimJ_{d+6}={tot_JD} ci={ci_dim_J(d+6)})  {rec["seconds"]}s', flush=True)
    nonzero = {c: s['ann_mod_J'] for c, s in sectors.items() if s['ann_mod_J']}
    if nonzero: print('   nonzero sectors:', nonzero, flush=True)
    log({'exp': 'e1_annihilator', 'mode': mode, 'prime': p, 'd': d, 'ann_mod_J': tot_ker - tot_J, 'dimJ_d': tot_J, 'ci_dimJ_d': ci_dim_J(d),
         'dimJ_dplus6_ok': tot_JD == ci_dim_J(d+6), 'nonzero_sectors': nonzero, 'seconds': rec['seconds']})
    receipt(f'e1_annihilator_{mode}_{dmin}_{dmax}', out)
