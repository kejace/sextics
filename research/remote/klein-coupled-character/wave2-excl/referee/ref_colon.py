"""Referee replay of T1 / E1-kernel tables at a new prime: dim (J:A)_d vs dim J_d per (parity, weight) sector, d = 4..10,
plus (J : Z_i)_6 = J_6 for the six coordinate linear forms (coordinate points of the NZD charts) and two generic forms."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
p = L.P1
F, Fi, H = L.Fpolys(p)
A = L.red(L.A, p)
out = {'prime': p, 'colon_A': {}, 'kernel_mod_J': {}, 'dimJ_total': {}, 'ci_total': {}, 'colon_linear': {}}
t0 = time.time()
for d in range(4, 11):
    ker = {0: [], 1: []}; totJ = 0
    for tp in (0, 1):
        for tk in range(7):
            r = L.colon_dim(d, tp, tk, A, Fi, p)
            ker[tp].append(r['kernel_mod_J']); totJ += r['dimJ']
    out['kernel_mod_J'][d] = {'even': ker[0], 'odd': ker[1], 'total': sum(ker[0]) + sum(ker[1])}
    out['dimJ_total'][d] = totJ; out['ci_total'][d] = L.comb(d + 5, 5) - L.ci_hilbert(d)
    print(f'd={d}: kernel of [A] mod J even {ker[0]} odd {ker[1]} total {sum(ker[0])+sum(ker[1])} | dim J_d {totJ} CI {out["ci_total"][d]} [{time.time()-t0:.0f}s]', flush=True)
# linear forms: (J : l)_6 = J_6 ?  (sector-wise for sigma-eigen l; full space otherwise)
def colon_full(d, mult):
    """dim (J:mult)_d - dim J_d over the full space (no sector split), mod p"""
    mons = L.monomials(d); dm = L.degree(mult)
    mons2 = L.monomials(d + dm); idx2 = {e: r for r, e in enumerate(mons2)}
    gens2 = [L.shift(Fi[i], m, 1, p) for i in range(6) for m in L.monomials(d + dm - 5)]
    ent = {}
    for c, g in enumerate(gens2):
        for e, v in g.items(): ent[(idx2[e], c)] = v
    nJ = len(gens2)
    for c, m in enumerate(mons):
        for e, v in L.shift(mult, m, 1, p).items(): ent[(idx2[e], nJ + c)] = v
    M = L.mat(len(mons2), nJ + len(mons), ent, p)
    rAll = M.rank(); rJ2 = L.mat(len(mons2), nJ, {k: v for k, v in ent.items() if k[1] < nJ}, p).rank()
    colon = len(mons) - (rAll - rJ2)
    dimJ = len(mons) - L.ci_hilbert(d)
    return colon, dimJ
forms = {'X0': {(1,0,0,0,0,0): 1}, 'X1': {(0,1,0,0,0,0): 1}, 'X2': {(0,0,1,0,0,0): 1}, 'Y0': {(0,0,0,1,0,0): 1}, 'Y1': {(0,0,0,0,1,0): 1}, 'Y2': {(0,0,0,0,0,1): 1},
         'X0+3X1-5X2': {(1,0,0,0,0,0): 1, (0,1,0,0,0,0): 3, (0,0,1,0,0,0): p-5}, 'Y0-2Y1+7Y2': {(0,0,0,1,0,0): 1, (0,0,0,0,1,0): p-2, (0,0,0,0,0,1): 7}}
for name, l in forms.items():
    c, j = colon_full(6, l)
    out['colon_linear'][name] = {'dim_colon_6': c, 'dim_J_6': j, 'NZD_on_R6': c == j}
    print(f'(J:{name})_6 = {c}, dim J_6 = {j} -> NZD on R_6: {c == j}', flush=True)
json.dump(out, open(f'referee_colon_p{p}.json', 'w'), indent=1)
print('COLON_DONE', f'[{time.time()-t0:.0f}s]')
