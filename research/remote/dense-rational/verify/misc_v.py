"""Cheap membership checks (rigorous non-membership when rank_p J_d = dim J_d): subset sums, blocks not in J_6; Hessian minors not in J_8.
Also R9 rank premises: restricted gradients on X0=0 and Y0=0 have no syzygies of multiplier degree 4 and only Koszul (15) in degree 5."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
import itertools
p = int(sys.argv[1])
F, A, q = parse_rec()
Fp = reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
S, R = hilbert()
def Jmat(d, par):
    gens = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(d - 5, (par - pi) % 2):
            gens.append(pshift(Fi[i], m, 1, p))
    rows = monomials(d, par)
    return gens, rows
def not_in_J(D, d, par):
    gens, rows = Jmat(d, par)
    M = build_mat(gens, rows, p); r0 = M.rank()
    M2 = build_mat(gens + [D], rows, p); r1 = M2.rank()
    return r0, r1, S[d][par] - R[d][par]
blocks = {0: {}, 2: {}, 4: {}, 6: {}}
for e, c in Fp.items(): blocks[e[3] + e[4] + e[5]][e] = c
cands6 = {'F60': blocks[0], 'F06': blocks[6], 'F42': blocks[2], 'F24': blocks[4], 'F60+F06': padd(blocks[0], blocks[6], p), 'sum Z_i^6': {tuple(6 if k == i else 0 for k in range(6)): 1 for i in range(6)}}
for I in itertools.combinations(range(6), 3):
    if (0 in I) != (3 in I):
        cands6['subset ' + ''.join(VARS[i] for i in I)] = {tuple(6 if k == i else 0 for k in range(6)): 1 for i in I}
res = {}
for nm, D in cands6.items():
    r0, r1, dj = not_in_J(D, 6, 0)
    res[nm] = (r0, r1, dj)
print('deg 6 (r0=rank J_6^+, r1=rank with D, exact dim J_6^+):', res)
ok6 = all(v[0] == v[2] and v[1] == v[0] + 1 for v in res.values())
log('12 Fermat subset sums, sum Z_i^6, F60, F06, F42, F24, F60+F06 are NOT in J_6 (hence fail E1 since (J:A)_6 = J_6)', 'CONFIRMED' if ok6 else 'REFUTED',
    f'p={p}: rank J_6^+ = {res["F60"][0]} = exact dim J_6^+ = {res["F60"][2]}; every candidate raises the rank by 1', 'non-membership rigorous: rank_p J_6 attains exact cap and augmented rank exceeds it')
H = [[pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
cands8 = {'F00F33-F03^2': padd(pmul(H[0][0], H[3][3], p), pmul(H[0][3], H[0][3], p), p, -1), 'F00F11-F01^2': padd(pmul(H[0][0], H[1][1], p), pmul(H[0][1], H[0][1], p), p, -1),
          'X0^8+X1^8+X2^8': {(8,0,0,0,0,0): 1, (0,8,0,0,0,0): 1, (0,0,8,0,0,0): 1}, 'F60*(X0^2+X1^2+X2^2)': pmul(blocks[0], {(2,0,0,0,0,0): 1, (0,2,0,0,0,0): 1, (0,0,2,0,0,0): 1}, p),
          'F60*(Y0^2+Y1^2+Y2^2)': pmul(blocks[0], {(0,0,0,2,0,0): 1, (0,0,0,0,2,0): 1, (0,0,0,0,0,2): 1}, p)}
res8 = {}
for nm, D in cands8.items():
    r0, r1, dj = not_in_J(D, 8, 0)
    res8[nm] = (r0, r1, dj)
print('deg 8:', res8)
ok8 = all(v[0] == v[2] and v[1] == v[0] + 1 for v in res8.values())
log('Hessian 2x2 minors F_X0X0F_Y0Y0-F_X0Y0^2, F_X0X0F_X1X1-F_X0X1^2, X0^8+X1^8+X2^8, F60*(sum X^2), F60*(sum Y^2) are NOT in J_8', 'CONFIRMED' if ok8 else 'REFUTED',
    f'p={p}: rank J_8^+ = {res8["X0^8+X1^8+X2^8"][0]} = exact {res8["X0^8+X1^8+X2^8"][2]}; all raise rank by 1', 'rigorous non-membership (cap attained)')
# --- R9 rank premises: restricted gradients ---
def restricted(var):
    G = [{e: c for e, c in Fi[i].items() if e[var] == 0} for i in range(6)]
    out = {}
    for mdeg, expect_null in ((4, 0), (5, 15)):
        cols = []
        for i in range(6):
            for m in monomials(mdeg):
                if m[var] == 0: cols.append(pshift(G[i], m, 1, p))
        rows = [m for m in monomials(mdeg + 5) if m[var] == 0]
        M = build_mat(cols, rows, p)
        rk = M.rank()
        out[mdeg] = (M.nrows(), M.ncols(), rk, M.ncols() - rk, expect_null)
    return out
r9 = {'X0=0': restricted(0), 'Y0=0': restricted(3)}
print('R9 premises:', r9)
ok9 = all(v[4][3] == 0 and v[5][3] == 15 for v in r9.values())
log('R9 rank premises: gradients restricted to X0=0 / Y0=0 have no syzygies with multiplier degree 4 (715x420 rank 420) and exactly the 15 Koszul ones in degree 5 (1001x756 rank 741)', 'CONFIRMED' if ok9 else 'REFUTED',
    f'p={p}: {r9}', 'rank lower bounds are rigorous; Koszul syzygies supply the 15 exact kernel vectors. The divisibility-lemma part of R9 (polynomial Piola atoms excluded) is packet-only and NOT re-derived here')
