import sys, time, itertools
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc import *
F, A, q = load()
p = P0
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]

# (b) diagonal symmetry group: gamma_i = zeta_n^{a_i}; brute force over n=42 roots (lcm of 6,7,14)
n = 42
cnt = 0; elems = []
exps = list(F.keys())
for a in itertools.product(range(n), repeat=5):
    aa = (0,) + a  # normalize gamma_0 = 1 projectively (X0 appears in F so fine)
    vals = {sum(e[i]*aa[i] for i in range(NV)) % n for e in exps}
    if len(vals) == 1:
        cnt += 1; elems.append(aa)
print('diagonal symmetries with gamma_0 = 1 (over mu_42):', cnt)
# projective: also allow overall scalar; F terms must all scale equally. Count orbits by weight of A under them
for aa in elems:
    print('  ', aa, 'F-char', (sum(e[i]*aa[i] for i in range(NV)) % n for e in [exps[0]]).__next__(), 'A-char', (3*aa[0]+3*aa[3]) % n)

# (a) sigma-weight sectors of R_6^- (Y-odd degree 6) and R_6^+ 
def jac_dim(deg, par, k):
    mons = wmonomials(deg, par, k)
    idx = {m: i for i, m in enumerate(mons)}
    cols = []
    for i in range(NV):
        for m in monomials(deg - 5):
            col = shift(Fi[i], m, 1, p)
            col = {e: v for e, v in col.items() if (par is None or yparity(e) == par) and (k is None or wt(e) == k)}
            if col: cols.append(col)
    if not cols: return len(mons), 0
    R, C = len(mons), len(cols)
    flat = [0]*(R*C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[idx[e]*C + j] = v
    M = nmod_mat(R, C, flat, p)
    return len(mons), M.rank()
tot = 0
print('deg 6, Y-odd (q=1 targets): weight k: #monomials, dim J_6, dim R_6')
for k in range(7):
    nm, rk = jac_dim(6, 1, k)
    print(f'  k={k}: {nm:4d} {rk:4d} {nm-rk:4d}'); tot += nm-rk
print('  total dim R_6^- =', tot)
tot = 0
print('deg 12, Y-odd (q=2 targets):')
for k in range(7):
    nm, rk = jac_dim(12, 1, k)
    print(f'  k={k}: {nm:4d} {rk:4d} {nm-rk:4d}'); tot += nm-rk
print('  total dim R_12^- =', tot)
# is A nonzero in R_6? (already: yes) -- compute the weight-3 sector basis containing A
