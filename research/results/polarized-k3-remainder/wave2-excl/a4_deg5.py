"""(a)/(b): degree <= 5 exclusion, replicated exactly over Q with the FULL Koszul relaxation (no L-trick):
unknown products c_i c_j (21), c_i k_{jl,m} (540) and G*m (126, deg-4 m) treated as independent; the claim is that
the 21 columns c_i c_j are independent modulo the span of the other 666 columns, forcing c = 0, i.e. D = 0.
Degrees 0..4: [A] injective on R_d (a3) and J(F)_d = 0 => no E1-admissible D at all."""
import sys, time, json
sys.path.insert(0, '.')
from relax import *
from logutil import log
t0 = time.time()
Rx = Relaxation(5, with_koszul=True)
keys = Rx.quad_keys()
cc = [k for k in keys if Rx.params[k[0]][0] == 'a' and Rx.params[k[1]][0] == 'a']
ck = [k for k in keys if (Rx.params[k[0]][0] == 'a') != (Rx.params[k[1]][0] == 'a')]
kk = [k for k in keys if Rx.params[k[0]][0] == 'k' and Rx.params[k[1]][0] == 'k']
print('quadratic keys: cc', len(cc), 'ck', len(ck), 'kk', len(kk))
assert not kk
rows = {}
for e in Rx.E2: rows.setdefault(e, len(rows))
Gcols = []
for m in monomials(4):
    Gm = imul(Rx.G, {m: 1})
    for e in Gm: rows.setdefault(e, len(rows))
    Gcols.append(Gm)
def colvec(key):
    return {rows[e]: v[key] for e, v in Rx.E2.items() if key in v}
other = [colvec(k) for k in ck] + [{rows[e]: c for e, c in Gm.items()} for Gm in Gcols]
ccv = [colvec(k) for k in cc]
nr = len(rows)
M_other = sparse_cols_to_fmpz(other, nr); r_other = M_other.rank()
M_all = sparse_cols_to_fmpz(other + ccv, nr); r_all = M_all.rank()
print(f'rows={nr} rank(other)={r_other} rank(other+cc)={r_all} (need +{len(cc)})')
ok = (r_all == r_other + len(cc))
print('DEGREE 5 EXCLUDED over Q (all c_i c_j forced to 0):', ok, f'{time.time()-t0:.1f}s')
log(exp='a4_deg5_exclusion', rigor='exact over Q (fmpz_mat.rank); full Koszul relaxation, parity-free, UV frame',
    rows=nr, rank_other=r_other, rank_all=r_all, n_cc=len(cc), excluded=ok, seconds=round(time.time() - t0, 1))
