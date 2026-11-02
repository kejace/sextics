"""(a) Exact Farkas certificates of inconsistency for fixed denominators (unsplit cleared identity, q = 1).

System (unknowns N_i of degree d+1 with Y-parity of e_i + 1 (A odd, D even), Q of degree 2d-6 Y-odd):
    L(N, Q) := -D sum_i N_i F_i + F (D div N - N(D)) - F^2 Q  =  -A D^2.
A vector y on the monomials of degree 2d+6 with  y . L(col) = 0 for every unknown column and  y . (-A D^2) = 1
proves that no (N, Q) exists over any field of characteristic 0 (and the check is exact integer arithmetic).
y is found on the G-sector of A D^2 (its support), and verified against ALL columns of the full unrestricted system.
Usage: python3 a2_farkas.py name 'D expr'
"""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from lib import *

name, Dsrc = sys.argv[1], sys.argv[2]
F, A, q, Fz, Az = load(); assert q == 1
D = core.parse(Dsrc); Dz = integral(D); d = core.degree(Dz)
t0 = time.time()
dF = [ideriv(Fz, i) for i in range(6)]; dD = [ideriv(Dz, i) for i in range(6)]
F2 = imul(Fz, Fz)
par = lambda m: (m[3] + m[4] + m[5]) % 2
parD = par(next(iter(Dz))); assert all(par(m) == parD for m in Dz), 'D must have pure parity'
# full column set (all parities allowed for the identity: N_i parity = parity(A) + parity(D) + parity(e_i) = 1 + parD + e_i-parity)
cols = []; meta = []
for i in range(6):
    pi = (1 + parD + (1 if i >= 3 else 0)) % 2
    for m in mons(d + 1):
        if par(m) != pi: continue
        col = {}
        t = imul({m: 1}, dF[i])
        for e, c in imul(Dz, t).items(): col[e] = col.get(e, 0) - c
        t2 = imul(Dz, ideriv({m: 1}, i))
        for e, c in imul({m: 1}, dD[i]).items(): t2[e] = t2.get(e, 0) - c
        for e, c in imul(Fz, t2).items(): col[e] = col.get(e, 0) + c
        col = {e: c for e, c in col.items() if c}
        cols.append(col); meta.append(('N', i, m))
pQ = (1 + 2*parD) % 2
for m in mons(2*d - 6):
    if par(m) != pQ: continue
    col = {e: -c for e, c in imul({m: 1}, F2).items()}
    cols.append(col); meta.append(('Q', -1, m))
rhs = {e: -c for e, c in imul(Az, imul(Dz, Dz)).items()}
print(f'{name}: d={d} full system: {len(cols)} columns, rows deg {2*d+6}; built {round(time.time()-t0,1)}s', flush=True)
# sector of the rhs under the stabilizer of (F, A, D)
G = semi_invariant_subgroup(group_elements(Fz), [Az, Dz]); gens = generators(G); chi = make_chi(gens)
cR = chi_of(rhs, chi)
sec_cols = [j for j, col in enumerate(cols) if col and chi_of(col, chi) == cR]
rows = {}
for j in sec_cols:
    for e in cols[j]: rows.setdefault(e, len(rows))
for e in rhs: rows.setdefault(e, len(rows))
R, C = len(rows), len(sec_cols)
print(f'  stabilizer order {len(G)}; sector system {R} rows x {C} cols', flush=True)
def build(p=None):
    flat = [0] * (R * (C + 1))
    for jj, j in enumerate(sec_cols):
        for e, v in cols[j].items(): flat[rows[e]*(C+1) + jj] = v if p is None else v % p
    for e, v in rhs.items(): flat[rows[e]*(C+1) + C] = v if p is None else v % p
    return flat
flat = build()
Mz = fmpz_mat(R, C + 1, flat)
t1 = time.time()
rank_aug = Mz.rank()
M0 = fmpz_mat(R, C, [flat[r*(C+1) + j] for r in range(R) for j in range(C)])
rank0 = M0.rank()
print(f'  exact ranks: system {rank0}, augmented {rank_aug} ({round(time.time()-t1,1)}s)', flush=True)
rec = {'name': name, 'D': core.to_str(D), 'd': d, 'full_columns': len(cols), 'sector_rows': R, 'sector_cols': C, 'stabilizer_order': len(G),
       'rank_exact': rank0, 'rank_aug_exact': rank_aug, 'consistent': rank0 == rank_aug}
if rank0 == rank_aug:
    print('  CONSISTENT (no Farkas certificate)'); rec['farkas'] = None
else:
    # Farkas vector: y with y^T M0 = 0, y^T rhs = 1.  Work with Nt = [M0 | rhs]^T  ((C+1) x R): Nt y = e_last.
    # choose pivot columns (of Nt, i.e. rows of M) and independent rows mod p, solve the square system exactly.
    p = P3
    flatp = build(p)
    Mp = nmod_mat(R, C + 1, flatp, p)
    NtP = Mp.transpose()                      # (C+1) x R
    rr, rk = NtP.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < R and rr[i, j] == 0: j += 1
        piv.append(j)
    # independent rows of Nt: rref of Nt^T = Mp gives pivot columns of Mp = independent rows of Nt
    rr2, rk2 = Mp.rref()
    indep = []
    for i in range(rk2):
        j = 0
        while j < C + 1 and rr2[i, j] == 0: j += 1
        indep.append(j)
    assert rk == rk2 == rank_aug, (rk, rk2, rank_aug)
    assert C in indep, 'rhs row of Nt not among independent rows?'
    # square system: Nt[indep][:, piv] y_piv = e (1 at the rhs row)
    sq = fmpz_mat(rk, rk, [flat[r*(C+1) + j] for j in indep for r in piv])
    e = fmpz_mat(rk, 1, [1 if j == C else 0 for j in indep])
    t2 = time.time()
    sol = sq.solve(e)      # fmpq_mat
    y = {}
    for k, r in enumerate(piv):
        v = sol[k, 0]
        if v != 0: y[r] = Fraction(int(v.p), int(v.q))
    # exact verification against ALL columns of the full system
    rowkey = {v: k for k, v in rows.items()}
    ysparse = {rowkey[r]: v for r, v in y.items()}
    bad = 0
    for col in cols:
        s = sum(v * ysparse[e] for e, v in col.items() if e in ysparse)
        if s != 0: bad += 1
    val = sum(v * ysparse[e] for e, v in rhs.items() if e in ysparse)
    ok = (bad == 0 and val == 1)
    print(f'  Farkas vector: {len(y)} nonzeros, max denominator bits {max(v.denominator.bit_length() for v in y.values())}; '
          f'verified against all {len(cols)} columns: {ok} (bad={bad}, y.rhs={val}) ({round(time.time()-t2,1)}s)', flush=True)
    rec['farkas'] = {'verified': ok, 'nonzeros': len(y), 'prime_used_for_pivots': p,
                     'vector': {str(rowkey[r]): str(v) for r, v in y.items()}}
rec['seconds'] = round(time.time() - t0, 1)
receipt(f'farkas_{name}', rec)
log({'exp': 'farkas_fixed_D', 'name': name, 'd': d, 'rank': rank0, 'rank_aug': rank_aug, 'consistent': rec['consistent'],
     'farkas_verified': (rec['farkas'] or {}).get('verified'), 'full_columns': len(cols), 'seconds': rec['seconds']})
