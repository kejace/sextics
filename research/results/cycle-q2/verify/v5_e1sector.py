"""Per-G-sector E1 kernel at degree d (default 6): K_sec = {D in S_d^sec : A D in J_{d+12}}, exact over Z (fmpz nullspace),
plus mod p.  Checks: sum of nontrivial dims = 26 (d=6); trivial sector kernel = (J_6)_0 = span{z_i F_i}; lists the kernel
elements (exact integer vectors) for the fixed-D sweep (A.6 replication).  Output JSON: e1sector_d{d}.json"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify')
import vcheck as V
import numpy as np
from flint import nmod_mat, fmpz_mat
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
p = 999983
F, A, q = V.load_instance()
F2 = {m: int(v) for m, v in V.scale(F, 2).items()}
Ai = {m: int(v) for m, v in A.items()}
G = V.diag_group(F)
def gens_of(Gm):
    gens = []; S = {(0,)*6}
    for w in Gm:
        tw = tuple(int(x) for x in w)
        if tw in S: continue
        gens.append(tw)
        frontier = list(S)
        while frontier:
            nxt = []
            for s in frontier:
                t = tuple((a+b) % 6 for a, b in zip(s, tw))
                if t not in S: S.add(t); nxt.append(t)
            frontier = nxt
        if len(S) == len(Gm): break
    return np.array(gens, dtype=np.int64)
gens = gens_of(G)
def key(m): return tuple(int(x) for x in (np.array(m) @ gens.T) % 6)
kA = key(list(Ai.keys())[0])
F2i = [V.pderiv(F2, i) for i in range(6)]
md = V.monomials(d); bykey = {}
for m in md: bykey.setdefault(key(m), []).append(m)
# J_{d+12} generators by key
mk5 = V.monomials(d + 12 - 5); gcols = {}
for i in range(6):
    kFi = key(list(F2i[i].keys())[0])
    for m in mk5:
        kk = tuple((a + b) % 6 for a, b in zip(key(m), kFi))
        gcols.setdefault(kk, []).append(V.mul(F2i[i], {m: 1}))
# J_d generators by key
gJd = {}
if d >= 5:
    for i in range(6):
        kFi = key(list(F2i[i].keys())[0])
        for m in V.monomials(d - 5):
            kk = tuple((a + b) % 6 for a, b in zip(key(m), kFi))
            gJd.setdefault(kk, []).append(V.mul(F2i[i], {m: 1}))
t0 = time.time()
results = {}; total_nontriv = 0; total_nontriv_p = 0
for kk, ms in sorted(bykey.items()):
    kE = tuple((a + b) % 6 for a, b in zip(kk, kA))
    Acols = [V.mul(Ai, {m: 1}) for m in ms]
    Jcols = gcols.get(kE, [])
    cols = Acols + Jcols
    rows = {}
    for c in cols:
        for m in c: rows.setdefault(m, len(rows))
    M = [[0] * len(cols) for _ in rows]
    for j, c in enumerate(cols):
        for m, v in c.items(): M[rows[m]][j] = v
    # exact nullspace over Z, project to the D-coordinates
    Mz = fmpz_mat(M)
    NS, nul = Mz.nullspace()
    # kernel vectors: columns 0..nul-1 of NS (rows = cols of M)
    Dvecs = []
    for c in range(nul):
        vec = [int(NS[j, c]) for j in range(len(ms))]
        if any(vec): Dvecs.append(vec)
    rkD = fmpz_mat(Dvecs).rank() if Dvecs else 0
    # mod p version
    flat = [x % p for r in M for x in r]
    Mp = nmod_mat(len(rows), len(cols), flat, p)
    rk_p = Mp.rank(); rkJ_p = nmod_mat(len(rows), len(Jcols), [x % p for r in M for x in r[len(ms):]], p).rank() if Jcols else 0
    dimK_p = len(ms) - (rk_p - rkJ_p)
    # J_d in this sector
    Jd = gJd.get(kk, [])
    dimJd = V.rank_of_cols(Jd)[0] if Jd else 0
    # exact basis of K (reduced): row-reduce Dvecs over Q
    basis = []
    if Dvecs:
        R = fmpz_mat(Dvecs)
        Rr, den = R.rref() if hasattr(R, 'rref') else (R, 1)
        for i in range(rkD):
            row = [int(Rr[i, j]) for j in range(len(ms))]
            if any(row): basis.append(row)
    results[str(kk)] = {'dimS': len(ms), 'dimJd': dimJd, 'dimK_exact': rkD, 'dimK_modp': dimK_p, 'nontrivial_exact': rkD - dimJd,
                        'monomials': [list(m) for m in ms], 'K_basis_exact': basis, 'Yparity': int(sum(ms[0][3:]) % 2)}
    total_nontriv += rkD - dimJd; total_nontriv_p += dimK_p - dimJd
    print(kk, 'dimS', len(ms), 'J_d', dimJd, 'K exact', rkD, 'K mod p', dimK_p, 'nontrivial', rkD - dimJd, round(time.time() - t0, 1), 's', flush=True)
summary = {'d': d, 'p': p, 'sum_nontrivial_exact': total_nontriv, 'sum_nontrivial_modp': total_nontriv_p, 'n_sectors': len(bykey), 'seconds': round(time.time() - t0, 1),
           'trivial_sector': results.get(str((0,) * len(gens)))}
json.dump({'summary': summary, 'sectors': results, 'gens': gens.tolist()}, open(f'/tmp/hodge/results/cycle-q2/verify/e1sector_d{d}.json', 'w'))
print(json.dumps({k: v for k, v in summary.items() if k != 'trivial_sector'}))
if d == 6:
    ts = results[str((0,) * len(gens))]
    V.log(claim='E1 kernel at d=6 per G-sector: total nontrivial = 26; trivial sector has 12 monomials, 6-dim kernel, which equals (J_6)_0 = span{z_i F_i} (nontrivial part 0)',
          source='dossier A.6 / attack1 exp_a', verdict='CONFIRMED' if (total_nontriv == 26 and ts['dimS'] == 12 and ts['dimK_exact'] == 6 and ts['nontrivial_exact'] == 0) else 'CHECK',
          evidence=summary, prime=p, notes='exact over Z (fmpz nullspace) and mod 999983; the "cones" of A.6 in the trivial sector are elements of J_6 (pure-power combinations of z_i F_i)')
