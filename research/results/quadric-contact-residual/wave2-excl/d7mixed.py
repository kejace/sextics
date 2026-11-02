"""Job: excl-d7mixed.  Degree-7 denominators with ALL G_aff sectors mixed: D = sum_{m in S_2, j} c_{m,j} m F_j (126 unknowns,
gauge (Z_k Z_5, j=5) dropped -> 120).  Necessary condition (Hessian relaxation): A*C^T Hess(F) C in J_14, i.e.
C^T Hess(F) C in (J:A)_8 = J_8 + <D8*>.  Quadrics = coordinates of NF_{J_8}(C^T Hess F C) modulo the D8* direction.
Steps: build (per S_8 character block), stats, monomials of I_2 and I_3 (per character block), export, msolve."""
import sys, time, json, os, subprocess
from w2lib import *
OUT = PK + '/d7mixed'
os.makedirs(OUT, exist_ok=True)
threads = int(os.environ.get('MSOLVE_THREADS', '60'))
msolve_timeout = int(os.environ.get('MSOLVE_TIMEOUT', '20000'))

t0 = time.time()
names = [(m, j) for m in monomials(2) for j in range(NV) if not (j == 5 and m[5] >= 1)]
n = len(names); assert n == 120
chi_u = [chi_add(chi(m), chi_neg(CHI_Z[j])) for m, j in names]
by_char = {}
for a in range(n):
    for b in range(a, n):
        by_char.setdefault(chi_add(chi_u[a], chi_u[b]), []).append((a, b))
quads = []
D8p = core.reduce_mod(D8, p)
sizes = {}
for cc, pairs in sorted(by_char.items()):
    S8 = sector(8, cc)
    vecs = []
    for (a, b) in pairs:
        ma, ia = names[a]; mb, ib = names[b]
        P = shift(Hess[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 2 if a != b else 1, p)
        vecs.append(S8.vec(P))
    NF = S8.nf_matrix(vecs)
    m = NF.ncols()
    cols = [[int(NF[r, j]) for r in range(len(pairs))] for j in range(m)]
    if cc == CHI0:
        v8 = S8.nf(D8p)
        assert any(v8), 'D8* in J_8?!'
        j0 = max(range(m), key=lambda j: 1 if v8[j] else 0)
        inv = pow(v8[j0], p - 2, p)
        newcols = []
        for j in range(m):
            if j == j0: continue
            f = v8[j] * inv % p
            newcols.append([(cols[j][r] - f * cols[j0][r]) % p for r in range(len(pairs))])
        cols = newcols
    for col in cols:
        cf = {pairs[r]: v for r, v in enumerate(col) if v}
        if cf: quads.append(cf)
    sizes[str(cc)] = {'pairs': len(pairs), 'dimR8': m}
print(f'built {len(quads)} quadrics in {n} unknowns ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
json.dump({'names': [[list(m), j] for m, j in names], 'chi': [list(c) for c in chi_u], 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads], 'sizes': sizes},
          open(OUT + '/d7mixed.quads.json', 'w'))

# ---- I_2 monomials (rref singleton rows) and rank
pairs_all = [(a, b) for a in range(n) for b in range(a, n)]
pidx = {pr: i for i, pr in enumerate(pairs_all)}
flat = [0] * (len(quads) * len(pairs_all))
for r, cf in enumerate(quads):
    for pr, v in cf.items(): flat[r * len(pairs_all) + pidx[pr]] = v
QM = nmod_mat(len(quads), len(pairs_all), flat, p)
R2, rk2 = QM.rref()
mono2 = []
for i in range(rk2):
    nz = [j for j in range(len(pairs_all)) if int(R2[i, j])]
    if len(nz) == 1: mono2.append(pairs_all[nz[0]])
sq2 = [a for (a, b) in mono2 if a == b]
nterms = [len(cf) for cf in quads]
log('d7mixed-build', d=7, n_unknowns=n, n_quadrics=len(quads), quadric_rank=rk2, n_pairs=len(pairs_all), max_terms=max(nterms), mean_terms=round(sum(nterms) / len(nterms), 1),
    n_monomials_in_I2=len(mono2), squares_in_I2=sq2, products_in_I2=[list(x) for x in mono2 if x[0] != x[1]][:50], seconds=round(time.time() - t0, 1),
    note='degree-8 Hessian formulation: NF_{J_8}(C^T Hess F C) mod <D8*>; dossier (deg-14 formulation, p=1e9+9): rank 950, no squares')

# ---- I_3 per character block: rows = c_v * q, columns = degree-3 monomials of that character
t1 = time.time()
chi_q = [chi_add(chi_u[k], chi_u[l]) for cf in quads for (k, l) in [next(iter(cf))]]
blocks = {}
for r, cq in enumerate(chi_q):
    for v in range(n):
        blocks.setdefault(chi_add(cq, chi_u[v]), []).append((v, r))
mono3, hilb3 = [], {}
for cc, lst in sorted(blocks.items()):
    cols3 = {}
    rows = []
    for (v, r) in lst:
        row = {}
        for (k, l), val in quads[r].items():
            key = tuple(sorted((v, k, l)))
            if key not in cols3: cols3[key] = len(cols3)
            row[cols3[key]] = (row.get(cols3[key], 0) + val) % p
        rows.append(row)
    # all degree-3 monomials of this character (some may not appear in any row)
    allm = [tuple(sorted((a, b, c))) for a in range(n) for b in range(a, n) for c in range(b, n) if chi_add(chi_add(chi_u[a], chi_u[b]), chi_u[c]) == cc]
    for key in allm:
        if key not in cols3: cols3[key] = len(cols3)
    inv = {i: k for k, i in cols3.items()}
    nc = len(cols3)
    flat = [0] * (len(rows) * nc)
    for r, row in enumerate(rows):
        for j, val in row.items(): flat[r * nc + j] = val
    M3 = nmod_mat(len(rows), nc, flat, p)
    R3, rk3 = M3.rref()
    found = []
    for i in range(rk3):
        nz = [j for j in range(nc) if int(R3[i, j])]
        if len(nz) == 1: found.append(inv[nz[0]])
    mono3.extend(found)
    hilb3[str(cc)] = {'rows': len(rows), 'cols': nc, 'rank': rk3, 'quotient': nc - rk3, 'monomials_in_I3': len(found)}
    print(f'  I_3 block {cc}: rows={len(rows)} cols={nc} rank={rk3} quot={nc-rk3} monos={len(found)} ({time.time()-t1:.0f}s)', file=sys.stderr, flush=True)
cubes = [m[0] for m in mono3 if m[0] == m[1] == m[2]]
sqlin = [m for m in mono3 if (m[0] == m[1] or m[1] == m[2]) and not (m[0] == m[1] == m[2])]
json.dump({'hilb3': hilb3, 'mono3': [list(m) for m in mono3], 'mono2': [list(m) for m in mono2]}, open(OUT + '/d7mixed.I3.json', 'w'))
log('d7mixed-I3', d=7, dim_S3=sum(v['cols'] for v in hilb3.values()), dim_I3=sum(v['rank'] for v in hilb3.values()), dim_quotient3=sum(v['quotient'] for v in hilb3.values()),
    n_monomials_in_I3=len(mono3), cubes_in_I3=cubes, n_square_times_linear_in_I3=len(sqlin), seconds=round(time.time() - t1, 1),
    note='cube c_a^3 in I_3 => c_a = 0 on V; c_a^2 c_b in I_3 => V subset {c_a=0} u {c_b=0}')

# ---- export and msolve
neq, nv = write_systems(quads, n, OUT + '/d7mixed')
log('d7mixed-export', n_eqs=neq, n_vars=nv, files=[OUT + '/d7mixed.ms', OUT + '/d7mixed.sing'])
t2 = time.time()
from msverdict import verdict
v = verdict(OUT + '/d7mixed.ms', threads=threads, timeout=msolve_timeout, verbose=2, logfile=OUT + '/d7mixed.msolve.log')
log('d7mixed-msolve', threads=threads, **{k: val for k, val in v.items() if k != 'log_tail'},
    note='msolve -g 1 (leading ideal, grevlex) on the full 120-unknown quadric system; ONLY-ORIGIN <=> every variable has a pure power among the leading monomials; rc=124 means timeout')
print('DONE', file=sys.stderr)
