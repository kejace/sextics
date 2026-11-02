"""EXACT stage-2 relaxation for a sigma-eigen sector (d, par, k) -- see ref_stage2.py for the construction.  All duals are exact
integer left-kernel vectors of the Koszul block matrix (fmpz_mat nullspace), so the quadrics y^T b(u) are exact necessary
conditions over Q; together with the exact Hessian relaxation they are reduced mod p1 and tested for m-primariness (msolve -g 1 /
Singular).  Projective emptiness of the reduced system implies emptiness over Q-bar.   usage: python3 ref_stage2_exact.py d par k [threads]"""
import sys, json, time, subprocess, os
from fractions import Fraction
from math import gcd
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from flint import fmpz_mat, fmpq_mat, nmod_mat
d, par, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); thr = int(sys.argv[4]) if len(sys.argv) > 4 else 4
p = L.P1; t0 = time.time()
F, Fi, H = L.Fpolys(None); A = L.red(L.A, None); pA, kA = 1, 3
Frest = {e: Fraction(-v, 7) for e, v in F.items() if e != (0,0,0,6,0,0)}
memo = {}
def nf_mono(e):
    if e[3] < 6: return {e: Fraction(1)}
    if e in memo: return memo[e]
    base = list(e); base[3] -= 6; base = tuple(base); out = {}
    for e2, v in Frest.items():
        ee = tuple(a + b for a, b in zip(base, e2))
        for e3, v3 in nf_mono(ee).items(): out[e3] = out.get(e3, 0) + v * v3
    memo[e] = {kk: v for kk, v in out.items() if v}; return memo[e]
def nf(P):
    out = {}
    for e, v in P.items():
        for e2, v2 in nf_mono(e).items(): out[e2] = out.get(e2, 0) + v * v2
    return {kk: v for kk, v in out.items() if v}
relax = L.hess_relaxation(d, [par], [k], Fi, H, None)
vars_ = relax['vars']; n = relax['n']
Dgen = [L.shift(Fi[i], m, 1, None) for (i, m) in vars_]
Dgen_d = [[L.pderiv(Dg, c, None) for c in range(6)] for Dg in Dgen]
AE = [L.shift(A, m, 1, None) for (i, m) in vars_]
def Bpart(a, b):
    ib = vars_[b][0]
    return L.add(L.mul(Dgen[a], L.pderiv(AE[b], ib, None), None), L.scale(L.mul(AE[b], Dgen_d[a][ib], None), -1, None), None)
pairs = [(a, b) for a in range(n) for b in range(a, n)]
B = {}
for (a, b) in pairs: B[(a, b)] = nf(Bpart(a, a)) if a == b else nf(L.add(Bpart(a, b), Bpart(b, a), None))
print(f'sector d={d} par={par} k={k}: n={n} relaxation eqs {relax["n_eqs"]}; b_ab built [{time.time()-t0:.0f}s]', flush=True)
cols = []
for a in range(n):
    for c in range(6):
        for dd in range(c + 1, 6):
            pp = (pA + par + L.PI[c] + L.PI[dd]) % 2; kp = (kA + k + L.W[c] + L.W[dd]) % 7
            for m in L.monomials(d - 4, pp, kp):
                dm_c, dm_d = L.pderiv({m: 1}, c, None), L.pderiv({m: 1}, dd, None)
                col = L.add(L.mul(Dgen[a], L.add(L.mul(dm_c, Fi[dd], None), L.scale(L.mul(dm_d, Fi[c], None), -1, None), None), None),
                            L.scale(L.shift(L.add(L.mul(Fi[dd], Dgen_d[a][c], None), L.scale(L.mul(Fi[c], Dgen_d[a][dd], None), -1, None), None), m, 1, None), -1, None), None)
                cols.append(nf(col))
print('Koszul columns', len(cols), f'[{time.time()-t0:.0f}s]', flush=True)
rows = {}
for c in cols:
    for e in c: rows.setdefault(e, len(rows))
for b in B.values():
    for e in b: rows.setdefault(e, len(rows))
nr, nc, nb = len(rows), len(cols), len(pairs)
# integral matrix: all denominators divide 49 (at most two Y0^6-reductions in degree 2d <= 16)
den = 1
for col in cols:
    for v in col.values(): den = den * v.denominator // gcd(den, v.denominator)
for b in B.values():
    for v in b.values(): den = den * v.denominator // gcd(den, v.denominator)
MT = fmpz_mat(nc, nr)          # M^T (columns of M become rows)
for c, col in enumerate(cols):
    for e, v in col.items(): MT[c, rows[e]] = int(v * den)
Y, nul = MT.nullspace()        # exact integer left kernel of M: columns of Y (first nul)
rkM = nr - nul
print(f'M: {nr} x {nc}, exact rank {rkM}, exact left-kernel dim {nul} [{time.time()-t0:.0f}s]', flush=True)
# quadrics: for each kernel vector y, q_y = (y . b_ab)_{a<=b}  -> exact integer matrix (nul x nb); its row space is Q_Q
Bm = fmpz_mat(nr, nb)
for j, pr in enumerate(pairs):
    for e, v in B[pr].items(): Bm[rows[e], j] = int(v * den)
YT = fmpz_mat(nul, nr)
for v in range(nul):
    for r in range(nr): YT[v, r] = Y[r, v]
QM = YT * Bm                   # nul x nb
# row-reduce exactly to a basis of Q_Q
Qq = fmpq_mat(QM)
R, rk = Qq.rref()
Qbasis = []
for i in range(rk):
    row = {pairs[j]: Fraction(int(R[i, j].p), int(R[i, j].q)) for j in range(nb) if R[i, j] != 0}
    Qbasis.append(row)
print(f'exact stage-2 necessary quadrics: {rk} (relaxation: {relax["n_eqs"]}) [{time.time()-t0:.0f}s]', flush=True)
tag = f'stage2exact_d{d}_p{par}_k{k}'
json.dump({'d': d, 'par': par, 'k': k, 'vars': vars_, 'n': n, 'rows': nr, 'cols_M': nc, 'rank_M': rkM, 'leftkernel_dim': nul, 'n_stage2': rk,
           'n_relax': relax['n_eqs'], 'pairs': pairs, 'quadrics': [{f'{a},{b}': str(v) for (a, b), v in q.items()} for q in Qbasis],
           'relaxation': [{f'{a},{b}': str(v) for (a, b), v in q.items()} for q in relax['eqs']]}, open(f'{tag}.json', 'w'))
# reduce everything mod p1 and test m-primariness
def redq(q): return {pr: (v.numerator * pow(v.denominator, -1, p)) % p for pr, v in q.items()}
allq = [redq(q) for q in Qbasis + relax['eqs']]
sysb = {'n': n, 'eqs': allq}
L.write_msolve(sysb, f'{tag}_both_p{p}.ms', p)
names = [f'u{i}' for i in range(n)]
if n <= 12:
    L.write_singular(sysb, f'{tag}_both_p{p}.sing', p, 'ideal SI = std(I);\n"gb_size", size(SI), "dim", dim(SI), "vdim", vdim(SI);\nquit;\n')
    r = subprocess.run(['Singular', '-q', f'{tag}_both_p{p}.sing'], capture_output=True, text=True, timeout=3000)
    print('Singular mod p1 (exact relaxation + exact stage-2, reduced):', r.stdout.strip(), flush=True)
else:
    t = time.time()
    r = subprocess.run(['timeout', '2700', 'msolve', '-f', f'{tag}_both_p{p}.ms', '-o', f'{tag}_both_p{p}.lead', '-g', '1', '-t', str(thr)], capture_output=True, text=True)
    txt = open(f'{tag}_both_p{p}.lead').read() if os.path.exists(f'{tag}_both_p{p}.lead') else ''
    if '[' in txt:
        body = txt[txt.index('[') + 1: txt.rindex(']')]; items = [x.strip() for x in body.replace('\n', '').split(',') if x.strip()]
        pure = {it.split('^')[0] for it in items if '*' not in it}
        missing = [v for v in names if v not in pure]
        print(f'msolve -g 1 rc {r.returncode}: {len(items)} leading monomials; variables without a pure power: {missing} => m-primary: {not missing} [{time.time()-t:.0f}s]', flush=True)
    else:
        print('msolve produced no leading ideal (timeout or error)', r.returncode, r.stderr[-300:], f'[{time.time()-t:.0f}s]', flush=True)
print('STAGE2_EXACT_DONE', tag, f'[{time.time()-t0:.0f}s]')
