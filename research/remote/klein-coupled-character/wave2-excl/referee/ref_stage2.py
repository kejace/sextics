"""Stage-2 relaxation for a sigma-eigen sector (d, par, k): constant left-kernel duals of the Koszul-form E2 system over the
whole sector family D(u) = sum_a u_a m_a F_{i_a} (u = gauge-fixed relaxation unknowns).  M(u) p = b(u) with M linear, b quadratic
in u; every y with y^T M_a = 0 for all blocks a gives the necessary quadric y^T b(u) = 0.  These quadrics contain the Hessian
relaxation (duals vanishing on J).  If (relaxation + stage-2 quadrics) is m-primary mod p, the sector has no certificate
denominator over Q-bar (projective emptiness transfers).   usage: python3 ref_stage2.py <d> <par> <k> [msolve_threads]"""
import sys, json, time, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from flint import nmod_mat
d, par, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); thr = int(sys.argv[4]) if len(sys.argv) > 4 else 4
p = L.P1; t0 = time.time()
F, Fi, H = L.Fpolys(p); A = L.red(L.A, p); pA, kA = 1, 3
inv7 = pow(7, -1, p)
Frest = {e: (-v * inv7) % p for e, v in F.items() if e != (0,0,0,6,0,0)}
memo = {}
def nf_mono(e):
    if e[3] < 6: return {e: 1}
    if e in memo: return memo[e]
    base = list(e); base[3] -= 6; base = tuple(base); out = {}
    for e2, v in Frest.items():
        ee = tuple(a + b for a, b in zip(base, e2))
        for e3, v3 in nf_mono(ee).items(): out[e3] = (out.get(e3, 0) + v * v3) % p
    memo[e] = {kk: v for kk, v in out.items() if v}; return memo[e]
def nf(P):
    out = {}
    for e, v in P.items():
        for e2, v2 in nf_mono(e).items(): out[e2] = (out.get(e2, 0) + v * v2) % p
    return {kk: v for kk, v in out.items() if v}
relax = L.hess_relaxation(d, [par], [k], Fi, H, p)
vars_ = relax['vars']; n = relax['n']
Dgen = [L.shift(Fi[i], m, 1, p) for (i, m) in vars_]                  # D^{(a)} = m_a F_{i_a}
Dgen_d = [[L.pderiv(Dg, c, p) for c in range(6)] for Dg in Dgen]
AE = [L.shift(A, m, 1, p) for (i, m) in vars_]                        # (A e_a)_{i_a} = A m_a
def Bpart(a, b):   # D^{(a)} div(A e_b) - (A e_b)(D^{(a)})
    ib = vars_[b][0]
    return L.add(L.mul(Dgen[a], L.pderiv(AE[b], ib, p), p), L.scale(L.mul(AE[b], Dgen_d[a][ib], p), p - 1, p), p)
pairs = [(a, b) for a in range(n) for b in range(a, n)]
B = {}
for (a, b) in pairs: B[(a, b)] = nf(Bpart(a, a)) if a == b else nf(L.add(Bpart(a, b), Bpart(b, a), p))
print(f'sector d={d} par={par} k={k}: n={n} relaxation eqs {relax["n_eqs"]}; b_ab built [{time.time()-t0:.0f}s]', flush=True)
cols = []
for a in range(n):
    for c in range(6):
        for dd in range(c + 1, 6):
            pp = (pA + par + L.PI[c] + L.PI[dd]) % 2; kp = (kA + k + L.W[c] + L.W[dd]) % 7
            for m in L.monomials(d - 4, pp, kp):
                dm_c, dm_d = L.pderiv({m: 1}, c, p), L.pderiv({m: 1}, dd, p)
                col = L.add(L.mul(Dgen[a], L.add(L.mul(dm_c, Fi[dd], p), L.scale(L.mul(dm_d, Fi[c], p), p - 1, p), p), p),
                            L.scale(L.shift(L.add(L.mul(Fi[dd], Dgen_d[a][c], p), L.scale(L.mul(Fi[c], Dgen_d[a][dd], p), p - 1, p), p), m, 1, p), p - 1, p), p)
                cols.append(nf(col))
print('Koszul columns', len(cols), f'[{time.time()-t0:.0f}s]', flush=True)
rows = {}
for c in cols:
    for e in c: rows.setdefault(e, len(rows))
for b in B.values():
    for e in b: rows.setdefault(e, len(rows))
nr, nc, nb = len(rows), len(cols), len(pairs)
flat = [0] * (nr * (nc + nb))
for c, col in enumerate(cols):
    for e, v in col.items(): flat[rows[e] * (nc + nb) + c] = v
for j, pr in enumerate(pairs):
    for e, v in B[pr].items(): flat[rows[e] * (nc + nb) + nc + j] = v
M = nmod_mat(nr, nc + nb, flat, p)
X, nul = M.nullspace()
K = nmod_mat(nul, nb, [int(X[nc + j, v]) for v in range(nul) for j in range(nb)], p)
rK = K.rank(); X2, nul2 = K.nullspace()
Q = [{pairs[j]: int(X2[j, v]) for j in range(nb) if int(X2[j, v])} for v in range(nul2)]
print(f'matrix {nr} x {nc + nb}; nullity {nul}; dim Q^perp {rK}; stage-2 necessary quadrics: {len(Q)} (relaxation had {relax["n_eqs"]}) [{time.time()-t0:.0f}s]', flush=True)
# is the relaxation contained in Q (as expected)?
names = [f'u{i}' for i in range(n)]
tag = f'stage2_d{d}_p{par}_k{k}'
sysd = {'n': n, 'eqs': Q}
L.write_msolve(sysd, f'{tag}_p{p}.ms', p)
sysb = {'n': n, 'eqs': Q + relax['eqs']}
L.write_msolve(sysb, f'{tag}_both_p{p}.ms', p)
json.dump({'prime': p, 'd': d, 'par': par, 'k': k, 'vars': vars_, 'n': n, 'n_relax': relax['n_eqs'], 'n_stage2': len(Q), 'rows': nr, 'cols_M': nc, 'nullity': nul}, open(f'{tag}_p{p}.json', 'w'))
if n <= 12:
    L.write_singular(sysb, f'{tag}_both_p{p}.sing', p, 'ideal SI = std(I);\n"gb_size", size(SI), "dim", dim(SI), "vdim", vdim(SI);\nquit;\n')
    r = subprocess.run(['Singular', '-q', f'{tag}_both_p{p}.sing'], capture_output=True, text=True, timeout=3000)
    print('Singular (relaxation + stage-2):', r.stdout.strip(), flush=True)
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
print('STAGE2_DONE', tag, f'[{time.time()-t0:.0f}s]')
