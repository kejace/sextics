"""Independent E1-kernel sweep: K_full(d,par) = {D in S_d^par : A D in J}, kernel mod J = dim K_full - dim J_d^par.
Also records rank of the N-block vs exact dim J_{d+6}^{par'} (needed for the 'rigorous' direction).
Usage: python3 e1sweep_v.py p mode d par   (mode direct|dual)"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
p = int(sys.argv[1]); mode = sys.argv[2]; d = int(sys.argv[3]); par = int(sys.argv[4])
F, A, q = parse_rec()
S, R = hilbert()
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
t0 = time.time()
if mode == 'direct':
    parT = (par + 1) % 2
    Dmons = monomials(d, par)
    cols = [pshift(Ap, m, 1, p) for m in Dmons]
    nD = len(cols)
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (parT - pi) % 2):
            cols.append(pshift(Fi[i], m, 1, p))
    rows = monomials(d + 6, parT)
    M = build_mat(cols, rows, p)
    Rn, C = M.nrows(), M.ncols()
    # rank of N-block alone
    MN = build_mat(cols[nD:], rows, p)
    rN = MN.rank()
    dimJ_t = S[d + 6][parT] - R[d + 6][parT]
    X, nullity = M.nullspace()
    rank = C - nullity
    proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p)
    _, kf = proj.rref()
    dimJd = S[d][par] - R[d][par]
    rec = {'p': p, 'mode': 'direct', 'd': d, 'par': par, 'rows': Rn, 'cols': C, 'nD': nD, 'rank': rank, 'nullity': nullity,
           'rankN': rN, 'dimJ_target_exact': dimJ_t, 'rankN_attains_cap': rN == dimJ_t,
           'Kfull': kf, 'dimJ_d': dimJd, 'kernel_mod_J': kf - dimJd, 'seconds': round(time.time() - t0, 1)}
else:
    e = 18 - d
    # rank of [A]: R_e^{1-par} -> R_{e+6}^{par}  (adjoint of [A]: R_d^par -> R_{d+6}^{1-par} under the even socle pairing)
    src = (par + 1) % 2; tgt = par
    Dmons = monomials(e, src)
    cols = [pshift(Ap, m, 1, p) for m in Dmons]
    nD = len(cols)
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(e + 1, (tgt - pi) % 2):
            cols.append(pshift(Fi[i], m, 1, p))
    rows = monomials(e + 6, tgt)
    M = build_mat(cols, rows, p)
    rank = M.rank()
    MN = build_mat(cols[nD:], rows, p)
    rN = MN.rank()
    dimJ_t = S[e + 6][tgt] - R[e + 6][tgt]
    rankA = rank - dimJ_t
    rec = {'p': p, 'mode': 'dual', 'd': d, 'par': par, 'e': e, 'rows': M.nrows(), 'cols': M.ncols(), 'rank': rank, 'rankN': rN,
           'dimJ_target_exact': dimJ_t, 'rankN_attains_cap': rN == dimJ_t, 'rankA': rankA, 'dimR_src': R[e][src],
           'A_injective_on_src': rankA == R[e][src], 'dimR_d': R[d][par], 'kernel_mod_J': R[d][par] - rankA,
           'seconds': round(time.time() - t0, 1)}
print(json.dumps(rec), flush=True)
with open('/tmp/hodge/results/dense-rational/verify/e1sweep_v.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
