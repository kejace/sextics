"""Degree-6 E2 reduction (packet R3): kernel of (S, N) -> A * sum_ij S_ij F_ij - sum_k N_k F_k  in S_12^{out}.
S symmetric 6x6 with quadratic entries of parity pattern 'sector'; N in S_7^6 with parities.
Checks: kernel projected to S-block equals the Euler span {S = z u^T + u z^T}.  Then the rank-one / mixed-parity
argument (diagonal entries) forces w_i = lambda z_i, all lambda equal => D = 6 lambda F illegal.
Usage: python3 hesskernel_v.py p"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
p = int(sys.argv[1])
F, A, q = parse_rec()
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
H = [[pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
pz = [0, 0, 0, 1, 1, 1]
out = {}
for sector in (0, 1):   # sector 0: w of type +, S_ij parity p_i+p_j (output of A*S*H odd); sector 1: cross terms, S_ij parity p_i+p_j+1 (output even)
    t0 = time.time()
    Spar = lambda i, j: (pz[i] + pz[j] + sector) % 2
    outpar = (1 + sector) % 2   # A odd * (S_ij F_ij: parity Spar + pz_i + pz_j = sector) -> 1 + sector
    cols = []; names = []
    for i in range(6):
        for j in range(i, 6):
            for m in monomials(2, Spar(i, j)):
                mult = 1 if i == j else 2
                cols.append(pscale(pmul(Ap, pshift(H[i][j], m, 1, p), p), mult, p)); names.append(('S', i, j, m))
    nS = len(cols)
    for k in range(6):
        for m in monomials(7, (outpar - pz[k]) % 2):
            cols.append(pscale(pshift(Fi[k], m, 1, p), p - 1, p)); names.append(('N', k, m))
    rows = monomials(12, outpar)
    M = build_mat(cols, rows, p)
    X, nul = M.nullspace()
    rank = M.ncols() - nul
    proj = nmod_mat(nul, nS, [int(X[j, k]) for k in range(nul) for j in range(nS)], p)
    _, rkS = proj.rref()
    # Euler span: S = z u^T + u z^T, u_i linear of parity (pz_i + sector)
    eul = []
    for i in range(6):
        for lm in monomials(1, (pz[i] + sector) % 2):
            v = {}
            # u = lm * e_i ; S_ab = z_a u_b + z_b u_a : nonzero for b = i (z_a * lm) and a = i
            for a in range(6):
                za = [0]*6; za[a] = 1
                mon = tuple(x + y for x, y in zip(za, lm))
                key = ('S', min(a, i), max(a, i), mon)
                v[key] = (v.get(key, 0) + (2 if a == i else 1)) % p
            eul.append(v)
    sidx = {nm: j for j, nm in enumerate(names[:nS])}
    E = nmod_mat(len(eul), nS, p)
    for r, v in enumerate(eul):
        for nm, c in v.items(): E[r, sidx[nm]] = c
    rkE = E.rank()
    both = nmod_mat(nul + len(eul), nS, p)
    for r in range(nul):
        for j in range(nS): both[r, j] = int(X[j, r])
    for r in range(len(eul)):
        for j in range(nS): both[nul + r, j] = int(E[r, j])
    rk_both = both.rank()
    rec = {'p': p, 'sector': sector, 'rows': M.nrows(), 'cols': M.ncols(), 'nS': nS, 'rank': rank, 'nullity': nul,
           'S_block_kernel_dim': rkS, 'euler_span_dim': rkE, 'S_kernel_equals_euler_span': (rkS == rkE == rk_both),
           'seconds': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True)
    out[sector] = rec
    with open('/tmp/hodge/results/dense-rational/verify/hesskernel_v.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
ok = all(out[s]['S_kernel_equals_euler_span'] for s in (0, 1))
log('[rerun with corrected Euler vectors; supersedes earlier REFUTED line caused by a verifier bug in the diagonal Euler entries] R3 (deg D >= 7): degree-6 Hessian map kernel has S-block exactly the Euler span in both output-parity sectors (so E2 forces w_i = lambda z_i, D = 6 lambda F)',
    'CONFIRMED' if ok else 'REFUTED', f"p={p}: sector0 rank {out[0]['rank']} S-kernel {out[0]['S_block_kernel_dim']} = Euler {out[0]['euler_span_dim']}; sector1 rank {out[1]['rank']} S-kernel {out[1]['S_block_kernel_dim']} = Euler {out[1]['euler_span_dim']}",
    'rigor direction: modular nullity is an upper bound on the Q-nullity, and the Euler vectors are exact kernel elements, so the S-block of the Q-kernel is exactly the Euler span. Mixed-parity w handled: diagonal entries of both sectors give z_i | w_i^+ w_i^- and z_i | (w_i^+)^2+(w_i^-)^2, parity kills w^-; then off-diagonals force equal lambdas.')
