"""E1-kernel sweep: dim of ker( [A] : R_d^+ -> R_{d+6}^- ) for the Jacobian ring R = S/J(F).
direct(d): Macaulay matrix with columns A*m (m in S_d^+) and m*F_i (m F_i odd, deg d+6); rows S_{d+6}^-.
           nullspace projected to the D-block gives K_full = {D in S_d^+ : A D in J}; kernel mod J = K_full - J_d^+.
dual(d):   ker on R_d^+ has dim = dim R_d^+ - rank([A]: R_{18-d}^- -> R_{24-d}^+).
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
import numpy as np

def macaulay_rank_and_null(F, A, d, parD, p, want_null=False):
    """Columns: A*m for m in S_d^{parD}; then m*F_i with total parity = parD+1 (A is odd).  Rows: S_{d+6}^{parD+1}."""
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    parT = (parD + 1) % 2
    cols = []
    Dmons = monomials(d, parD)
    for m in Dmons:
        cols.append(shift(Ap, m, 1, p))
    nD = len(cols)
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (parT - pi) % 2):
            cols.append(shift(Fi[i], m, 1, p))
    rows = {}
    for e in monomials(d + 6, parT):
        rows[e] = len(rows)
    R, C = len(rows), len(cols)
    t0 = time.time()
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    del flat
    t1 = time.time()
    if want_null:
        X, nullity = M.nullspace()
        rank = C - nullity
        # project nullspace onto D-block and get its rank / basis
        proj = nmod_mat(nD, nullity, [int(X[j, k]) for j in range(nD) for k in range(nullity)], p)
        Tr, rk = proj.transpose().rref()
        basis = []
        for i in range(rk):
            poly = {}
            for j in range(nD):
                v = int(Tr[i, j])
                if v: poly[Dmons[j]] = v
            basis.append(poly)
        t2 = time.time()
        return {'rows': R, 'cols': C, 'nD': nD, 'rank': rank, 'nullity': nullity, 'Kfull': rk, 'basis': basis,
                't_build': round(t1 - t0, 1), 't_null': round(t2 - t1, 1)}
    else:
        rank = M.rank()
        t2 = time.time()
        return {'rows': R, 'cols': C, 'nD': nD, 'rank': rank, 't_build': round(t1 - t0, 1), 't_rank': round(t2 - t1, 1)}

if __name__ == '__main__':
    F, A, q = load()
    S, Rh = hilbert_counts()
    p = P
    mode = sys.argv[1]
    ds = [int(x) for x in sys.argv[2:]]
    out = []
    for d in ds:
        t0 = time.time()
        if mode == 'direct':
            r = macaulay_rank_and_null(F, A, d, 0, p, want_null=True)
            Jd = S[d][0] - Rh[d][0]
            kmodJ = r['Kfull'] - Jd
            print(f"[direct d={d}] rows={r['rows']} cols={r['cols']} rank={r['rank']} nullity={r['nullity']} Kfull={r['Kfull']} dimJ_d^+={Jd} kernel_mod_J={kmodJ} ({time.time()-t0:.1f}s)", flush=True)
            rec = {'mode': 'direct', 'd': d, 'Kfull': r['Kfull'], 'dimJ': Jd, 'kernel_mod_J': kmodJ, 'rows': r['rows'], 'cols': r['cols'], 'seconds': round(time.time()-t0,1)}
            if kmodJ > 0:
                # save basis of K_full (mod p, signed) for later
                with open(f'/tmp/hodge/results/dense-rational/e1kernel_d{d}_basis.json', 'w') as f:
                    json.dump([{str(k): v for k, v in b.items()} for b in r['basis']], f)
        else:
            e = 18 - d
            # rank of [A]: R_e^- -> R_{e+6}^+  (parD = 1 for the source)
            r = macaulay_rank_and_null(F, A, e, 1, p, want_null=False)
            dimJ_target = S[e + 6][0] - Rh[e + 6][0]
            rankA = r['rank'] - dimJ_target
            ker = Rh[d][0] - rankA
            print(f"[dual d={d} via e={e}] rows={r['rows']} cols={r['cols']} rank={r['rank']} dimJ_{e+6}^+={dimJ_target} rank[A]={rankA} dimR_d^+={Rh[d][0]} kernel_mod_J={ker} ({time.time()-t0:.1f}s)", flush=True)
            rec = {'mode': 'dual', 'd': d, 'e': e, 'rankA': rankA, 'dimR': Rh[d][0], 'kernel_mod_J': ker, 'rows': r['rows'], 'cols': r['cols'], 'seconds': round(time.time()-t0,1)}
        out.append(rec)
        with open('/tmp/hodge/results/dense-rational/e1kernel_sweep.jsonl', 'a') as f:
            f.write(json.dumps(rec) + '\n')
