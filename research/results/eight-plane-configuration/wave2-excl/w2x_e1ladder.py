"""Independent recomputation of Ann_R(A)_d (E1 kernel modulo J) by the normal-form method, new prime.

Method (differs from wave-1 e1kernel.py, which took the nullspace of [A*m | J-columns]):
  1. rref of the generator matrix of J_{d+6}^{par'} (rows = m*F_i) -> pivots, standard monomials, rank_p J_{d+6};
  2. normal forms NF(A*m) for all D-monomials m of degree d, parity par, as NF = V_std - V_piv * Rref_std;
  3. Ann = left kernel of NF (nullspace of NF^T); kerdim = dim{D : A D in J};
  4. quotient by J_d^par via a second reducer at degree d -> kermodJ and a basis of Ann_R(A)_d.
Rigor: rank_p J_k <= dim_Q J_k always; equality is checked against the exact Koszul dimension (w2x_lib.J_dim).
Usage: python3 w2x_e1ladder.py d1 d2 ... [--p P] [--tag TAG]
Output: appends to wave2-excl/e1ladder.<p>.jsonl ; saves bases (mod p) to wave2-excl/annA_d<d>_par<par>.p<p>.pkl
"""
import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *

def submatrix(M, rows, cols, p):
    R, C = len(rows), len(cols)
    flat = [int(M[i, j]) for i in rows for j in cols]
    return nmod_mat(R, C, flat, p)

def e1_ladder(d, par, p, want_basis=True, verbose=True):
    t0 = time.time()
    Fp, Fi = partials(p)
    Ap = reduce_mod(A, p)
    parp = (parity_of(A) + par) % 2
    # 1. J_{d+6}^{parp} rref
    mons = monomials(d + 6, parp); idx = {m: i for i, m in enumerate(mons)}
    gens = J_generators(Fi, d + 6, parp, p)
    R, C = len(gens), len(mons)
    flat = [0] * (R * C)
    for r, g in enumerate(gens):
        for e, v in g.items(): flat[r * C + idx[e]] = v
    M = nmod_mat(R, C, flat, p); del flat
    Mr, rankJ = M.rref(); del M
    piv = []
    for i in range(rankJ):
        j = 0
        while j < C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    pivset = set(piv); std = [j for j in range(C) if j not in pivset]
    t1 = time.time()
    # 2. normal forms of A*m
    Dm = monomials(d, par); nD = len(Dm)
    Vflat = [0] * (nD * C)
    for r, m in enumerate(Dm):
        for e, v in shift(Ap, m, 1, p).items(): Vflat[r * C + idx[e]] = v
    V = nmod_mat(nD, C, Vflat, p); del Vflat
    Vstd = submatrix(V, range(nD), std, p)
    Vpiv = submatrix(V, range(nD), piv, p); del V
    Rstd = submatrix(Mr, range(rankJ), std, p); del Mr
    NF = Vstd - Vpiv * Rstd
    t2 = time.time()
    # 3. left kernel
    X, kerdim = NF.transpose().nullspace()   # columns of X: vectors x with NF^T x = 0, i.e. x^T NF = 0
    kernel = []
    for k in range(kerdim):
        v = {Dm[j]: int(X[j, k]) for j in range(nD) if int(X[j, k])}
        kernel.append(v)
    t3 = time.time()
    # 4. quotient by J_d^par
    rankJd = 0; kermodJ = kerdim; basis = []
    if d >= 5:
        red = Reducer(Fi, d, par, p)
        rankJd = red.rank
        # kernel vectors in standard coordinates of degree d
        if kerdim:
            Cd = len(red.mons)
            KF = [0] * (kerdim * Cd)
            for r, v in enumerate(kernel):
                for e, c in v.items(): KF[r * Cd + red.idx[e]] = c
            K = nmod_mat(kerdim, Cd, KF, p); del KF
            Kstd = submatrix(K, range(kerdim), red.standard, p)
            Kpiv = submatrix(K, range(kerdim), red.pivots, p)
            # rref rows restricted to standard cols
            RS = nmod_mat(red.rank, len(red.standard), [c for row in red.rows_std for c in row], p)
            NFk = Kstd - Kpiv * RS
            NFr, kermodJ = NFk.rref()
            if want_basis:
                for i in range(kermodJ):
                    basis.append({red.mons[red.standard[j]]: int(NFr[i, j]) for j in range(len(red.standard)) if int(NFr[i, j])})
    t4 = time.time()
    info = {'d': d, 'parD': par, 'p': p, 'method': 'normal-form', 'nD': nD, 'deg_target': d + 6, 'par_target': parp,
            'rank_p_J_target': rankJ, 'exact_J_target': J_dim(d + 6, parp), 'rank_p_J_d': rankJd, 'exact_J_d': J_dim(d, par),
            'kerdim': kerdim, 'kermodJ': kermodJ, 'exact_R_d': R_dim(d, par), 'exact_R_target': R_dim(d + 6, parp),
            'forced_lower_bound': max(0, R_dim(d, par) - R_dim(d + 6, parp)),
            't_rrefJ': round(t1 - t0, 1), 't_nf': round(t2 - t1, 1), 't_kernel': round(t3 - t2, 1), 't_quot': round(t4 - t3, 1)}
    info['hilbert_ok'] = (rankJ == info['exact_J_target']) and (rankJd == info['exact_J_d'])
    return info, basis

if __name__ == '__main__':
    args = sys.argv[1:]
    p = P_NEW; tag = ''
    if '--p' in args:
        i = args.index('--p'); p = int(args[i + 1]); del args[i:i + 2]
    if '--tag' in args:
        i = args.index('--tag'); tag = args[i + 1]; del args[i:i + 2]
    ds = [int(a) for a in args]
    for d in ds:
        for par in (0, 1):
            info, basis = e1_ladder(d, par, p)
            print(json.dumps(info), flush=True)
            with open(f'{PK}/e1ladder.{p}.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
            if basis:
                with open(f'{PK}/annA_d{d}_par{par}.p{p}.pkl', 'wb') as f: pickle.dump(basis, f)
            wlog({'exp': 'a2-e1ladder-nf', 'd': d, 'parD': par, 'prime': p, 'result': f"kermodJ={info['kermodJ']} (forced>={info['forced_lower_bound']}) hilbert_ok={info['hilbert_ok']}", 'rank': info['rank_p_J_target'], 'seconds': round(sum(info[k] for k in ('t_rrefJ','t_nf','t_kernel','t_quot')),1), 'job': tag or 'local', 'note': 'normal-form method, independent of wave-1 nullspace code'}, fn=f'{PK}/wave2-excl-remote-{tag}.jsonl' if tag else None)
