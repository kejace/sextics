"""Memory-lean variant of fixed_d.solve_fixed_D_fast: same FixedDSystem, but the random row
compression is done in chunks so peak memory stays ~1-2 GB at d=8 (the engine's _compress
allocates three dense k x nrows int64 arrays, ~7 GB at d=8, which OOM-killed the earlier batch)."""
import sys, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from flint import nmod_mat
from hodge import core, fixed_d
from hodge.fixed_d import FixedDSystem, _sparse_matrix

def compress_chunked(M, p, k, seed=1, chunk=512):
    rng = np.random.default_rng(seed)
    nrows, ncols = M.shape
    MT = M.T.tocsr()                       # ncols x nrows
    out = np.empty((k, ncols), dtype=np.int64)
    for s in range(0, k, chunk):
        kk = min(chunk, k - s)
        R = rng.integers(0, p, size=(kk, nrows), dtype=np.int64)
        lo = (MT @ (R & 0xFFFF).T) % p     # ncols x kk
        hi = (MT @ (R >> 16).T) % p
        out[s:s+kk, :] = ((lo + (hi * (1 << 16)) % p) % p).T
    return out

def solve_lean(F, A, q, D, p=None, verbose=True, k_extra=40, seed=1, want_exact=True):
    p = p or core.PRIMES[0]
    t0 = time.time()
    S = FixedDSystem(F, A, q, D, p, verbose=verbose)
    M = _sparse_matrix(S)
    C = S.ncols
    k = min(S.nrows, C + 1 + k_extra)
    t1 = time.time()
    Cm = compress_chunked(M, p, k, seed)
    del M
    t2 = time.time()
    Mn = nmod_mat(k, C + 1, Cm.ravel().tolist(), p)
    del Cm
    Mr, rank = Mn.rref()
    t3 = time.time()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    out = {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'rows': S.nrows, 'cols': C, 'k': k,
           't_build': round(t1-t0,1), 't_compress': round(t2-t1,1), 't_rref': round(t3-t2,1)}
    if verbose:
        print(f'[lean p={p}] k={k} consistent={consistent} rank={rank} nullity={C-rank} build {t1-t0:.1f}s compress {t2-t1:.1f}s rref {t3-t2:.1f}s', flush=True)
    if consistent and want_exact:
        for sd in range(3):
            sol = fixed_d.exact_solution_from_pivots(F, A, q, D, pivots, verbose=verbose, seed=sd + 7)
            if sol is not None:
                out['N'], out['Q'] = sol; out['exact'] = True; return out
        out['exact'] = False
    return out
