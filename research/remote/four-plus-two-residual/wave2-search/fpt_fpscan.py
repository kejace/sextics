"""Exact exhaustive scan over P^(K-1)(F_p) of a small E1-admissible sector family: for every F_p-point a, test the
Koszul-form certificate system  M(a) x_C - MQ x_Q = -v(a)  for consistency mod p (flint rank).  A rational witness D in
the family would reduce to an F_p-point whose system is consistent mod p (unless p is a bad prime for that point,
probability ~1/p), so 'no consistent non-spurious point' is a rigorous mod-p exclusion of the whole family.
usage: python3 fpt_fpscan.py <d> <h,r2,s|auto> <p> [maxK=4] [procs=8] [family=kernelJ|vars:...|prod:...]
Output: wave2-search/fpscan.jsonl (one line per (d, sector, p))."""
import os, sys, json, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
import fpt_msolve as fm
import numpy as np, scipy.sparse as sp
from flint import nmod_mat
import multiprocessing as mp

d = int(sys.argv[1]); labarg = sys.argv[2]; p = int(sys.argv[3])
maxK = int(sys.argv[4]) if len(sys.argv) > 4 else 4
PROCS = int(sys.argv[5]) if len(sys.argv) > 5 else 8
OUT = f'{LANE}/wave2-search/fpscan.jsonl'
inst = Instance()
G = {}


def projective_points(K, p, chunk, nchunks):
    """Points of P^(K-1)(F_p) with first nonzero coordinate 1, split into nchunks by a hash of the point."""
    import itertools
    idx = 0
    for lead in range(K):
        for rest in itertools.product(range(p), repeat=K - 1 - lead):
            if idx % nchunks == chunk:
                yield (0,) * lead + (1,) + rest
            idx += 1


def setup(pieces):
    K0, K, nN, nC, nQ, R = pieces['K0'], pieces['K'], pieces['nN'], pieces['nC'], pieces['nQ'], pieces['R']
    B = np.array(pieces['B'], dtype=np.int64)            # K0 x K
    N = np.array(pieces['Nmat'], dtype=np.int64)         # nN x K
    T = sp.csr_matrix((np.array(pieces['tvals'], dtype=np.int64), (np.array(pieces['tri']), np.array(pieces['tci']))), shape=(R, K0 * nN))
    Mj = np.zeros((K0, R, nC), dtype=np.int64)
    for j in range(K0):
        for r, c, v in pieces['Mj_entries'][j]:
            Mj[j, r, c] = (Mj[j, r, c] + v) % p
    MQ = np.zeros((R, nQ), dtype=np.int64)
    for r, c, v in pieces['MQ_entries']:
        MQ[r, c] = (MQ[r, c] + v) % p
    # F-multiples subspace in Dm coordinates (to flag spurious points)
    Dm = pieces['Dm']; didx = {m: j for j, m in enumerate(Dm)}
    Gm = inst.sector_monomials(d - 6, G['lab']) if d >= 6 else []
    FM = None
    if Gm:
        Fp = reduce_mod(inst.F, p)
        FM = [[0] * len(Gm) for _ in range(K0)]
        for jj, g in enumerate(Gm):
            for e, v in shift(Fp, g, 1, p).items():
                FM[didx[e]][jj] = int(v) % p
    G.update({'B': B, 'N': N, 'T': T, 'Mj': Mj, 'MQ': MQ, 'K0': K0, 'K': K, 'nN': nN, 'nC': nC, 'nQ': nQ, 'R': R, 'FM': FM})


def is_fmultiple(c):
    FM = G['FM']
    if FM is None:
        return False
    K0 = G['K0']; ng = len(FM[0])
    aug = fm.nm([FM[j] + [int(c[j])] for j in range(K0)], p, R=K0, C=ng + 1)
    base = fm.nm(FM, p, R=K0, C=ng)
    return aug.rank() == base.rank()


def test_point(a):
    a = np.array(a, dtype=np.int64)
    c = (G['B'] @ a) % p
    if not c.any():
        return None
    n = (G['N'] @ a) % p
    v = (G['T'] @ np.outer(c, n).ravel()) % p
    M = np.tensordot(c, G['Mj'], axes=(0, 0)) % p           # R x nC
    R = G['R']
    Aug = np.concatenate([M, G['MQ'], ((-v) % p)[:, None]], axis=1)
    C = Aug.shape[1]
    flat = Aug.ravel().tolist()
    Mr, rank = nmod_mat(R, C, flat, p).rref()
    # consistent iff no pivot in the last column
    pivot_last = False
    for i in range(rank):
        j = 0
        while j < C and int(Mr[i, j]) == 0:
            j += 1
        if j == C - 1:
            pivot_last = True; break
    return (not pivot_last), rank, c


def worker(args):
    chunk, nchunks = args
    K = G['K']
    n = 0; cons = []; t0 = time.time()
    for a in projective_points(K, p, chunk, nchunks):
        res = test_point(a)
        if res is None:
            continue
        n += 1
        ok, rank, c = res
        if ok:
            spur = is_fmultiple(c)
            cons.append({'a': [int(x) for x in a], 'rank': rank, 'spurious': bool(spur)})
    return chunk, n, cons, time.time() - t0


if __name__ == '__main__':
    labs = sorted(inst.labels_in_degree(d)) if labarg == 'auto' else [tuple(int(x) for x in labarg.split(','))]
    for lab in labs:
        G.clear(); G['lab'] = lab
        text, info = fm.build_system(inst, d, lab, p, verbose=False)
        K = info.get('K', 0)
        if K == 0 or text is None:
            continue
        if K > maxK:
            log_line(OUT, {'kind': 'fpscan', 'd': d, 'sector': list(lab), 'p': p, 'K': K, 'status': f'skipped (K > {maxK})'})
            print(f'd={d} lab={lab}: K={K} > maxK, skipped', flush=True); continue
        setup(fm.PIECES)
        npts = sum(p ** (K - 1 - l) for l in range(K))
        print(f'd={d} lab={lab}: K={K} K0={G["K0"]} system {G["R"]}x{G["nC"]+G["nQ"]} -> scanning {npts} points of P^{K-1}(F_{p}) with {PROCS} procs', flush=True)
        t0 = time.time()
        nchunks = PROCS * 4 if npts > 2000 else 1
        total = 0; consistent = []
        with mp.get_context('fork').Pool(min(PROCS, nchunks)) as pool:
            for chunk, n, cons, dt in pool.imap_unordered(worker, [(i, nchunks) for i in range(nchunks)]):
                total += n; consistent += cons
        nonspur = [c for c in consistent if not c['spurious']]
        rec = {'kind': 'fpscan', 'd': d, 'sector': list(lab), 'p': p, 'K': K, 'K0': G['K0'], 'rows': G['R'], 'cols': G['nC'] + G['nQ'],
               'points': total, 'consistent': len(consistent), 'consistent_nonspurious': len(nonspur), 'examples': nonspur[:5],
               'seconds': round(time.time() - t0, 1), 'engine': ENGINE}
        log_line(OUT, rec)
        print(json.dumps({k: rec[k] for k in rec if k != 'kind'}), flush=True)
        if nonspur:
            print('!!!!! CONSISTENT NON-SPURIOUS POINT(S) mod p:', nonspur[:3], flush=True)
