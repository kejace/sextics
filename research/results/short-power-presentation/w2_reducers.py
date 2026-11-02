"""Build and save Reducer objects (rref of J(F)_deg^par mod p) as .npz for reuse.  Usage: python3 w2_reducers.py p deg:par [deg:par ...]"""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

def save_reducer(red, path):
    n = len(red.mons); r = red.rank
    Rm = np.zeros((r, n), dtype=np.uint32)
    if r:
        Rm[:] = np.array([int(red.R[i, j]) for i in range(r) for j in range(n)], dtype=np.uint32).reshape(r, n)
    np.savez_compressed(path, R=Rm, pivots=np.array(red.pivots, dtype=np.int64), mons=np.array(red.mons, dtype=np.int64).reshape(n, 6),
                        p=red.p, deg=red.deg, par=red.par)

def load_reducer(path, p=None):
    z = np.load(path)
    class R: pass
    red = R()
    red.p = int(z['p']); red.deg = int(z['deg']); red.par = int(z['par'])
    red.mons = [tuple(int(x) for x in row) for row in z['mons']]
    red.idx = {m: i for i, m in enumerate(red.mons)}
    red.pivots = [int(x) for x in z['pivots']]
    red.rank = len(red.pivots)
    Rm = z['R']
    red.Rnp = Rm.astype(np.int64)
    pset = set(red.pivots)
    red.nonpiv = [j for j in range(len(red.mons)) if j not in pset]
    red.standard = [red.mons[j] for j in red.nonpiv]
    return red

def nf_rows_np(red, polys):
    """normal forms (numpy int64 rows mod p) of polys using the saved reducer; uses blocked modular matmul."""
    p = red.p; n = len(red.mons); k = len(polys)
    U = np.zeros((k, n), dtype=np.int64)
    for i, P in enumerate(polys):
        for e, v in P.items(): U[i, red.idx[e]] = int(v) % p
    if red.rank == 0: return U
    Up = U[:, red.pivots]                     # k x r
    # U - Up @ R  (mod p), blocked to avoid overflow: entries < 2^31, products < 2^62; sum over r terms overflows -> block r
    out = U.copy()
    r = red.rank
    bs = max(1, (2**63 - 1) // (p * p) // 2)   # ~ 2^63 / 2^62 = 2 ... too small; use float-free splitting instead
    # split R into 16-bit halves: R = R1*2^16 + R0; Up < 2^31. Up*R0 < 2^47, sum over up to 2^15 terms < 2^62.
    R0 = red.Rnp & 0xFFFF; R1 = red.Rnp >> 16
    Up = Up % p
    for s in range(0, r, 32768):
        e = min(r, s + 32768)
        A0 = (Up[:, s:e] @ R0[s:e]) % p
        A1 = (Up[:, s:e] @ R1[s:e]) % p
        out = (out - A0 - (A1 * 65536) % p) % p
    return out % p

if __name__ == '__main__':
    p = int(sys.argv[1])
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p)
    for spec in sys.argv[2:]:
        deg, par = map(int, spec.split(':'))
        path = f'{RES}/w2/reducer_deg{deg}_par{par}_mod{p}.npz'
        if os.path.exists(path): print('exists', path); continue
        t0 = time.time()
        red = Reducer(Fp, deg, par, p, verbose=True)
        save_reducer(red, path)
        wlog(batch='reducer', degree=deg, parity=par, p=p, dim_S=len(red.mons), rank_J=red.rank, dim_R=len(red.mons)-red.rank, seconds=round(time.time()-t0,1))
        print('saved', path, f'{time.time()-t0:.0f}s', flush=True)
        del red
