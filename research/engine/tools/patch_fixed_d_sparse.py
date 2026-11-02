"""Idempotent: add set_threads / chunked _compress / solve_sparse to hodge/fixed_d.py."""
import sys
P = '/home/kejace/devel/sextics-work/hodge/fixed_d.py'
src = open(P).read()
if 'def solve_sparse' in src:
    print('fixed_d sparse: already patched'); sys.exit(0)
if 'import os as _os' not in src:
    src = src.replace('from __future__ import annotations\nimport sys, time\n', 'from __future__ import annotations\nimport os as _os\nimport sys, time\n', 1)
old = '''try:
    import flint
    from flint import nmod_mat, fmpq_mat, fmpz_mat
except ImportError:  # pragma: no cover
    flint = None
'''
new = old + '''

def set_threads(n=None):
    """Set the number of FLINT threads used by dense nmod_mat rref/nullspace (4x at 8 threads for n=4000).

    Default: FLINT_NUM_THREADS from the environment (rjob submit --cpus N sets it), else 1.  Returns the value."""
    if n is None:
        try:
            n = int(_os.environ.get('FLINT_NUM_THREADS', '1'))
        except ValueError:
            n = 1
    n = max(1, int(n))
    if flint is not None:
        flint.ctx.threads = n
    return n


THREADS = set_threads()
'''
assert old in src, 'flint import block not found'; src = src.replace(old, new)
old = '''def _compress(M, p, k, seed=1):
    """Return (k x ncols) int64 array = (R @ M) mod p for random R with entries in [0,p)."""
    rng = np.random.default_rng(seed)
    nrows = M.shape[0]
    R = rng.integers(0, p, size=(k, nrows), dtype=np.int64)
    R_lo = R & 0xFFFF
    R_hi = R >> 16
    MT = M.T.tocsr()  # (ncols x nrows)
    # (MT @ R^T)^T = R @ M ; scipy sparse @ dense -> dense
    lo = (MT @ R_lo.T) % p
    hi = (MT @ R_hi.T) % p
    out = (lo + (hi * (1 << 16)) % p) % p
    return np.ascontiguousarray(out.T)  # k x ncols
'''
new = '''def _compress(M, p, k, seed=1, chunk=4096):
    """Return (k x ncols) int64 array = (R @ M) mod p for random R with entries in [0,p).

    R is generated in chunks of `chunk` rows of M so memory stays ~ k*chunk*8 bytes (the old
    one-shot version needed k*nrows*8 bytes: 54 GB at d=12).  Exact: 16-bit split, needs column
    weights < 2^16 (checked)."""
    from .sge import compress_rows, _check_colweights
    _check_colweights(M)
    return compress_rows(M, p, k, seed=seed, chunk=chunk)
'''
assert old in src, '_compress not found'; src = src.replace(old, new)
old = '''def exact_solution_from_pivots(F, A, q, D, pivots, verbose=False, seed=7):'''
new = '''def solve_sparse(S, stop_cols=400, verbose=False, want_nullspace=False, seed=1, **kw):
    """Consistency + rank + particular solution mod p via sparse elimination (hodge.sge) with a dense finish.

    Same result schema as solve_compressed (consistent, rank, nullity, x, pivots, rows, cols, timings) plus
    the sparse-phase statistics.  `x` is verified by an exact sparse matvec (key 'verified'); on failure
    consistent=None.  With want_nullspace the full nullspace basis (uint32 array nullity x cols) is returned.
    Extra keyword arguments go to sge.solve_sparse_system (need_rank, max_nnz, cands, extra, ...)."""
    from . import sge
    t0 = time.time()
    M = _sparse_matrix(S)
    r = sge.solve_sparse_system(M, S.p, with_rhs=True, stop_cols=stop_cols, seed=seed, verbose=verbose,
                                want_nullspace=want_nullspace, **kw)
    r.update({'rows': S.nrows, 'cols': S.ncols, 't_matrix': round(time.time() - t0 - r.get('t_total', 0), 2)})
    if verbose:
        print(f'[sparse p={S.p}] consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} '
              f'(sparse {r["t_sparse"]}s, dense {r.get("t_dense")}s, total {r["t_total"]}s)', file=sys.stderr)
    return r


''' + old
assert old in src, 'exact_solution_from_pivots not found'; src = src.replace(old, new)
open(P, 'w').write(src); print('fixed_d sparse: patched')
