"""Sparse modular Gaussian elimination (ctypes wrapper around csrc/sge.c).

The C core does structured Gaussian elimination mod p with Markowitz-style pivoting on a CSR matrix
whose last column may hold the right-hand side.  It stops when the active submatrix is small or starts
to fill in; the Python side finishes that remainder densely with python-flint (random row compression
first when it is much taller than wide), and back-substitutes through the sparse pivots.

Public API
----------
build()                                  -> path of the compiled shared object (compiled on first use)
SparseElim(M_csr, p, rhs_col=-1)         -> handle; .eliminate(...), .export_active(), .backsub(x), ...
solve_sparse_system(M_csr, p, with_rhs=True, ...) -> dict(consistent, rank, nullity, x, pivots, ...)
nullspace_sparse(M_csr, p, ...)          -> (nullity, basis as uint32 array nullity x ncols) (basis optional)

M_csr: scipy.sparse.csr_matrix with int64/uint32 data already reduced mod p (zeros allowed, removed here).
"""
from __future__ import annotations
import ctypes, hashlib, os, subprocess, sys, tempfile, time
import numpy as np
import scipy.sparse as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, 'csrc', 'sge.c')
_BUILD = os.path.join(_HERE, '_build')
_lib = None


def build(force=False):
    """Compile csrc/sge.c into hodge/_build/sge-<srchash>.so (atomic; safe under concurrent imports)."""
    src = open(_SRC, 'rb').read()
    h = hashlib.sha1(src).hexdigest()[:12]
    os.makedirs(_BUILD, exist_ok=True)
    so = os.path.join(_BUILD, f'sge-{h}.so')
    if os.path.exists(so) and not force:
        return so
    fd, tmp = tempfile.mkstemp(prefix='sge-', suffix='.so', dir=_BUILD)
    os.close(fd)
    cmd = ['gcc', '-O3', '-march=native', '-shared', '-fPIC', '-o', tmp, _SRC]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        os.unlink(tmp)
        raise RuntimeError('sge build failed: ' + r.stderr[-2000:])
    os.replace(tmp, so)
    return so


def _load():
    global _lib
    if _lib is not None:
        return _lib
    lib = ctypes.CDLL(build())
    P = ctypes.c_void_p
    i32 = ctypes.c_int32; i64 = ctypes.c_int64; u32 = ctypes.c_uint32
    lib.sge_create.restype = P
    lib.sge_create.argtypes = [i32, i32, P, P, P, u32, i32]
    lib.sge_free.argtypes = [P]
    lib.sge_eliminate.restype = ctypes.c_int
    lib.sge_eliminate.argtypes = [P, i32, i32, i64, i64, ctypes.c_int, ctypes.c_int]
    for name in ('sge_npiv', 'sge_active_rows', 'sge_active_cols'):
        getattr(lib, name).restype = i32; getattr(lib, name).argtypes = [P]
    for name in ('sge_active_nnz', 'sge_fill_events'):
        getattr(lib, name).restype = i64; getattr(lib, name).argtypes = [P]
    lib.sge_inconsistent.restype = ctypes.c_int; lib.sge_inconsistent.argtypes = [P]
    lib.sge_col_states.argtypes = [P, P, P]
    lib.sge_export_active.argtypes = [P, P, P, P, P]
    lib.sge_backsub.argtypes = [P, P, ctypes.c_int]
    lib.sge_backsub_many.argtypes = [P, P, i32, ctypes.c_int]
    _lib = lib
    return lib


def _ptr(a):
    return a.ctypes.data_as(ctypes.c_void_p)


class SparseElim:
    def __init__(self, M, p, rhs_col=-1):
        lib = _load()
        M = sp.csr_matrix(M)
        M.sum_duplicates(); M.sort_indices()
        data = np.asarray(M.data, dtype=np.int64) % p
        M = sp.csr_matrix((data, M.indices, M.indptr), shape=M.shape)
        M.eliminate_zeros(); M.sort_indices()
        self.shape = M.shape
        self.p = int(p); self.rhs_col = int(rhs_col)
        self._indptr = np.ascontiguousarray(M.indptr, dtype=np.int64)
        self._indices = np.ascontiguousarray(M.indices, dtype=np.int32)
        self._vals = np.ascontiguousarray(M.data, dtype=np.uint32)
        self.h = lib.sge_create(M.shape[0], M.shape[1], _ptr(self._indptr), _ptr(self._indices), _ptr(self._vals),
                                self.p, self.rhs_col)
        self.lib = lib
        self.nnz0 = int(M.nnz)

    def __del__(self):
        try:
            if getattr(self, 'h', None):
                self.lib.sge_free(self.h); self.h = None
        except Exception:
            pass

    def eliminate(self, max_pivots=1 << 30, stop_cols=0, max_cost=1 << 40, max_nnz=1 << 40, cands=8,
                  stop_on_inconsistent=False):
        """Returns status: -1 inconsistent (only with stop_on_inconsistent), 0 finished, 1 stop_cols reached,
        2 fill limit, 3 nnz limit, 4 max_pivots.  Without stop_on_inconsistent an inconsistency only sets the
        sticky .inconsistent flag; the rhs-only rows stay active and show up in the dense remainder."""
        return self.lib.sge_eliminate(self.h, int(max_pivots), int(stop_cols), int(max_cost), int(max_nnz), int(cands),
                                      1 if stop_on_inconsistent else 0)

    @property
    def npiv(self): return self.lib.sge_npiv(self.h)
    @property
    def active_rows(self): return self.lib.sge_active_rows(self.h)
    @property
    def active_cols(self): return self.lib.sge_active_cols(self.h)
    @property
    def active_nnz(self): return self.lib.sge_active_nnz(self.h)
    @property
    def fill_events(self): return self.lib.sge_fill_events(self.h)
    @property
    def inconsistent(self): return bool(self.lib.sge_inconsistent(self.h))

    def col_states(self):
        """(states uint8[ncols]: 0 active, 1 pivot, 2 free, 3 rhs;  pivot columns in elimination order)"""
        st = np.zeros(self.shape[1], dtype=np.uint8)
        piv = np.zeros(max(self.npiv, 1), dtype=np.int32)
        self.lib.sge_col_states(self.h, _ptr(st), _ptr(piv))
        return st, piv[:self.npiv]

    def export_active(self):
        """Active submatrix as CSR with original column ids, plus the original row ids."""
        nr = self.active_rows; nnz = self.active_nnz
        indptr = np.zeros(nr + 1, dtype=np.int64)
        indices = np.zeros(max(nnz, 1), dtype=np.int32); vals = np.zeros(max(nnz, 1), dtype=np.uint32)
        rowids = np.zeros(max(nr, 1), dtype=np.int32)
        self.lib.sge_export_active(self.h, _ptr(indptr), _ptr(indices), _ptr(vals), _ptr(rowids))
        M = sp.csr_matrix((vals[:nnz].astype(np.int64), indices[:nnz], indptr), shape=(nr, self.shape[1]))
        return M, rowids[:nr]

    def backsub(self, x, use_rhs=True):
        x = np.ascontiguousarray(x, dtype=np.uint32)
        self.lib.sge_backsub(self.h, _ptr(x), 1 if use_rhs else 0)
        return x

    def backsub_many(self, X, use_rhs=False):
        X = np.ascontiguousarray(X, dtype=np.uint32)
        self.lib.sge_backsub_many(self.h, _ptr(X), X.shape[0], 1 if use_rhs else 0)
        return X


# ----------------------------------------------------------------------------
# modular helpers (numpy, exact)
# ----------------------------------------------------------------------------

def matvec_mod(M, x, p):
    """(M @ x) mod p for CSR int64 M (entries in [0,p)) and x in [0,p); exact via 16-bit splitting."""
    x = np.asarray(x, dtype=np.int64) % p
    lo = x & 0xFFFF; hi = x >> 16
    out = (M @ lo) % p
    out = (out + ((M @ hi) % p) * (1 << 16)) % p
    return out


def matmat_mod(M, X, p):
    """(M @ X) mod p, X dense int64 (n x k) with entries in [0,p).  Row weights must be < 2^16."""
    X = np.asarray(X, dtype=np.int64) % p
    lo = X & 0xFFFF; hi = X >> 16
    out = (M @ lo) % p
    out = (out + ((M @ hi) % p) * (1 << 16)) % p
    return out


def compress_rows(M, p, k, seed=1, chunk=4096):
    """(k x ncols) int64 = (R @ M) mod p for a dense random R (entries in [0,p)), generated in row
    chunks of M so memory stays ~ k*chunk*8 bytes.  Exact (16-bit splitting; column weights < 2^16)."""
    rng = np.random.default_rng(seed)
    nrows, ncols = M.shape
    out = np.zeros((ncols, k), dtype=np.int64)
    MT = M.T.tocsc()  # ncols x nrows, column-slicing by original row ranges is cheap on CSC
    for a in range(0, nrows, chunk):
        b = min(nrows, a + chunk)
        R = rng.integers(0, p, size=(b - a, k), dtype=np.int64)
        blk = MT[:, a:b].tocsr()
        lo = (blk @ (R & 0xFFFF)) % p
        hi = (blk @ (R >> 16)) % p
        out = (out + lo + (hi * (1 << 16)) % p) % p
    return np.ascontiguousarray(out.T)


def _check_colweights(M):
    w = np.diff(M.tocsc().indptr)
    if w.size and w.max() >= (1 << 16):
        raise ValueError(f'column weight {int(w.max())} too large for the 16-bit split (max 65535)')


def dense_finish(M, p, with_rhs=True, extra=40, seed=1, want_nullspace=False, verbose=False):
    """Rank/consistency/particular solution of a (moderate) CSR system densely with flint.

    M: nr x nc CSR with entries in [0,p); if with_rhs the last column is the rhs.  When nr > nc+extra the
    rows are first compressed by a random dense combination (rank preserved with prob >= 1 - nc/p^extra).
    Returns dict(consistent, rank (of the coefficient part, i.e. without the rhs column), rank_aug,
    x (length nc-1 if with_rhs else nc; free vars 0), pivots, nullspace)."""
    from flint import nmod_mat
    nr, nc = M.shape
    ncol_unk = nc - 1 if with_rhs else nc
    if nr == 0:
        x = np.zeros(ncol_unk, dtype=np.int64)
        ns = np.eye(ncol_unk, dtype=np.int64) if want_nullspace else None
        return {'consistent': True, 'rank': 0, 'x': x, 'pivots': [], 'nullspace': ns, 'compressed': False, 'k': 0}
    t0 = time.time()
    if nr > ncol_unk + extra:
        _check_colweights(M)
        k = ncol_unk + extra
        Dm = compress_rows(M, p, k, seed)
        compressed = True
    else:
        k = nr
        Dm = M.toarray().astype(np.int64) % p
        compressed = False
    t1 = time.time()
    Mn = nmod_mat(k, nc, Dm.ravel().tolist(), p)
    Mr, rank = Mn.rref()
    t2 = time.time()
    # pivot columns: first nonzero of each of the `rank` rows (scan rows as python lists for speed)
    pivots = []
    rows_l = Mr.tolist() if hasattr(Mr, 'tolist') else None
    C = ncol_unk
    for i in range(rank):
        if rows_l is not None:
            row = rows_l[i]
            j = 0
            while j < nc and int(row[j]) == 0:
                j += 1
        else:
            j = 0
            while j < nc and int(Mr[i, j]) == 0:
                j += 1
        pivots.append(j)
    consistent = not (with_rhs and pivots and pivots[-1] == C)
    rank_aug = rank
    if not consistent:
        rank -= 1                    # rank of the coefficient part (the last pivot sits in the rhs column)
        pivots = pivots[:-1]
    x = np.zeros(C, dtype=np.int64)
    if consistent and with_rhs:
        for i, j in enumerate(pivots):
            x[j] = int(rows_l[i][C]) if rows_l is not None else int(Mr[i, C])
    ns = None
    if want_nullspace:
        pivset = set(pivots)
        free = [j for j in range(C) if j not in pivset]
        ns = np.zeros((len(free), C), dtype=np.int64)
        for t, f in enumerate(free):
            ns[t, f] = 1
            for i, j in enumerate(pivots):
                v = int(rows_l[i][f]) if rows_l is not None else int(Mr[i, f])
                if v:
                    ns[t, j] = (-v) % p
    if verbose:
        print(f'[dense_finish] {nr}x{nc} -> k={k} compressed={compressed} rank={rank} consistent={consistent} '
              f'(compress {t1-t0:.1f}s, rref {t2-t1:.1f}s)', file=sys.stderr)
    return {'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'x': x, 'pivots': pivots, 'nullspace': ns,
            'compressed': compressed, 'k': k, 't_compress': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}


def solve_sparse_system(M, p, with_rhs=True, stop_cols=400, max_cost=None, max_nnz=None, cands=8,
                        extra=40, seed=1, verbose=False, want_nullspace=False, verify=True, need_rank=True):
    """Solve M x = b (last column of M = b when with_rhs) mod p by sparse elimination + dense finish.

    Returns dict(consistent, rank (of the coefficient matrix), nullity, x (list, free vars 0), pivots (column
    basis), nullspace (uint32 array or None), stats...).  `verify` re-checks M x = b with an exact sparse
    matvec.  With need_rank=False an inconsistency stops the elimination at once (rank=None, only
    rank_lower_bound); otherwise the elimination continues and the true rank/nullity are reported."""
    t0 = time.time()
    M = sp.csr_matrix(M)
    nr, nc = M.shape
    C = nc - 1 if with_rhs else nc
    nnz0 = M.nnz
    if max_nnz is None:
        max_nnz = max(50 * nnz0, 2 * 10**8)
    if max_cost is None:
        max_cost = 1 << 40
    E = SparseElim(M, p, rhs_col=C if with_rhs else -1)
    status = E.eliminate(stop_cols=stop_cols, max_cost=max_cost, max_nnz=max_nnz, cands=cands,
                         stop_on_inconsistent=not need_rank)
    t1 = time.time()
    st, piv_sparse = E.col_states()
    npiv = int(E.npiv)
    stats = {'status_sparse': status, 'npiv_sparse': npiv, 'active_rows': int(E.active_rows),
             'active_cols': int(E.active_cols), 'active_nnz': int(E.active_nnz), 'fill_events': int(E.fill_events),
             'nnz0': int(nnz0), 't_sparse': round(t1 - t0, 2)}
    if verbose:
        print(f'[sparse] {nr}x{nc} nnz={nnz0}: status={status} pivots={npiv} active {stats["active_rows"]}x'
              f'{stats["active_cols"]} nnz={stats["active_nnz"]} fill={stats["fill_events"]} ({t1-t0:.1f}s)', file=sys.stderr)
    stats['inconsistent_sparse'] = bool(E.inconsistent)
    if status == -1:
        return {'consistent': False, 'rank': None, 'nullity': None, 'x': None, 'pivots': None, 'nullspace': None,
                'rank_lower_bound': npiv, **stats, 't_total': round(time.time() - t0, 2)}
    # dense remainder on the active columns (+ rhs)
    Ma, rowids = E.export_active()
    active = np.nonzero(st == 0)[0]
    amap = np.full(nc, -1, dtype=np.int64)
    amap[active] = np.arange(len(active))
    ncols_d = len(active) + (1 if with_rhs else 0)
    if with_rhs:
        amap[C] = len(active)
    Ma = Ma.tocoo()
    cols_d = amap[Ma.col]
    assert (cols_d >= 0).all(), 'active row references a non-active column'
    Md = sp.csr_matrix((Ma.data, (Ma.row, cols_d)), shape=(Ma.shape[0], ncols_d))
    dres = dense_finish(Md, p, with_rhs=with_rhs, extra=extra, seed=seed, want_nullspace=want_nullspace, verbose=verbose)
    t2 = time.time()
    rank = npiv + dres['rank']
    consistent = dres['consistent'] and not E.inconsistent
    assert not (E.inconsistent and dres['consistent']), 'sparse phase found an rhs-only row but the dense finish is consistent'
    out = {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'x': None, 'pivots': None, 'nullspace': None,
           **stats, 'rank_dense': dres['rank'], 'k_dense': dres['k'], 't_dense': round(t2 - t1, 2)}
    pivots = [int(c) for c in piv_sparse] + [int(active[j]) for j in dres['pivots'] if j < len(active)]
    out['pivots'] = pivots
    if consistent:
        x = np.zeros(nc, dtype=np.uint32)
        if with_rhs:
            x[active] = dres['x'].astype(np.uint32)
            x[C] = 0
        else:
            x[active] = dres['x'].astype(np.uint32)
        E.backsub(x, use_rhs=with_rhs)
        xs = x[:C].astype(np.int64)
        if verify and with_rhs:
            r = matvec_mod(M[:, :C], xs, p)
            b = np.zeros(nr, dtype=np.int64)
            bcol = M[:, C].toarray().ravel().astype(np.int64) % p
            ok = bool(np.all(r == bcol))
            out['verified'] = ok
            if not ok:
                out['consistent'] = None
                out['note'] = 'particular solution failed M x = b re-check (compression/elimination bug?)'
        out['x'] = xs.tolist()
    if want_nullspace and dres['nullspace'] is not None:
        free_retired = np.nonzero(st == 2)[0]
        nd = dres['nullspace'].shape[0]
        NS = np.zeros((nd + len(free_retired), nc), dtype=np.uint32)
        if nd:
            NS[:nd, active] = dres['nullspace'].astype(np.uint32)
        for t, f in enumerate(free_retired):
            NS[nd + t, f] = 1
        E.backsub_many(NS, use_rhs=False)
        out['nullspace'] = NS[:, :C]
        if verify:
            ok = True
            for t in range(min(NS.shape[0], 5)):
                ok &= bool(np.all(matvec_mod(M[:, :C], NS[t, :C].astype(np.int64), p) == 0))
            out['nullspace_verified'] = ok
    out['t_total'] = round(time.time() - t0, 2)
    del E
    return out
