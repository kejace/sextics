"""Blocked modular Gaussian elimination (rank-revealing LU) over F_p using float64 BLAS.
Exactness: block size b and prime p must satisfy b * p^2 < 2^53.  Default p < 2^22, b = 128.
mod_echelon(M, p) -> (rank, pivot_columns).  Consistency of [M | b]: last column is not a pivot.
"""
import numpy as np

P22 = 4194301   # prime < 2^22

def _is_prime(n):
    if n < 2: return False
    for q in range(2, int(n ** 0.5) + 1):
        if n % q == 0: return False
    return True
assert _is_prime(P22)

def _red(X, p, invp):
    """X - floor(X/p)*p, exact for |X| < 2^53."""
    return X - np.floor(X * invp) * p

def mod_echelon(M, p=P22, b=128):
    assert b * p * p < 2 ** 53, 'block size / prime too large for exact float64 arithmetic'
    A = np.array(M, dtype=np.float64, order='C')
    invp = 1.0 / p
    n, m = A.shape
    r = 0
    pivcols = []
    c0 = 0
    while c0 < m and r < n:
        c1 = min(c0 + b, m)
        w = c1 - c0
        nr = n - r
        panel = np.ascontiguousarray(A[r:, c0:c1])   # contiguous copy (nr x w)
        perm = np.arange(nr)
        k = 0
        L = np.zeros((nr, w))                          # multipliers, column per pivot
        for j in range(w):
            col = panel[k:, j]
            nz = np.flatnonzero(col)
            if len(nz) == 0:
                continue
            i = k + int(nz[0])
            if i != k:
                panel[[k, i]] = panel[[i, k]]
                perm[[k, i]] = perm[[i, k]]
                L[[k, i]] = L[[i, k]]
            inv = pow(int(panel[k, j]), -1, p)
            mult = _red(panel[k + 1:, j] * inv, p, invp)
            if j + 1 < w:
                panel[k + 1:, j + 1:] = _red(panel[k + 1:, j + 1:] - np.outer(mult, panel[k, j + 1:]), p, invp)
            panel[k + 1:, j] = 0
            L[k + 1:, k] = mult
            pivcols.append(c0 + j)
            k += 1
            if k == nr:
                break
        if k > 0:
            # write back the panel (permuted rows) -- needed only for completeness
            A[r:, c0:c1] = panel
            if c1 < m:
                T = A[r:, c1:][perm]                       # permuted copy (nr x (m-c1))
                Lk = L[:, :k]
                # inverse of unit lower triangular L11 (k x k) mod p
                L11 = Lk[:k]
                Linv = np.eye(k)
                for j in range(k):
                    if j + 1 < k:
                        Linv[j + 1:] = _red(Linv[j + 1:] - np.outer(L11[j + 1:, j], Linv[j]), p, invp)
                U12 = _red(Linv @ T[:k], p, invp)
                T[:k] = U12
                if k < nr:
                    T[k:] = _red(T[k:] - Lk[k:] @ U12, p, invp)
                A[r:, c1:] = T
            # rows of A[r:, :c0] are already zero
        r += k
        c0 = c1
    return r, pivcols

def consistent_rank(M_aug, p=P22, b=128):
    n, m1 = M_aug.shape
    rank_aug, piv = mod_echelon(M_aug, p, b)
    inconsistent = (len(piv) > 0 and piv[-1] == m1 - 1)
    rank = rank_aug - (1 if inconsistent else 0)
    return (not inconsistent), rank, rank_aug

if __name__ == '__main__':
    import time
    from flint import nmod_mat
    rng = np.random.default_rng(1)
    p = P22
    for (n, m, rk) in [(50, 30, 20), (300, 200, 150), (1000, 700, 600), (2500, 1200, 900)]:
        X = rng.integers(0, 1 << 10, size=(n, rk), dtype=np.int64); Y = rng.integers(0, 1 << 10, size=(rk, m), dtype=np.int64)
        M = (X @ Y) % p
        M[:, 5] = 0; M[:, 7] = M[:, 6]
        x = rng.integers(0, 1 << 10, size=(m,), dtype=np.int64)
        b_cons = (M @ x) % p
        b_inc = rng.integers(0, p, size=(n,), dtype=np.int64)
        for b, name in [(b_cons, 'consistent'), (b_inc, 'random')]:
            Maug = np.concatenate([M, b[:, None]], axis=1)
            t0 = time.time(); c, r, ra = consistent_rank(Maug, p); t1 = time.time()
            if n <= 1000:
                Mn = nmod_mat(n, m + 1, Maug.ravel().tolist(), p); ra_f = Mn.rank()
                M0 = nmod_mat(n, m, M.ravel().tolist(), p); r_f = M0.rank()
                ok = (c == (r_f == ra_f)) and r == r_f and ra == ra_f
            else:
                r_f = ra_f = None; ok = (c == (name == 'consistent')) and r == rk
            print(f'n={n} m={m} rk={rk}: {name}: ours consistent={c} rank={r} rank_aug={ra} | flint rank={r_f} rank_aug={ra_f} ({t1-t0:.2f}s)  OK={ok}', flush=True)
    n, m, rk = 10500, 4061, 3000
    X = rng.integers(0, 1 << 10, size=(n, rk), dtype=np.int64); Y = rng.integers(0, 1 << 10, size=(rk, m), dtype=np.int64)
    M = (X @ Y) % p
    t0 = time.time(); r, piv = mod_echelon(M, p); print(f'speed: {n}x{m} rank {r} (expected {rk}) in {time.time()-t0:.1f}s', flush=True)
