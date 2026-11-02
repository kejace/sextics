"""Koszul-reduced E2 system for Jacobian denominators D = w(F) = sum_i w_i F_i  (q = 1).

Derivation (exact, over Q; see README):
  E1 general solution: N = A*w + sum_{i<j} c_ij K_ij + (G/6) E,  K_ij = F_j d_i - F_i d_j, E = Euler field,
  deg c_ij = d-4, G in S_d (free).  The Euler part cancels identically in the certificate, so
  certificate exists  <=>  exists c, Q with
     P(w) + sum_{i<j} [ D*K_ij(c_ij) - c_ij*K_ij(D) ]  =  F * Q,     deg Q = 2d-6,
  where P(w) = D*div(Aw) - (Aw)(D) = D*(w(A) + A*div w) - A*w(D).
All arithmetic uses the integer form Fint = LAM*F (scale invariance of the certificate problem).
We build the integer matrix M (rows: monomials of degree 2d, optionally one parity) and rhs b = -P(w),
and decide consistency mod p AND report the exact structural data needed for rigor:
  rank_p(M), rank_p([M|b]), ncols; if rank_p(M) == ncols then rank_Q(M) == ncols and modular inconsistency
  is rigorous over Q.  Otherwise an exact rational kernel basis is attempted (ratrecon) and verified.
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl/vm')
from w2lib import *
from flint import nmod_mat, fmpz_mat, fmpq_mat
import numpy as np, scipy.sparse as sp

def field_apply(w, P):
    """w(P) = sum_i w_i * dP/dz_i  (exact ints/Fractions)."""
    out = {}
    for i in range(6):
        if w[i]:
            out = add(out, mul(w[i], pderiv(P, i)))
    return out

def divergence(w):
    out = {}
    for i in range(6):
        if w[i]:
            out = add(out, pderiv(w[i], i))
    return out

def BD(D, N):
    """D div N - N(D)."""
    return sub(mul(D, divergence(N)), field_apply(N, D))

def koszul_field(i, j, c):
    K = [{} for _ in range(6)]
    K[i] = mul(Fi_int[j], c)
    K[j] = scale(mul(Fi_int[i], c), -1)
    return K

class KoszulSystem:
    def __init__(self, w, parity_rows=None, verbose=False):
        """w: list of 6 polys (exact) of degree d-5; D = w(Fint)."""
        t0 = time.time()
        self.w = w
        D = field_apply(w, Fint)
        self.D = D
        d = core.degree(D)
        self.d = d
        pD = parity_of(D)
        self.pD = pD
        N0 = [mul(Aint, w[i]) for i in range(6)]
        self.N0 = N0
        degc = d - 4
        cols, colpolys = [], []
        # parity bookkeeping (only when D has pure parity): then N_i parity = pA + pD + (i>=3) ... we keep ALL
        # monomials unless parity_rows is given, in which case we filter c-monomials to the matching parity.
        for i in range(6):
            for j in range(i + 1, 6):
                for m in monomials(degc):
                    if parity_rows is not None and pD is not None:
                        # component i of c K_ij is c F_j ; needed parity of N_i: (pA - pD - (i>=3)) mod 2
                        need = (1 - pD - (0 if i < 3 else 1)) % 2
                        if (yparity(m) + (0 if j < 3 else 1)) % 2 != need:
                            continue
                    c = {m: 1}
                    cols.append(('c', i, j, m)); colpolys.append(BD(D, koszul_field(i, j, c)))
        Fneg = scale(Fint, -1)
        for m in monomials(2 * d - 6):
            if parity_rows is not None and yparity(m) != parity_rows:
                continue
            cols.append(('Q', m)); colpolys.append(shift(Fneg, m, 1))
        rhs = scale(BD(D, N0), -1)
        rows = {}
        for col in colpolys:
            for e in col:
                if parity_rows is not None and yparity(e) != parity_rows: continue
                rows.setdefault(e, len(rows))
        for e in rhs:
            if parity_rows is not None and yparity(e) != parity_rows: continue
            rows.setdefault(e, len(rows))
        self.cols, self.colpolys, self.rhs, self.rows = cols, colpolys, rhs, rows
        self.nrows, self.ncols = len(rows), len(cols)
        self.parity_rows = parity_rows
        self.t_build = time.time() - t0
        if verbose:
            print(f'[kred] d={d} rows={self.nrows} cols={self.ncols} build {self.t_build:.1f}s', file=sys.stderr, flush=True)

    def dense_mod(self, p, with_rhs=True):
        R, C = self.nrows, self.ncols
        W = C + (1 if with_rhs else 0)
        flat = [0] * (R * W)
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                if e in self.rows:
                    flat[self.rows[e] * W + j] = v % p
        if with_rhs:
            for e, v in self.rhs.items():
                if e in self.rows:
                    flat[self.rows[e] * W + C] = v % p
        return nmod_mat(R, W, flat, p)

    def dense_int(self, with_rhs=True):
        R, C = self.nrows, self.ncols
        W = C + (1 if with_rhs else 0)
        flat = [0] * (R * W)
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                if e in self.rows:
                    flat[self.rows[e] * W + j] = int(v)
        if with_rhs:
            for e, v in self.rhs.items():
                if e in self.rows:
                    flat[self.rows[e] * W + C] = int(v)
        return fmpz_mat(R, W, flat)

    def solve_mod(self, p, verbose=False):
        t0 = time.time()
        M = self.dense_mod(p, True)
        Mr, rank_aug = M.rref()
        C = self.ncols
        pivots = []
        for i in range(rank_aug):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        rank = rank_aug - (0 if consistent else 1)
        out = {'p': p, 'consistent': consistent, 'rank_p': rank, 'rank_aug_p': rank_aug, 'ncols': C, 'nrows': self.nrows,
               'nullity_p': C - rank, 'full_column_rank': rank == C, 't_solve': round(time.time() - t0, 1)}
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
            out['x'] = x
        if verbose: print(f'[kred p={p}] consistent={consistent} rank={rank} nullity={C-rank} ({out["t_solve"]}s)', file=sys.stderr, flush=True)
        return out

    def kernel_mod(self, p):
        M = self.dense_mod(p, False)
        X, nullity = M.nullspace()
        return X, nullity

    def exact_rank(self, with_rhs):
        """Exact rank over Q via fmpz_mat (fraction-free). Feasible for moderate sizes."""
        t0 = time.time()
        M = self.dense_int(with_rhs)
        r = M.rank()
        return r, round(time.time() - t0, 1)

def unit_field(i, m=None):
    w = [{} for _ in range(6)]
    w[i] = {m if m is not None else (0,) * 6: 1}
    return w

def rand_field(deg, par=None, rng=None, lo=-3, hi=3):
    """random integer vector field of degree deg; if par is given, make w(F) have Y-parity par."""
    import random
    rng = rng or random
    w = [{} for _ in range(6)]
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(deg):
            if par is not None and (yparity(m) + pi) % 2 != par: continue
            c = rng.randint(lo, hi)
            if c: w[i][m] = c
    return w
