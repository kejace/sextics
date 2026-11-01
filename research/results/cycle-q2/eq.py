"""Equivariant (character-graded) fixed-denominator solver for the hodge.computer
certificate equation, specialised to diagonal-symmetric sextics.

Conventions as in hodge.core:  E = A D^2 - q D sum N_i F_i + F (D sum d_i N_i - sum N_i d_i D) = F^(q+1) Q,
deg N_i = d + 6q - 5, deg Q = 2d - 6.

If every monomial of D has the same weight under a subgroup H of the diagonal symmetry
group of F (a subgroup of (Z/6)^6 acting by sixth roots of unity), then any solution (N,Q) can be
averaged (with the twisted action) to one where N_i has character chi_A chi_D chi_{x_i} and Q has
character chi_A chi_D^2; the rows (monomials of E) have character chi_A chi_D^2.
"""
import sys, time, json, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import numpy as np
from flint import nmod_mat

NV = 6
BASE = 64          # exponent code base: code = sum e_i * 64^i (exponents < 64)
POW = np.array([BASE**i for i in range(NV)], dtype=np.int64)
PRIME = 2147483647

def enc(e):
    return int(sum(int(e[i]) * BASE**i for i in range(NV)))

def dec(c):
    e = []
    for i in range(NV):
        e.append(c % BASE); c //= BASE
    return tuple(e)

def poly_to_arrays(P, p):
    """dict poly -> (codes int64 array, vals int64 array mod p)"""
    Pp = core.reduce_mod(P, p)
    if not Pp:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    codes = np.array([enc(e) for e in Pp], dtype=np.int64)
    vals = np.array([int(v) for v in Pp.values()], dtype=np.int64)
    return codes, vals

def monomial_exps(deg):
    return np.array(core.monomials(deg), dtype=np.int64).reshape(-1, NV)

# ----------------------------------------------------------------------------
# diagonal symmetry group
# ----------------------------------------------------------------------------

def symmetry_group(F):
    """All w in (Z/6)^6 with every monomial of F having the same weight (mod 6).  Returns int array (n x 6)."""
    exps = np.array(list(F.keys()), dtype=np.int64)
    allw = np.array(list(itertools.product(range(6), repeat=NV)), dtype=np.int64)
    wts = (exps @ allw.T) % 6          # (#monomials x 6^6)
    ok = np.all(wts == wts[0:1, :], axis=0)
    return allw[ok]

def stabilizer(Gw, D):
    """Subgroup of Gw under which all monomials of D have equal weight."""
    exps = np.array(list(D.keys()), dtype=np.int64)
    wts = (exps @ Gw.T) % 6
    ok = np.all(wts == wts[0:1, :], axis=0)
    return Gw[ok]

def generators(Hw):
    """Greedy generating set of the subgroup Hw (rows), as an int array (k x 6)."""
    Hset = {tuple(w) for w in Hw}
    gens = []
    span = {(0,) * NV}
    for w in Hw:
        tw = tuple(w)
        if tw in span:
            continue
        gens.append(tw)
        # recompute span
        new = set(span)
        frontier = list(span)
        while frontier:
            nxt = []
            for s in frontier:
                t = tuple((s[i] + tw[i]) % 6 for i in range(NV))
                if t not in new:
                    new.add(t); nxt.append(t)
            frontier = nxt
        span = new
        if len(span) == len(Hset):
            break
    return np.array(gens, dtype=np.int64).reshape(-1, NV)

def keys_of(exps, gens):
    """character keys (as int64 scalars encoding the weight tuple base 6) for exponent rows."""
    if gens.shape[0] == 0:
        return np.zeros(len(exps), dtype=np.int64)
    wts = (exps @ gens.T) % 6
    mult = np.array([6**j for j in range(gens.shape[0])], dtype=np.int64)
    return wts @ mult

def key_of_poly(P, gens):
    exps = np.array(list(P.keys()), dtype=np.int64)
    k = keys_of(exps, gens)
    assert np.all(k == k[0]), 'polynomial is not semi-invariant under given generators'
    return int(k[0])

def key_add(k1, k2, ngens):
    """add two keys (base-6 digit-wise mod 6)."""
    out = 0
    m = 1
    for j in range(ngens):
        d = ((k1 // m) % 6 + (k2 // m) % 6) % 6
        out += d * m
        m *= 6
    return out

def key_scale(k, c, ngens):
    out = 0; m = 1
    for j in range(ngens):
        out += (((k // m) % 6) * c % 6) * m
        m *= 6
    return out

# ----------------------------------------------------------------------------
# the fixed-D system in one character sector
# ----------------------------------------------------------------------------

class EqSystem:
    def __init__(self, F, A, q, D, p=PRIME, with_target=False, verbose=False, Hw=None):
        t0 = time.time()
        self.p = p; self.q = q
        d = core.degree(D); self.d = d
        Gw = symmetry_group(F) if Hw is None else Hw
        Hw = stabilizer(Gw, D)
        Hw = stabilizer(Hw, A)   # A must be semi-invariant too (it is, as a monomial; for general A restrict)
        self.Hw = Hw
        gens = generators(Hw); self.gens = gens; ng = gens.shape[0]
        self.order = len(Hw)
        kA = key_of_poly(A, gens); kD = key_of_poly(D, gens) if D else 0
        kx = [int(keys_of(np.eye(NV, dtype=np.int64)[i:i+1], gens)[0]) for i in range(NV)]
        self.kA, self.kD = kA, kD
        kE = key_add(kA, key_scale(kD, 2, ng), ng)      # character of A D^2 (rows, Q)
        kN = [key_add(key_add(kA, kD, ng), kx[i], ng) for i in range(NV)]
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(NV)]
        Di = [core.pderiv(Dp, i, p) for i in range(NV)]
        FD = core.mul(Fp, Dp, p)
        # G_i = -(q D F_i - F D_i)
        G = [core.scale(core.add(core.scale(core.mul(Dp, Fi[i], p), q, p), core.mul(Fp, Di[i], p), p), p - 1, p) for i in range(NV)]  # -(q D F_i + F D_i)
        Fq1 = core.scale(core.power(Fp, q + 1, p), p - 1, p)
        D2 = core.mul(Dp, Dp, p)
        rhs = core.scale(core.mul(Ap, D2, p), p - 1, p)
        # columns
        col_codes = []; col_vals = []; col_ids = []; col_meta = []
        ncol = 0
        allN = monomial_exps(degN)
        kN_all = keys_of(allN, gens)
        self.nN = 0
        FDc, FDv = poly_to_arrays(FD, p)
        for i in range(NV):
            sel = allN[kN_all == kN[i]]
            if len(sel) == 0: continue
            mcodes = sel @ POW
            Gc, Gv = poly_to_arrays(G[i], p)
            # part 1: G_i shifted by m
            if len(Gc):
                cc = (Gc[:, None] + mcodes[None, :]).ravel()
                vv = np.broadcast_to(Gv[:, None], (len(Gc), len(sel))).ravel()
                ii = np.broadcast_to(np.arange(ncol, ncol + len(sel))[None, :], (len(Gc), len(sel))).ravel()
                col_codes.append(cc); col_vals.append(vv.copy()); col_ids.append(ii.copy())
            # part 2: m_i * FD shifted by m - e_i
            mi = sel[:, i]
            has = mi > 0
            if np.any(has) and len(FDc):
                sel2 = sel[has]; mc2 = mcodes[has] - BASE**i
                cc = (FDc[:, None] + mc2[None, :]).ravel()
                vv = ((FDv[:, None] * mi[has][None, :]) % p).ravel()
                ii = np.broadcast_to((np.arange(ncol, ncol + len(sel))[has])[None, :], (len(FDc), int(has.sum()))).ravel()
                col_codes.append(cc); col_vals.append(vv); col_ids.append(ii.copy())
            for m in sel:
                col_meta.append(('N', i, tuple(int(x) for x in m)))
            ncol += len(sel); self.nN += len(sel)
        self.nQ = 0
        if degQ >= 0:
            allQ = monomial_exps(degQ)
            selQ = allQ[keys_of(allQ, gens) == kE]
            if len(selQ):
                mcodes = selQ @ POW
                Fc, Fv = poly_to_arrays(Fq1, p)
                cc = (Fc[:, None] + mcodes[None, :]).ravel()
                vv = np.broadcast_to(Fv[:, None], (len(Fc), len(selQ))).ravel()
                ii = np.broadcast_to(np.arange(ncol, ncol + len(selQ))[None, :], (len(Fc), len(selQ))).ravel()
                col_codes.append(cc); col_vals.append(vv.copy()); col_ids.append(ii.copy())
                for m in selQ:
                    col_meta.append(('Q', -1, tuple(int(x) for x in m)))
                ncol += len(selQ); self.nQ = len(selQ)
        self.nT = 0
        if with_target:
            allT = monomial_exps(6 * q)
            selT = allT[keys_of(allT, gens) == kA]
            mcodes = selT @ POW
            Dc, Dv = poly_to_arrays(D2, p)
            cc = (Dc[:, None] + mcodes[None, :]).ravel()
            vv = np.broadcast_to(Dv[:, None], (len(Dc), len(selT))).ravel()
            ii = np.broadcast_to(np.arange(ncol, ncol + len(selT))[None, :], (len(Dc), len(selT))).ravel()
            col_codes.append(cc); col_vals.append(vv.copy()); col_ids.append(ii.copy())
            for m in selT:
                col_meta.append(('A', -1, tuple(int(x) for x in m)))
            ncol += len(selT); self.nT = len(selT)
        rc, rv = poly_to_arrays(rhs, p)
        col_codes.append(rc); col_vals.append(rv); col_ids.append(np.full(len(rc), ncol, dtype=np.int64))
        codes = np.concatenate(col_codes); vals = np.concatenate(col_vals); ids = np.concatenate(col_ids)
        # rows
        urows, rinv = np.unique(codes, return_inverse=True)
        self.nrows = len(urows); self.ncols = ncol
        self.row_codes = urows
        # accumulate duplicates: sort by (col,row)
        lin = rinv.astype(np.int64) * (ncol + 1) + ids
        order = np.argsort(lin, kind='stable')
        lin_s = lin[order]; vals_s = vals[order]
        uniq, start = np.unique(lin_s, return_index=True)
        sums = np.add.reduceat(vals_s % p, start) % p
        keep = sums != 0
        self.r_idx = (uniq[keep] // (ncol + 1)).astype(np.int64)
        self.c_idx = (uniq[keep] % (ncol + 1)).astype(np.int64)
        self.v = sums[keep].astype(np.int64)
        self.col_meta = col_meta
        self.build_time = time.time() - t0
        if verbose:
            print(f'[eq] d={d} q={q} |H|={self.order} gens={ng} rows={self.nrows} cols={ncol} (N {self.nN}, Q {self.nQ}, T {self.nT}) nnz={len(self.v)} build {self.build_time:.1f}s', file=sys.stderr)

    def dense(self, with_rhs=True):
        C = self.ncols + (1 if with_rhs else 0)
        M = np.zeros((self.nrows, C), dtype=np.int64)
        if with_rhs:
            M[self.r_idx, self.c_idx] = self.v
        else:
            sel = self.c_idx < self.ncols
            M[self.r_idx[sel], self.c_idx[sel]] = self.v[sel]
        return M

    def compressed(self, k_extra=40, seed=1, with_rhs=True):
        """random row compression to k = min(rows, cols+1+k_extra) rows, exact mod p."""
        import scipy.sparse as sp
        p = self.p
        C = self.ncols + (1 if with_rhs else 0)
        if with_rhs:
            M = sp.csr_matrix((self.v, (self.r_idx, self.c_idx)), shape=(self.nrows, C))
        else:
            sel = self.c_idx < self.ncols
            M = sp.csr_matrix((self.v[sel], (self.r_idx[sel], self.c_idx[sel])), shape=(self.nrows, C))
        k = min(self.nrows, C + k_extra)
        if k == self.nrows:
            return np.asarray(M.todense()) % p
        rng = np.random.default_rng(seed)
        R = rng.integers(0, p, size=(k, self.nrows), dtype=np.int64)
        MT = M.T.tocsr()
        out = np.zeros((k, C), dtype=np.int64)
        # split R into 16-bit chunks to avoid overflow: R = r0 + r1*2^16 (r1 < 2^15)
        r0 = R & 0xFFFF; r1 = R >> 16
        a = (MT @ r0.T) % p       # (C x k), entries < p * 2^16 * nnz_per_col... use chunking if needed
        b = (MT @ r1.T) % p
        out = (a + (b * (1 << 16)) % p) % p
        return np.ascontiguousarray(out.T)

    def solve(self, verbose=False, k_extra=40, seed=1):
        t0 = time.time()
        Mc = self.compressed(k_extra=k_extra, seed=seed)
        t1 = time.time()
        k, C1 = Mc.shape
        Mn = nmod_mat(k, C1, Mc.ravel().tolist(), self.p)
        Mr, rank = Mn.rref()
        t2 = time.time()
        C = self.ncols
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        res = {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'rows': self.nrows, 'cols': C,
               'k': k, 't_compress': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2), 'pivots': pivots, 'x': x}
        if verbose:
            print(f'[eq solve] consistent={consistent} rank={rank} nullity={C-rank} (compress {t1-t0:.1f}s rref {t2-t1:.1f}s)', file=sys.stderr)
        return res

    def target_image(self, verbose=False):
        """requires with_target=True.  Returns (dim served in sector, nT, basis list of polys mod p)."""
        t0 = time.time()
        Mc = self.compressed(with_rhs=False)
        k, C = Mc.shape
        Mn = nmod_mat(k, C, Mc.ravel().tolist(), self.p)
        X, nullity = Mn.nullspace()
        a_idx = [j for j, c in enumerate(self.col_meta) if c[0] == 'A']
        if nullity == 0:
            return 0, self.nT, []
        proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], self.p)
        T = proj.transpose()
        Tr, rk = T.rref()
        basis = []
        for i in range(rk):
            poly = {}
            for jj, j in enumerate(a_idx):
                v = int(Tr[i, jj])
                if v:
                    poly[self.col_meta[j][2]] = v
            basis.append(poly)
        if verbose:
            print(f'[target_image] served dim {rk} of sector dim {self.nT} ({time.time()-t0:.1f}s)', file=sys.stderr)
        return rk, self.nT, basis

    def unpack(self, x):
        N = [{} for _ in range(NV)]; Q = {}
        for (kind, i, m), v in zip(self.col_meta, x):
            if not v: continue
            if kind == 'N': N[i][m] = v
            elif kind == 'Q': Q[m] = v
        return N, Q


def e1_system(F, A, q, D, p=PRIME, verbose=False):
    """E1 only: q N(F) = A D in the sector. returns dict(consistent, rank, nullity, rows, cols)."""
    t0 = time.time()
    d = core.degree(D)
    Gw = symmetry_group(F)
    Hw = stabilizer(stabilizer(Gw, D), A)
    gens = generators(Hw); ng = gens.shape[0]
    kA = key_of_poly(A, gens); kD = key_of_poly(D, gens)
    kx = [int(keys_of(np.eye(NV, dtype=np.int64)[i:i+1], gens)[0]) for i in range(NV)]
    kN = [key_add(key_add(kA, kD, ng), kx[i], ng) for i in range(NV)]
    degN = d + 6 * q - 5
    Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(NV)]
    allN = monomial_exps(degN); kall = keys_of(allN, gens)
    codes = []; vals = []; ids = []; ncol = 0
    for i in range(NV):
        sel = allN[kall == kN[i]]
        if not len(sel): continue
        mcodes = sel @ POW
        Fc, Fv = poly_to_arrays(core.scale(Fi[i], q, p), p)
        cc = (Fc[:, None] + mcodes[None, :]).ravel()
        vv = np.broadcast_to(Fv[:, None], (len(Fc), len(sel))).ravel().copy()
        ii = np.broadcast_to(np.arange(ncol, ncol + len(sel))[None, :], (len(Fc), len(sel))).ravel().copy()
        codes.append(cc); vals.append(vv); ids.append(ii); ncol += len(sel)
    rc, rv = poly_to_arrays(core.mul(Ap, Dp, p), p)
    codes.append(rc); vals.append(rv); ids.append(np.full(len(rc), ncol, dtype=np.int64))
    codes = np.concatenate(codes); vals = np.concatenate(vals); ids = np.concatenate(ids)
    urows, rinv = np.unique(codes, return_inverse=True)
    M = np.zeros((len(urows), ncol + 1), dtype=np.int64)
    np.add.at(M, (rinv, ids), vals)
    M %= p
    R, C1 = M.shape
    Mn = nmod_mat(R, C1, M.ravel().tolist(), p)
    rank_aug = Mn.rank()
    M0 = nmod_mat(R, ncol, M[:, :ncol].ravel().tolist(), p)
    rank0 = M0.rank()
    return {'consistent': rank_aug == rank0, 'rank': rank0, 'nullity': ncol - rank0, 'rows': R, 'cols': ncol,
            'order': len(Hw), 'seconds': round(time.time() - t0, 2)}


LOG = '/tmp/hodge/results/cycle-q2/log.jsonl'

def log_entry(**kw):
    with open(LOG, 'a') as f:
        f.write(json.dumps(kw) + '\n')


def test_D(F, A, q, D, family, note='', verbose=True, do_e1=True, do_full=True, target=False, Dstr=None):
    """Run E1 and the full sector system for D; log; return dict."""
    Dstr = Dstr or core.to_str(D)
    d = core.degree(D)
    out = {'D': Dstr, 'd': d, 'family': family}
    t0 = time.time()
    e1 = None
    if do_e1:
        r1 = e1_system(F, A, q, D)
        e1 = r1['consistent']
        out['e1'] = e1
        if verbose:
            print(f'  E1: {e1} rank={r1["rank"]} rows={r1["rows"]} cols={r1["cols"]} |H|={r1["order"]} ({r1["seconds"]}s)', file=sys.stderr)
    res = None
    if do_full and (e1 is None or e1):
        S = EqSystem(F, A, q, D, with_target=target, verbose=verbose)
        if target:
            dim, nT, basis = S.target_image(verbose=verbose)
            out['served_dim'] = dim; out['sector_dim'] = nT
            # A in served space?  (A is a single monomial in the sector; check membership)
            # do the plain solve as well
        S2 = S if not target else EqSystem(F, A, q, D, with_target=False, verbose=False)
        res = S2.solve(verbose=verbose)
        out.update({'consistent': res['consistent'], 'rank': res['rank'], 'nullity': res['nullity'],
                    'rows': res['rows'], 'cols': res['cols'], 'order': S2.order})
        out['_system'] = S2; out['_res'] = res
    else:
        out.update({'consistent': False if e1 is False else None, 'rank': -1, 'nullity': -1})
    out['seconds'] = round(time.time() - t0, 2)
    out['note'] = note
    log_entry(**{k: v for k, v in out.items() if not k.startswith('_')})
    return out


def solve_exact(F, A, q, D, max_primes=10, verbose=False):
    """CRT + rational reconstruction over several primes of the sector system; exact verification."""
    acc = None; piv0 = None; meta = None
    for k, p in enumerate(core.PRIMES[:max_primes]):
        S = EqSystem(F, A, q, D, p=p)
        r = S.solve()
        if not r['consistent']:
            return 'inconsistent', r
        if piv0 is None or piv0 != r['pivots']:
            acc = None; piv0 = r['pivots']; meta = S.col_meta
        x = r['x']
        if acc is None:
            acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, i, mon), v in zip(meta, rat):
                if v == 0: continue
                if kind == 'N': N[i][mon] = v
                elif kind == 'Q': Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                return 'solved', {'N': N, 'Q': Q, 'primes': k + 1, 'rank': r['rank'], 'nullity': r['nullity']}
            elif verbose:
                print(f'[exact] reconstruction failed verification after {k+1} primes', file=sys.stderr)
    return 'lift-failed', {}


# ----------------------------------------------------------------------------
# fast path: sparse row compression + blocked float64 modular LU (modlu)
# ----------------------------------------------------------------------------
import modlu

def _compressed_int(S, k_extra=40, seed=1, with_rhs=True):
    """(k x C[+1]) int64 compressed matrix with small random multipliers (< 2^16); exact mod S.p
    provided S.p < 2^22 and columns have < 2^14 nonzeros."""
    import scipy.sparse as sp
    p = S.p
    C = S.ncols + (1 if with_rhs else 0)
    if with_rhs:
        M = sp.csr_matrix((S.v, (S.r_idx, S.c_idx)), shape=(S.nrows, C), dtype=np.int64)
    else:
        sel = S.c_idx < S.ncols
        M = sp.csr_matrix((S.v[sel], (S.r_idx[sel], S.c_idx[sel])), shape=(S.nrows, C), dtype=np.int64)
    k = min(S.nrows, C + k_extra)
    if k == S.nrows:
        return np.asarray(M.todense(), dtype=np.int64) % p
    rng = np.random.default_rng(seed)
    R = rng.integers(0, 1 << 16, size=(k, S.nrows), dtype=np.int64)
    MT = M.T.tocsr()
    out = (MT @ R.T) % p          # (C x k)
    return np.ascontiguousarray(out.T)

def solve_fast(S, verbose=False, k_extra=40, seed=1):
    """consistency + rank via modlu on the compressed system. S must be built with p < 2^22."""
    assert S.p * S.p * 128 < 2 ** 53
    t0 = time.time()
    Mc = _compressed_int(S, k_extra=k_extra, seed=seed)
    t1 = time.time()
    consistent, rank, rank_aug = modlu.consistent_rank(Mc, S.p)
    t2 = time.time()
    C = S.ncols
    res = {'consistent': bool(consistent), 'rank': int(rank), 'nullity': int(C - rank), 'rows': S.nrows, 'cols': C,
           'k': Mc.shape[0], 't_compress': round(t1 - t0, 2), 't_lu': round(t2 - t1, 2)}
    if verbose:
        print(f'[fast p={S.p}] consistent={consistent} rank={rank} nullity={C-rank} (compress {t1-t0:.1f}s LU {t2-t1:.1f}s)', file=sys.stderr)
    return res

def test_D_fast(F, A, q, D, family, note='', verbose=True, Dstr=None, do_e1=False, Hw=None):
    Dstr = Dstr or core.to_str(D)
    d = core.degree(D)
    out = {'D': Dstr, 'd': d, 'family': family}
    t0 = time.time()
    if do_e1:
        r1 = e1_system(F, A, q, D)
        out['e1'] = r1['consistent']
        if not r1['consistent']:
            out.update({'consistent': False, 'rank': -1, 'nullity': -1, 'order': r1['order'], 'seconds': round(time.time() - t0, 2), 'note': note + ' (E1 fails)'})
            log_entry(**out)
            return out
    S = EqSystem(F, A, q, D, p=modlu.P22, verbose=verbose, Hw=Hw)
    res = solve_fast(S, verbose=verbose)
    out.update({'consistent': res['consistent'], 'rank': res['rank'], 'nullity': res['nullity'],
                'rows': res['rows'], 'cols': res['cols'], 'order': S.order, 'e1': out.get('e1', None)})
    out['seconds'] = round(time.time() - t0, 2)
    out['note'] = note
    log_entry(**out)
    out['_system'] = S
    return out
