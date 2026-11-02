"""Exact scan of 2-dimensional (pencil) E1-kernel sectors: all F_p-rational points of the pencil D(t) = D1 + t D2 (and t = inf)
that pass the full E2 system, via interpolation of a maximal minor of the augmented sector system.

For the sector system M(t) x = b(t) (full fixed-D system of eq.EqSystem: M linear in t, b quadratic in t), let r be the
generic rank of M(t).  Choose a random compression S ((r+1) x rows) and the pivot columns C (r+1 of them, including b) of
the augmented matrix at a random t0; g(t) = det(S [M(t)|b(t)]_{:,C}) is a polynomial of degree <= r+2 that vanishes at every
consistent t.  Interpolate g from r+3 evaluations, factor over F_p, test every F_p-root with the full solve, and record the
degrees of the irreducible factors (a point defined over an extension of degree k shows up as a degree-k factor).
Usage: python3 w2_pencil.py d [p]"""
import os as _os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    _os.environ[_v] = "1"
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
from flint import nmod_mat, nmod_poly, nmod
import w2lib as W
import eq
from hodge import core


def sector_system_modp(F, A, q, D, p):
    S = eq.EqSystem(F, A, q, D, p=p)
    return S


def _compress_blocks(S_, row_codes, mats, p, m, seed=7):
    """random (m x rows) compression applied to sparse int64 matrices given as (r_idx, c_idx, v, ncols) over row_codes"""
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    R = rng.integers(0, p, size=(m, len(row_codes)), dtype=np.int64)
    r0 = R & 0xFFFF; r1 = R >> 16
    outs = []
    for (ri, ci, v, ncols) in mats:
        M = sp.csr_matrix((v, (ri, ci)), shape=(len(row_codes), ncols), dtype=np.int64)
        MT = M.T.tocsr()
        a = (MT @ r0.T) % p
        b = (MT @ r1.T) % p
        outs.append(np.ascontiguousarray(((a + (b * (1 << 16)) % p) % p).T))   # m x ncols
    return outs


def pencil_scan(F, A, q, D1, D2, p, verbose=False):
    """all F_p-rational t with D1 + t D2 consistent (plus t = inf), via interpolation of a maximal minor polynomial.
    M(t) = [MN1 + t MN2 | MQ],  b(t) = b11 + t b12 + t^2 b22  (all systems share the sector of kappa, since D1, D2 are
    G-semi-invariant with the same character)."""
    import os
    t0 = time.time()
    D1p, D2p = core.reduce_mod(D1, p), core.reduce_mod(D2, p)
    S1 = eq.EqSystem(F, A, q, D1p, p=p)
    S2 = eq.EqSystem(F, A, q, D2p, p=p)
    assert S1.ncols == S2.ncols and S1.col_meta == S2.col_meta, 'column sets differ (sector mismatch)'
    ncols, nN, nQ = S1.ncols, S1.nN, S1.nQ
    Ap = core.reduce_mod(A, p)
    b12 = core.scale(core.mul(Ap, core.mul(D1p, D2p, p), p), (p - 2) % p, p)   # -2 A D1 D2
    codes12 = np.array([eq.enc(e) for e in b12], dtype=np.int64)
    row_codes = np.unique(np.concatenate([S1.row_codes, S2.row_codes, codes12]))
    def remap(S_):
        pos = np.searchsorted(row_codes, S_.row_codes)
        return pos[S_.r_idx]
    r1 = remap(S1); r2 = remap(S2)
    selN1 = S1.c_idx < nN; selQ1 = (S1.c_idx >= nN) & (S1.c_idx < ncols); selb1 = S1.c_idx == ncols
    selN2 = S2.c_idx < nN; selb2 = S2.c_idx == ncols
    # generic rank / augmented rank at a random t
    rng = np.random.default_rng(1)
    tr = int(rng.integers(2, p - 1))
    Dt = core.add(D1p, core.scale(D2p, tr, p), p)
    St = eq.EqSystem(F, A, q, Dt, p=p)
    Mc = St.compressed(k_extra=40, seed=3)
    k, C1 = Mc.shape
    Mn = nmod_mat(k, C1, Mc.ravel().tolist(), p)
    Mr, rank_aug = Mn.rref()
    piv = []
    for i in range(rank_aug):
        j = 0
        while j < C1 and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    generic_consistent = not (piv and piv[-1] == ncols)
    info = {'rank_aug_generic': rank_aug, 'rank_generic': rank_aug - (0 if generic_consistent else 1), 'cols': ncols, 'rows': len(row_codes),
            'generic_consistent': generic_consistent}
    if generic_consistent:
        info['result'] = 'GENERICALLY CONSISTENT pencil (!)'
        info['seconds'] = round(time.time() - t0, 1)
        return info
    m = rank_aug
    piv = np.array(piv)
    # compressed blocks
    mats = [(r1[selN1], S1.c_idx[selN1], S1.v[selN1], nN), (r2[selN2], S2.c_idx[selN2], S2.v[selN2], nN),
            (r1[selQ1], S1.c_idx[selQ1] - nN, S1.v[selQ1], nQ),
            (r1[selb1], np.zeros(int(selb1.sum()), dtype=np.int64), S1.v[selb1], 1),
            (np.searchsorted(row_codes, codes12), np.zeros(len(codes12), dtype=np.int64), np.array([b12[e] for e in b12], dtype=np.int64) if b12 else np.zeros(0, dtype=np.int64), 1),
            (r2[selb2], np.zeros(int(selb2.sum()), dtype=np.int64), S2.v[selb2], 1)]
    CN1, CN2, CQ, c11, c12, c22 = _compress_blocks(S1, row_codes, mats, p, m)
    pivN = piv[piv < nN]; pivQ = piv[(piv >= nN) & (piv < ncols)] - nN
    has_b = bool(np.any(piv == ncols))
    assert has_b, 'b column not a pivot although inconsistent?'
    CN1p, CN2p, CQp = CN1[:, pivN], CN2[:, pivN], CQ[:, pivQ]

    def g_at(t):
        t %= p
        MN = (CN1p + (t * CN2p) % p) % p
        bb = (c11[:, 0] + (t * c12[:, 0]) % p + (((t * t) % p) * c22[:, 0]) % p) % p
        Mt = np.hstack([MN, CQp, bb[:, None]])
        assert Mt.shape == (m, m), Mt.shape
        return int(nmod_mat(m, m, Mt.ravel().tolist(), p).det())
    deg = m + 1
    if os.environ.get('W2_PENCIL_TEST'):
        deg = 4
    ts = list(range(2, 2 + deg + 1))
    vals = [g_at(t) for t in ts]
    n = len(ts)
    coef = list(vals)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) * pow((ts[i] - ts[i - j]) % p, -1, p) % p
    poly_newton = [coef[n - 1]]
    for i in range(n - 2, -1, -1):
        new = [0] * (len(poly_newton) + 1)
        for k2, c in enumerate(poly_newton):
            new[k2 + 1] = (new[k2 + 1] + c) % p
            new[k2] = (new[k2] - c * ts[i]) % p
        new[0] = (new[0] + coef[i]) % p
        poly_newton = new
    while len(poly_newton) > 1 and poly_newton[-1] == 0:
        poly_newton.pop()
    g = nmod_poly(poly_newton, p)
    if g.degree() <= 0 and int(g[0]) == 0:
        info['result'] = 'g identically zero (degenerate compression)'
        info['seconds'] = round(time.time() - t0, 1)
        return info
    info['deg_g'] = g.degree()
    # sanity: g at a fresh random point equals the direct evaluation
    tchk = int(rng.integers(2, p - 1))
    info['interp_check_ok'] = (int(g(nmod(tchk, p))) == g_at(tchk))
    unit, factors = g.factor()
    fdegs = sorted([f.degree() for f, e in factors])
    info['factor_degrees'] = fdegs[:40]
    info['n_factors'] = len(fdegs)
    roots = [(-int(f[0]) * pow(int(f[1]), -1, p)) % p for f, e in factors if f.degree() == 1]
    info['n_roots'] = len(roots)
    consistent_roots = []
    for rt in roots:
        Dt = core.add(D1p, core.scale(D2p, rt, p), p)
        Sx = eq.EqSystem(F, A, q, Dt, p=p)
        rs = Sx.solve()
        rec = {'t': rt, 'consistent': rs['consistent'], 'rank_aug': rs['rank'], 'ratrecon': str(core.ratrecon(rt, p))}
        if rs['consistent']:
            consistent_roots.append(rec)
    info['roots_checked'] = len(roots)
    info['consistent_roots'] = consistent_roots
    r2s = S2.solve()
    info['t_inf_consistent'] = r2s['consistent']
    info['seconds'] = round(time.time() - t0, 1)
    info['result'] = ('CONSISTENT POINT(S) FOUND' if consistent_roots or r2s['consistent'] else
                      f'no F_p-rational consistent point in the pencil (g degree {g.degree()}, {len(roots)} roots, irreducible factor degrees {fdegs[:12]})')
    return info


if __name__ == '__main__':
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    p = int(sys.argv[2]) if len(sys.argv) > 2 else W.PRIME
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    out = []
    for kappa in S.all_keys(d):
        r = W.e1_kernel_modp(S, d, kappa, q, p=p)
        if r['dimK'] != 2:
            continue
        D1, D2 = r['basis']
        info = pencil_scan(F, A, q, D1, D2, p)
        info.update({'d': d, 'key': kappa, 'p': p, 'dimJ': r['dimJ'], 'D1': core.to_str(D1)[:150], 'D2': core.to_str(D2)[:150]})
        print(json.dumps(info), flush=True)
        W.log(step='pencil-scan', **info, note='exact scan of a 2-dim E1-kernel sector family: F_p-roots of a maximal minor polynomial, each verified by the full sector solve')
        out.append(info)
    print(f'{len(out)} pencils scanned')
