"""Exact scan of 2-dimensional (pencil) E1-kernel sectors: all F_p-rational points of the pencil D(t) = D1 + t D2 (and t = inf)
that pass the full E2 system, via interpolation of a maximal minor of the augmented sector system.

For the sector system M(t) x = b(t) (full fixed-D system of eq.EqSystem: M linear in t, b quadratic in t), let r be the
generic rank of M(t).  Choose a random compression S ((r+1) x rows) and the pivot columns C (r+1 of them, including b) of
the augmented matrix at a random t0; g(t) = det(S [M(t)|b(t)]_{:,C}) is a polynomial of degree <= r+2 that vanishes at every
consistent t.  Interpolate g from r+3 evaluations, factor over F_p, test every F_p-root with the full solve, and record the
degrees of the irreducible factors (a point defined over an extension of degree k shows up as a degree-k factor).
Usage: python3 w2_pencil.py d [p]"""
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


def pencil_scan(F, A, q, D1, D2, p, verbose=False):
    t0 = time.time()
    # generic rank from a random t
    rng = np.random.default_rng(1)
    tr = int(rng.integers(1, p - 1))
    Dt = core.add(core.reduce_mod(D1, p), core.scale(core.reduce_mod(D2, p), tr, p), p)
    S0 = eq.EqSystem(F, A, q, Dt, p=p)
    res0 = S0.solve()
    r = res0['rank'] if not res0['consistent'] else res0['rank']
    # augmented rank
    Mc = S0.compressed(k_extra=40, seed=3)
    k, C1 = Mc.shape
    Mn = nmod_mat(k, C1, Mc.ravel().tolist(), p)
    Mr, rank_aug = Mn.rref()
    piv = []
    for i in range(rank_aug):
        j = 0
        while j < C1 and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    cols = S0.ncols
    generic_consistent = not (piv and piv[-1] == cols)
    info = {'rank': res0['rank'], 'rank_aug': rank_aug, 'cols': cols, 'rows': S0.nrows, 'generic_consistent': generic_consistent}
    if generic_consistent:
        info['result'] = 'GENERICALLY CONSISTENT pencil (!)'
        return info
    m = rank_aug  # = r + 1
    # fixed compression S (m x rows) and column set piv (includes b column)
    rng2 = np.random.default_rng(7)
    Srows = rng2.integers(0, p, size=(m, S0.nrows), dtype=np.int64)
    # the row index set must be consistent across t: EqSystem builds rows from the union of supports, which may vary with t
    # (zero coefficients at special t drop rows).  Use the monomial codes of degree 2d+12 in the sector as a fixed row basis.
    row_codes0 = S0.row_codes

    def g_at(t):
        Dt = core.add(core.reduce_mod(D1, p), core.scale(core.reduce_mod(D2, p), t % p, p), p)
        St = eq.EqSystem(F, A, q, Dt, p=p)
        # map St rows into the fixed row set
        pos = np.searchsorted(row_codes0, St.row_codes)
        ok = (pos < len(row_codes0)) & (row_codes0[np.minimum(pos, len(row_codes0) - 1)] == St.row_codes)
        if not np.all(ok):
            # new rows appeared (should not happen: supports shrink or stay)
            return None
        import scipy.sparse as sp
        M = sp.csr_matrix((St.v, (pos[St.r_idx], St.c_idx)), shape=(len(row_codes0), St.ncols + 1), dtype=np.int64)
        sub = M[:, piv]
        r0 = Srows & 0xFFFF; r1 = Srows >> 16
        subT = sub.T.tocsr()
        a = (subT @ r0.T) % p
        b = (subT @ r1.T) % p
        comp = ((a + (b * (1 << 16)) % p) % p).T   # m x m
        Mm = nmod_mat(m, m, np.ascontiguousarray(comp).ravel().tolist(), p)
        return int(Mm.det())
    deg = m + 2
    ts = []
    vals = []
    t = 2
    while len(ts) < deg + 1:
        v = g_at(t)
        if v is not None:
            ts.append(t); vals.append(v)
        t += 1
    # Newton interpolation mod p
    n = len(ts)
    coef = list(vals)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) * pow((ts[i] - ts[i - j]) % p, -1, p) % p
    # convert Newton form to monomial coefficients
    poly = [0] * n
    poly_newton = [coef[n - 1]]
    for i in range(n - 2, -1, -1):
        # poly = poly * (x - ts[i]) + coef[i]
        new = [0] * (len(poly_newton) + 1)
        for k2, c in enumerate(poly_newton):
            new[k2 + 1] = (new[k2 + 1] + c) % p
            new[k2] = (new[k2] - c * ts[i]) % p
        new[0] = (new[0] + coef[i]) % p
        poly_newton = new
    g = nmod_poly(poly_newton, p)
    # trim
    while g.degree() >= 0 and int(g[g.degree()]) == 0:
        g = nmod_poly(list(g.coeffs())[:-1], p)
    if g.degree() < 0:
        info['result'] = 'g identically zero (compression degenerate); retry with another seed'
        return info
    info['deg_g'] = g.degree()
    fac = g.factor()
    unit, factors = fac
    fdegs = sorted([f.degree() for f, e in factors])
    info['factor_degrees'] = fdegs
    roots = [(-int(f[0]) * pow(int(f[1]), -1, p)) % p for f, e in factors if f.degree() == 1]
    info['n_roots'] = len(roots)
    consistent_roots = []
    for rt in roots:
        Dt = core.add(core.reduce_mod(D1, p), core.scale(core.reduce_mod(D2, p), rt, p), p)
        St = eq.EqSystem(F, A, q, Dt, p=p)
        rs = St.solve()
        rec = {'t': rt, 'consistent': rs['consistent'], 'rank': rs['rank'], 'ratrecon': str(core.ratrecon(rt, p))}
        if rs['consistent']:
            consistent_roots.append(rec)
    info['consistent_roots'] = consistent_roots
    # t = infinity: D2 alone
    S2 = eq.EqSystem(F, A, q, core.reduce_mod(D2, p), p=p)
    r2 = S2.solve()
    info['t_inf_consistent'] = r2['consistent']
    info['seconds'] = round(time.time() - t0, 1)
    info['result'] = ('CONSISTENT POINT(S) FOUND' if consistent_roots or r2['consistent'] else
                      f'no F_p-rational consistent point in the pencil (g has {len(roots)} roots, factor degrees {fdegs[:12]}{"..." if len(fdegs) > 12 else ""})')
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
