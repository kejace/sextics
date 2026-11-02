"""Exact Jacobian relaxation of E2 at degree 6 (and optionally 7) for ALL E1-admissible denominators (no symmetry assumption).

K_d = {D in S_d : A D in J(F)} = ⊕_kappa K_kappa (sector decomposition, dim 62 at d = 6).  For D in K_d with any E1 lift N0(D)
(linear in D), E2 forces  Psi(D) := [B_D(N0(D))] = 0 in R_{2d+6} = (S/J)_{2d+6}  (im T_D, F^2 Q lie in J; the class is
independent of the lift).  Psi is a quadratic map K_d -> R_{2d+6}; V(Psi) is a cone containing the line C.F (F in K_6).
If, mod p, the ideal I_p of the quadrics equals (after saturation) the ideal of that line, then by flatness of the closure
of V(Psi)_Q over Z and the Hilbert-polynomial comparison, V(Psi)(C) = C.F, and no degree-d certificate exists with any D.

Usage: python3 w2_relax6.py [d] [p] [singular_timeout]
Writes relax{d}_{p}.sing and runs Singular: dim, mult, saturation vs. the F-line ideal."""
import os as _os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    _os.environ[_v] = "1"
import sys, time, json, os, subprocess
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
import w2lib as W
from flint import nmod_mat
from hodge import core

LANE = '/tmp/hodge/results/cycle-q2'
NV = 6
BASE = W.BASE
POW = W.POW


def nmod_rows(mat):
    return [[int(mat[i, j]) for j in range(mat.ncols())] for i in range(mat.nrows())]


class SectorReducer:
    """normal forms modulo (J_deg)_key over F_p; class coordinates = non-pivot entries"""
    def __init__(self, S, deg, key, p):
        self.mons = S.mons(deg, key)
        self.codes = np.sort(self.mons @ POW)
        self.n = len(self.codes)
        jb = W.jac_basis(S, deg, key)
        rows = []
        for c in jb:
            cp = core.reduce_mod(c, p)
            row = [0] * self.n
            for e, v in cp.items():
                j = int(np.searchsorted(self.codes, W.enc(e)))
                row[j] = v
            rows.append(row)
        self.p = p
        if rows:
            Jm = nmod_mat(len(rows), self.n, [x for r in rows for x in r], p)
            Jr, rk = Jm.rref()
            self.rank = rk
            self.Jr = Jr
            piv = []
            for i in range(rk):
                j = 0
                while j < self.n and int(Jr[i, j]) == 0:
                    j += 1
                piv.append(j)
            self.piv = piv
        else:
            self.rank = 0; self.piv = []; self.Jr = None
        self.nonpiv = [j for j in range(self.n) if j not in set(self.piv)]
        self.dimR = self.n - self.rank

    def coords(self, vecs):
        """vecs: list of dict polys mod p (degree deg, key) -> nmod_mat (len x dimR) of class coordinates"""
        m = len(vecs)
        flat = [0] * (m * self.n)
        for r, P in enumerate(vecs):
            for e, v in P.items():
                j = int(np.searchsorted(self.codes, W.enc(e)))
                assert self.codes[j] == W.enc(e)
                flat[r * self.n + j] = v % self.p
        V = nmod_mat(m, self.n, flat, self.p)
        if self.rank:
            # NF = V - V[:, piv] * Jr[:rank]
            Vp = nmod_mat(m, self.rank, [int(V[r, j]) for r in range(m) for j in self.piv], self.p)
            Jtop = nmod_mat(self.rank, self.n, [int(self.Jr[i, j]) for i in range(self.rank) for j in range(self.n)], self.p)
            NF = V - Vp * Jtop
        else:
            NF = V
        return [[int(NF[r, j]) for j in self.nonpiv] for r in range(m)]


def e1_lift_modp(S, D, kappa, p, q=2):
    """particular solution N (list of 6 dict polys mod p) of q N(F) = A D in the sector (D has key kappa)"""
    F, A = S.F, S.A
    d = core.degree(D)
    degN = d + 6 * q - 5
    kN = S.kN(kappa)
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(NV)]
    cols = []; meta = []
    for i in range(NV):
        for n in S.mons(degN, kN[i]):
            cols.append(core.shift(Fi[i], tuple(int(x) for x in n), q % p, p)); meta.append((i, tuple(int(x) for x in n)))
    rhs = core.mul(Ap, Dp, p)
    rows = {}
    for c in cols + [rhs]:
        for e in c:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, c in enumerate(cols):
        for e, v in c.items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items():
        flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    if piv and piv[-1] == C:
        raise ValueError('E1 inconsistent for a kernel element?!')
    N = [{} for _ in range(NV)]
    for i, j in enumerate(piv):
        v = int(Mr[i, C])
        if v:
            var, mon = meta[j]
            N[var][mon] = v
    return N


def _parr_mod(P, p):
    if not P:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return np.array([W.enc(e) for e in P], dtype=np.int64), np.array([int(v) % p for v in P.values()], dtype=np.int64)


def _acc_mod(codes, vals, p):
    if len(codes) == 0:
        return codes, vals
    u, inv = np.unique(codes, return_inverse=True)
    out = np.zeros(len(u), dtype=np.int64)
    np.add.at(out, inv, vals % p)
    out %= p
    keep = out != 0
    return u[keep], out[keep]


def _mul_mod(c1, v1, c2, v2, p):
    if len(c1) == 0 or len(c2) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    cc = (c1[:, None] + c2[None, :]).ravel()
    vv = ((v1[:, None] * v2[None, :]) % p).ravel()
    return _acc_mod(cc, vv, p)


def _deriv_mod(codes, vals, i, p):
    ei = (codes // (BASE ** i)) % BASE
    has = ei > 0
    return codes[has] - BASE ** i, (vals[has] * ei[has]) % p


class LiftArr:
    """E1 lift as arrays: components N_i and div N"""
    def __init__(self, N, p):
        self.N = [_parr_mod(N[i], p) for i in range(NV)]
        dc, dv = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
        for i in range(NV):
            c, v = _deriv_mod(*self.N[i], i, p)
            dc = np.concatenate([dc, c]); dv = np.concatenate([dv, v])
        self.div = _acc_mod(dc, dv, p)


def B_D_arr(Dc, Dv, Dder, L, p):
    """B_D(N) = D div N - sum_i N_i d_i D (arrays mod p); Dder = list of derivative arrays of D"""
    cc, vv = _mul_mod(Dc, Dv, L.div[0], L.div[1], p)
    parts_c = [cc]; parts_v = [vv]
    for i in range(NV):
        c, v = _mul_mod(L.N[i][0], L.N[i][1], Dder[i][0], Dder[i][1], p)
        parts_c.append(c); parts_v.append((p - v) % p)
    return _acc_mod(np.concatenate(parts_c), np.concatenate(parts_v), p)


def arr_to_dict_mod(codes, vals):
    return {W.dec(int(c)): int(v) for c, v in zip(codes, vals)}


def main(d=6, p=1000003, sing_timeout=600):
    t0 = time.time()
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    # 1. kernel bases per sector (mod p)
    basis = []  # (kappa, Dpoly mod p)
    dims = {}
    for kappa in S.all_keys(d):
        r = W.e1_kernel_modp(S, d, kappa, q, p=p)
        dims[kappa] = (r['dimS'], r['dimK'], r['dimJ'])
        for b in r['basis']:
            basis.append((kappa, b))
    K = len(basis)
    print(f'[relax{d}] p={p}: dim K_{d} = {K} (sum over sectors), J part {sum(v[2] for v in dims.values())} ({time.time()-t0:.1f}s)', flush=True)
    # locate F in the basis span (sanity: F in K_d for d=6)
    # 2. E1 lifts
    lifts = []
    for kappa, Dk in basis:
        lifts.append(LiftArr(e1_lift_modp(S, Dk, kappa, p, q), p))
    Darr = [_parr_mod(Dk, p) for (_, Dk) in basis]
    Dder = [[_deriv_mod(*Darr[k], i, p) for i in range(NV)] for k in range(len(basis))]
    print(f'[relax{d}] lifts done ({time.time()-t0:.1f}s)', flush=True)
    # 3. pair polys and class coordinates per target sector
    degE = 2 * d + 6
    reducers = {}
    pair_entries = {}  # key -> list of (k, l, poly)
    for k in range(K):
        kk, Dk = basis[k]
        for l in range(k, K):
            kl, Dl = basis[l]
            pc, pv = B_D_arr(Darr[k][0], Darr[k][1], Dder[k], lifts[l], p)
            if l != k:
                qc, qv = B_D_arr(Darr[l][0], Darr[l][1], Dder[l], lifts[k], p)
                pc, pv = _acc_mod(np.concatenate([pc, qc]), np.concatenate([pv, qv]), p)
            if len(pc) == 0:
                continue
            P = arr_to_dict_mod(pc, pv)
            key = S.add(S.add(kk, kl), S.kA)
            pair_entries.setdefault(key, []).append((k, l, P))
    print(f'[relax{d}] {sum(len(v) for v in pair_entries.values())} nonzero pair polys in {len(pair_entries)} target sectors ({time.time()-t0:.1f}s)', flush=True)
    quadrics = []  # list of dict {(k,l): coeff}
    n_coords = 0
    for key, ents in pair_entries.items():
        red = SectorReducer(S, degE, key, p)
        n_coords += red.dimR
        if red.dimR == 0:
            continue
        co = red.coords([P for (_, _, P) in ents])
        for c in range(red.dimR):
            qd = {}
            for (k, l, _), row in zip(ents, co):
                if row[c]:
                    qd[(k, l)] = row[c]
            if qd:
                quadrics.append(qd)
    print(f'[relax{d}] {len(quadrics)} nonzero quadrics (of {n_coords} class coordinates) in {K} unknowns ({time.time()-t0:.1f}s)', flush=True)
    # sanity: Psi(F) = 0 -- express F mod p in the basis (F lies in the trivial sector)
    kF = S.key(F)
    idxF = [i for i, (kap, _) in enumerate(basis) if kap == kF]
    Fp = core.reduce_mod(F, p)
    Dm = S.mons(d, kF)
    codes = np.sort(Dm @ POW)
    rowsF = []
    for i in idxF:
        v = [0] * len(codes)
        for e, c in basis[i][1].items():
            v[int(np.searchsorted(codes, W.enc(e)))] = c
        rowsF.append(v)
    vF = [0] * len(codes)
    for e, c in Fp.items():
        vF[int(np.searchsorted(codes, W.enc(e)))] = c
    Msol = nmod_mat(len(idxF) + 1, len(codes), [x for r in rowsF + [vF] for x in r], p)
    rkF = Msol.rank(); rk0 = nmod_mat(len(idxF), len(codes), [x for r in rowsF for x in r], p).rank()
    F_in_span = (rkF == rk0) if d == 6 else None
    # coefficients of F in the basis: solve transpose system
    aF = None
    if d == 6 and F_in_span:
        T = nmod_mat(len(codes), len(idxF) + 1, [rowsF[i][j] if i < len(idxF) else vF[j] for j in range(len(codes)) for i in range(len(idxF) + 1)], p)
        Tr, rk = T.rref()
        aF = [0] * K
        for i in range(rk):
            j = 0
            while j <= len(idxF) and int(Tr[i, j]) == 0:
                j += 1
            if j < len(idxF):
                aF[idxF[j]] = int(Tr[i, len(idxF)])
        # check Psi(F) = 0
        bad = 0
        for qd in quadrics:
            s = 0
            for (k, l), c in qd.items():
                s = (s + c * aF[k] * aF[l]) % p
            if s:
                bad += 1
        print(f'[relax{d}] F in K_6: {F_in_span}; quadrics not vanishing at F: {bad} (must be 0)', flush=True)
    # 4. Singular input
    sing = os.path.join(LANE, f'relax{d}_{p}.sing')
    with open(sing, 'w') as f:
        f.write(f'ring r = {p}, (a(1..{K})), dp;\n')
        f.write('ideal I;\n')
        for qd in quadrics:
            terms = []
            for (k, l), c in qd.items():
                terms.append(f'{c}*a({k+1})*a({l+1})')
            f.write('I[size(I)+1] = ' + '+'.join(terms) + ';\n')
        if aF is not None:
            # ideal of the F-line: all 2x2 minors a_i aF_j - a_j aF_i, with aF having a nonzero entry j0 -> linear forms a_i*aF_j0 - a_j0*aF_i
            j0 = max(range(K), key=lambda i: 1 if aF[i] else 0)
            f.write('ideal L;\n')
            for i in range(K):
                if i == j0:
                    continue
                f.write(f'L[size(L)+1] = {aF[j0]}*a({i+1}) - {aF[i]}*a({j0+1});\n')
        f.write('option(redSB);\nint t0 = timer;\nideal G = std(I);\n')
        f.write('"dim(std(I)) ="; dim(G); "mult ="; mult(G); "size(G) ="; size(G); "maxdeg of GB ="; deg(G[size(G)]);\n')
        if aF is not None:
            f.write('LIB "elim.lib";\nideal Gs = sat(G, maxideal(1))[1];\n"dim(sat) ="; dim(std(Gs)); "sat equals line ideal (1=yes) ="; (std(Gs) == std(L));\n')
            f.write('"line ideal contained in sat (1=yes) ="; size(reduce(L, std(Gs))) == 0;\n')
        f.write('"singular seconds ="; timer - t0;\nquit;\n')
    print(f'[relax{d}] wrote {sing}; running Singular (timeout {sing_timeout}s)', flush=True)
    t1 = time.time()
    try:
        r = subprocess.run(['timeout', str(sing_timeout), 'Singular', '-q', sing], capture_output=True, text=True)
        out = r.stdout
    except Exception as ex:
        out = f'error {ex}'
    print(out, flush=True)
    res = {'step': f'jacobian-relaxation-d{d}', 'd': d, 'p': p, 'dimK': K, 'sector_dims': {str(k): v for k, v in dims.items()},
           'n_quadrics': len(quadrics), 'n_class_coords': n_coords, 'F_in_span': F_in_span, 'singular_output': out[-1500:],
           'seconds': round(time.time() - t0, 1), 'singular_seconds': round(time.time() - t1, 1),
           'note': f'exact Jacobian relaxation of E2 at d={d}: quadrics Psi(D) = [B_D(N0)] in R_{degE} on the full E1 kernel (no symmetry assumption); Singular dim/mult/saturation vs the F-line ideal'}
    W.log(**res)
    with open(os.path.join(LANE, f'relax{d}_{p}.json'), 'w') as f:
        json.dump(res, f)


if __name__ == '__main__':
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    p = int(sys.argv[2]) if len(sys.argv) > 2 else 1000003
    to = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    main(d, p, to)
