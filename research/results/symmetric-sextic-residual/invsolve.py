"""Symmetry-reduced fixed-D solver for G = S3(X) x S3(Y)-invariant denominators.

If D, F, A are G-invariant, averaging a certificate field over G (v -> phi^* v, which maps
certificates to certificates) yields a G-invariant certificate over the SAME D.  Invariant
vector fields N satisfy N_{sigma(j)} = sigma . N_j, so N is determined by N_0 (invariant under
Stab(0) = S2(X1,X2) x S3(Y)) and N_3 (invariant under S3(X) x S2(Y1,Y2)); Q is G-invariant.
The cleared identity E = 0 is G-invariant, so one equation per G-orbit of monomials suffices.
Numpy-vectorised column construction, flint nmod_mat for linear algebra.
"""
import sys, time, itertools, json
import numpy as np
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from flint import nmod_mat

P0 = core.PRIMES[0]
BASE = 64

def to_arr(Pol, p):
    """dict poly -> (E int64 n x 6, C int64 n) mod p"""
    if not Pol:
        return np.zeros((0, 6), dtype=np.int64), np.zeros(0, dtype=np.int64)
    items = list(core.reduce_mod(Pol, p).items())
    E = np.array([e for e, _ in items], dtype=np.int64)
    C = np.array([c for _, c in items], dtype=np.int64)
    return E, C

def canon_codes(E):
    """G-orbit canonical code: sort X exps desc, Y exps desc, encode base 64."""
    Xs = -np.sort(-E[:, :3], axis=1)
    Ys = -np.sort(-E[:, 3:], axis=1)
    code = np.zeros(E.shape[0], dtype=np.int64)
    for i in range(3):
        code = code * BASE + Xs[:, i]
    for i in range(3):
        code = code * BASE + Ys[:, i]
    return code

def canon_H(E, fixed):
    """orbit code under stabiliser of index `fixed` (0 or 3)."""
    if fixed == 0:
        a = E[:, 0:1]; b = -np.sort(-E[:, 1:3], axis=1); c = -np.sort(-E[:, 3:], axis=1)
        parts = [a, b, c]
    else:
        a = -np.sort(-E[:, :3], axis=1); b = E[:, 3:4]; c = -np.sort(-E[:, 4:], axis=1)
        parts = [a, b, c]
    M = np.concatenate(parts, axis=1)
    code = np.zeros(E.shape[0], dtype=np.int64)
    for i in range(6):
        code = code * BASE + M[:, i]
    return code

def orbit_reps(deg, par, canon):
    mons = core.monomials(deg, par)
    if not mons:
        return []
    E = np.array(mons, dtype=np.int64)
    codes = canon(E)
    _, idx = np.unique(codes, return_index=True)
    return [tuple(int(x) for x in E[i]) for i in sorted(idx)]

def is_invariant(Pol):
    for e, c in Pol.items():
        for perm in [(1,0,2,3,4,5), (1,2,0,3,4,5), (0,1,2,4,3,5), (0,1,2,5,3,4)]:
            e2 = [0]*6
            for i in range(6):
                e2[perm[i]] = e[i]
            if Pol.get(tuple(e2), 0) != c:
                return False
    return True

def act(Pol, perm):
    """(sigma . P)(z) = P(z_{sigma(0)}, ..., z_{sigma(5)}); perm[i] = sigma(i).  Monomial z^e -> prod z_{sigma(i)}^{e_i}."""
    out = {}
    for e, c in Pol.items():
        e2 = [0]*6
        for i in range(6):
            e2[perm[i]] += e[i]
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}

class InvSystem:
    def __init__(self, F, A, q, D, p=P0, targets=False, verbose=False):
        t0 = time.time()
        assert q == 1, 'invariant solver written for q=1'
        assert is_invariant(D), 'D must be S3xS3-invariant'
        self.p = p; self.F = F; self.A = A; self.D = D; self.q = q
        d = core.degree(D); self.d = d
        pD = core.parity_of(D); pA = core.parity_of(A)
        self.pN = [(pA - pD) % 2, (pA - pD - 1) % 2]   # slot 0 (X) and slot 3 (Y)
        degN = d + 6*q - 5; degQ = 2*d - 6; degE = 2*d + 6*q
        Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        Di = [core.pderiv(Dp, i, p) for i in range(6)]
        FD = core.mul(Fp, Dp, p)
        # rows: G-orbits of degree-degE monomials of parity pA
        rowmons = np.array(core.monomials(degE, pA), dtype=np.int64)
        rcodes = canon_codes(rowmons)
        self.orbit_codes = np.unique(rcodes)
        self.nrows = len(self.orbit_codes)
        cols = []; colinfo = []
        def fold(E, C):
            codes = canon_codes(E)
            idx = np.searchsorted(self.orbit_codes, codes)
            v = np.bincount(idx, weights=C.astype(np.float64), minlength=self.nrows)
            return (np.rint(v).astype(np.int64)) % p
        for slot in (0, 3):
            G_s = core.sub(core.scale(core.mul(Dp, Fi[slot], p), p - 1, p), core.mul(Fp, Di[slot], p), p)  # -D F_s - F D_s
            GE, GC = to_arr(G_s, p)
            FDE, FDC = to_arr(FD, p)
            reps = orbit_reps(degN, self.pN[0 if slot == 0 else 1], lambda E: canon_H(E, slot))
            for m in reps:
                mv = np.array(m, dtype=np.int64)
                E1 = GE + mv; C1 = GC
                if m[slot]:
                    m2 = mv.copy(); m2[slot] -= 1
                    E2 = FDE + m2; C2 = (FDC * m[slot]) % p
                    E1 = np.concatenate([E1, E2]); C1 = np.concatenate([C1, C2])
                cols.append(fold(E1, C1)); colinfo.append(('N', slot, m))
        self.nN = len(cols)
        if degQ >= 0:
            F2 = core.mul(Fp, Fp, p)
            F2E, F2C = to_arr(core.scale(F2, p - 1, p), p)
            for m in orbit_reps(degQ, pA, canon_codes):
                cols.append(fold(F2E + np.array(m, dtype=np.int64), F2C)); colinfo.append(('Q', -1, m))
        self.nQ = len(cols) - self.nN
        self.nT = 0
        if targets:
            D2 = core.mul(Dp, Dp, p)
            D2E, D2C = to_arr(D2, p)
            for m in orbit_reps(6*q, 1, canon_codes):
                cols.append(fold(D2E + np.array(m, dtype=np.int64), D2C)); colinfo.append(('A', -1, m))
            self.nT = len(cols) - self.nN - self.nQ
        AD2 = core.mul(Ap, core.mul(Dp, Dp, p), p)
        rE, rC = to_arr(core.scale(AD2, p - 1, p), p)
        self.rhs = fold(rE, rC)
        self.cols = cols; self.colinfo = colinfo; self.ncols = len(cols)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[inv] d={d} rows(orbits)={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}, T {self.nT}) built {self.build_time:.1f}s', file=sys.stderr)

    def matrix(self, with_rhs=True):
        M = np.stack(self.cols + ([self.rhs] if with_rhs else []), axis=1)  # rows x cols
        return nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), self.p)

    def solve(self, verbose=False):
        t0 = time.time()
        M = self.matrix(True)
        Mr, rank = M.rref()
        C = self.ncols
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0]*C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        # rank of coefficient part
        res = {'consistent': consistent, 'rank': rank if consistent else rank - 1, 'nullity': C - (rank if consistent else rank - 1),
               'rows': self.nrows, 'cols': C, 'x': x, 'pivots': pivots, 't_rref': round(time.time() - t0, 2)}
        if verbose:
            print(f"[inv solve] consistent={consistent} rank={res['rank']} nullity={res['nullity']} ({res['t_rref']}s)", file=sys.stderr)
        return res

    def unfold(self, xs):
        """xs: list of exact Fractions (or ints mod p) per column -> (N list, Q) exact polys.
        Column ('N', slot, m) stands for the field Sym(m d_slot) = sum_sigma (sigma.m) d_{sigma(slot)}."""
        N = [{} for _ in range(6)]; Q = {}
        perms = [(a + b) for a in itertools.permutations((0, 1, 2)) for b in itertools.permutations((3, 4, 5))]
        for (kind, slot, m), v in zip(self.colinfo, xs):
            if not v: continue
            mono = {m: (Fraction(v) if not isinstance(v, Fraction) else v) / 36}
            if kind == 'N':
                for perm in perms:
                    N[perm[slot]] = core.add(N[perm[slot]], act(mono, perm))
            elif kind == 'Q':
                for perm in perms:
                    Q = core.add(Q, act(mono, perm))
        return N, Q

def e1_inv(F, A, q, D, p=P0):
    """Invariant E1 test: A D in J(F) via invariant N (equivalent to the unrestricted test)."""
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A)
    pN = [(pA - pD) % 2, (pA - pD - 1) % 2]
    degN = d + 6*q - 5
    Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    rowmons = np.array(core.monomials(d + 6*q, (pA + pD) % 2), dtype=np.int64)
    orbit_codes = np.unique(canon_codes(rowmons)); nrows = len(orbit_codes)
    def fold(E, C):
        idx = np.searchsorted(orbit_codes, canon_codes(E))
        v = np.bincount(idx, weights=C.astype(np.float64), minlength=nrows)
        return (np.rint(v).astype(np.int64)) % p
    cols = []
    for slot in (0, 3):
        FE, FC = to_arr(core.scale(Fi[slot], q, p), p)
        for m in orbit_reps(degN, pN[0 if slot == 0 else 1], lambda E: canon_H(E, slot)):
            cols.append(fold(FE + np.array(m, dtype=np.int64), FC))
    rE, rC = to_arr(core.mul(Ap, Dp, p), p)
    rhs = fold(rE, rC)
    M0 = np.stack(cols, axis=1); M1 = np.concatenate([M0, rhs[:, None]], axis=1)
    r0 = nmod_mat(M0.shape[0], M0.shape[1], M0.ravel().tolist(), p).rank()
    r1 = nmod_mat(M1.shape[0], M1.shape[1], M1.ravel().tolist(), p).rank()
    return r0 == r1, {'rank': r0, 'rank_aug': r1, 'rows': nrows, 'cols': M0.shape[1]}

def solve_inv(F, A, q, D, verbose=False, max_primes=6, lift=True):
    """Full invariant solve with rational lift.  Returns (status, data)."""
    acc = None; piv0 = None
    for k, p in enumerate(core.PRIMES[:max_primes]):
        S = InvSystem(F, A, q, D, p=p, verbose=verbose)
        r = S.solve(verbose=verbose)
        if not r['consistent']:
            return 'inconsistent', {'prime': p, 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols,
                                    't_build': round(S.build_time, 1), 't_rref': r['t_rref']}
        if not lift:
            return 'consistent-modp', {'prime': p, 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols,
                                        't_build': round(S.build_time, 1), 't_rref': r['t_rref']}
        if piv0 is None: piv0 = r['pivots']
        elif piv0 != r['pivots']:
            acc = None; piv0 = r['pivots']
        x = r['x']
        if acc is None: acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N, Q = S.unfold(rat)
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                return 'solved', {'N': N, 'Q': Q, 'primes_used': k + 1, 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols}
            elif verbose:
                print('[inv lift] reconstruction failed verification; continuing', file=sys.stderr)
    return 'lift-failed', {'primes_used': max_primes}

def served_invariant_targets(F, q, D, p=P0, verbose=False):
    """Dimension of the space of G-invariant Y-odd sextic targets served by D (mod p), and a basis."""
    A0 = {core.monomials(6*q, 1)[0]: 1}
    # need an invariant dummy A for parity bookkeeping only
    S = InvSystem(F, A0, q, D, p=p, targets=True, verbose=verbose)
    M = np.stack(S.cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    X, nullity = Mn.nullspace()
    tidx = [j for j, c in enumerate(S.colinfo) if c[0] == 'A']
    if nullity == 0:
        return 0, [], S
    proj = nmod_mat(len(tidx), nullity, [int(X[j, k]) for j in tidx for k in range(nullity)], p)
    T = proj.transpose(); Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        vec = {}
        for jj, j in enumerate(tidx):
            v = int(Tr[i, jj])
            if v: vec[S.colinfo[j][2]] = v if v < p // 2 else v - p
        basis.append(vec)
    return rk, basis, S

if __name__ == '__main__':
    # self-test: Fermat F0 with D = p6(X) must be consistent (subset-Euler certificate), D = QX inconsistent
    import sys
    sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
    from h import *
    F, A, q = load()
    F0 = core.add(p_k(X, 6), p_k(Y, 6))
    D = p_k(X, 6)
    t0 = time.time()
    st, data = solve_inv(F0, A, q, D, verbose=True)
    print('Fermat, D=p6X:', st, {k: v for k, v in data.items() if k not in ('N', 'Q')}, f'{time.time()-t0:.1f}s')
    if st == 'solved':
        print('   N term counts', [len(n) for n in data['N']], 'Q terms', len(data['Q']))
        print('   N0 =', core.to_str(data['N'][0])[:200])
    # cross-check against the full solver on our F with D = QX and D = p3Y + e3Y
    for D, lab in [(e_k(X, 2), 'QX'), (core.add(p_k(Y, 3), e_k(Y, 3)), 'p3Y+e3Y')]:
        t0 = time.time(); st, data = solve_inv(F, A, q, D, verbose=True); t1 = time.time()
        st2, data2 = fixed_d.solve_fixed_D(F, A, q, D, max_primes=1); t2 = time.time()
        print(f'{lab}: inv -> {st} ({t1-t0:.1f}s) ; full -> {st2} ({t2-t1:.1f}s)')
        ok, info = e1_inv(F, A, q, D); ok2, info2 = fixed_d.e1_solvable(F, A, q, D)
        print(f'   e1 inv={ok} {info} ; e1 full={ok2} {info2}')


# ---------------------------------------------------------------------------
# invariant E1 kernel: { invariant D of degree d : A D in J(F) } modulo J(F)_d^G
# ---------------------------------------------------------------------------
def orbit_sum(m):
    perms = [(a + b) for a in itertools.permutations((0, 1, 2)) for b in itertools.permutations((3, 4, 5))]
    s = set()
    for perm in perms:
        e2 = [0]*6
        for i in range(6):
            e2[perm[i]] = m[i]
        s.add(tuple(e2))
    return {e: Fraction(1) for e in s}

def sym(Pol):
    perms = [(a + b) for a in itertools.permutations((0, 1, 2)) for b in itertools.permutations((3, 4, 5))]
    out = {}
    for perm in perms:
        out = core.add(out, act(Pol, perm))
    return out

def e1_kernel_inv(F, A, q, d, pD, p=P0):
    pA = core.parity_of(A)
    pN = [(pA - pD) % 2, (pA - pD - 1) % 2]
    degN = d + 6*q - 5
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    rowmons = np.array(core.monomials(d + 6*q, (pA + pD) % 2), dtype=np.int64)
    orbit_codes = np.unique(canon_codes(rowmons)); nrows = len(orbit_codes)
    def fold(E, C):
        idx = np.searchsorted(orbit_codes, canon_codes(E))
        v = np.bincount(idx, weights=C.astype(np.float64), minlength=nrows)
        return (np.rint(v).astype(np.int64)) % p
    cols = []
    for slot in (0, 3):
        FE, FC = to_arr(core.scale(Fi[slot], q, p), p)
        for m in orbit_reps(degN, pN[0 if slot == 0 else 1], lambda E: canon_H(E, slot)):
            cols.append(fold(FE + np.array(m, dtype=np.int64), FC))
    nN = len(cols)
    Dreps = orbit_reps(d, pD, canon_codes)
    Dbasis = [orbit_sum(m) for m in Dreps]
    for b in Dbasis:
        E_, C_ = to_arr(core.mul(Ap, core.reduce_mod(b, p), p), p)
        cols.append(fold(E_, C_))
    nD = len(Dbasis)
    M = np.stack(cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    if nullity == 0:
        return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': nD, 'dimK': 0, 'dimJ': None, 'ker': 0}, []
    proj = nmod_mat(nD, nullity, [int(Xn[nN + j, k]) for j in range(nD) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    Kvecs = []
    for i in range(rk):
        v = [int(Tr[i, j]) for j in range(nD)]
        Kvecs.append(v)
    # invariant part of J_d: Sym(m F_slot), m of degree d-5
    Jvecs = []
    if d >= 5:
        for slot in (0, 3):
            par_m = (pD - (0 if slot == 0 else 1)) % 2
            for m in orbit_reps(d - 5, par_m, lambda E: canon_H(E, slot)):
                Pm = sym(core.shift(Fi[slot], m, 1, p))
                Pm = core.reduce_mod(Pm, p)
                Jvecs.append([Pm.get(r, 0) for r in Dreps])
    # kernel modulo J: rref [J | K]
    allv = Jvecs + Kvecs
    if not allv:
        return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': nD, 'dimK': rk, 'dimJ': 0, 'ker': rk}, [sum((core.scale(Dbasis[j], v[j]) for j in range(nD)), {}) for v in Kvecs]
    MJ = nmod_mat(nD, len(allv), [allv[c][r] for r in range(nD) for c in range(len(allv))], p)
    Mr, rk2 = MJ.rref()
    pivots = []
    for i in range(rk2):
        j = 0
        while j < MJ.ncols() and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    dimJ = sum(1 for j in pivots if j < len(Jvecs))
    extra_idx = [j - len(Jvecs) for j in pivots if j >= len(Jvecs)]
    polys = []
    for k in extra_idx:
        v = Kvecs[k]
        Pk = {}
        for j in range(nD):
            if v[j]:
                c = v[j] if v[j] < p // 2 else v[j] - p
                Pk = core.add(Pk, core.scale(Dbasis[j], Fraction(c)))
        polys.append(Pk)
    return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': nD, 'dimK': rk, 'dimJ': dimJ, 'ker': len(extra_idx)}, polys
