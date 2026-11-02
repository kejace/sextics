"""Character-twisted reduced fixed-D solver.

H = prod of symmetric groups on index blocks (e.g. [[0,1,2],[3,4,5]] = S3xS3, [[1,2],[4,5]] = S2xS2),
chi = sign character per block (+1/-1).  D must be chi-semi-invariant: sigma.D = chi(sigma) D.
Averaging a certificate over H (N -> chi(sigma) sigma.N, Q -> sigma.Q) shows WLOG N is chi-equivariant
as a vector field and Q is H-invariant.  The cleared identity E is then H-invariant, so one equation per
H-orbit of monomials suffices (plain fold).  The E1 subsystem N(F) = A D alone is chi-semi-invariant and
uses the sign-twisted fold (dead orbits = those whose stabiliser meets a chi=-1 block nontrivially).
"""
import sys, time, itertools
import numpy as np
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
from flint import nmod_mat

P0 = core.PRIMES[0]
BASE = 64

def to_arr(Pol, p):
    if not Pol:
        return np.zeros((0, 6), dtype=np.int64), np.zeros(0, dtype=np.int64)
    items = list(core.reduce_mod(Pol, p).items())
    return np.array([e for e, _ in items], dtype=np.int64), np.array([c for _, c in items], dtype=np.int64)

def act(Pol, perm):
    out = {}
    for e, c in Pol.items():
        e2 = [0]*6
        for i in range(6):
            e2[perm[i]] += e[i]
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}

class Group:
    def __init__(self, blocks, chi):
        self.blocks = [list(b) for b in blocks]
        self.chi = list(chi)
        self.perms = []   # list of (perm tuple, chi value)
        for choice in itertools.product(*[itertools.permutations(b) for b in self.blocks]):
            perm = list(range(6)); s = 1
            for b, pb, c in zip(self.blocks, choice, self.chi):
                for i, j in zip(b, pb):
                    perm[i] = j
                if c == -1:
                    s *= perm_sign(b, pb)
            self.perms.append((tuple(perm), s))
        self.order = len(self.perms)
        inb = set(i for b in self.blocks for i in b)
        self.slots = [b[0] for b in self.blocks] + [i for i in range(6) if i not in inb]
    def stab_blocks(self, s):
        return [[i for i in b if i != s] for b in self.blocks]
    def canon(self, E, blocks=None):
        blocks = self.blocks if blocks is None else blocks
        M = E.copy()
        for b in blocks:
            if len(b) > 1:
                M[:, b] = -np.sort(-E[:, b], axis=1)
        code = np.zeros(E.shape[0], dtype=np.int64)
        for i in range(6):
            code = code * BASE + M[:, i]
        return code
    def tsign(self, E):
        """sign chi(g_m) of the sorting permutation for chi=-1 blocks; 0 if dead (ties)."""
        s = np.ones(E.shape[0], dtype=np.int64)
        for b, c in zip(self.blocks, self.chi):
            if c != -1 or len(b) < 2: continue
            if len(b) == 2:
                a, bb = E[:, b[0]], E[:, b[1]]
                s = s * np.where(a == bb, 0, np.where(a >= bb, 1, -1))
            elif len(b) == 3:
                a, bb, cc = E[:, b[0]], E[:, b[1]], E[:, b[2]]
                inv = (a < bb).astype(np.int64) + (a < cc) + (bb < cc)
                tie = (a == bb) | (a == cc) | (bb == cc)
                s = s * np.where(tie, 0, np.where(inv % 2 == 0, 1, -1))
            else:
                raise NotImplementedError
        return s
    def is_semi_invariant(self, Pol):
        for perm, s in self.perms:
            if act(Pol, perm) != {e: s*c for e, c in Pol.items()}:
                return False
        return True
    def sym_chi(self, Pol):
        out = {}
        for perm, s in self.perms:
            out = core.add(out, core.scale(act(Pol, perm), s))
        return out
    def sym(self, Pol):
        out = {}
        for perm, s in self.perms:
            out = core.add(out, act(Pol, perm))
        return out

def perm_sign(b, pb):
    # sign of permutation mapping b[k] -> pb[k]
    pos = {v: k for k, v in enumerate(b)}
    arr = [pos[v] for v in pb]
    inv = sum(1 for i in range(len(arr)) for j in range(i+1, len(arr)) if arr[i] > arr[j])
    return -1 if inv % 2 else 1

def orbit_reps(deg, par, canon):
    mons = core.monomials(deg, par)
    if not mons: return []
    E = np.array(mons, dtype=np.int64)
    codes = canon(E)
    _, idx = np.unique(codes, return_index=True)
    return [tuple(int(x) for x in E[i]) for i in sorted(idx)]

def live_reps(deg, par, G):
    """reps r of H-orbits with Sym_chi(r) != 0, plus the signed orbit sum polynomial of each."""
    reps = orbit_reps(deg, par, G.canon)
    E = np.array(reps, dtype=np.int64) if reps else np.zeros((0, 6), dtype=np.int64)
    s = G.tsign(E) if reps else np.zeros(0)
    out = []
    for r, sg in zip(reps, s):
        if sg == 0: continue
        # signed orbit sum: sum over orbit of chi(g_m) m
        P = {}
        for perm, c in G.perms:
            m = act({r: 1}, perm)
            (e, _), = m.items()
            P[e] = c
        out.append((r, P))
    return out

def pN_of(pA, pD, s):
    return (pA + pD + (1 if s >= 3 else 0)) % 2

class CharSystem:
    def __init__(self, F, A, q, D, G, p=P0, targets=False, verbose=False):
        t0 = time.time()
        assert q == 1
        assert G.is_semi_invariant(D), 'D must be chi-semi-invariant under H'
        self.p = p; self.G = G; self.F = F; self.A = A; self.D = D; self.q = q
        d = core.degree(D); self.d = d
        pD = core.parity_of(D); pA = core.parity_of(A)
        degN = d + 6*q - 5; degQ = 2*d - 6; degE = 2*d + 6*q
        Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
        Fi = [core.pderiv(Fp, i, p) for i in range(6)]
        Di = [core.pderiv(Dp, i, p) for i in range(6)]
        FD = core.mul(Fp, Dp, p)
        rowmons = np.array(core.monomials(degE, pA), dtype=np.int64)
        self.orbit_codes = np.unique(G.canon(rowmons)); self.nrows = len(self.orbit_codes)
        def fold(E, C):
            idx = np.searchsorted(self.orbit_codes, G.canon(E))
            v = np.bincount(idx, weights=C.astype(np.float64), minlength=self.nrows)
            return (np.rint(v).astype(np.int64)) % p
        cols = []; colinfo = []
        FDE, FDC = to_arr(FD, p)
        for s in G.slots:
            G_s = core.sub(core.scale(core.mul(Dp, Fi[s], p), p - 1, p), core.mul(Fp, Di[s], p), p)
            GE, GC = to_arr(G_s, p)
            sb = G.stab_blocks(s)
            for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
                mv = np.array(m, dtype=np.int64)
                E1 = GE + mv; C1 = GC
                if m[s]:
                    m2 = mv.copy(); m2[s] -= 1
                    E1 = np.concatenate([E1, FDE + m2]); C1 = np.concatenate([C1, (FDC * m[s]) % p])
                col = fold(E1, C1)
                cols.append(col); colinfo.append(('N', s, m))
        self.nN = len(cols)
        if degQ >= 0:
            F2E, F2C = to_arr(core.scale(core.mul(Fp, Fp, p), p - 1, p), p)
            for m in orbit_reps(degQ, pA, G.canon):
                cols.append(fold(F2E + np.array(m, dtype=np.int64), F2C)); colinfo.append(('Q', -1, m))
        self.nQ = len(cols) - self.nN
        self.nT = 0
        if targets:
            D2E, D2C = to_arr(core.mul(Dp, Dp, p), p)
            for m in orbit_reps(6*q, 1, G.canon):
                cols.append(fold(D2E + np.array(m, dtype=np.int64), D2C)); colinfo.append(('A', -1, m))
            self.nT = len(cols) - self.nN - self.nQ
        rE, rC = to_arr(core.scale(core.mul(Ap, core.mul(Dp, Dp, p), p), p - 1, p), p)
        self.rhs = fold(rE, rC)
        self.cols = cols; self.colinfo = colinfo; self.ncols = len(cols)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[char] d={d} rows={self.nrows} cols={self.ncols} (N {self.nN}, Q {self.nQ}, T {self.nT}) built {self.build_time:.1f}s', file=sys.stderr)

    def solve(self, verbose=False):
        t0 = time.time()
        M = np.stack(self.cols + [self.rhs], axis=1)
        Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), self.p)
        Mr, rank = Mn.rref()
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
            for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
        rk = rank if consistent else rank - 1
        res = {'consistent': consistent, 'rank': rk, 'nullity': C - rk, 'rows': self.nrows, 'cols': C, 'x': x,
               'pivots': pivots, 't_rref': round(time.time() - t0, 2)}
        if verbose: print(f"[char solve] consistent={consistent} rank={rk} nullity={C-rk} ({res['t_rref']}s)", file=sys.stderr)
        return res

    def unfold(self, xs):
        G = self.G
        N = [{} for _ in range(6)]; Q = {}; Ap = {}
        for (kind, s, m), v in zip(self.colinfo, xs):
            if not v: continue
            mono = {m: Fraction(v) / G.order}
            for perm, c in G.perms:
                if kind == 'N': N[perm[s]] = core.add(N[perm[s]], core.scale(act(mono, perm), c))
                elif kind == 'Q': Q = core.add(Q, act(mono, perm))
                else: Ap = core.add(Ap, act(mono, perm))
        return N, Q, Ap

    def served_targets(self):
        """nullspace projected on A'-coords -> list of H-invariant target polys (exact Sym normalisation) mod p."""
        M = np.stack(self.cols, axis=1)
        Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), self.p)
        Xn, nullity = Mn.nullspace()
        tidx = [j for j, c in enumerate(self.colinfo) if c[0] == 'A']
        if nullity == 0: return []
        proj = nmod_mat(len(tidx), nullity, [int(Xn[j, k]) for j in tidx for k in range(nullity)], self.p)
        Tr, rk = proj.transpose().rref()
        out = []
        for i in range(rk):
            P = {}
            for jj, j in enumerate(tidx):
                v = int(Tr[i, jj])
                if v: P = core.add(P, core.scale(self.G.sym({self.colinfo[j][2]: 1}), v), self.p)
            out.append(core.reduce_mod(P, self.p))
        return out

def solve_char(F, A, q, D, G, verbose=False, max_primes=6, lift=True):
    acc = None; piv0 = None
    for k, p in enumerate(core.PRIMES[:max_primes]):
        S = CharSystem(F, A, q, D, G, p=p, verbose=verbose)
        r = S.solve(verbose=verbose)
        base = {'prime': p, 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols,
                't_build': round(S.build_time, 1), 't_rref': r['t_rref']}
        if not r['consistent']: return 'inconsistent', base
        if not lift: return 'consistent-modp', base
        if piv0 is None: piv0 = r['pivots']
        elif piv0 != r['pivots']: acc = None; piv0 = r['pivots']
        x = r['x']
        if acc is None: acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m * p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            N, Q, _ = S.unfold(rat)
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok: return 'solved', dict(base, N=N, Q=Q, primes_used=k + 1)
    return 'lift-failed', {'primes_used': max_primes}

# ---------------------------------------------------------------- twisted E1
def _tfold_factory(G, deg, par, p):
    rowmons = np.array(core.monomials(deg, par), dtype=np.int64)
    codes = G.canon(rowmons); signs = G.tsign(rowmons)
    live = np.unique(codes[signs != 0]); n = len(live)
    def tfold(E, C):
        s = G.tsign(E); keep = s != 0
        idx = np.searchsorted(live, G.canon(E[keep]))
        v = np.bincount(idx, weights=(C[keep] * s[keep]).astype(np.float64), minlength=n)
        return (np.rint(v).astype(np.int64)) % p
    return tfold, n

def e1_char(F, A, q, D, G, p=P0):
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A)
    degN = d + 6*q - 5
    Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    tfold, nrows = _tfold_factory(G, d + 6*q, (pA + pD) % 2, p)
    cols = []
    for s in G.slots:
        FE, FC = to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            cols.append(tfold(FE + np.array(m, dtype=np.int64), FC))
    rE, rC = to_arr(core.mul(Ap, Dp, p), p)
    rhs = tfold(rE, rC)
    M0 = np.stack(cols, axis=1); M1 = np.concatenate([M0, rhs[:, None]], axis=1)
    r0 = nmod_mat(M0.shape[0], M0.shape[1], M0.ravel().tolist(), p).rank()
    r1 = nmod_mat(M1.shape[0], M1.shape[1], M1.ravel().tolist(), p).rank()
    return r0 == r1, {'rank': r0, 'rank_aug': r1, 'rows': nrows, 'cols': M0.shape[1]}

def J_char_basis(F, d, pD, G, p=None):
    """spanning set of J_d^chi: Sym_chi(m dF_s), m of degree d-5."""
    out = []
    if d < 5: return out
    Fi = [core.pderiv(F, i, p) for i in range(6)]
    for s in G.slots:
        sb = G.stab_blocks(s)
        for m in orbit_reps(d - 5, (pD + (1 if s >= 3 else 0)) % 2, lambda E: G.canon(E, sb)):
            P = G.sym_chi(core.shift(Fi[s], m, 1, p))
            if P: out.append(core.reduce_mod(P, p) if p else P)
    return out

def e1_kernel_char(F, A, q, d, pD, G, p=P0):
    """{chi-type D of degree d : A D in J} modulo J_d^chi. Returns info, list of polys (signed ints mod p)."""
    pA = core.parity_of(A); degN = d + 6*q - 5
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    tfold, nrows = _tfold_factory(G, d + 6*q, (pA + pD) % 2, p)
    cols = []
    for s in G.slots:
        FE, FC = to_arr(core.scale(Fi[s], q, p), p); sb = G.stab_blocks(s)
        for m in orbit_reps(degN, pN_of(pA, pD, s), lambda E: G.canon(E, sb)):
            cols.append(tfold(FE + np.array(m, dtype=np.int64), FC))
    nN = len(cols)
    Dl = live_reps(d, pD, G)
    Dreps = [r for r, _ in Dl]; Dbasis = [P for _, P in Dl]
    for P in Dbasis:
        E_, C_ = to_arr(core.mul(Ap, core.reduce_mod(P, p), p), p)
        cols.append(tfold(E_, C_))
    nD = len(Dbasis)
    if nD == 0:
        return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': 0, 'dimK': 0, 'dimJ': 0, 'ker': 0}, []
    M = np.stack(cols, axis=1)
    Mn = nmod_mat(M.shape[0], M.shape[1], M.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    if nullity == 0:
        return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': nD, 'dimK': 0, 'dimJ': None, 'ker': 0}, []
    proj = nmod_mat(nD, nullity, [int(Xn[nN + j, k]) for j in range(nD) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    Kvecs = [[int(Tr[i, j]) for j in range(nD)] for i in range(rk)]
    Jpolys = J_char_basis(Fp, d, pD, G, p)
    Jvecs = [[Pj.get(r, 0) for r in Dreps] for Pj in Jpolys]   # coefficient on rep = coordinate in signed-orbit-sum basis
    allv = Jvecs + Kvecs
    MJ = nmod_mat(nD, len(allv), [allv[c][r] for r in range(nD) for c in range(len(allv))], p)
    Mr, rk2 = MJ.rref()
    pivots = []
    for i in range(rk2):
        j = 0
        while j < MJ.ncols() and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    dimJ = sum(1 for j in pivots if j < len(Jvecs))
    extra = [j - len(Jvecs) for j in pivots if j >= len(Jvecs)]
    polys = []
    for k in extra:
        v = Kvecs[k]; Pk = {}
        for j in range(nD):
            if v[j]:
                c = v[j] if v[j] < p // 2 else v[j] - p
                Pk = core.add(Pk, core.scale(Dbasis[j], Fraction(c)))
        polys.append(Pk)
    return {'d': d, 'pD': pD, 'rows': nrows, 'nN': nN, 'nD': nD, 'dimK': rk, 'dimJ': dimJ, 'ker': len(extra)}, polys

# ---------------------------------------------------------------- full (unreduced) sector kernel
def sector_kernel(F, q, D, sector, p=P0):
    """dim and basis of {A' in span(sector) : A' D in J(F)} (sector = list of Y-odd degree-6q polys)."""
    d = core.degree(D); pD = core.parity_of(D); pA = 1
    degN = d + 6*q - 5
    Fp = core.reduce_mod(F, p); Dp = core.reduce_mod(D, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    rowmons = core.monomials(d + 6*q, (pA + pD) % 2)
    ridx = {m: i for i, m in enumerate(rowmons)}; R = len(rowmons)
    cols = []
    for s in range(6):
        FE, FC = to_arr(core.scale(Fi[s], q, p), p)
        for m in core.monomials(degN, pN_of(pA, pD, s)):
            E = FE + np.array(m, dtype=np.int64)
            cols.append((E, FC))
    nN = len(cols)
    for Ap in sector:
        cols.append(to_arr(core.mul(core.reduce_mod(Ap, p), Dp, p), p))
    C = len(cols)
    flat = np.zeros((R, C), dtype=np.int64)
    rcode = {m: i for i, m in enumerate(rowmons)}
    # vectorised row lookup via codes
    allcodes = np.array([sum(e * BASE**(5-i) for i, e in enumerate(m)) for m in rowmons], dtype=np.int64)
    order = np.argsort(allcodes); sorted_codes = allcodes[order]
    for j, (E, Cc) in enumerate(cols):
        if E.shape[0] == 0: continue
        codes = np.zeros(E.shape[0], dtype=np.int64)
        for i in range(6): codes = codes * BASE + E[:, i]
        pos = order[np.searchsorted(sorted_codes, codes)]
        np.add.at(flat[:, j], pos, Cc)
    flat %= p
    Mn = nmod_mat(R, C, flat.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    nS = len(sector)
    if nullity == 0: return 0, []
    proj = nmod_mat(nS, nullity, [int(Xn[nN + j, k]) for j in range(nS) for k in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    return rk, [[int(Tr[i, j]) for j in range(nS)] for i in range(rk)]

if __name__ == '__main__':
    sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
    from h import *
    import invsolve
    F, A, q = load()
    GG = Group([[0,1,2],[3,4,5]], [1, 1])
    GsX = Group([[0,1,2],[3,4,5]], [-1, 1])
    GsY = Group([[0,1,2],[3,4,5]], [1, -1])
    Gss = Group([[0,1,2],[3,4,5]], [-1, -1])
    H22 = Group([[1,2],[4,5]], [1, 1])
    mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
    def radial(F, idx):
        out = {}
        for i in idx: out = add(out, mul(var(i), core.pderiv(F, i)))
        return out
    EXF = radial(F, (0,1,2))
    def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
    DX, DY = dlt(X), dlt(Y)
    QX = e_k(X, 2); e1X = e_k(X, 1); e3X = e_k(X, 3)
    e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2))
    D7 = mul(DX, e2sqX)
    # (1) invariant reproduction
    t0 = time.time(); S = CharSystem(F, A, q, EXF, GG); r = S.solve()
    print('EX(F) invariant via CharSystem:', r['consistent'], r['rank'], r['nullity'], r['rows'], r['cols'], f'{time.time()-t0:.1f}s  (invsolve: rank 140 null 21 rows 838 cols 161)')
    # (1b) same D with the S2xS2 subgroup
    t0 = time.time(); S = CharSystem(F, A, q, EXF, H22); r = S.solve()
    print('EX(F) via S2xS2:', r['consistent'], r['rank'], r['nullity'], r['rows'], r['cols'], f'{time.time()-t0:.1f}s')
    # (2) twisted E1 on D7 (sgn_X type) must be consistent; invariant E1 on EX(F) consistent; DX (sgn_X, deg 3) inconsistent
    print('D7 semi-inv sgnX:', GsX.is_semi_invariant(D7), ' e1_char:', e1_char(F, A, q, D7, GsX))
    print('DX e1_char (expect False):', e1_char(F, A, q, DX, GsX))
    print('DX*DY sgn⊗sgn e1_char (expect False):', e1_char(F, A, q, mul(DX, DY), Gss))
    print('EX(F) e1_char inv:', e1_char(F, A, q, EXF, GG))
    # (3) twisted kernel at d=7 sgnX even: expect 1 (D7)
    for G_, lab in [(GsX, 'sgnX'), (GsY, 'sgnY'), (Gss, 'sgnsgn'), (GG, 'triv')]:
        for pD in (0, 1):
            t0 = time.time(); info, polys = e1_kernel_char(F, A, q, 7, pD, G_)
            print(f'kernel d=7 pD={pD} {lab}:', info, f'{time.time()-t0:.1f}s')
    # (4) full solve of D7 with sgnX reduction (expect inconsistent)
    t0 = time.time(); S = CharSystem(F, A, q, D7, GsX, verbose=True); r = S.solve(verbose=True)
    print('D7 reduced sgnX solve:', r['consistent'], r['rank'], r['nullity'], r['rows'], r['cols'], f'{time.time()-t0:.1f}s')
    # (5) sector kernel for QX (expect 2 = J part) and D7
    sector = [invsolve.orbit_sum(m) for m in invsolve.orbit_reps(6, 1, invsolve.canon_codes)]
    t0 = time.time(); rk, _ = sector_kernel(F, q, QX, sector); print('sector_kernel QX:', rk, f'{time.time()-t0:.1f}s')
    t0 = time.time(); rk, _ = sector_kernel(F, q, D7, sector); print('sector_kernel D7:', rk, f'{time.time()-t0:.1f}s')
