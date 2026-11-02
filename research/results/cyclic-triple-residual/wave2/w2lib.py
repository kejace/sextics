"""Wave-2 search library for cyclic-triple-residual (q = 1, A = X0^3 Y0^3).

Character-sector-reduced, E1-eliminated variable projection over an E1-admissible linear family of denominators.

Group: projective diagonal group G of F (order 18, generators mu3 = (0,0,2,0,0,0) and g = (1,5,1,2,2,0)),
characters chi = (a, b) with a in {0,2,4} (mu3), b in Z/6 (g); Y-parity = b mod 2.  chi(A) = (0,3).

Split certificate system (q = 1), D semi-invariant of character c:
   E1:  sum_i N_i F_i = A D            rows: degree d+6, character chiA + c
   E2:  D div N - N(D) = F Q           rows: degree 2d,  character chiA + 2c
   chi(N_i) = chiA + c + e_i,  chi(Q) = chiA + 2c,  deg N_i = d+1, deg Q = 2d-6.

For a family D(a) = sum_k a_k M_k with every M_k E1-admissible (A M_k in J(F)) we eliminate E1 numerically:
   N(a, z) = N0(a) + K z,  N0(a) = sum a_k N0_k (particular lifts), K = null space of the E1 map,
and minimize   f(a) = min_{z,Q} || E2_{D(a)}(N0(a) + K z) - F Q ||^2   on the unit sphere of the complement
of {a : D(a) in F S_{d-6}} (spurious trivial solutions).  Q is eliminated by projecting onto the orthogonal
complement of the column space of (F * S_{2d-6}).  Gradient by Golub-Pereyra (exact at the optimum z).

Hessian relaxation: z = 0 (N = A w for Jacobian elements M_k = w . grad F), Q free.
"""
from __future__ import annotations
import sys, os, time, json, itertools, math
from fractions import Fraction
from math import lcm
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import minimize
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from hodge import core
from flint import nmod_mat, fmpz_mat, fmpq_mat
import sect

LANE = 'cyclic-triple-residual'
WD = '/tmp/hodge/results/cyclic-triple-residual'
W2 = WD + '/wave2'
LOG = WD + '/wave2-search.log.jsonl'
PRIMES = core.PRIMES
B = 32  # exponent code base (degree <= 31)
E6 = [tuple(int(j == i) for j in range(6)) for i in range(6)]


def log_batch(entry):
    entry = dict(entry)
    entry.setdefault('engine', 2)
    entry.setdefault('lane', LANE)
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')


# ---------------------------------------------------------------- instance / group
class Inst:
    def __init__(self, iid=LANE, F=None, A=None):
        if F is None:
            F, A, q, rec = core.load_instance(iid)
        else:
            q = 1
        self.F, self.A, self.q = F, A, 1
        self.Fz, self.Az = sect.integral(F), sect.integral(A)
        G = sect.semi_invariant_subgroup(sect.elems(self.Fz), [self.Az])
        self.G = G
        self.gens = sect.generators(G)
        self.chi = sect.make_chi(self.gens)
        self.chiA = sect.chi_of(self.Az, self.chi)
        self.Fi = [core.pderiv(F, i) for i in range(6)]
        self.chars = sorted({self.chi(m) for m in sect.mons(8)})

    def cadd(self, *cs):
        return sect.cadd(*cs)

    def cneg(self, c):
        return tuple((-x) % 6 for x in c)

    def sector_mons(self, deg, c):
        if deg < 0:
            return []
        return [m for m in sect.mons(deg) if self.chi(m) == c]

    def chi_ei(self, i):
        return self.chi(E6[i])

    # ------------------------------------------------------------ E1-admissible space in one sector
    def jac_generators(self, d, c):
        """(i, m) with chi(m F_i) = c, deg m = d-5 (chi(F_i) = -e_i)."""
        out = []
        if d < 5:
            return out
        for i in range(6):
            ci = self.cadd(c, self.chi_ei(i))
            for m in sect.mons(d - 5):
                if self.chi(m) == ci:
                    out.append((i, m))
        return out

    def jac_poly(self, i, m):
        return core.mul({m: Fraction(1)}, self.Fi[i])

    def e1_sector(self, d, c, p=PRIMES[0], want_vectors=True):
        """mod p: dims and kernel of [A]: S_d^c -> R_{d+6}, reduced mod J_d^c (canonical rref rows),
        plus an independent subset of the J_d^c generators (indices into jac_generators)."""
        t0 = time.time()
        Dm = self.sector_mons(d, c)
        idx = {m: j for j, m in enumerate(Dm)}
        Az = self.Az
        dF = [sect.ideriv(self.Fz, i) for i in range(6)]
        target = self.cadd(self.chiA, c)
        cols = []
        for m in Dm:
            cols.append({tuple(a + b for a, b in zip(m, e)): v % p for e, v in Az.items()})
        for i in range(6):
            ci = self.cadd(target, self.chi_ei(i))
            for m in sect.mons(d + 1):
                if self.chi(m) != ci:
                    continue
                col = {}
                for e, v in dF[i].items():
                    k = tuple(a + b for a, b in zip(m, e)); col[k] = (col.get(k, 0) + v) % p
                cols.append(col)
        rows = {}
        for col in cols:
            for e in col:
                rows.setdefault(e, len(rows))
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, col in enumerate(cols):
            for e, v in col.items():
                flat[rows[e] * C + j] = v
        M = nmod_mat(R, C, flat, p)
        X, nullity = M.nullspace()
        nA = len(Dm)
        kerA = 0
        vecs = []
        if nullity and nA:
            proj = nmod_mat(nA, nullity, [int(X[j, k]) for j in range(nA) for k in range(nullity)], p)
            Pr, kerA = proj.transpose().rref()
            kvecs = [[int(Pr[i, j]) for j in range(nA)] for i in range(kerA)]
        else:
            kvecs = []
        # J_d^c generators and an independent subset
        gens = self.jac_generators(d, c)
        Jrows = []
        for (i, m) in gens:
            row = [0] * nA
            for e, v in dF[i].items():
                k = tuple(a + b for a, b in zip(m, e)); row[idx[k]] = (row[idx[k]] + v) % p
            Jrows.append(row)
        dimJ = 0
        indep = []
        piv_rows = []
        if Jrows:
            J = nmod_mat(len(Jrows), nA, [x for r in Jrows for x in r], p)
            dimJ = J.rank()
            # independent subset: rref of transpose -> pivot columns = independent generators
            Jt, rk = J.transpose().rref()
            for i in range(rk):
                j = 0
                while j < len(Jrows) and int(Jt[i, j]) == 0:
                    j += 1
                indep.append(j)
            Jr, rk2 = J.rref()
            piv_rows = [[int(Jr[i, j]) for j in range(nA)] for i in range(rk2)]
        # reduce kernel vectors modulo J (row space), keep independent remainder, canonical rref
        red = []
        for v in kvecs:
            v = list(v)
            for row in piv_rows:
                j0 = next(j for j in range(nA) if row[j])
                if v[j0]:
                    f = v[j0] * pow(row[j0], -1, p) % p
                    v = [(x - f * y) % p for x, y in zip(v, row)]
            if any(v):
                red.append(v)
        kmodJ = 0
        reps = []
        if red:
            Mo = nmod_mat(len(red), nA, [x for r in red for x in r], p)
            Mr, kmodJ = Mo.rref()
            reps = [[int(Mr[i, j]) for j in range(nA)] for i in range(kmodJ)]
        return {'d': d, 'c': c, 'p': p, 'dim_S': nA, 'kerA_in_S': kerA, 'dimJ': dimJ, 'ker_mod_J': kmodJ,
                'rows': R, 'cols': C, 'gens': gens, 'indep': indep, 'Dm': Dm,
                'reps': reps if want_vectors else None, 'seconds': round(time.time() - t0, 2)}

    def e1_space(self, d, c, nprimes=4, verbose=False):
        """Exact E1-admissible family in sector c, degree d: list of exact polys = independent J generators
        + rationally reconstructed kernel representatives (verified mod an extra prime).  Returns dict."""
        info = self.e1_sector(d, c, PRIMES[0])
        Dm = info['Dm']
        Jb = [self.jac_poly(i, m) for (i, m) in (info['gens'][j] for j in info['indep'])]
        reps = []
        if info['ker_mod_J']:
            res = [info]
            for p in PRIMES[1:nprimes]:
                res.append(self.e1_sector(d, c, p))
            k = info['ker_mod_J']
            assert all(r['ker_mod_J'] == k for r in res), 'kernel dim differs between primes'
            nA = len(Dm)
            for r in range(k):
                vals = []
                for j in range(nA):
                    a, m = res[0]['reps'][r][j], res[0]['p']
                    for rr in res[1:]:
                        a, m = core.crt_pair(a, m, rr['reps'][r][j], rr['p'])
                    x = core.ratrecon(a, m)
                    if x is None:
                        raise RuntimeError(f'ratrecon failed d={d} c={c} rep {r} coeff {j}')
                    vals.append(x)
                D = {Dm[j]: vals[j] for j in range(nA) if vals[j]}
                # verify with an extra prime: A D in J_{d+6} mod p
                pv = PRIMES[nprimes]
                ok, _ = self.e1_check_modp(D, pv)
                if not ok:
                    raise RuntimeError(f'kernel rep failed verification d={d} c={c} rep {r}')
                reps.append(D)
        out = {'d': d, 'c': list(c), 'dim_S': info['dim_S'], 'dimJ': info['dimJ'], 'kerA_in_S': info['kerA_in_S'],
               'ker_mod_J': info['ker_mod_J'], 'J_indep': [info['gens'][j] for j in info['indep']],
               'reps': [core.to_str(D) for D in reps], 'seconds': info['seconds']}
        return out, Jb, reps

    def e1_check_modp(self, D, p):
        """Sector-wise check that A D lies in J(F)_{d+6} mod p (rank test on the sector matrix)."""
        c = self.chi_of(D)
        d = core.degree(D)
        c1 = self.cadd(self.chiA, c)
        dF = [sect.ideriv(self.Fz, i) for i in range(6)]
        AD = sect.imul(self.Az, sect.integral(D))
        cols = []
        for i in range(6):
            ci = self.cadd(c1, self.chi_ei(i))
            for m in sect.mons(d + 1):
                if self.chi(m) != ci:
                    continue
                col = {}
                for e, v in dF[i].items():
                    k = tuple(a + b for a, b in zip(m, e)); col[k] = (col.get(k, 0) + v) % p
                cols.append(col)
        rows = {}
        for col in cols:
            for e in col:
                rows.setdefault(e, len(rows))
        for e in AD:
            rows.setdefault(e, len(rows))
        R, C = len(rows), len(cols)
        flat = [0] * (R * (C + 1))
        for j, col in enumerate(cols):
            for e, v in col.items():
                flat[rows[e] * (C + 1) + j] = v
        M0 = nmod_mat(R, C + 1, flat, p)
        r0 = nmod_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)], p).rank()
        for e, v in AD.items():
            flat[rows[e] * (C + 1) + C] = v % p
        r1 = nmod_mat(R, C + 1, flat, p).rank()
        return r0 == r1, {'rank': r0, 'rank_aug': r1}

    def f_multiples(self, d, c):
        """Basis of F * S_{d-6} inside sector c (exact polys)."""
        if d < 6:
            return []
        cF = self.chi_of(self.F)
        cG = self.cadd(c, self.cneg(cF))
        return [core.mul({m: Fraction(1)}, self.F) for m in self.sector_mons(d - 6, cG)]

    def chi_of(self, P):
        return sect.chi_of(sect.integral(P), self.chi)


# ---------------------------------------------------------------- numerics helpers
def code(m):
    return sum(int(e) * B ** j for j, e in enumerate(m))


def decode(x):
    out = []
    for _ in range(6):
        out.append(x % B); x //= B
    return tuple(out)


def poly_float(P):
    return {e: float(v) for e, v in P.items()}


def poly_norm(P):
    return math.sqrt(sum(float(v) ** 2 for v in P.values()))


class SectorVP:
    """E1-eliminated sector variable projection for family basis (exact polys of one character c, one degree)."""

    def __init__(self, inst: Inst, basis, c, relax=False, exclude_F=True, verbose=False, cond=1e-10,
                 w_of_basis=None, exact_rank_check=True):
        t0 = time.time()
        self.inst = inst
        self.c = tuple(c)
        self.relax = relax
        F, A = inst.F, inst.A
        self.basis = [dict(Bk) for Bk in basis]
        K = len(basis); self.K = K
        degs = {core.degree(Bk) for Bk in basis}
        assert len(degs) == 1, degs
        d = degs.pop(); self.d = d
        chiA = inst.chiA
        c1 = inst.cadd(chiA, self.c); c2 = inst.cadd(chiA, self.c, self.c)
        self.c1, self.c2 = c1, c2
        # unknown columns
        Ncols = []
        for i in range(6):
            ci = inst.cadd(chiA, self.c, inst.chi_ei(i))
            for m in sect.mons(d + 1):
                if inst.chi(m) == ci:
                    Ncols.append((i, m))
        self.Ncols = Ncols; nN = len(Ncols); self.nN = nN
        Qm = inst.sector_mons(2 * d - 6, c2) if 2 * d - 6 >= 0 else []
        self.Qm = Qm; nQ = len(Qm); self.nQ = nQ
        r1 = inst.sector_mons(d + 6, c1); r2 = inst.sector_mons(2 * d, c2)
        self.r1, self.r2 = r1, r2
        ri1 = {code(m): k for k, m in enumerate(r1)}; ri2 = {code(m): k for k, m in enumerate(r2)}
        self.n1, self.n2 = len(r1), len(r2)
        # E1 matrix (dense) : column (i,m) = m F_i
        M1 = np.zeros((self.n1, nN))
        Fi_f = [poly_float(inst.Fi[i]) for i in range(6)]
        for j, (i, m) in enumerate(Ncols):
            cm = code(m)
            for e, v in Fi_f[i].items():
                M1[ri1[cm + code(e)], j] += v
        self.M1 = M1
        # E2_k sparse matrices: column (i,m): sum over terms (e,v) of M_k: row e+m-eps_i, val v (m_i - e_i)
        col_codes = np.array([code(m) for (i, m) in Ncols], dtype=np.int64)
        col_i = np.array([i for (i, m) in Ncols], dtype=np.int64)
        col_mi = np.array([m[i] for (i, m) in Ncols], dtype=np.int64)
        self.E2 = []
        r2codes = np.array(sorted(ri2), dtype=np.int64)
        r2index = np.array([ri2[x] for x in r2codes], dtype=np.int64)
        for Bk in self.basis:
            es = list(Bk.items())
            ecodes = np.array([code(e) for e, v in es], dtype=np.int64)
            evals = np.array([float(v) for e, v in es])
            rows_all = []; cols_all = []; vals_all = []
            for t_i in range(6):
                sel = np.where(col_i == t_i)[0]
                if sel.size == 0:
                    continue
                ei = np.array([e[t_i] for e, v in es], dtype=np.int64)
                # outer: columns sel x terms
                coef = (col_mi[sel][:, None] - ei[None, :]) * evals[None, :]
                rc = col_codes[sel][:, None] + ecodes[None, :] - B ** t_i
                mask = (coef != 0) & ((col_mi[sel][:, None] + ei[None, :]) > 0)
                rr = rc[mask]; cc = np.broadcast_to(sel[:, None], rc.shape)[mask]; vv = coef[mask]
                pos = np.searchsorted(r2codes, rr)
                assert np.all(r2codes[pos] == rr), 'E2 row outside sector'
                rows_all.append(r2index[pos]); cols_all.append(cc); vals_all.append(vv)
            if rows_all:
                rows_c = np.concatenate(rows_all); cols_c = np.concatenate(cols_all); vals_c = np.concatenate(vals_all)
            else:
                rows_c = np.zeros(0, dtype=np.int64); cols_c = rows_c.copy(); vals_c = np.zeros(0)
            self.E2.append(sp.csr_matrix((vals_c, (rows_c, cols_c)), shape=(self.n2, nN)))
        # Q columns: -F m
        Ff = poly_float(F)
        MQ = np.zeros((self.n2, nQ))
        for j, m in enumerate(Qm):
            cm = code(m)
            for e, v in Ff.items():
                MQ[ri2[cm + code(e)], j] -= v
        self.MQ = MQ
        # orthonormal complement of col(MQ)
        if nQ:
            Qf, Rf = np.linalg.qr(MQ, mode='complete')
            s = np.abs(np.diag(Rf)) if Rf.shape[0] and Rf.shape[1] else np.zeros(0)
            rQ = int(np.sum(s > 1e-10 * max(1.0, s.max() if s.size else 1.0)))
            self.rQ = rQ
            self.Q2 = Qf[:, rQ:]
        else:
            self.rQ = 0
            self.Q2 = np.eye(self.n2)
        # rhs E1: A M_k ; particular lifts N0_k and null space K
        Af = poly_float(A)
        self.AM = np.zeros((self.n1, K))
        for k, Bk in enumerate(self.basis):
            prod = core.mul(A, Bk)
            for e, v in prod.items():
                self.AM[ri1[code(e)], k] += float(v)
        U, sv, Vt = np.linalg.svd(M1, full_matrices=True)
        tol = 1e-11 * (sv[0] if sv.size else 1.0)
        rank1 = int(np.sum(sv > tol))
        self.rank1 = rank1
        self.Knull = Vt[rank1:].T  # nN x dimK
        self.dimK = self.Knull.shape[1]
        # min-norm particular solutions
        Sinv = np.zeros(nN); Sinv[:rank1] = 1.0 / sv[:rank1]
        self.N0 = Vt.T[:, :rank1] @ ((U.T[:rank1] @ self.AM) / sv[:rank1, None])
        self.e1_resid = float(np.linalg.norm(M1 @ self.N0 - self.AM) / max(np.linalg.norm(self.AM), 1e-300))
        if relax:
            # Hessian relaxation: N = A w for Jacobian elements (w given), z absent
            if w_of_basis is not None:
                colidx = {(i, m): j for j, (i, m) in enumerate(Ncols)}
                for k, w in enumerate(w_of_basis):
                    if w is None:
                        continue
                    (i, m) = w
                    vec = np.zeros(nN)
                    for e, v in core.mul({m: Fraction(1)}, A).items():
                        vec[colidx[(i, e)]] += float(v)
                    self.N0[:, k] = vec
            self.Knull = np.zeros((nN, 0)); self.dimK = 0
        # exclusion of F-multiples: subspace V = {a : D(a) in F S_{d-6}}
        self.P = np.eye(K); self.excluded = 0
        if exclude_F and d >= 6:
            Fm = inst.f_multiples(d, self.c)
            if Fm:
                mons = {}
                for Pp in self.basis + Fm:
                    for e in Pp:
                        mons.setdefault(e, len(mons))
                Mb = np.zeros((len(mons), K + len(Fm)))
                for k, Pp in enumerate(self.basis):
                    for e, v in Pp.items():
                        Mb[mons[e], k] = float(v)
                for j, Pp in enumerate(Fm):
                    for e, v in Pp.items():
                        Mb[mons[e], K + j] = -float(v)
                ns = sla.null_space(Mb, rcond=1e-11)
                S = ns[:K, :] if ns.size else np.zeros((K, 0))
                if S.shape[1]:
                    Us, ss, _ = np.linalg.svd(S, full_matrices=True)
                    r = int(np.sum(ss > 1e-9 * max(1.0, ss[0])))
                    self.P = Us[:, r:]; self.excluded = r
        self.dimP = self.P.shape[1]
        # family independence check (basis as polynomials)
        mons = {}
        for Pp in self.basis:
            for e in Pp:
                mons.setdefault(e, len(mons))
        Mb = np.zeros((len(mons), K))
        for k, Pp in enumerate(self.basis):
            for e, v in Pp.items():
                Mb[mons[e], k] = float(v)
        self.family_rank = int(np.linalg.matrix_rank(Mb, tol=1e-9 * max(1.0, np.abs(Mb).max())))
        self.build_time = time.time() - t0
        self.nfev = 0
        if verbose:
            print(f'[svp] d={d} c={self.c} K={K} (rank {self.family_rank}) relax={relax} n1={self.n1} n2={self.n2} '
                  f'nN={nN} nQ={nQ} rank1={rank1} dimK={self.dimK} e1_resid={self.e1_resid:.1e} '
                  f'excludedF={self.excluded} dimP={self.dimP} built {self.build_time:.1f}s', file=sys.stderr, flush=True)

    # ------------------------------------------------------------
    def E2_of(self, a):
        M = self.E2[0] * a[0]
        for k in range(1, self.K):
            if a[k] != 0.0:
                M = M + self.E2[k] * a[k]
        return M

    def value_full(self, a):
        """returns f, z, residual (full E2 rows, projected), b-norm, N vector at normalized a."""
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        E = self.E2_of(a)
        N0a = self.N0 @ a
        b = -(E @ N0a)
        bp = self.Q2.T @ b
        if self.dimK:
            CK = E @ self.Knull
            CKp = self.Q2.T @ CK
            z, *_ = sla.lstsq(CKp, bp, cond=1e-10, lapack_driver='gelsd', check_finite=False)
            rp = CKp @ z - bp
        else:
            z = np.zeros(0); rp = -bp
        f = float(rp @ rp)
        self.nfev += 1
        return f, z, rp, float(np.linalg.norm(bp)), N0a + (self.Knull @ z if self.dimK else 0.0), E

    def value_and_grad(self, u):
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        a = a_full / na
        f, z, rp, bn, Nfull, E = self.value_full(a)
        r = self.Q2 @ rp  # projected residual in full E2 coordinates
        # d r / d a_k = E2_k Nfull + E2_{D(a)} N0_k
        g = np.zeros(self.K)
        Et_r = E.T @ r
        for k in range(self.K):
            g[k] = 2.0 * (float(r @ (self.E2[k] @ Nfull)) + float(Et_r @ self.N0[:, k]))
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def rel(self, a):
        f, z, rp, bn, Nfull, E = self.value_full(a)
        return math.sqrt(f) / max(bn, 1e-300), f, bn

    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, init=None, polish=True,
               stop_rel=1e-13, time_limit=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            if init is not None and it == 0:
                u0 = self.P.T @ np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.dimP)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            u = res.x
            f = res.fun
            if polish and f < 1e-6:
                res2 = minimize(self.value_and_grad, u, jac=True, method='L-BFGS-B',
                                options={'maxiter': 4 * maxiter, 'gtol': 1e-16, 'ftol': 1e-20, 'maxcor': 30})
                u = res2.x; f = res2.fun
            a = self.P @ u; a /= np.linalg.norm(a)
            rel, f2, bn = self.rel(a)
            out.append({'rel': rel, 'f': f2, 'bnorm': bn, 'a': a, 'iters': int(res.nit), 'restart': it,
                        'seed': seed})
            if verbose:
                print(f'[svp] restart {it}: f={f2:.3e} rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)',
                      file=sys.stderr, flush=True)
            if rel < stop_rel:
                break
            if time_limit and time.time() - t0 > time_limit:
                break
        out.sort(key=lambda c: c['rel'])
        return out

    # ------------------------------------------------------------ exact side
    def D_exact(self, frac):
        D = {}
        for k, Bk in enumerate(self.basis):
            if frac[k]:
                D = core.add(D, core.scale(Bk, Fraction(frac[k])))
        return D

    def reconstruct(self, a, max_dens=(10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6), tol=1e-7):
        """Rational candidates (list of tuples of Fractions) for the direction a."""
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a)))
        a = a / a[j]
        cands = []
        for md in max_dens:
            fr = tuple(Fraction(float(v)).limit_denominator(md) for v in a)
            err = max(abs(float(f) - v) for f, v in zip(fr, a))
            if err < tol * max(1.0, md ** 0.5):
                if fr not in cands:
                    cands.append(fr)
        # common-denominator integer relation (small)
        for den in range(1, 2001):
            v = a * den
            if np.max(np.abs(v - np.round(v))) < 1e-6 * den:
                fr = tuple(Fraction(int(round(x)), den) for x in v)
                if fr not in cands:
                    cands.append(fr)
                break
        return cands

    def confirm(self, frac, verbose=False):
        """Exact sector-reduced solve for the rational candidate (sect.test), then core.verify_witness."""
        D = self.D_exact(frac)
        if not D:
            return 'inconsistent', {'note': 'D = 0'}
        from hodge import jacobian
        divF, _ = jacobian.divisible_by(self.inst.F, D)
        if divF:
            return 'spurious', {'note': 'F divides D', 'D': core.to_str(D)}
        out = sect.test(self.inst.F, self.inst.A, 1, D, exact=True, verbose=verbose)
        out['D'] = core.to_str(D)
        if out.get('verified'):
            return 'solved', out
        return ('consistent-unlifted' if out.get('consistent') else 'inconsistent'), out


def save_witness_candidate(D, N, Q, d, path=WD + '/witness-candidate.json', family=None):
    w = {'instanceId': LANE, 'd': d, 'D': D if isinstance(D, str) else core.to_str(D),
         'N': [n if isinstance(n, str) else core.to_str(n) for n in N],
         'Q': Q if isinstance(Q, str) else core.to_str(Q), 'family': family,
         'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    with open(path, 'w') as f:
        json.dump(w, f, indent=1)
    return path


def load_space(d, c):
    """Load the saved exact E1 space for (d, c) from wave2/e1space_d{d}_c{a}{b}.json -> (basis polys, w_of_basis, meta)."""
    fn = f'{W2}/e1space_d{d}_c{c[0]}{c[1]}.json'
    meta = json.load(open(fn))
    inst_F = None
    return meta


def basis_from_meta(inst, meta):
    Jb = [inst.jac_poly(i, tuple(m)) for (i, m) in meta['J_indep']]
    w = [(i, tuple(m)) for (i, m) in meta['J_indep']]
    reps = [core.parse(s) for s in meta['reps']]
    basis = Jb + reps
    w_of = w + [None] * len(reps)
    # normalize basis polys to unit float norm for conditioning (keeps exactness: rational scale)
    out = []
    for Pp in basis:
        nrm = poly_norm(Pp)
        s = Fraction(1)
        if nrm > 0:
            # rational scale close to 1/nrm with small denominator
            s = Fraction(1.0 / nrm).limit_denominator(1000)
            if s == 0:
                s = Fraction(1)
        out.append(core.scale(Pp, s))
    return out, w_of
