"""Reduced variable-projection search over an E1-admissible linear family (q = 1).

D(a) = sum_k a_k M_k (pure Y-parity parD, degree d).  For each family member a particular E1 lift N0_k
(N0_k(F) = A M_k - F G0_k) is fixed: for Jacobian generators M = m F_j: N0 = A m e_j, G0 = 0; for kernel reps:
float least squares.  The complete E2 residual (book ch. 5.2 'full correction map') is
    r(a, c, Q) = D G + D div N - N(D) - F Q,   N = N0(a) + K(c) grad F,  K antisymmetric, deg c_ij = d - 4,
linear in (c, Q), quadratic in a.  f(a) = min_{c,Q} |r|^2 (Tikhonov-regularized) with analytic gradient
(Golub-Pereyra), minimized by L-BFGS over the unit sphere of a restricted to the complement of F-multiples.
"""
import sys, os, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
import numpy as np
import scipy.sparse as sps
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction

PAR = [0, 0, 0, 1, 1, 1]

class Space:
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = core.monomials(deg, par) if deg >= 0 else []
        self.idx = {m: i for i, m in enumerate(self.mons)}
        self.n = len(self.mons)
        self.E = np.array(self.mons, dtype=np.int64).reshape(self.n, 6)
        self.keys = self.E @ (31 ** np.arange(6, dtype=np.int64))
        order = np.argsort(self.keys)
        self.skeys, self.sorder = self.keys[order], order
    def lookup(self, keys):
        pos = np.searchsorted(self.skeys, keys)
        pos = np.clip(pos, 0, self.n - 1)
        ok = self.skeys[pos] == keys
        return self.sorder[pos], ok
    def vec(self, P):
        v = np.zeros(self.n)
        for e, c in P.items():
            v[self.idx[e]] = float(c)
        return v
    def poly(self, v, tol=0.0):
        return {m: float(v[i]) for i, m in enumerate(self.mons) if abs(v[i]) > tol}

_space_cache = {}
def space(deg, par):
    k = (deg, par)
    if k not in _space_cache: _space_cache[k] = Space(deg, par)
    return _space_cache[k]

class PairMap:
    """All products Sa x Sb -> Sc (parities must match)."""
    def __init__(self, Sa, Sb, Sc):
        assert (Sa.par + Sb.par) % 2 == Sc.par and Sa.deg + Sb.deg == Sc.deg
        ia = np.repeat(np.arange(Sa.n), Sb.n); ib = np.tile(np.arange(Sb.n), Sa.n)
        keys = Sa.keys[ia] + Sb.keys[ib]
        ic, ok = Sc.lookup(keys)
        assert ok.all()
        self.ia, self.ib, self.ic = ia, ib, ic
        self.Sa, self.Sb, self.Sc = Sa, Sb, Sc
    def prod(self, x, y):
        return np.bincount(self.ic, weights=x[self.ia] * y[self.ib], minlength=self.Sc.n)
    def corr_a(self, r, y):
        return np.bincount(self.ia, weights=r[self.ic] * y[self.ib], minlength=self.Sa.n)
    def corr_b(self, r, x):
        return np.bincount(self.ib, weights=r[self.ic] * x[self.ia], minlength=self.Sb.n)
    def shift_matrix(self, y):
        """dense matrix (Sc.n x Sa.n) with column m = m * y  (y in Sb)."""
        return sps.coo_matrix((y[self.ib], (self.ic, self.ia)), shape=(self.Sc.n, self.Sa.n)).toarray()

def mult_matrix(P, Sa, Sc):
    """sparse (Sc.n x Sa.n) matrix of multiplication by the poly P (float dict)."""
    rows, cols, vals = [], [], []
    PE = np.array(list(P.keys()), dtype=np.int64).reshape(-1, 6)
    Pv = np.array([float(v) for v in P.values()])
    Pk = PE @ (31 ** np.arange(6, dtype=np.int64))
    for t in range(len(Pk)):
        ic, ok = Sc.lookup(Sa.keys + Pk[t])
        assert ok.all()
        rows.append(ic); cols.append(np.arange(Sa.n)); vals.append(np.full(Sa.n, Pv[t]))
    if not rows: return sps.csr_matrix((Sc.n, Sa.n))
    return sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Sc.n, Sa.n))

def deriv_matrix(i, Sa, Sc):
    e = Sa.E[:, i]
    nz = np.nonzero(e)[0]
    keys = Sa.keys[nz] - 31 ** i
    ic, ok = Sc.lookup(keys)
    assert ok.all()
    return sps.csr_matrix((e[nz].astype(float), (ic, nz)), shape=(Sc.n, Sa.n))

def fpoly(P):
    return {e: float(v) for e, v in P.items()}

class ReducedVarPro:
    def __init__(self, F, A, members, lam=1e-12, verbose=True):
        """members: list of dicts {'M': exact poly, 'tag': str, 'kind': 'J'|'kernel', 'N0': [6 float polys] or None, 'G0': float poly or None}"""
        t0 = time.time()
        self.F, self.A = F, A
        self.members = members
        K = len(members); self.K = K
        d = core.degree(members[0]['M']); self.d = d
        parD = core.parity_of(members[0]['M']); self.parD = parD
        parA = core.parity_of(A); self.parA = parA
        self.parN = [(parA + parD + PAR[i]) % 2 for i in range(6)]
        self.lam = lam
        Ff = fpoly(F); Af = fpoly(A)
        Fi = [fpoly(core.pderiv(F, i)) for i in range(6)]
        # spaces
        self.SD = space(d, parD)
        self.SDi = [space(d - 1, (parD + PAR[i]) % 2) for i in range(6)]
        self.SN = [space(d + 1, self.parN[i]) for i in range(6)]
        self.SG = space(d, (parA + parD) % 2)          # G, div N
        self.SR = space(2 * d, parA)                    # residual
        self.SQ = space(2 * d - 6, parA)
        self.Sc = {}                                     # c-spaces by parity
        for pi in (0, 1): self.Sc[pi] = space(d - 4, pi)
        self.SP = {}                                     # P_i spaces: S_{2d-5}^{pi}
        for pi in (0, 1): self.SP[pi] = space(2 * d - 5, pi)
        # product maps
        self.pm1 = PairMap(self.SD, self.SG, self.SR)
        self.pm2 = [PairMap(self.SN[i], self.SDi[i], self.SR) for i in range(6)]
        # c-block structure: pairs (i<j), parity pi_ij of c_ij
        self.cblocks = []
        for i in range(6):
            for j in range(i + 1, 6):
                pij = (parA + parD + PAR[i] + PAR[j]) % 2
                self.cblocks.append((i, j, pij))
        self.nK = sum(self.Sc[pij].n for (_, _, pij) in self.cblocks)
        self.nQ = self.SQ.n
        # shift maps for P_i^{(pi)}: S_{d-4}^{pi} x S_d^{parD} -> S_{2d-4}^{pi+parD} then d_i -> S_{2d-5}^{pi+parD+par(i)};
        # and S_{d-4}^{pi} x S_{d-1}^{parD+par(i)} -> S_{2d-5}^{pi+parD+par(i)}
        self.Pmaps = {}
        for i in range(6):
            for pi in (0, 1):
                if not any((ii == i or jj == i) and pij == pi for (ii, jj, pij) in self.cblocks): continue
                Sm = self.Sc[pi]
                S_mD = space(2 * d - 4, (pi + parD) % 2)
                S_out = self.SP[(pi + parD + PAR[i]) % 2]
                self.Pmaps[(i, pi)] = (PairMap(Sm, self.SD, S_mD), deriv_matrix(i, S_mD, S_out),
                                       PairMap(Sm, self.SDi[i], S_out))
        # multiplication by F_j: S_{2d-5}^{pi} -> S_{2d}^{parA} (pi = parA + par(j)) and S_{d-4}^{pi} -> S_{d+1}^{parN_i}
        self.MF_big = {}
        self.MF_small = {}
        for j in range(6):
            pi = (parA + PAR[j]) % 2
            self.MF_big[j] = mult_matrix(Fi[j], self.SP[pi], self.SR)
        for (i, j, pij) in self.cblocks:
            # V_i += c_ij F_j ; V_j -= c_ij F_i
            self.MF_small[(i, j)] = mult_matrix(Fi[j], self.Sc[pij], self.SN[i])
            self.MF_small[(j, i)] = mult_matrix(Fi[i], self.Sc[pij], self.SN[j])
        self.Dmaps_N = [deriv_matrix(i, self.SN[i], self.SG) for i in range(6)]   # div of N
        # Q columns: F * S_{2d-6}^{parA} -> S_{2d}^{parA}
        self.FQ = mult_matrix(Ff, self.SQ, self.SR).tocsc()
        G = (self.FQ.T @ self.FQ).toarray() + lam * np.eye(self.nQ)
        self.cho_Q = sla.cho_factor(G)
        # family matrices
        self.UD = np.zeros((self.SD.n, K))
        self.UDi = [np.zeros((self.SDi[i].n, K)) for i in range(6)]
        self.UN = [np.zeros((self.SN[i].n, K)) for i in range(6)]
        self.UdivN = np.zeros((self.SG.n, K))
        self.UG = np.zeros((self.SG.n, K))
        for k, mem in enumerate(members):
            M = mem['M']
            self.UD[:, k] = self.SD.vec(M)
            for i in range(6): self.UDi[i][:, k] = self.SDi[i].vec(core.pderiv(M, i))
            if mem['kind'] == 'J':
                j, m = mem['j'], mem['m']
                N0 = [{} for _ in range(6)]; N0[j] = core.shift(A, m, 1)
                G0 = {}
            else:
                N0, G0 = mem['N0'], mem['G0']
            for i in range(6):
                self.UN[i][:, k] = self.SN[i].vec(N0[i])
                self.UdivN[:, k] += self.Dmaps_N[i] @ self.UN[i][:, k]
            self.UG[:, k] = self.SG.vec(G0)
        self.UGdiv = self.UG + self.UdivN
        # complement of F-multiples in a-space
        SFm = space(d - 6, parD)
        if SFm.n > 0:
            MFm = mult_matrix(Ff, SFm, self.SD).toarray()
            Mbig = np.hstack([self.UD, -MFm])
            # nullspace
            u, s, vt = np.linalg.svd(Mbig, full_matrices=True)
            tol = s[0] * 1e-10
            rank = int((s > tol).sum())
            null = vt[rank:].T          # (K + nFm) x nullity
            aF = null[:K]                # a-components of F-multiples inside family
            # orthonormal basis of a-directions to exclude
            if aF.shape[1] > 0:
                qF, _ = np.linalg.qr(aF)
                self.nF = qF.shape[1]
                # complement
                u2, s2, vt2 = np.linalg.svd(qF.T, full_matrices=True)
                self.P = vt2[self.nF:].T    # K x (K - nF)
            else:
                self.nF = 0; self.P = np.eye(K)
            self.MFm = MFm
        else:
            self.nF = 0; self.P = np.eye(K); self.MFm = None
        self.dimP = self.P.shape[1]
        self.rows = self.SR.n
        self.nevals = 0
        if verbose:
            print(f'[rvp] d={d} parD={parD} K={K} (F-multiple directions excluded: {self.nF}, dimP={self.dimP}) '
                  f'rows={self.rows} nK={self.nK} nQ={self.nQ} setup {time.time()-t0:.1f}s', file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    def build_CK(self, D, Di):
        blocks = []
        Pcache = {}
        for (i, j, pij) in self.cblocks:
            for ii in (i, j):
                key = (ii, pij)
                if key not in Pcache:
                    pmD, Der, pmDi = self.Pmaps[key]
                    Pcache[key] = Der @ pmD.shift_matrix(D) - 2.0 * pmDi.shift_matrix(Di[ii])
            blocks.append(self.MF_big[j] @ Pcache[(i, pij)] - self.MF_big[i] @ Pcache[(j, pij)])
        return np.hstack(blocks)

    def V_from_c(self, c):
        V = [np.zeros(self.SN[i].n) for i in range(6)]
        off = 0
        for (i, j, pij) in self.cblocks:
            n = self.Sc[pij].n
            cij = c[off:off + n]; off += n
            V[i] += self.MF_small[(i, j)] @ cij
            V[j] -= self.MF_small[(j, i)] @ cij
        return V

    def evaluate(self, a, want_grad=True):
        self.nevals += 1
        D = self.UD @ a
        Di = [self.UDi[i] @ a for i in range(6)]
        N0 = [self.UN[i] @ a for i in range(6)]
        Gd = self.UGdiv @ a                      # G + div N0
        b = self.pm1.prod(D, Gd)
        for i in range(6): b -= self.pm2[i].prod(N0[i], Di[i])
        CK = self.build_CK(D, Di)
        # project out Q-range: W = P_perp CK, bp = P_perp b
        def Pperp(X):
            Y = self.FQ.T @ X
            Y = sla.cho_solve(self.cho_Q, Y)
            return X - self.FQ @ Y
        W = Pperp(CK); bp = Pperp(b)
        Gram = W.T @ W + self.lam * np.eye(self.nK)
        c = -np.linalg.solve(Gram, W.T @ bp)
        qv = -sla.cho_solve(self.cho_Q, self.FQ.T @ (b + CK @ c))
        r = b + CK @ c + self.FQ @ qv
        f = float(r @ r)
        bn2 = float(b @ b)
        V = self.V_from_c(c)
        divV = np.zeros(self.SG.n)
        for i in range(6): divV += self.Dmaps_N[i] @ V[i]
        Gd_tot = Gd + divV
        N = [N0[i] + V[i] for i in range(6)]
        scale = float(D @ D) * (float(Gd_tot @ Gd_tot) + sum(float(N[i] @ N[i]) for i in range(6)))
        out = {'f': f, 'b2': bn2, 'rel': f / max(scale, 1e-300), 'relb': f / max(bn2, 1e-300), 'c': c, 'q': qv, 'r': r, 'D': D}
        if want_grad:
            u1 = self.pm1.corr_a(r, Gd_tot)          # in SD
            u2 = self.pm1.corr_b(r, D)               # in SG
            g = self.UD.T @ u1 + self.UGdiv.T @ u2
            for i in range(6):
                u3 = self.pm2[i].corr_a(r, Di[i])    # in SN_i
                u4 = self.pm2[i].corr_b(r, N[i])     # in SDi
                g -= self.UN[i].T @ u3 + self.UDi[i].T @ u4
            out['grad'] = 2.0 * g
        return out

    def value_and_grad(self, u):
        nu = np.linalg.norm(u)
        a = self.P @ (u / nu)
        ev = self.evaluate(a)
        ga = ev['grad']
        # chain rule: a = P u/|u|
        gu = self.P.T @ ga
        gu = (gu - (u / nu) * float(gu @ (u / nu))) / nu
        return ev['f'], gu

    def is_spurious(self, a, tol=1e-8):
        D = self.UD @ a
        if self.MFm is None or self.MFm.shape[1] == 0: return False, 1.0
        g, res, rk, s = np.linalg.lstsq(self.MFm, D, rcond=None)
        rel = np.linalg.norm(D - self.MFm @ g) / max(np.linalg.norm(D), 1e-300)
        return rel < tol, float(rel)

    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, log_path=None, tag=''):
        rng = np.random.default_rng(seed)
        out = []
        for it in range(n_restarts):
            t0 = time.time()
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 30})
            a = self.P @ (res.x / np.linalg.norm(res.x))
            ev = self.evaluate(a, want_grad=False)
            spur, srel = self.is_spurious(a)
            cand = {'restart': it, 'seed': seed, 'rel': ev['rel'], 'relb': ev['relb'], 'f': ev['f'], 'iters': int(res.nit),
                    'spurious': spur, 'div_rel': srel, 'seconds': round(time.time() - t0, 1), 'a': a.tolist(), 'tag': tag}
            out.append(cand)
            if verbose:
                print(f'[rvp] restart {it}: rel={ev["rel"]:.3e} f={ev["f"]:.3e} relb={ev["relb"]:.2e} iters={res.nit} spur={spur} ({cand["seconds"]}s)', file=sys.stderr, flush=True)
            if log_path:
                with open(log_path, 'a') as fh:
                    fh.write(json.dumps({k: v for k, v in cand.items()}) + '\n')
        out.sort(key=lambda c: (c['spurious'], c['rel']))
        return out

    # ------------------------------------------------------------------
    def D_exact(self, frac):
        D = {}
        for k, mem in enumerate(self.members):
            if frac[k] == 0: continue
            D = core.add(D, core.scale(mem['M'], Fraction(frac[k])))
        return D

    def reconstruct(self, a, max_den=1000, tol=1e-7):
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a))); a = a / a[j]
        cands = []
        for md in (10, 100, max_den):
            fr = [Fraction(float(v)).limit_denominator(md) for v in a]
            if max(abs(float(f) - v) for f, v in zip(fr, a)) < tol:
                cands.append(tuple(fr)); break
        for den in range(1, max_den + 1):
            v = a * den
            if np.max(np.abs(v - np.round(v))) < tol * den:
                t = tuple(Fraction(int(round(x)), den) for x in v)
                if t not in cands: cands.append(t)
                break
        return cands

    def confirm(self, frac, verbose=False):
        D = self.D_exact(frac)
        if not D: return 'inconsistent', {'note': 'D=0'}
        ok, info = fixed_d.e1_solvable(self.F, self.A, 1, D)
        if not ok: return 'inconsistent', {'note': 'e1 fails', 'e1': info, 'D': D}
        out = fixed_d.solve_fixed_D_fast(self.F, self.A, 1, D, verbose=verbose)
        if not out['consistent']:
            return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity'], 'D': D}
        if out.get('exact'):
            return 'solved', {'N': out['N'], 'Q': out['Q'], 'D': D}
        return 'lift-failed', {'D': D}

# ----------------------------------------------------------------------
def e1_particular(F, A, D0, d, parD):
    """float least-squares E1 lift N0 and quotient G0 for a kernel rep D0 (exact or float poly)."""
    parA = core.parity_of(A)
    parN = [(parA + parD + PAR[i]) % 2 for i in range(6)]
    St = space(d + 6, (parA + parD) % 2)
    cols, meta = [], []
    Fi = [fpoly(core.pderiv(F, i)) for i in range(6)]
    blocks = []
    for i in range(6):
        Sn = space(d + 1, parN[i])
        blocks.append(mult_matrix(Fi[i], Sn, St))
    Mat = sps.hstack(blocks).toarray()
    rhs = St.vec(core.mul(A, D0))
    x, res, rk, s = np.linalg.lstsq(Mat, rhs, rcond=1e-13)
    resid = np.linalg.norm(Mat @ x - rhs) / np.linalg.norm(rhs)
    N0 = []
    off = 0
    for i in range(6):
        Sn = space(d + 1, parN[i])
        N0.append(Sn.poly(x[off:off + Sn.n])); off += Sn.n
    # G0 = (A D0 - N0(F)) / F  -- exact A D0 - N0(F) is F*G0 only if N0 exact; use lstsq
    SG = space(d, (parA + parD) % 2)
    MF = mult_matrix(fpoly(F), SG, St).toarray()
    rhs2 = rhs - Mat @ x
    g, _, _, _ = np.linalg.lstsq(MF, rhs2, rcond=1e-13)
    resid2 = np.linalg.norm(MF @ g - rhs2) / max(np.linalg.norm(rhs2), 1e-300)
    return N0, SG.poly(g), float(resid), float(resid2), float(np.linalg.norm(rhs2) / np.linalg.norm(rhs))

def J_members(F, A, d, parD):
    mems = []
    for j in range(6):
        for m in core.monomials(d - 5, (parD + PAR[j]) % 2):
            mems.append({'M': core.shift(core.pderiv(F, j), m, 1), 'kind': 'J', 'j': j, 'm': m, 'tag': f'm*F_{j}'})
    return mems

def kernel_members(F, A, d, parD, reps, verbose=True):
    mems = []
    for k, D0 in enumerate(reps):
        N0, G0, r1, r2, gn = e1_particular(F, A, D0, d, parD)
        if verbose: print(f'[rvp] kernel rep {k}: E1 lstsq resid {r1:.2e}, G resid {r2:.2e}, |FG|/|AD| {gn:.2e}', file=sys.stderr, flush=True)
        mems.append({'M': D0, 'kind': 'kernel', 'N0': N0, 'G0': G0, 'tag': f'kernel{k}'})
    return mems

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--instance', default='short-power-presentation')
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--par', type=int, default=0)
    ap.add_argument('--family', default='J+kernel', help='J | kernel | J+kernel | file:<json list of D strings>')
    ap.add_argument('--kernel-file', default=None)
    ap.add_argument('--restarts', type=int, default=10)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--maxiter', type=int, default=400)
    ap.add_argument('--out', default=None)
    ap.add_argument('--tag', default='')
    ap.add_argument('--confirm-rel', type=float, default=1e-9)
    ap.add_argument('--lam', type=float, default=1e-12)
    args = ap.parse_args()
    F, A, q, rec = core.load_instance(args.instance)
    assert q == 1
    d, parD = args.d, args.par
    mems = []
    if 'J' in args.family.split('+'):
        mems += J_members(F, A, d, parD)
    if 'kernel' in args.family.split('+'):
        path = args.kernel_file or f'{RES}/w2/kernel_d{d}_p{parD}_rational.json'
        reps = [core.parse(it['D']) for it in json.load(open(path))]
        mems += kernel_members(F, A, d, parD, reps)
    if args.family.startswith('file:'):
        items = json.load(open(args.family[5:]))
        reps = [core.parse(it['D'] if isinstance(it, dict) else it) for it in items]
        mems += kernel_members(F, A, d, parD, reps)
    fam = ReducedVarPro(F, A, mems, lam=args.lam)
    out = args.out or f'{RES}/w2/vp_d{d}_p{parD}_{args.tag or "run"}_s{args.seed}.jsonl'
    t0 = time.time()
    if args.procs > 1:
        import multiprocessing as mp
        per = args.restarts // args.procs
        def worker(w):
            fam.search(n_restarts=per, seed=args.seed * 1000 + w, maxiter=args.maxiter, log_path=out, tag=args.tag)
        ps = [mp.Process(target=worker, args=(w,)) for w in range(args.procs)]
        for p_ in ps: p_.start()
        for p_ in ps: p_.join()
        cands = [json.loads(l) for l in open(out)]
        cands.sort(key=lambda c: (c['spurious'], c['rel']))
    else:
        cands = fam.search(n_restarts=args.restarts, seed=args.seed, maxiter=args.maxiter, verbose=True, log_path=out, tag=args.tag)
    rels = sorted(c['rel'] for c in cands if not c['spurious'])
    hist = {}
    for r in rels:
        b = int(np.floor(np.log10(max(r, 1e-300))))
        hist[b] = hist.get(b, 0) + 1
    summary = {'batch': 'varpro', 'family': args.family, 'degree': d, 'parity': parD, 'K': fam.K, 'dimP': fam.dimP,
               'nK': fam.nK, 'nQ': fam.nQ, 'rows': fam.rows, 'restarts': len(cands), 'seed': args.seed,
               'best_rel': rels[0] if rels else None, 'median_rel': rels[len(rels)//2] if rels else None,
               'hist_log10': {str(k): v for k, v in sorted(hist.items())}, 'n_spurious': sum(c['spurious'] for c in cands),
               'seconds': round(time.time() - t0, 1), 'out': out, 'tag': args.tag}
    print(json.dumps(summary), flush=True)
    # confirmation attempts
    confirmed = []
    for c in cands[:20]:
        if c['spurious'] or c['rel'] > args.confirm_rel: continue
        for fr in fam.reconstruct(c['a']):
            status, data = fam.confirm(fr)
            print(f'[rvp] confirm restart {c["restart"]} rel {c["rel"]:.2e}: {status}', flush=True)
            wlog(batch='varpro-confirm', degree=d, parity=parD, rel=c['rel'], status=status, a=[str(x) for x in fr])
            if status == 'solved':
                W = {'d': d, 'D': core.to_str(data['D']), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
                json.dump(W, open(f'{RES}/witness-candidate.json', 'w'))
                print('!!! WITNESS FOUND (exact) — saved witness-candidate.json', flush=True)
                confirmed.append(W)
                break
        if confirmed: break
    summary['confirmed'] = len(confirmed)
    wlog(**summary)

if __name__ == '__main__':
    main()
