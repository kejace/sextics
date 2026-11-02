"""Reduced variable-projection search over an E1-admissible family (kernel reps + Jacobian generators), q = 1.

For D(a) = sum_k a_k G_k with particular numerators P_k (P_k(F) = A G_k), every E1 solution is
N = sum_k a_k P_k + K_c (K_c Koszul: (K_c)_i = sum_j c_ij F_j, c antisymmetric), and the certificate condition is
    B(a) + L_c(a) c + F*Q = 0,   B(a) = D div(P(a)) - P(a)(D) = sum_{k<=l} a_k a_l v_kl   (degree 2d).
For fixed a the unknowns (c, Q) enter linearly; f(a) = min_{c,Q} |B(a) + L_c(a) c + F Q|^2 on the unit sphere of the
gauge-fixed a-space (F-multiple directions removed).  Gradient by the Golub-Pereyra formula.
Usage: python3 w2_varpro.py <d> <par> --restarts N --procs P --seed S [--maxiter 400] [--tag T] [--no-kernel]
       [--confirm-top k] [--confirm-tol 1e-9]
Rational kernel reps + numerators are read from kernel_d{d}_p{par}.json (numerators_rational required)."""
import sys, os, json, time, argparse, math
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
import numpy as np
import scipy.sparse as sps
import scipy.linalg as sla
from scipy.optimize import minimize

def mon_index(deg, par):
    ms = monomials(deg, par)
    return ms, {m: i for i, m in enumerate(ms)}

def poly_vec(P, idx, n):
    v = np.zeros(n)
    for e, c in P.items():
        v[idx[e]] += float(c)
    return v

def mult_matrix(P, deg_in, par_in, deg_out, par_out):
    """Sparse matrix of multiplication by P: S_{deg_in}^{par_in} -> S_{deg_out}^{par_out}."""
    mi, ii = mon_index(deg_in, par_in); mo, io = mon_index(deg_out, par_out)
    ri, ci, vals = [], [], []
    for e, c in P.items():
        fc = float(c)
        for k, m in enumerate(mi):
            ri.append(io[mon_add(m, e)]); ci.append(k); vals.append(fc)
    return sps.csr_matrix((vals, (ri, ci)), shape=(len(mo), len(mi)))

def deriv_matrix(i, deg_in, par_in):
    mi, ii = mon_index(deg_in, par_in)
    par_out = (par_in - (1 if i >= 3 else 0)) % 2
    mo, io = mon_index(deg_in - 1, par_out)
    ri, ci, vals = [], [], []
    for k, m in enumerate(mi):
        if m[i]:
            e = list(m); e[i] -= 1
            ri.append(io[tuple(e)]); ci.append(k); vals.append(float(m[i]))
    return sps.csr_matrix((vals, (ri, ci)), shape=(len(mo), len(mi))), par_out

def shift_index(deg_in, par_in, m, mo_idx):
    mi = monomials(deg_in, par_in)
    return np.array([mo_idx[mon_add(e, m)] for e in mi], dtype=np.int64)

class ReducedVarPro:
    def __init__(self, F, A, d, parD, kernel_reps=(), kernel_nums=(), verbose=True, include_jacobian=True):
        t0 = time.time()
        self.F, self.A, self.d, self.parD = F, A, d, parD
        pA = core.parity_of(A)
        fam = Family(F, A, d, kernel_reps=kernel_reps, kernel_nums=kernel_nums, parities=(parD,), include_jacobian=include_jacobian)
        self.fam = fam; K = fam.K; self.K = K
        # normalize generators
        self.scale = []
        for mem in fam.members:
            nrm = math.sqrt(sum(float(c) ** 2 for c in mem['G'].values()))
            self.scale.append(nrm)
        self.rows_par = pA  # parity of B = A D^2
        e = 2 * d
        self.R_mons, self.R_idx = mon_index(e, pA); R = len(self.R_mons); self.R = R
        # ---- V tensor: v_kl for k<=l (symmetrized T_kl) exactly, then float
        terms = relaxation_terms(fam, None, verbose=False) if False else None
        # exact rational T_kl via the same structure as relaxation_terms but with p=None
        self.pairs = [(k, l) for k in range(K) for l in range(k, K)]
        Vt = np.zeros((len(self.pairs), R))
        mem = fam.members
        Fi = [pderiv(F, i) for i in range(NV)]
        Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
        Ai = [pderiv(A, i) for i in range(NV)]
        FA = [[mul(Fi[j], Ai[jp]) for jp in range(NV)] for j in range(NV)]
        AF = [mul(A, Fi[j]) for j in range(NV)]
        AH = [[mul(A, Fij[j][jp]) for jp in range(NV)] for j in range(NV)]
        def jac_T(k, l):
            j, m = mem[k]['j'], mem[k]['m']; jp, mp = mem[l]['j'], mem[l]['m']
            mm = mon_add(m, mp)
            out = shift(FA[j][jp], mm, 1)
            if mp[jp]:
                ee = list(mp); ee[jp] -= 1
                out = add(out, shift(AF[j], mon_add(m, tuple(ee)), mp[jp]))
            if m[jp]:
                ee = list(m); ee[jp] -= 1
                out = sub(out, shift(AF[j], mon_add(mp, tuple(ee)), m[jp]))
            out = sub(out, shift(AH[j][jp], mm, 1))
            return out
        for pi, (k, l) in enumerate(self.pairs):
            if mem[k]['tag'] == 'J' and mem[l]['tag'] == 'J':
                t = jac_T(k, l)
                if l != k: t = add(t, jac_T(l, k))
            else:
                t = T_kl(mem[k]['G'], mem[l]['P'])
                if l != k: t = add(t, T_kl(mem[l]['G'], mem[k]['P']))
            s = 1.0 / (self.scale[k] * self.scale[l])
            for ee, c in t.items():
                Vt[pi, self.R_idx[ee]] += float(c) * s
        self.Vt = Vt
        # pair index helpers for gradient: for each k, list of (pair index, l, factor)
        self.pair_pos = {pr: i for i, pr in enumerate(self.pairs)}
        t1 = time.time()
        # ---- generator vectors in S_d
        self.D_mons, self.D_idx = mon_index(d, parD); nD = len(self.D_mons)
        Gmat = np.zeros((K, nD))
        for k, mm in enumerate(mem):
            for ee, c in mm['G'].items():
                Gmat[k, self.D_idx[ee]] += float(c) / self.scale[k]
        self.Gmat = Gmat
        # ---- Koszul column structure
        # c_ij parity = 1 + parD + par_i + par_j ; mu of degree d-4
        self.cpairs = []   # (i, j, mu)
        self.cmons = {}
        for i in range(NV):
            for j in range(i + 1, NV):
                pc = (pA + parD + (i >= 3) + (j >= 3)) % 2
                for mu in monomials(d - 4, pc):
                    self.cpairs.append((i, j, mu))
        self.nc = len(self.cpairs)
        # multiplication maps: F_j: S_d -> S_{d+5};  F_j: S_{d-1}^{par'} -> S_{d+4}; derivative maps on S_d
        self.MF_d = [mult_matrix(Fi[j], d, parD, d + 5, (parD + (j >= 3)) % 2) for j in range(NV)]
        self.dD = []; self.dD_par = []
        for i in range(NV):
            Mi, po = deriv_matrix(i, d, parD); self.dD.append(Mi); self.dD_par.append(po)
        self.MF_dm1 = [[None] * NV for _ in range(NV)]  # [j][i]: F_j * (S_{d-1}^{par of d_i D}) -> S_{d+4}
        for i in range(NV):
            for j in range(NV):
                po = self.dD_par[i]
                self.MF_dm1[j][i] = mult_matrix(Fi[j], d - 1, po, d + 4, (po + (j >= 3)) % 2)
        # index maps for scatters: S_{d+5}^{par} shifted by nu (deg d-5) -> rows ; S_{d+4}^{par} shifted by mu (deg d-4) -> rows
        self.shift_d5 = {}; self.shift_d4 = {}
        for i, j, mu in self.cpairs:
            for t in (i, j):
                if mu[t]:
                    nu = list(mu); nu[t] -= 1; nu = tuple(nu)
                    pn = (self.rows_par - 0) % 2
                    key = (nu, (parD + (j >= 3)) % 2) if t == i else (nu, (parD + (i >= 3)) % 2)
                    if key not in self.shift_d5:
                        self.shift_d5[key] = shift_index(d + 5, key[1], nu, self.R_idx)
            po = (self.dD_par[i] + (j >= 3)) % 2  # parity of F_j d_i D
            key = (mu, po)
            if key not in self.shift_d4:
                self.shift_d4[key] = shift_index(d + 4, po, mu, self.R_idx)
        # ---- Q block: columns -F * nu, nu of degree 2d-6 parity pA
        self.Q_mons = monomials(2 * d - 6, pA); nQ = len(self.Q_mons); self.nQ = nQ
        LQ = np.zeros((R, nQ))
        for qi, nu in enumerate(self.Q_mons):
            for ee, c in F.items():
                LQ[self.R_idx[mon_add(ee, nu)], qi] += float(c)
        self.LQ = LQ
        # ---- gradient helpers: Koszul field maps.  K_c,i = sum_j c_ij F_j (S_{d-4} -> S_{d+1});
        self.c_blocks = {}   # (i,j) -> (col slice indices, mult matrix of F_j on S_{d-4}^{pc} -> S_{d+1}^{par N_i}, parity)
        self.parN = [(pA + parD + (i >= 3)) % 2 for i in range(NV)]
        self.N_mons = [mon_index(d + 1, self.parN[i]) for i in range(NV)]
        col = 0
        self.cslices = []
        for i in range(NV):
            for j in range(i + 1, NV):
                pc = (pA + parD + (i >= 3) + (j >= 3)) % 2
                n_mu = len(monomials(d - 4, pc))
                self.cslices.append((i, j, pc, col, col + n_mu))
                col += n_mu
        self.MF_c = {}
        for (i, j, pc, c0, c1) in self.cslices:
            self.MF_c[(i, j)] = (mult_matrix(Fi[j], d - 4, pc, d + 1, self.parN[i]),   # contributes to K_i
                                 mult_matrix(Fi[i], d - 4, pc, d + 1, self.parN[j]))   # contributes to K_j (minus)
        self.dN = [deriv_matrix(i, d + 1, self.parN[i])[0] for i in range(NV)]     # d_i K_i : S_{d+1} -> S_d^{parD}
        # products needed in the gradient: F_j * divK (S_d^{parD} -> S_{d+5}) [= MF_d], F_ji * K_i (S_{d+1}^{parN_i} -> S_{d+5}^{..}),
        # F_j * K_i (S_{d+1}^{parN_i} -> S_{d+6})
        self.pK = (pA + parD) % 2
        self.MF_dK = [mult_matrix(Fi[j], d, self.pK, d + 5, (self.pK + (j >= 3)) % 2) for j in range(NV)]
        self.MFij_N = [[mult_matrix(Fij[j][i], d + 1, self.parN[i], d + 5, (self.pK + (j >= 3)) % 2) for i in range(NV)] for j in range(NV)]
        self.MFj_N = [[mult_matrix(Fi[j], d + 1, self.parN[i], d + 6, (self.parN[i] + (j >= 3)) % 2) for i in range(NV)] for j in range(NV)]
        # index maps for Jacobian members: shift of S_{d+5}^{par(F_j)+parD} by m (deg d-5) and of S_{d+6}^{..} by d_i m (deg d-6)
        self.jac_idx = []
        for k, mm in enumerate(mem):
            if mm['tag'] != 'J':
                self.jac_idx.append(None); continue
            j, m = mm['j'], mm['m']
            pj = (self.pK + (j >= 3)) % 2
            idx_m = shift_index(d + 5, pj, m, self.R_idx)
            dlist = []
            for i in range(NV):
                if m[i]:
                    nu = list(m); nu[i] -= 1; nu = tuple(nu)
                    pz = (self.parN[i] + (j >= 3)) % 2
                    dlist.append((i, float(m[i]), shift_index(d + 6, pz, nu, self.R_idx)))
            self.jac_idx.append((j, idx_m, dlist))
        # kernel members: dense multiplication maps for G_k and d_i G_k
        self.ker_maps = {}
        for k, mm in enumerate(mem):
            if mm['tag'] == 'K':
                MG = mult_matrix(mm['G'], d, (pA + parD) % 2, 2 * d, pA)
                MdG = [mult_matrix(pderiv(mm['G'], i), d + 1, self.parN[i], 2 * d, pA) for i in range(NV)]
                self.ker_maps[k] = (MG, MdG)
        # ---- gauge: complement of the spurious subspace
        sv = fam.spurious_vectors()
        if sv:
            S = np.array([[float(x) * self.scale[k] for k, x in enumerate(v)] for v in sv]).T  # K x s (in scaled coords)
            U, s, _ = np.linalg.svd(S, full_matrices=True)
            r = int(np.sum(s > 1e-10 * max(1.0, s[0])))
            self.P = U[:, r:]
        else:
            self.P = np.eye(K)
        self.dimP = self.P.shape[1]
        self.nfev = 0
        rng0 = np.random.default_rng(12345)
        bb = []
        for _ in range(8):
            v = rng0.standard_normal(K); v /= np.linalg.norm(v); bv = self.bvec(v); bb.append(float(bv @ bv))
        self.scale0 = float(np.mean(bb))
        self.build_time = time.time() - t0
        if verbose:
            print(f'[rvp] d={d} par={parD} K={K} (kernel {fam.nK}) rows={R} nc={self.nc} nQ={nQ} dimP={self.dimP} scale0={self.scale0:.3e} '
                  f'(V {t1-t0:.0f}s, total {self.build_time:.0f}s)', file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    def Dvec(self, a):
        return a @ self.Gmat

    def koszul_matrix(self, a):
        """R x nc matrix L_c(a)."""
        Dv = self.Dvec(a)
        DF = [self.MF_d[j] @ Dv for j in range(NV)]
        dDv = [self.dD[i] @ Dv for i in range(NV)]
        L = np.zeros((self.R, self.nc))
        col = 0
        for (i, j, pc, c0, c1) in self.cslices:
            Psi = self.MF_dm1[j][i] @ dDv[i] - self.MF_dm1[i][j] @ dDv[j]   # F_j d_i D - F_i d_j D  (deg d+4)
            po = (self.dD_par[i] + (j >= 3)) % 2
            for mu in monomials(self.d - 4, pc):
                colv = np.zeros(self.R)
                if mu[i]:
                    nu = list(mu); nu[i] -= 1; nu = tuple(nu)
                    idx = self.shift_d5[(nu, (self.parD + (j >= 3)) % 2)]
                    np.add.at(colv, idx, mu[i] * DF[j])
                if mu[j]:
                    nu = list(mu); nu[j] -= 1; nu = tuple(nu)
                    idx = self.shift_d5[(nu, (self.parD + (i >= 3)) % 2)]
                    np.add.at(colv, idx, -mu[j] * DF[i])
                idx = self.shift_d4[(mu, po)]
                np.add.at(colv, idx, -Psi)
                L[:, col] = colv; col += 1
        return L

    def bvec(self, a):
        pa = np.array([a[k] * a[l] for (k, l) in self.pairs])
        return pa @ self.Vt

    def dbvec(self, a):
        """K x R: d b / d a_k."""
        out = np.zeros((self.K, self.R))
        for pi, (k, l) in enumerate(self.pairs):
            if k == l:
                out[k] += 2 * a[k] * self.Vt[pi]
            else:
                out[k] += a[l] * self.Vt[pi]; out[l] += a[k] * self.Vt[pi]
        return out

    def solve(self, a):
        L = np.hstack([self.koszul_matrix(a), self.LQ])
        b = self.bvec(a)
        x, *_ = sla.lstsq(L, -b, cond=1e-10, lapack_driver='gelsy', check_finite=False)
        r = L @ x + b
        self.nfev += 1
        return float(r @ r), x, r, b

    def grad(self, a, x, r):
        c = x[:self.nc]
        # Koszul field K_c
        Kc = [np.zeros(len(self.N_mons[i][0])) for i in range(NV)]
        for (i, j, pc, c0, c1) in self.cslices:
            Mi, Mj = self.MF_c[(i, j)]
            Kc[i] += Mi @ c[c0:c1]; Kc[j] -= Mj @ c[c0:c1]
        divK = sum(self.dN[i] @ Kc[i] for i in range(NV))
        Y = [self.MF_dK[j] @ divK - sum(self.MFij_N[j][i] @ Kc[i] for i in range(NV)) for j in range(NV)]
        Z = [[self.MFj_N[j][i] @ Kc[i] for i in range(NV)] for j in range(NV)]
        db = self.dbvec(a)
        g = np.zeros(self.K)
        for k in range(self.K):
            if self.jac_idx[k] is not None:
                j, idx_m, dlist = self.jac_idx[k]
                val = float(Y[j] @ r[idx_m])
                for (i, coef, idx) in dlist:
                    val -= coef * float(Z[j][i] @ r[idx])
                val /= self.scale[k]
            else:
                MG, MdG = self.ker_maps[k]
                vk = MG @ divK - sum(MdG[i] @ Kc[i] for i in range(NV))
                val = float(vk @ r) / self.scale[k]
            g[k] = 2.0 * (val + float(db[k] @ r))
        return g

    def value_and_grad(self, u):
        af = self.P @ u; na = np.linalg.norm(af); a = af / na
        f, x, r, b = self.solve(a)
        g = self.grad(a, x, r)
        g_full = (g - float(g @ a) * a) / na
        return f, self.P.T @ g_full

    def rel(self, a):
        f, x, r, b = self.solve(a)
        return f / max(self.scale0, 1e-300), f, float(b @ b)

    def search(self, n_restarts, seed, maxiter=400, verbose=False, stop_rel=1e-14):
        rng = np.random.default_rng(seed)
        out = []
        for it in range(n_restarts):
            t0 = time.time()
            u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': 1e-12, 'ftol': 1e-18})
            a = self.P @ res.x; a /= np.linalg.norm(a)
            rel, f, bb = self.rel(a)
            cand = {'rel': rel, 'f': f, 'bnorm2': bb, 'a': a.tolist(), 'iters': int(res.nit), 'restart': it,
                    'seed': seed, 'seconds': round(time.time() - t0, 1)}
            out.append(cand)
            if verbose:
                print(f'[rvp] seed={seed} restart {it}: rel={rel:.3e} f={f:.3e} iters={res.nit} ({time.time()-t0:.0f}s)',
                      file=sys.stderr, flush=True)
            if rel < stop_rel:
                break
        return out

    # ------------------------------------------------------------------
    def D_exact(self, frac):
        D = {}
        for k, mem in enumerate(self.fam.members):
            if frac[k]:
                D = add(D, scale(mem['G'], Fraction(frac[k]) / Fraction(self.scale[k]).limit_denominator(10**12)))
        return D

    def reconstruct(self, a, max_dens=(10, 100, 1000, 10000, 100000), tol=1e-7):
        a = np.asarray(a, float)
        j = int(np.argmax(np.abs(a))); a = a / a[j]
        cands = []
        for md in max_dens:
            fr = [Fraction(float(v)).limit_denominator(md) for v in a]
            err = max(abs(float(f) - v) for f, v in zip(fr, a))
            cands.append((err, md, tuple(fr)))
        return cands

def worker(args):
    (d, par, kfile, seed, n_restarts, maxiter, no_kernel) = args
    F, A, q, rec = load()
    reps, nums = [], []
    if not no_kernel and kfile and os.path.exists(kfile):
        J = json.load(open(kfile))
        if J.get('basis_rational') and J.get('numerators_rational'):
            reps = [core.parse(s) for s in J['basis_rational']]
            nums = [[core.parse(c) for c in P] for P in J['numerators_rational']]
    rvp = ReducedVarPro(F, A, d, par, reps, nums, verbose=(seed % 1000 == 0))
    res = rvp.search(n_restarts, seed, maxiter=maxiter, verbose=True)
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('par', type=int)
    ap.add_argument('--restarts', type=int, default=8); ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--maxiter', type=int, default=400)
    ap.add_argument('--tag', default=''); ap.add_argument('--no-kernel', action='store_true')
    ap.add_argument('--confirm-top', type=int, default=2); ap.add_argument('--confirm-tol', type=float, default=1e-9)
    a = ap.parse_args()
    t0 = time.time()
    kfile = os.path.join(RES, f'kernel_d{a.d}_p{a.par}.json')
    per = max(1, a.restarts // a.procs)
    jobs = [(a.d, a.par, kfile, a.seed + 1000 * i, per, a.maxiter, a.no_kernel) for i in range(a.procs)]
    if a.procs > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(a.procs) as pool:
            results = pool.map(worker, jobs)
    else:
        results = [worker(jobs[0])]
    cands = [c for r in results for c in r]
    cands.sort(key=lambda c: c['rel'])
    rels = np.array([c['rel'] for c in cands])
    hist = {}
    for r in rels:
        b = int(math.floor(math.log10(max(r, 1e-300))))
        hist[b] = hist.get(b, 0) + 1
    out_path = os.path.join(RES, f'rvp_d{a.d}_p{a.par}_s{a.seed}' + (f'_{a.tag}' if a.tag else '') + '.json')
    with open(out_path, 'w') as f:
        json.dump({'d': a.d, 'par': a.par, 'seed': a.seed, 'restarts': len(cands), 'maxiter': a.maxiter,
                   'candidates': cands[:50]}, f)
    logrec = {'kind': 'varpro-reduced', 'd': a.d, 'parity': a.par, 'family': 'K_d + J(F)_d (gauge-fixed)',
              'restarts': len(cands), 'procs': a.procs, 'seed': a.seed, 'maxiter': a.maxiter,
              'best_rel': float(rels[0]) if len(rels) else None, 'median_rel': float(np.median(rels)) if len(rels) else None,
              'hist_log10': {str(k): v for k, v in sorted(hist.items())}, 'out': out_path, 'tag': a.tag,
              'seconds': round(time.time() - t0, 1)}
    # exact confirmation of the best few
    F, A, q, rec = load()
    J = json.load(open(kfile)) if os.path.exists(kfile) and not a.no_kernel else {}
    reps = [core.parse(s) for s in J['basis_rational']] if J.get('basis_rational') and J.get('numerators_rational') else []
    nums = [[core.parse(c) for c in P] for P in J['numerators_rational']] if reps else []
    rvp = ReducedVarPro(F, A, a.d, a.par, reps, nums, verbose=False)
    confirmed = []
    for c in cands[:a.confirm_top]:
        if c['rel'] > a.confirm_tol:
            break
        for err, md, fr in rvp.reconstruct(c['a']):
            D = rvp.D_exact(fr)
            if not D: continue
            ok, info = fixed_d.e1_solvable(F, A, 1, D)
            res = {'rel': c['rel'], 'max_den': md, 'recon_err': err, 'e1': bool(ok)}
            if ok:
                out = fixed_d.solve_fixed_D_fast(F, A, 1, D, verbose=False, want_exact=True)
                res.update({'consistent': bool(out['consistent']), 'rank': out['rank'], 'nullity': out['nullity']})
                if out.get('exact'):
                    okv, _ = core.verify_witness(F, A, 1, D, out['N'], out['Q'])
                    res['WITNESS'] = bool(okv)
                    if okv:
                        with open(os.path.join(RES, 'witness-candidate.json'), 'w') as f:
                            json.dump({'instance': IID, 'd': a.d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']],
                                       'Q': core.to_str(out['Q']), 'family': 'varpro-reduced'}, f, indent=1)
                        print('!!! WITNESS FOUND', file=sys.stderr, flush=True)
            confirmed.append(res)
            if res.get('consistent'):
                break
    logrec['confirm'] = confirmed
    w2log(logrec)
    print(json.dumps(logrec, default=str), flush=True)

if __name__ == '__main__':
    main()
