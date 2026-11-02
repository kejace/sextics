"""w2s_nvp.py -- normalized-E2 variable projection over a linear family of E1-admissible denominators.

Family: D(w) = sum_k w_k D_k with E1 lifts N0(w) = sum_k w_k N0_k (N0_k . grad F = A D_k).  Every E1 lift is
N = N0 + K grad F with K skew (entries of degree d-4, prescribed parity).  E2 (normalized):
    B_D(N0) + T_D(K) = F C,    B_D(N) = D div N - N(D),
    T_D(K) = sum_{i<j} ( D [(d_i K_ij) F_j - (d_j K_ij) F_i] - K_ij (D_i F_j - D_j F_i) ).
We eliminate C by the normal form modulo F (division by F with respect to X0^6: kernel exactly F*S_{2d-6}^-), so the
inner problem is the linear least squares  min_K | NF(T_D(w)(K)) - NF(-B_D(w)(N0(w))) |^2  and the outer functional
f(w) = min_K ... is minimized on the unit sphere of a gauge-free parameter space (variable projection; Golub-Pereyra
gradient).  f(w) = 0 iff D(w) admits a certificate (modulo floating point)."""
import sys, os, time, json, math
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2s_common import *
from scipy.optimize import minimize


class NormE2:
    def __init__(self, F, A, d, pD, Dmat, N0mat, P=None, cond=1e-10, verbose=True, label=''):
        t0 = time.time()
        self.F, self.A, self.d, self.pD = F, A, d, pD
        self.Fv = vec(F, 6); self.Av = vec(A, 6)
        self.Fi = [vec(pderiv(F, i), 5) for i in range(NV)]
        self.Dmat = np.asarray(Dmat, float); self.N0 = np.asarray(N0mat, float)
        self.K = self.Dmat.shape[0]
        self.P = np.eye(self.K) if P is None else np.asarray(P, float)
        self.dimP = self.P.shape[1]
        self.cond = cond; self.label = label
        n2d = nmons(2 * d)
        # rows: S_{2d}^- ; top = X0^6-divisible
        self.rows = sel(2 * d, 1); self.nrows = len(self.rows)
        rowmap = -np.ones(n2d, dtype=np.int64); rowmap[self.rows] = np.arange(self.nrows); self.rowmap = rowmap
        E = E_of(2 * d)[self.rows]
        self.top = np.nonzero(E[:, 0] >= 6)[0]; self.bot = np.nonzero(E[:, 0] < 6)[0]
        k2 = 2 * d - 6
        Etop = E[self.top].copy(); Etop[:, 0] -= 6
        gidx = lookup(k2, (Etop * (BASE ** np.arange(NV))).sum(1))
        Fsup = np.nonzero(self.Fv)[0]
        T6 = mult_table(6, k2)[np.ix_(Fsup, gidx)]
        R = rowmap[T6]
        assert np.all(R >= 0)
        toppos = -np.ones(self.nrows, dtype=np.int64); toppos[self.top] = np.arange(len(self.top))
        botpos = -np.ones(self.nrows, dtype=np.int64); botpos[self.bot] = np.arange(len(self.bot))
        tp = toppos[R]; bp = botpos[R]
        vals = np.broadcast_to(self.Fv[Fsup][:, None], R.shape)
        cols = np.broadcast_to(np.arange(len(gidx))[None, :], R.shape)
        mt = tp >= 0
        L = np.zeros((len(self.top), len(self.top)))
        L[tp[mt], cols[mt]] = vals[mt]
        self.LU = sla.lu_factor(L, check_finite=False)
        mb = bp >= 0
        self.Fbot = sp.csr_matrix((np.asarray(vals[mb]), (bp[mb], cols[mb])), shape=(len(self.bot), len(self.top)))
        self.FbotT = self.Fbot.T.tocsr()
        self.nbot = len(self.bot)
        # K unknowns
        self.kcols = []
        pd4 = par_of(d - 4)
        for a in range(NV):
            for b in range(a + 1, NV):
                par = (1 - pD - PVEC[a] - PVEC[b]) % 2
                for c in np.nonzero(pd4 == par)[0]:
                    self.kcols.append((a, b, c))
        self.nK = len(self.kcols)
        # P_c = (d_a m_c) F_b - (d_b m_c) F_a  (degree d)
        self.cP = sel(d, (1 - pD) % 2)
        Pfull = np.zeros((nmons(d), self.nK))
        Ed4 = E_of(d - 4)
        T5 = mult_table(5, d - 5)
        for col, (a, b, c) in enumerate(self.kcols):
            e = Ed4[c]
            if e[a]:
                e2 = e.copy(); e2[a] -= 1
                ci = lookup(d - 5, np.array([(e2 * (BASE ** np.arange(NV))).sum()]))[0]
                nz = np.nonzero(self.Fi[b])[0]
                Pfull[T5[nz, ci], col] += e[a] * self.Fi[b][nz]
            if e[b]:
                e2 = e.copy(); e2[b] -= 1
                ci = lookup(d - 5, np.array([(e2 * (BASE ** np.arange(NV))).sum()]))[0]
                nz = np.nonzero(self.Fi[a])[0]
                Pfull[T5[nz, ci], col] -= e[b] * self.Fi[a][nz]
        self.Pfull = Pfull
        self.Pr = Pfull[self.cP, :]
        # MultD index template: rows for e_D in S_d^pD times e_c in cP
        self.sD = sel(d, pD)
        self.RIDX = rowmap[mult_table(d, d)[np.ix_(self.sD, self.cP)]]
        assert np.all(self.RIDX >= 0)
        self.CIDX = np.broadcast_to(np.arange(len(self.cP))[None, :], self.RIDX.shape)
        # M2 templates per pair
        self.pairs = []
        colidx = {}
        for col, (a, b, c) in enumerate(self.kcols):
            colidx.setdefault((a, b), []).append((col, c))
        Tg = mult_table(d + 4, d - 4)
        for (a, b), lst in colidx.items():
            parG = (pD + PVEC[a] + PVEC[b]) % 2
            sG = sel(d + 4, parG)
            cc = np.array([c for _, c in lst]); gc = np.array([col for col, _ in lst])
            RI = rowmap[Tg[np.ix_(sG, cc)]]
            assert np.all(RI >= 0)
            CI = np.broadcast_to(gc[None, :], RI.shape)
            self.pairs.append((a, b, sG, RI, CI, gc, cc))
        # derivative matrices
        self.dD = [deriv_matrix(d, i) for i in range(NV)]           # n_{d-1} x n_d
        self.dN = [deriv_matrix(d + 1, i) for i in range(NV)]       # n_d x n_{d+1}
        # precomputed family derivatives
        self.dDmat = [self.Dmat @ self.dD[i].T for i in range(NV)]  # K x n_{d-1}
        self.divN0mat = sum(self.N0[:, i, :] @ self.dN[i].T for i in range(NV))   # K x n_d
        # F-multiple test
        if d >= 6:
            FS = np.stack([polymul(self.Fv, 6, np.eye(nmons(d - 6))[c], d - 6) for c in sel(d - 6, pD)], axis=1)
            self.FS_Q, _ = np.linalg.qr(FS)
        else:
            self.FS_Q = None
        self._nfev = 0
        self.build_time = time.time() - t0
        if verbose:
            print(f'[NormE2 {label}] d={d} pD={pD} K={self.K} dimP={self.dimP} rows={self.nrows} (bot {self.nbot}, top {len(self.top)}) '
                  f'nK={self.nK} built in {self.build_time:.1f}s', file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    def NF(self, V):
        g = sla.lu_solve(self.LU, V[self.top], check_finite=False)
        return V[self.bot] - self.Fbot @ g

    def NFT(self, y):
        z = np.zeros(self.nrows)
        z[self.bot] = y
        z[self.top] = -sla.lu_solve(self.LU, self.FbotT @ y, trans=1, check_finite=False)
        return z

    def assemble(self, w):
        d = self.d
        Dv = self.Dmat.T @ w
        M = np.zeros((self.nrows, len(self.cP)))
        M[self.RIDX, self.CIDX] = np.broadcast_to(Dv[self.sD][:, None], self.RIDX.shape)
        MK = M @ self.Pr
        dDv = [self.dD[i] @ Dv for i in range(NV)]
        for (a, b, sG, RI, CI, gc, cc) in self.pairs:
            G = polymul(dDv[a], d - 1, self.Fi[b], 5) - polymul(dDv[b], d - 1, self.Fi[a], 5)
            MK[RI, CI] -= np.broadcast_to(G[sG][:, None], RI.shape)
        N0w = np.einsum('k,kin->in', w, self.N0)
        divN0 = sum(self.dN[i] @ N0w[i] for i in range(NV))
        Bf = polymul(Dv, d, divN0, d) - sum(polymul(N0w[i], d + 1, dDv[i], d - 1) for i in range(NV))
        b = -Bf[self.rows]
        return MK, b, Dv, dDv, N0w, divN0

    def value_full(self, w):
        """f, x, rtilde, parts at the given w (no normalization)."""
        MK, b, Dv, dDv, N0w, divN0 = self.assemble(w)
        Mt = self.NF(MK); bt = self.NF(b)
        x, res, rk, sv = sla.lstsq(Mt, bt, cond=self.cond, lapack_driver='gelsy', check_finite=False)
        rt = Mt @ x - bt
        self._nfev += 1
        return float(rt @ rt), x, rt, (Dv, dDv, N0w, divN0, b, bt)

    def value(self, w):
        w = np.asarray(w, float); w = w / np.linalg.norm(w)
        f, x, rt, parts = self.value_full(w)
        bt = parts[5]
        rel = f / max(float(bt @ bt), 1e-300)
        return f, rel, x, rt

    def grad_w(self, w, x, rt, parts):
        d = self.d
        Dv, dDv, N0w, divN0, b, bt = parts
        R = np.zeros(nmons(2 * d)); R[self.rows] = self.NFT(rt)
        # K part: <R, D_k (P x)> - sum_pairs <R, K_ab(x) G_k^ab>
        u = self.Pfull @ x
        g1 = self.Dmat @ mult_adjoint(u, d, R, 2 * d)
        V = [np.zeros(nmons(d - 1)) for _ in range(NV)]
        for (a, b_, sG, RI, CI, gc, cc) in self.pairs:
            Kab = np.zeros(nmons(d - 4)); Kab[cc] = x[gc]
            z = mult_adjoint(Kab, d - 4, R, 2 * d)                    # degree d+4
            V[a] += mult_adjoint(self.Fi[b_], 5, z, d + 4)
            V[b_] -= mult_adjoint(self.Fi[a], 5, z, d + 4)
        g2 = sum(self.dDmat[i] @ V[i] for i in range(NV))
        gK = g1 - g2
        # b part: b = -[D divN0 - sum_i N0_i d_i D]
        t1 = self.Dmat @ mult_adjoint(divN0, d, R, 2 * d)
        t2 = self.divN0mat @ mult_adjoint(Dv, d, R, 2 * d)
        t3 = sum(self.N0[:, i, :] @ mult_adjoint(dDv[i], d - 1, R, 2 * d) for i in range(NV))
        t4 = sum(self.dDmat[i] @ mult_adjoint(N0w[i], d + 1, R, 2 * d) for i in range(NV))
        gb = -(t1 + t2 - t3 - t4)
        return 2.0 * (gK - gb)

    def value_and_grad(self, u):
        u = np.asarray(u, float)
        a_full = self.P @ u
        na = np.linalg.norm(a_full)
        w = a_full / na
        f, x, rt, parts = self.value_full(w)
        g = self.grad_w(w, x, rt, parts)
        g_full = (g - float(g @ w) * w) / na
        return f, self.P.T @ g_full

    def is_spurious(self, w, tol=1e-8):
        if self.FS_Q is None:
            return False, 1.0
        Dv = self.Dmat.T @ w
        resid = Dv - self.FS_Q @ (self.FS_Q.T @ Dv)
        rel = float(np.linalg.norm(resid) / max(np.linalg.norm(Dv), 1e-300))
        return rel < tol, rel

    def search_one(self, seed, maxiter=300, gtol=1e-12, init=None, verbose=False):
        rng = np.random.default_rng(seed)
        if init is not None:
            u0 = self.P.T @ np.asarray(init, float)
        else:
            u0 = rng.standard_normal(self.dimP)
        u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
        w = self.P @ res.x; w /= np.linalg.norm(w)
        f, rel, x, rt = self.value(w)
        spur, srel = self.is_spurious(w)
        return {'seed': int(seed), 'value': f, 'rel': rel, 'iters': int(res.nit), 'nfev': int(res.nfev),
                'spurious': bool(spur), 'div_rel': srel, 'seconds': round(time.time() - t0, 2),
                'msg': str(res.message)[:60], 'w': w}

    def gradient_check(self, seed=0, eps=1e-6, n=5):
        rng = np.random.default_rng(seed)
        u = rng.standard_normal(self.dimP); u /= np.linalg.norm(u)
        f0, g = self.value_and_grad(u)
        out = []
        for _ in range(n):
            v = rng.standard_normal(self.dimP); v /= np.linalg.norm(v)
            fp, _ = self.value_and_grad(u + eps * v); fm, _ = self.value_and_grad(u - eps * v)
            fd = (fp - fm) / (2 * eps)
            out.append((float(g @ v), fd))
        return f0, out


def family(F, A, d, pD, kind, verbose=True):
    """Build (Dmat, N0, P, label) for kind in {'jac', 'kfull', 'ker'}."""
    Fi = [vec(pderiv(F, i), 5) for i in range(NV)]
    if kind == 'jac':
        Dmat, N0, labels, gauge = jac_family(F, A, Fi, d, pD)
        P = complement_projector(Dmat.shape[0], gauge)
        return Dmat, N0, P, f'jac d={d} pD={pD}'
    cp = cache_path(f'family_{kind}_d{d}_par{pD}')
    if os.path.exists(cp):
        z = np.load(cp)
        return z['Dmat'], z['N0'], None, str(z['label'])
    Kb, rankA, cols = K_full(F, A, Fi, d, pD, verbose=verbose)          # n_d x dimK (orthonormal columns)
    Fv = vec(F, 6)
    FS = np.stack([polymul(Fv, 6, np.eye(nmons(d - 6))[c], d - 6) for c in sel(d - 6, pD)], axis=1)  # n_d x g
    # coordinates of F-multiples inside K_full (they lie in J subset K_full)
    C = Kb.T @ FS                                                            # dimK x g
    resid = np.linalg.norm(Kb @ C - FS) / np.linalg.norm(FS)
    if verbose:
        print(f'[family {kind} d={d} pD={pD}] F-multiples inside K_full: residual {resid:.2e}', file=sys.stderr, flush=True)
    U, s, Vt = np.linalg.svd(C, full_matrices=True)
    r = int(np.sum(s > 1e-9 * s[0]))
    Pc = U[:, r:]                                                            # dimK x (dimK - g)
    if kind == 'ker':
        # remove the J_d part as well: orthogonal complement of J inside K_full
        Jb, rJ = J_basis(Fi, d, pD)
        CJ = Kb.T @ Jb
        UJ, sJ, _ = np.linalg.svd(CJ, full_matrices=True)
        rJ2 = int(np.sum(sJ > 1e-9 * sJ[0]))
        Pc = UJ[:, rJ2:]
        if verbose:
            print(f'[family ker d={d} pD={pD}] dim J_d={rJ} (num {rJ2}), kernel-only dim {Pc.shape[1]}', file=sys.stderr, flush=True)
    Dmat = (Kb @ Pc).T                                                       # K' x n_d
    N0, e1res = e1_lifts(F, A, Fi, d, pD, Dmat, verbose=verbose)
    label = f'{kind} d={d} pD={pD} dim={Dmat.shape[0]} (E1 lift max resid {e1res.max():.1e})'
    np.savez(cp, Dmat=Dmat, N0=N0, label=label)
    return Dmat, N0, None, label


# ----------------------------------------------------------------------------------------------
# multiprocessing driver
# ----------------------------------------------------------------------------------------------
_G = {}


def _init(obj):
    _G['obj'] = obj


def _run(args):
    seed, maxiter = args
    return _G['obj'].search_one(seed, maxiter=maxiter)


def main(argv=None):
    import argparse, multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--pD', type=int, default=0)
    ap.add_argument('--kind', default='jac', choices=['jac', 'kfull', 'ker'])
    ap.add_argument('--restarts', type=int, default=100)
    ap.add_argument('--seed0', type=int, default=0)
    ap.add_argument('--maxiter', type=int, default=300)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--max-minutes', type=float, default=1e9)
    ap.add_argument('--name', required=True)
    ap.add_argument('--gradcheck', action='store_true')
    ap.add_argument('--rel-save', type=float, default=1e-6)
    ap.add_argument('--confirm-rel', type=float, default=1e-10)
    a = ap.parse_args(argv)
    F, A, q = load()
    t0 = time.time()
    Dmat, N0, P, lab = family(F, A, a.d, a.pD, a.kind)
    obj = NormE2(F, A, a.d, a.pD, Dmat, N0, P=P, label=lab)
    out_path = os.path.join(LANE, f'{a.name}.restarts.jsonl')
    if a.gradcheck:
        f0, chk = obj.gradient_check()
        print('gradcheck f0', f0, chk, file=sys.stderr, flush=True)
    # timing probe
    tt = time.time(); f0, rel0, _, _ = obj.value(np.random.default_rng(123).standard_normal(obj.K) @ np.eye(obj.K)); t_eval = time.time() - tt
    log_event(event='nvp-build', name=a.name, d=a.d, pD=a.pD, kind=a.kind, K=obj.K, dimP=obj.dimP, rows=obj.nrows,
              rows_modF=obj.nbot, nK=obj.nK, build_seconds=round(obj.build_time, 1), eval_seconds=round(t_eval, 3),
              rel_at_random=rel0, label=lab)
    print(f'[nvp {a.name}] one evaluation {t_eval:.3f}s; rel at random point {rel0:.3e}', file=sys.stderr, flush=True)
    seeds = list(range(a.seed0, a.seed0 + a.restarts))
    best = None; nrun = 0; rels = []
    tstart = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(a.procs, initializer=_init, initargs=(obj,)) as pool:
        for rec in pool.imap_unordered(_run, [(s, a.maxiter) for s in seeds], chunksize=1):
            nrun += 1
            w = rec.pop('w')
            rels.append(rec['rel'])
            rec.update({'name': a.name, 'd': a.d, 'pD': a.pD, 'kind': a.kind})
            if rec['rel'] < a.rel_save and not rec['spurious']:
                np.save(os.path.join(LANE, f'{a.name}.seed{rec["seed"]}.w.npy'), w)
                np.save(os.path.join(LANE, f'{a.name}.seed{rec["seed"]}.D.npy'), Dmat.T @ w)
                rec['saved'] = True
            with open(out_path, 'a') as f:
                f.write(json.dumps(rec) + '\n')
            if best is None or (not rec['spurious'] and rec['rel'] < best['rel']):
                best = dict(rec); best['w'] = w
            if nrun % 20 == 0 or nrun == len(seeds):
                print(f'[nvp {a.name}] {nrun} restarts, best rel {best["rel"]:.3e} (seed {best["seed"]}), '
                      f'median {np.median(rels):.3e} ({time.time()-tstart:.0f}s)', file=sys.stderr, flush=True)
            if (time.time() - tstart) / 60 > a.max_minutes:
                pool.terminate(); break
    rels = np.array(rels)
    summary = dict(event='nvp-batch', name=a.name, d=a.d, pD=a.pD, kind=a.kind, family=lab, restarts=int(nrun),
                   maxiter=a.maxiter, outer_dim=obj.dimP, inner_dim=obj.nK, rows_modF=obj.nbot,
                   best_rel=float(best['rel']) if best else None, best_seed=best['seed'] if best else None,
                   median_rel=float(np.median(rels)) if len(rels) else None,
                   quantiles=[float(np.quantile(rels, t)) for t in (0.0, 0.1, 0.5, 0.9, 1.0)] if len(rels) else None,
                   n_below_1e6=int(np.sum(rels < 1e-6)), n_below_1e10=int(np.sum(rels < 1e-10)),
                   seconds=round(time.time() - tstart, 1))
    log_event(**summary)
    print(json.dumps(summary), flush=True)
    # confirmation of the best candidates
    if best is not None and best['rel'] < a.confirm_rel and not best['spurious']:
        Dv = Dmat.T @ best['w']
        for cand in rational_candidates(Dv, a.d):
            status, data = exact_confirm(F, A, q, cand['poly'])
            log_event(event='nvp-confirm', name=a.name, d=a.d, den=cand['den'], err=cand['err'], status=status,
                      note=str({k: v for k, v in data.items() if k not in ('N', 'Q')})[:300])
            if status == 'solved':
                rec = save_witness_candidate(F, A, q, cand['poly'], data['N'], data['Q'])
                log_event(event='WITNESS', name=a.name, d=a.d, verified=rec['verified'], F_divides_D=rec['F_divides_D'])
                print('WITNESS FOUND', rec['verified'], file=sys.stderr, flush=True)
                break
    return 0


if __name__ == '__main__':
    sys.exit(main())
