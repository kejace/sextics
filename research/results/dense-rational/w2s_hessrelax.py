"""w2s_hessrelax.py -- Hessian relaxation for Jacobian-type denominators D = w(F) of degree d:
necessary condition  w^T Hess(F) w  in  (J : A)_{2d-6}   (i.e. A w^T H w in J_{2d}).  Equations = projection onto the
orthogonal complement of K_full(2d-6) (numerical, from w2s_common.K_full).  Gauge w -> w + l*Euler removed.
Sectors: 0 / 1 (tau-pure, parity(w_i) = pD + p_i) or 'mixed' (both parities, both equation parities).
Modes: real (L-BFGS on the unit sphere, many restarts) and complex (Levenberg-Marquardt with a random affine slice).
Relaxation solutions (pure sectors) are promoted to the full normalized-E2 functional (w2s_nvp.NormE2)."""
import sys, os, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2s_common import *
from scipy.optimize import minimize, least_squares


def build(F, A, d, sector, verbose=True):
    Fi = [vec(pderiv(F, i), 5) for i in range(NV)]
    H = [[vec(pderiv(pderiv(F, i), j), 4) for j in range(NV)] for i in range(NV)]
    e = d - 5
    labels = []
    for i in range(NV):
        pars = [0, 1] if sector == 'mixed' else [(int(sector) + PVEC[i]) % 2]
        for par in pars:
            for c in sel(e, par):
                labels.append((i, c))
    nw = len(labels)
    # gauge
    gauge = []
    Ed6 = E_of(d - 6)
    lpars = [0, 1] if sector == 'mixed' else [int(sector)]
    idx = {lab: k for k, lab in enumerate(labels)}
    for lp in lpars:
        for c in sel(d - 6, lp):
            v = np.zeros(nw)
            for i in range(NV):
                ee = Ed6[c].copy(); ee[i] += 1
                cc = lookup(e, np.array([(ee * (BASE ** np.arange(NV))).sum()]))[0]
                v[idx[(i, cc)]] = 1.0
            gauge.append(v)
    P = complement_projector(nw, np.array(gauge))
    k2 = 2 * d - 6
    # equation projector(s)
    Cs = []
    eq_pars = [0, 1] if sector == 'mixed' else [0]
    for par in eq_pars:
        Kb, rankA, cols = K_full(F, A, Fi, k2, par, verbose=verbose)
        Kr = Kb[cols, :]                                   # n_par x dimK, orthonormal columns
        if Kr.shape[1] == 0:
            r = 0; C = np.eye(Kr.shape[0])
        else:
            U, s, Vt = np.linalg.svd(Kr, full_matrices=True)
            r = int(np.sum(s > 1e-9 * s[0]))
            C = U[:, r:]                                   # n_par x n_eq
        Cs.append((cols, C))
        if verbose:
            print(f'[hessrelax d={d} sector={sector}] eq parity {par}: dim K_full({k2})={Kb.shape[1]} (num rank {r}) n_eq={C.shape[1]}', file=sys.stderr, flush=True)
    neq = sum(C.shape[1] for _, C in Cs)
    T = np.zeros((neq, nw, nw))
    Tm = mult_table(e, e)
    for k in range(nw):
        ik, ck = labels[k]
        for l in range(k, nw):
            il, cl = labels[l]
            h = H[ik][il]
            nz = np.nonzero(h)[0]
            if len(nz) == 0:
                continue
            mi = Tm[ck, cl]
            v = np.zeros(nmons(k2))
            Tmm = mult_table(4, 2 * e)
            v[Tmm[nz, mi]] = h[nz]
            off = 0
            for cols, C in Cs:
                q = C.T @ v[cols]
                T[off:off + C.shape[1], k, l] = q
                if l != k:
                    T[off:off + C.shape[1], l, k] = q
                off += C.shape[1]
    return dict(T=T, P=P, labels=labels, nw=nw, neq=neq, gauge=np.array(gauge), dimP=P.shape[1])


class Relax:
    def __init__(self, data):
        self.T, self.P = data['T'], data['P']
        self.nw, self.neq, self.dimP = data['nw'], data['neq'], data['dimP']
        self.T2 = self.T.reshape(self.neq, -1)

    def q(self, w):
        return self.T2 @ np.outer(w, w).ravel()

    def fg(self, u):
        a = self.P @ u; na = np.linalg.norm(a); w = a / na
        Tw = np.einsum('mkl,l->mk', self.T, w)
        q = Tw @ w
        f = float(q @ q)
        g = 4.0 * (q @ Tw)
        g_full = (g - float(g @ w) * w) / na
        return f, self.P.T @ g_full

    def real_restart(self, seed, maxiter=500):
        rng = np.random.default_rng(seed)
        u0 = rng.standard_normal(self.dimP); u0 /= np.linalg.norm(u0)
        t0 = time.time()
        res = minimize(self.fg, u0, jac=True, method='L-BFGS-B', options={'maxiter': maxiter, 'gtol': 1e-14, 'ftol': 1e-20, 'maxcor': 30})
        w = self.P @ res.x; w /= np.linalg.norm(w)
        f = float(self.q(w) @ self.q(w))
        return {'seed': int(seed), 'mode': 'real', 'f': f, 'iters': int(res.nit), 'seconds': round(time.time() - t0, 2), 'w': w}

    def complex_restart(self, seed, max_nfev=400):
        rng = np.random.default_rng(seed)
        L = rng.standard_normal(self.dimP) + 1j * rng.standard_normal(self.dimP)
        n = self.dimP
        P = self.P
        def resid(x):
            u = x[:n] + 1j * x[n:]
            w = P @ u
            r = self.T2 @ np.outer(w, w).ravel()
            lin = L @ u - 1.0
            return np.concatenate([r.real, r.imag, [lin.real, lin.imag]])
        def jac(x):
            u = x[:n] + 1j * x[n:]
            w = P @ u
            Jc = 2.0 * np.einsum('mkl,l->mk', self.T, w) @ P
            top = np.hstack([Jc.real, -Jc.imag]); bot = np.hstack([Jc.imag, Jc.real])
            lt = np.hstack([L.real, -L.imag]); lb = np.hstack([L.imag, L.real])
            return np.vstack([top, bot, lt[None, :], lb[None, :]])
        x0 = rng.standard_normal(2 * n)
        u0 = x0[:n] + 1j * x0[n:]; u0 /= (L @ u0); x0 = np.concatenate([u0.real, u0.imag])
        t0 = time.time()
        res = least_squares(resid, x0, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=max_nfev)
        u = res.x[:n] + 1j * res.x[n:]
        w = P @ u
        r = self.T2 @ np.outer(w, w).ravel()
        f = float(np.linalg.norm(r) ** 2 / max(np.linalg.norm(w) ** 4, 1e-300))
        return {'seed': int(seed), 'mode': 'complex', 'f': f, 'iters': int(res.nfev), 'seconds': round(time.time() - t0, 2), 'w': w}


_G = {}


def _init(obj):
    _G['obj'] = obj


def _run(args):
    seed, mode, maxiter = args
    return _G['obj'].real_restart(seed, maxiter) if mode == 'real' else _G['obj'].complex_restart(seed, maxiter)


def main(argv=None):
    import argparse, multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--sector', default='0')
    ap.add_argument('--mode', default='real', choices=['real', 'complex'])
    ap.add_argument('--restarts', type=int, default=1000)
    ap.add_argument('--seed0', type=int, default=0)
    ap.add_argument('--maxiter', type=int, default=500)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--name', required=True)
    ap.add_argument('--max-minutes', type=float, default=1e9)
    ap.add_argument('--promote', type=float, default=1e-9, help='relaxation f below which pure-sector w is promoted to full E2')
    a = ap.parse_args(argv)
    F, A, q = load()
    t0 = time.time()
    data = build(F, A, a.d, a.sector)
    obj = Relax(data)
    # scale reference: f at random points
    rng = np.random.default_rng(1)
    ref = []
    for _ in range(20):
        u = rng.standard_normal(obj.dimP); w = obj.P @ u; w /= np.linalg.norm(w); ref.append(float(obj.q(w) @ obj.q(w)))
    log_event(event='hessrelax-build', name=a.name, d=a.d, sector=a.sector, nw=obj.nw, dimP=obj.dimP, neq=obj.neq,
              f_random_median=float(np.median(ref)), build_seconds=round(time.time() - t0, 1))
    print(f'[hessrelax {a.name}] nw={obj.nw} dimP={obj.dimP} neq={obj.neq} f(random) median {np.median(ref):.3e} build {time.time()-t0:.0f}s', file=sys.stderr, flush=True)
    seeds = list(range(a.seed0, a.seed0 + a.restarts))
    out_path = os.path.join(LANE, f'{a.name}.restarts.jsonl')
    fs = []; best = None; nrun = 0; tstart = time.time()
    ctx = mp.get_context('fork')
    promoted = []
    with ctx.Pool(a.procs, initializer=_init, initargs=(obj,)) as pool:
        for rec in pool.imap_unordered(_run, [(s, a.mode, a.maxiter) for s in seeds], chunksize=4):
            nrun += 1
            w = rec.pop('w'); fs.append(rec['f'])
            rec.update({'name': a.name, 'd': a.d, 'sector': a.sector})
            if rec['f'] < a.promote:
                np.save(os.path.join(LANE, f'{a.name}.seed{rec["seed"]}.w.npy'), w)
                promoted.append((rec['seed'], w)); rec['saved'] = True
            with open(out_path, 'a') as f:
                f.write(json.dumps(rec) + '\n')
            if best is None or rec['f'] < best['f']:
                best = dict(rec)
            if nrun % 100 == 0 or nrun == len(seeds):
                print(f'[hessrelax {a.name}] {nrun} restarts best f {best["f"]:.3e} (seed {best["seed"]}) median {np.median(fs):.3e} ({time.time()-tstart:.0f}s)', file=sys.stderr, flush=True)
            if (time.time() - tstart) / 60 > a.max_minutes:
                pool.terminate(); break
    fs = np.array(fs)
    summary = dict(event='hessrelax-batch', name=a.name, d=a.d, sector=a.sector, mode=a.mode, restarts=int(nrun), nw=obj.nw,
                   dimP=obj.dimP, neq=obj.neq, best_f=float(best['f']) if best else None, best_seed=best['seed'] if best else None,
                   quantiles=[float(np.quantile(fs, t)) for t in (0.0, 0.1, 0.5, 0.9, 1.0)] if len(fs) else None,
                   f_random_median=float(np.median(ref)), n_below_1e6=int(np.sum(fs < 1e-6)), n_below_1e9=int(np.sum(fs < 1e-9)),
                   n_promoted=len(promoted), seconds=round(time.time() - tstart, 1))
    log_event(**summary)
    print(json.dumps(summary), flush=True)
    # promotion to the full E2 functional (pure sectors, real w only)
    if promoted and a.sector in ('0', '1') and a.mode == 'real':
        from w2s_nvp import family, NormE2
        Dmat, N0, P, lab = family(F, A, a.d, int(a.sector), 'jac')
        E2 = NormE2(F, A, a.d, int(a.sector), Dmat, N0, P=P, label=lab)
        for seed, w in promoted[:50]:
            f, rel, x, rt = E2.value(np.asarray(w, float))
            log_event(event='hessrelax-promote', name=a.name, d=a.d, sector=a.sector, seed=int(seed), e2_rel=rel, e2_f=f)
            print(f'[promote] seed {seed}: full E2 rel {rel:.3e}', file=sys.stderr, flush=True)
            if rel < 1e-10:
                Dv = Dmat.T @ np.asarray(w, float)
                for cand in rational_candidates(Dv, a.d):
                    status, data = exact_confirm(F, A, q, cand['poly'])
                    log_event(event='hessrelax-confirm', name=a.name, d=a.d, seed=int(seed), den=cand['den'], err=cand['err'], status=status)
                    if status == 'solved':
                        rec = save_witness_candidate(F, A, q, cand['poly'], data['N'], data['Q'])
                        log_event(event='WITNESS', name=a.name, d=a.d, verified=rec['verified'], F_divides_D=rec['F_divides_D'])
                        print('WITNESS FOUND', file=sys.stderr, flush=True)
                        return 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
