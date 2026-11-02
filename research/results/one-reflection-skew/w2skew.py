"""Wave-2 joint-search library for one-reflection-skew (q = 1, A = X0^3 Y0^3).

Five-variable reduction (F = f(x) + Y2^6, x = (X0,X1,X2,Y0,Y1)):
    A D^2 - D N(f) + (6/5) f (D div N - N(D)) = 0,   deg D = d, deg N_i = d+1 (5 components)
A solution lifts to a strict 6-variable certificate  (D^2, (N D, Y2 (D div N - N(D))/5), Q = 0)  of degree 2d
(checked below with core.verify_witness).  Sectors: rho = (Y0,Y1) sign, sigma = (X1,Y1) sign.

Variable projection over a linear family D = sum_k a_k B_k (B_k given as rows of a float matrix over the
sector monomials): f(a) = min_N || M(D(a)) N + A D(a)^2 ||^2 with ridge-regularized normal equations;
gradient by the Golub-Pereyra formula.  Everything here is scratch code: the engine is not modified.
"""
import sys, os, time, json, math, fcntl, itertools, pickle
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
from fractions import Fraction
from flint import nmod_mat
from hodge import core, fixed_d, jacobian
from hodge.core import mul, add, sub, scale, pderiv, shift
import five
from five import f, A, A_plus, A_minus, X, P1, P2, rpar, spar

LANE = 'one-reflection-skew'
WD = f'/tmp/hodge/results/{LANE}'
LOG = f'{WD}/wave2-search.log.jsonl'
REPORT = f'/tmp/hodge/notes/wave2-{LANE}-search.md'
F6, A6, q6, _rec = core.load_instance(LANE)
assert A6 == A and q6 == 1
C65 = Fraction(6, 5)

# ----------------------------------------------------------------------------- logging
def logbatch(**rec):
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    rec.setdefault('engine', 2)
    with open(LOG, 'a') as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.write(json.dumps(rec, default=str) + '\n')
        fcntl.flock(fh, fcntl.LOCK_UN)
    return rec

# ----------------------------------------------------------------------------- monomials / sectors (5-var)
def mons5(deg, rD=None, sD=None):
    out = []
    if deg < 0:
        return out
    for c in itertools.combinations(range(deg + 4), 4):
        e = []; prev = -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + 3 - prev)
        e = tuple(e) + (0,)
        if (rD is None or rpar(e) == rD) and (sD is None or spar(e) == sD):
            out.append(e)
    return out

RX = [0, 0, 0, 1, 1]      # rho-parity of x_i
SX = [0, 1, 0, 0, 1]      # sigma-parity of x_i

def npar5(rD, sD, rT=1, sT=0):
    """(rho, sigma) parities of N_i for D of parity (rD,sD) and target of parity (rT,sT)."""
    return [((rT - rD - RX[i]) % 2, (sT - sD - SX[i]) % 2) for i in range(5)]

def e1_kernel5(d, rD, sD, p=P1, want_float=True, verbose=False):
    """E1 kernel {D in S5_d[rD,sD] : A D in J_f} : exact dims mod p and a float orthonormal basis.

    Returns dict(mons, nD, kerD, dimJ, nonJ, Kbasis (kerD x nD), Jbasis (dimJ x nD), Cbasis (nonJ x nD)).
    Cbasis = kernel vectors orthogonal to the span of J_d (float)."""
    t0 = time.time()
    Dm = mons5(d, rD, sD)
    nD = len(Dm)
    rT, sT = (1 + rD) % 2, sD           # parity of A D
    # columns: A*m (D monomials), f_i * m' (m' of degree d+1, parity so that product has (rT,sT))
    fp = core.reduce_mod(f, p); Ap = core.reduce_mod(A, p)
    fi_p = [pderiv(fp, i, p) for i in range(5)]
    fi_q = [pderiv(f, i) for i in range(5)]
    Acols_p = [shift(Ap, m, 1, p) for m in Dm]
    Jm = [(i, m) for i in range(5) for m in mons5(d + 1, (rT - RX[i]) % 2, (sT - SX[i]) % 2)]
    Jcols_p = [shift(fi_p[i], m, 1, p) for i, m in Jm]
    rows = {}
    for col in Acols_p + Jcols_p:
        for e in col:
            rows.setdefault(e, len(rows))
    R = len(rows); C = nD + len(Jm)
    flat = [0] * (R * C)
    for j, col in enumerate(Acols_p + Jcols_p):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    Xn, nullity = M.nullspace()
    if nullity:
        proj = nmod_mat(nD, nullity, [int(Xn[j, k]) for j in range(nD) for k in range(nullity)], p)
        Tr, kerD = proj.transpose().rref()
    else:
        kerD = 0
    # J_d in the sector
    Jd = [(i, m) for i in range(5) for m in mons5(d - 5, (rD - RX[i]) % 2, (sD - SX[i]) % 2)] if d >= 5 else []
    idxD = {m: k for k, m in enumerate(Dm)}
    if Jd:
        flatJ = [0] * (nD * len(Jd))
        colsJ_p = [shift(fi_p[i], m, 1, p) for i, m in Jd]
        for j, col in enumerate(colsJ_p):
            for e, v in col.items():
                flatJ[idxD[e] * len(Jd) + j] = v
        dimJ = nmod_mat(nD, len(Jd), flatJ, p).rank()
    else:
        dimJ = 0
    out = {'d': d, 'rD': rD, 'sD': sD, 'mons': Dm, 'nD': nD, 'kerD': kerD, 'dimJ': dimJ, 'nonJ': kerD - dimJ,
           'rows': R, 'cols': C, 'p': p}
    if want_float and kerD:
        # float: projector onto complement of J_{d+6}, kernel of projected A-columns
        Jmat = np.zeros((R, len(Jm)))
        for j, (i, m) in enumerate(Jm):
            for e, v in shift(fi_q[i], m, 1).items():
                Jmat[rows[e], j] = float(v)
        Amat = np.zeros((R, nD))
        for j, m in enumerate(Dm):
            for e, v in shift(A, m, 1).items():
                Amat[rows[e], j] = float(v)
        Qj, Rj, piv = sla.qr(Jmat, mode='economic', pivoting=True)
        diag = np.abs(np.diag(Rj))
        rk = int(np.sum(diag > 1e-10 * diag[0])) if diag.size else 0
        Qj = Qj[:, :rk]
        T = Amat - Qj @ (Qj.T @ Amat)
        U, s, Vt = np.linalg.svd(T, full_matrices=True)
        # kernel = last kerD right singular vectors
        Kb = Vt[nD - kerD:, :]
        s_in = s[nD - kerD - 1] if nD - kerD - 1 >= 0 else np.inf
        s_out = s[nD - kerD] if nD - kerD < len(s) else 0.0
        out['sv_gap'] = (float(s_in), float(s_out))
        out['Kbasis'] = Kb
        # J_d float basis (orthonormal) inside the sector, and its complement in the kernel
        if Jd:
            JdM = np.zeros((nD, len(Jd)))
            for j, (i, m) in enumerate(Jd):
                for e, v in shift(fi_q[i], m, 1).items():
                    JdM[idxD[e], j] = float(v)
            Qd, Rd, pv = sla.qr(JdM, mode='economic', pivoting=True)
            dg = np.abs(np.diag(Rd)); rkd = int(np.sum(dg > 1e-10 * dg[0]))
            Jb = Qd[:, :rkd].T                      # rkd x nD
            assert rkd == dimJ, (rkd, dimJ)
            # complement: kernel vectors orthogonal to Jb
            Cb = Kb - (Kb @ Jb.T) @ Jb
            Uc, sc, Vtc = np.linalg.svd(Cb, full_matrices=False)
            rc = int(np.sum(sc > 1e-8))
            Cb = Vtc[:rc, :]
            out['Jbasis'] = Jb; out['Cbasis'] = Cb
            out['nonJ_float'] = rc
        else:
            out['Jbasis'] = np.zeros((0, nD)); out['Cbasis'] = Kb; out['nonJ_float'] = kerD
    out['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k not in ('mons', 'Kbasis', 'Jbasis', 'Cbasis')}), file=sys.stderr)
    return out

def poly_from_vec(vec, mons, tol=0.0):
    return {m: float(v) for m, v in zip(mons, vec) if abs(v) > tol}

def vec_from_poly(P, mons):
    idx = {m: k for k, m in enumerate(mons)}
    v = np.zeros(len(mons))
    for e, c in P.items():
        v[idx[e]] = float(c)
    return v

# ----------------------------------------------------------------------------- VarPro, 5 variables
class VP5:
    """Variable projection for the 5-variable equation over D = sum_k a_k B_k (rows of Basis over mons)."""
    def __init__(self, d, rD, sD, mons, Basis, target=None, sigma=True, ridge=1e-10, verbose=False):
        t0 = time.time()
        self.d, self.rD, self.sD = d, rD, sD
        self.mons = list(mons); self.nD = len(mons)
        self.Basis = np.asarray(Basis, float)            # K x nD
        self.K = self.Basis.shape[0]
        self.target = A if target is None else target
        rT, sT = five.rparity_of(self.target), five.sparity_of(self.target)
        assert rT is not None
        if sT is None:
            sigma = False; sT = 0
        self.sigma = sigma
        pars = npar5(rD, sD, rT, sT)
        self.pars = pars
        Ncols = []
        for i in range(5):
            for m in mons5(d + 1, pars[i][0], pars[i][1] if sigma else None):
                Ncols.append((i, m))
        self.Ncols = Ncols; self.nN = len(Ncols)
        fi = [pderiv(f, i) for i in range(5)]
        rows = {}
        def ridx(e):
            if e not in rows:
                rows[e] = len(rows)
            return rows[e]
        # tensor entries (e-index, row, col, value) for monomial D = x^e
        EE, RR, CC, VV = [], [], [], []
        for k, e in enumerate(self.mons):
            De = {e: Fraction(1)}
            fD = scale(mul(f, De), C65)
            G = []
            for i in range(5):
                Di = pderiv(De, i)
                G.append(scale(add(mul(De, fi[i]), scale(mul(f, Di), C65)), -1))
            for j, (i, m) in enumerate(Ncols):
                col = shift(G[i], m, 1)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(fD, tuple(m2), m[i]))
                for ee, v in col.items():
                    EE.append(k); RR.append(ridx(ee)); CC.append(j); VV.append(float(v))
        # rhs: A * x^{e1} * x^{e2} index table
        self.prod_idx = np.zeros((self.nD, self.nD), dtype=np.int64)
        Ae = next(iter(A.keys())); Ac = float(A[Ae])
        for k1, e1 in enumerate(self.mons):
            for k2, e2 in enumerate(self.mons):
                ee = tuple(a + b + c for a, b, c in zip(e1, e2, Ae))
                self.prod_idx[k1, k2] = ridx(ee)
        self.Ac = Ac
        self.rows = len(rows); self.row_index = rows
        self.EE = np.array(EE, dtype=np.int64); self.RR = np.array(RR, dtype=np.int64)
        self.CC = np.array(CC, dtype=np.int64); self.VV = np.array(VV)
        self.flat = self.RR * self.nN + self.CC
        self.ridge = ridge
        self.build_time = time.time() - t0
        self._nfev = 0
        if verbose:
            print(f'[vp5] d={d} sector=({rD},{sD}) K={self.K} nD={self.nD} rows={self.rows} cols={self.nN} '
                  f'nnz={len(VV)} built in {self.build_time:.1f}s', file=sys.stderr)

    def Dvec(self, a):
        return self.Basis.T @ np.asarray(a, float)

    def assemble(self, Dv):
        M = np.bincount(self.flat, weights=Dv[self.EE] * self.VV, minlength=self.rows * self.nN).reshape(self.rows, self.nN)
        b = -self.Ac * np.bincount(self.prod_idx.ravel(), weights=np.outer(Dv, Dv).ravel(), minlength=self.rows)
        return M, b

    def solve_inner(self, M, b):
        G = M.T @ M
        tr = np.trace(G) / max(G.shape[0], 1)
        G[np.diag_indices_from(G)] += self.ridge * tr
        c = M.T @ b
        try:
            cf = sla.cho_factor(G, lower=True, check_finite=False)
            x = sla.cho_solve(cf, c, check_finite=False)
        except Exception:
            x = sla.lstsq(M, b, cond=1e-10, check_finite=False)[0]
        return x

    def value(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        Dv = self.Dvec(a)
        M, b = self.assemble(Dv)
        x = self.solve_inner(M, b)
        r = M @ x - b
        self._nfev += 1
        return float(r @ r), x, r, b

    def value_and_grad(self, u):
        u = np.asarray(u, float); nu = np.linalg.norm(u); a = u / nu
        Dv = self.Dvec(a)
        M, b = self.assemble(Dv)
        x = self.solve_inner(M, b)
        r = M @ x - b
        fval = float(r @ r)
        # d/dD_e of r.(M x): z[e] = sum_{nnz with index e} r[row]*v*x[col]
        z = np.bincount(self.EE, weights=r[self.RR] * self.VV * x[self.CC], minlength=self.nD)
        # d/dD_e of r.b: b = -Ac * sum_{e1,e2} D[e1] D[e2] delta(idx) -> db/dD_e = -2 Ac * sum_{e1} D[e1] r[idx(e1,e)]
        y = -2.0 * self.Ac * (r[self.prod_idx] @ Dv)      # (nD x nD) @ nD -> y[e1]... symmetric table
        gD = 2.0 * (z - y)
        g = self.Basis @ gD
        g_u = (g - float(g @ a) * a) / nu
        self._nfev += 1
        return fval, g_u

    def bnorm2(self, a):
        a = np.asarray(a, float); a = a / np.linalg.norm(a)
        _, b = self.assemble(self.Dvec(a))
        return float(b @ b)

    def search(self, n_restarts=10, seed=0, maxiter=400, gtol=1e-12, verbose=False, stop_rel=1e-13, init=None,
               keep_rel=1e-6, on_restart=None):
        rng = np.random.default_rng(seed)
        out = []
        t0 = time.time()
        for it in range(n_restarts):
            t1 = time.time()
            if init is not None and it == 0:
                u0 = np.asarray(init, float)
            else:
                u0 = rng.standard_normal(self.K)
            u0 /= np.linalg.norm(u0)
            res = minimize(self.value_and_grad, u0, jac=True, method='L-BFGS-B',
                           options={'maxiter': maxiter, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 20})
            a = res.x / np.linalg.norm(res.x)
            fval, x, r, b = self.value(a)
            rel = fval / max(float(b @ b), 1e-300)
            cand = {'rel': rel, 'value': fval, 'iters': int(res.nit), 'restart': it, 'seed': seed,
                    'secs': round(time.time() - t1, 2), 'a': a if rel < keep_rel else None}
            out.append(cand)
            if on_restart:
                on_restart(cand)
            if verbose:
                print(f'[vp5] seed {seed} restart {it}: rel={rel:.3e} iters={res.nit} ({time.time()-t0:.1f}s)', file=sys.stderr)
            if rel < stop_rel:
                break
        return out

    def D_float(self, a):
        return poly_from_vec(self.Dvec(np.asarray(a, float) / np.linalg.norm(a)), self.mons)

# ----------------------------------------------------------------------------- reconstruction + exact confirmation
def reconstruct_D(Dfl, max_dens=(10, 100, 1000, 10000, 100000), tol=1e-7):
    """Rational candidates for a float polynomial (normalized so the largest |coeff| is 1).  Returns list of exact polys."""
    items = sorted(Dfl.items(), key=lambda kv: -abs(kv[1]))
    if not items:
        return []
    top = items[0][1]
    vals = {e: v / top for e, v in items}
    cands = []
    for md in max_dens:
        fr = {e: Fraction(v).limit_denominator(md) for e, v in vals.items()}
        err = max(abs(float(c) - v) for (e, c), v in zip(fr.items(), vals.values()))
        if err < tol * max(1, md ** 0.5):
            P = {e: c for e, c in fr.items() if c != 0}
            if P and P not in cands:
                cands.append(P)
    return cands

def confirm_exact5(D, target=None, verbose=False):
    """Exact 5-var fixed-D test (mod p served check, then rational lift and 6-var verification)."""
    target = A if target is None else target
    r = five.analyze(D, use_sigma=(five.sparity_of(D) is not None and five.sparity_of(target) == 0))
    key = 'serves_A' if target == A else ('serves_A+' if target == A_plus else 'serves_A-')
    if not r.get(key):
        return 'inconsistent', {'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'], 'rank_N': r['rank_N']}
    st, data = five.solve_exact(D, target, verbose=verbose)
    data['served_dim'] = r['served_dim']
    return st, data

def save_witness6(D6, N6, Q6, note=''):
    ok, info = core.verify_witness(F6, A6, 1, D6, N6, Q6)
    w = {'instance': LANE, 'q': 1, 'd': core.degree(D6), 'D': core.to_str(D6), 'N': [core.to_str(n) for n in N6],
         'Q': core.to_str(Q6) if Q6 else '0', 'verified': ok, 'info': str(info), 'note': note}
    with open(f'{WD}/witness-candidate.json', 'w') as fh:
        json.dump(w, fh, indent=1)
    return ok, info

# ----------------------------------------------------------------------------- report regeneration
def rewrite_report(status_line=None):
    recs = []
    try:
        with open(LOG) as fh:
            for line in fh:
                try:
                    recs.append(json.loads(line))
                except Exception:
                    pass
    except FileNotFoundError:
        pass
    head = ''
    try:
        with open(REPORT) as fh:
            txt = fh.read()
        head = txt.split('## Batches')[0]
    except FileNotFoundError:
        head = f'# Wave 2 — {LANE} — SEARCH lens\n\nSTATUS: running\n\n'
    if status_line:
        lines = head.split('\n')
        for i, l in enumerate(lines):
            if l.startswith('STATUS:'):
                lines[i] = 'STATUS: ' + status_line
                break
        head = '\n'.join(lines)
    rows = ['| ts | batch | family | d | sector | restarts | best rel | min/median iters | secs | note |', '|---|---|---|---|---|---|---|---|---|---|']
    for r in recs:
        br = r.get('best_rel')
        rows.append(f"| {r.get('ts','')[5:16]} | {r.get('batch','')} | {r.get('family','')} | {r.get('degree','')} | {r.get('sector','')} | "
                    f"{r.get('restarts','')} | {('%.2e' % br) if isinstance(br,(int,float)) else br} | {r.get('iters','')} | {r.get('seconds','')} | {str(r.get('note',''))[:160]} |")
    body = head + '## Batches\n' + '\n'.join(rows) + '\n'
    with open(REPORT + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        with open(REPORT, 'w') as fh:
            fh.write(body)
        fcntl.flock(lk, fcntl.LOCK_UN)

# ----------------------------------------------------------------------------- exact rational E1 kernel (CRT over primes)
def _kernel_rref_modp(d, rD, sD, p):
    """Row-reduced basis (pivots, rows of ints mod p) of {D in S5_d[rD,sD] : A D in J_f}; also rank J_d."""
    Dm = mons5(d, rD, sD); nD = len(Dm)
    rT, sT = (1 + rD) % 2, sD
    fp = core.reduce_mod(f, p); Ap = core.reduce_mod(A, p)
    fi_p = [pderiv(fp, i, p) for i in range(5)]
    Acols = [shift(Ap, m, 1, p) for m in Dm]
    Jm = [(i, m) for i in range(5) for m in mons5(d + 1, (rT - RX[i]) % 2, (sT - SX[i]) % 2)]
    Jcols = [shift(fi_p[i], m, 1, p) for i, m in Jm]
    rows = {}
    for col in Acols + Jcols:
        for e in col:
            rows.setdefault(e, len(rows))
    R = len(rows); C = nD + len(Jm)
    flat = [0] * (R * C)
    for j, col in enumerate(Acols + Jcols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    Xn, nullity = M.nullspace()
    if not nullity:
        return Dm, [], [], R, C
    proj = nmod_mat(nD, nullity, [int(Xn[j, k]) for j in range(nD) for k in range(nullity)], p)
    Tr, kerD = proj.transpose().rref()
    rowsK = []; pivs = []
    for i in range(kerD):
        row = [int(Tr[i, j]) for j in range(nD)]
        piv = next(j for j, v in enumerate(row) if v)
        pivs.append(piv); rowsK.append(row)
    return Dm, pivs, rowsK, R, C

def e1_kernel5_exact(d, rD, sD, primes=None, max_primes=6, verbose=False):
    """Exact rational rref basis of the E1 kernel (CRT + rational reconstruction, checked with one extra prime).

    Returns dict(mons, nD, kerD, dimJ, nonJ, Kexact (list of Fraction rows), Kbasis (orthonormal float, kerD x nD),
                 Jbasis (dimJ x nD float), Cbasis (complement of J_d inside the kernel, float)).
    """
    t0 = time.time()
    primes = list(primes or core.PRIMES[:max_primes])
    acc = None; mod = None; pivs0 = None; Dm = None
    result = None
    for k, p in enumerate(primes):
        Dm, pivs, rowsK, R, C = _kernel_rref_modp(d, rD, sD, p)
        if pivs0 is None:
            pivs0 = pivs; acc = rowsK; mod = p
        elif pivs != pivs0:
            # unlucky prime; restart accumulation with this one
            pivs0 = pivs; acc = rowsK; mod = p
            continue
        else:
            acc = [[core.crt_pair(a, mod, b, p)[0] for a, b in zip(ra, rb)] for ra, rb in zip(acc, rowsK)]
            mod *= p
        rat = [[core.ratrecon(v, mod) for v in row] for row in acc]
        if all(v is not None for row in rat for v in row):
            # verify with the next prime
            pn = primes[k + 1] if k + 1 < len(primes) else core.PRIMES[max_primes]
            Dm2, pivs2, rowsK2, _, _ = _kernel_rref_modp(d, rD, sD, pn)
            ok = pivs2 == pivs0 and all((v.numerator * pow(v.denominator, -1, pn) - w) % pn == 0
                                        for row, row2 in zip(rat, rowsK2) for v, w in zip(row, row2))
            if ok:
                result = rat; break
    if result is None:
        raise RuntimeError(f'rational reconstruction of the E1 kernel failed d={d} sector=({rD},{sD})')
    nD = len(Dm); kerD = len(result)
    out = {'d': d, 'rD': rD, 'sD': sD, 'mons': Dm, 'nD': nD, 'kerD': kerD, 'Kexact': result, 'primes_used': k + 1,
           'rows': R, 'cols': C}
    # J_d exact generators
    fi_q = [pderiv(f, i) for i in range(5)]
    idxD = {m: j for j, m in enumerate(Dm)}
    Jd = [(i, m) for i in range(5) for m in mons5(d - 5, (rD - RX[i]) % 2, (sD - SX[i]) % 2)] if d >= 5 else []
    JdM = np.zeros((nD, len(Jd)))
    for j, (i, m) in enumerate(Jd):
        for e, v in shift(fi_q[i], m, 1).items():
            JdM[idxD[e], j] = float(v)
    if Jd:
        p = core.PRIMES[0]
        fi_p = [pderiv(core.reduce_mod(f, p), i, p) for i in range(5)]
        flatJ = [0] * (nD * len(Jd))
        for j, (i, m) in enumerate(Jd):
            for e, v in shift(fi_p[i], m, 1, p).items():
                flatJ[idxD[e] * len(Jd) + j] = v
        dimJ = nmod_mat(nD, len(Jd), flatJ, p).rank()
        Qd, Rd, pv = sla.qr(JdM, mode='economic', pivoting=True)
        Jb = Qd[:, :dimJ].T
    else:
        dimJ = 0; Jb = np.zeros((0, nD))
    out['dimJ'] = dimJ; out['nonJ'] = kerD - dimJ
    Kf = np.array([[float(v) for v in row] for row in result]) if kerD else np.zeros((0, nD))
    if kerD:
        Qk, Rk = np.linalg.qr(Kf.T)
        Kb = Qk.T                                   # orthonormal rows spanning the kernel
        Cb = Kb - (Kb @ Jb.T) @ Jb if dimJ else Kb
        Uc, sc, Vtc = np.linalg.svd(Cb, full_matrices=False)
        out['comp_sv'] = (float(sc[kerD - dimJ - 1]) if kerD - dimJ > 0 else None, float(sc[kerD - dimJ]) if kerD - dimJ < len(sc) else None)
        Cb = Vtc[:kerD - dimJ, :]
    else:
        Kb = np.zeros((0, nD)); Cb = Kb
    out['Kbasis'] = Kb; out['Jbasis'] = Jb; out['Cbasis'] = Cb
    out['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(json.dumps({kk: v for kk, v in out.items() if kk not in ('mons', 'Kexact', 'Kbasis', 'Jbasis', 'Cbasis')}), file=sys.stderr)
    return out

def load_kernel(d, rD, sD, exact=True):
    fn = f'{WD}/w2ker5x_d{d}_r{rD}s{sD}.pkl' if exact else f'{WD}/w2ker5_d{d}_r{rD}s{sD}.pkl'
    if os.path.exists(fn):
        with open(fn, 'rb') as fh:
            return pickle.load(fh)
    k = e1_kernel5_exact(d, rD, sD) if exact else e1_kernel5(d, rD, sD)
    with open(fn, 'wb') as fh:
        pickle.dump(k, fh)
    return k
