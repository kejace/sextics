"""Exact (mod p) decision of the joint certificate problem on a whole sector family via msolve.

For one (d, sector) the Koszul-form E2 system is bilinear in (a, x):  v(a) + M(a) x_C - MQ x_Q = 0 with
v(a) = T (c kron n), c = B a (D's monomial coefficients), n = N a (linear E1 lift).  Over F_p, with Q eliminated by
the left null space L of MQ, the F-gauge fixed (s.a = 0 for the F-multiple directions), a normalisation lambda.a = 1
and n_g random linear forms on x_C (n_g = generic Koszul-gauge dimension), the solutions over F_p-bar are exactly the
(gauge-fixed) certificate denominators of the family over F_p-bar.  msolve output [-1] = no solution over F_p-bar.
usage: python3 fpt_msolve.py <instance|control:IID> <d> <h,r2,s|auto> [p] [threads]
Writes wave2-search/msolve/<tag>.ms (input), .out (msolve output), and a JSON line to wave2-search/msolve.jsonl."""
import os, sys, json, time, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
import numpy as np, scipy.sparse as sp
from flint import nmod_mat

MSDIR = f'{LANE}/wave2-search/msolve'
os.makedirs(MSDIR, exist_ok=True)
OUT = f'{LANE}/wave2-search/msolve.jsonl'


def nm(rows_list, p, R=None, C=None):
    R = len(rows_list) if R is None else R
    C = len(rows_list[0]) if (C is None and rows_list) else (C or 0)
    return nmod_mat(R, C, [int(v) % p for row in rows_list for v in row], p)


def table(M):
    try:
        return [[int(v) for v in row] for row in M.table()]
    except Exception:
        return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


def nullspace_cols(M):
    X, nul = M.nullspace()
    T = table(X)
    return [[T[i][k] for i in range(len(T))] for k in range(nul)]   # list of vectors (length ncols)


def modp_family(inst, d, lab, p):
    Dm = inst.sector_monomials(d, lab); K0 = len(Dm)
    tl = inst.addlab(inst.labA, lab)
    gens = inst.jac_generators_sector(d + 6, tl, p)
    Ap = reduce_mod(inst.A, p)
    AV = [mul(Ap, {m: 1}, p) for m in Dm]
    cols = AV + gens
    rows = {}
    for P in cols:
        for e in P:
            rows.setdefault(e, len(rows))
    M = nmod_mat(len(rows), len(cols), [0] * (len(rows) * len(cols)), p)
    flat = [0] * (len(rows) * len(cols))
    for j, P in enumerate(cols):
        for e, v in P.items():
            flat[rows[e] * len(cols) + j] = int(v) % p
    M = nmod_mat(len(rows), len(cols), flat, p)
    vecs = nullspace_cols(M)
    proj = [v[:K0] for v in vecs if any(v[:K0])]
    if not proj:
        return Dm, []
    Pr, rk = nm(proj, p).rref()
    Tt = table(Pr)[:rk]
    return Dm, Tt   # basis vectors (rows) in Dm coordinates


def particular_solutions(E, RHS, p):
    """E: list of rows (R x C), RHS: list of rows (R x K). Returns X (C x K) with E X = RHS (free vars 0) or None."""
    R, C = len(E), len(E[0]); K = len(RHS[0])
    aug = nm([E[i] + RHS[i] for i in range(R)], p)
    Mr, rank = aug.rref()
    Tt = table(Mr)
    pivots = []
    for i in range(rank):
        j = 0
        while j < C + K and Tt[i][j] == 0:
            j += 1
        pivots.append(j)
    if any(j >= C for j in pivots):
        return None
    X = [[0] * K for _ in range(C)]
    for i, j in enumerate(pivots):
        for k in range(K):
            X[j][k] = Tt[i][C + k]
    return X


def build_system(inst, d, lab, p, seed=1, verbose=True):
    t0 = time.time()
    F, A, q = inst.F, inst.A, inst.q
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Dm, Bt = modp_family(inst, d, lab, p)
    K0 = len(Dm); K = len(Bt)
    info = {'d': d, 'sector': list(lab), 'p': p, 'K0': K0, 'K': K}
    if K == 0:
        return None, info
    didx = {m: j for j, m in enumerate(Dm)}
    B = [[Bt[k][j] for k in range(K)] for j in range(K0)]       # K0 x K
    D0 = {Dm[0]: 1}
    allowed, subg = symmetry.sector_filter(F, A, q, D0, inst.Gsub)
    pars = fixed_d.numerator_parities(D0, A)
    degN = d + 1; degQ = 2 * d - 6; degC = d - 4
    Ncols = []
    for i in range(NV):
        ms = monomials(degN, pars[i])
        Ncols += [(i, m) for m in symmetry.filter_monomials(ms, 'N', i, allowed)]
    nN = len(Ncols)
    # E1 system mod p
    erows = {}
    E_entries = []
    for nu, (i, m) in enumerate(Ncols):
        for e, v in shift(Fi[i], m, q, p).items():
            erows.setdefault(e, len(erows)); E_entries.append((erows[e], nu, int(v)))
    Bpolys = []
    for k in range(K):
        P = {}
        for j in range(K0):
            if Bt[k][j]:
                P[Dm[j]] = Bt[k][j]
        Bpolys.append(P)
    RHSp = [mul(Ap, P, p) for P in Bpolys]
    for P in RHSp:
        for e in P:
            erows.setdefault(e, len(erows))
    R1 = len(erows)
    E = [[0] * nN for _ in range(R1)]
    for r, c, v in E_entries:
        E[r][c] = (E[r][c] + v) % p
    RHS = [[0] * K for _ in range(R1)]
    for k, P in enumerate(RHSp):
        for e, v in P.items():
            RHS[erows[e]][k] = int(v) % p
    Nmat = particular_solutions(E, RHS, p)      # nN x K
    if Nmat is None:
        info['error'] = 'E1 lift inconsistent mod p'; return None, info
    # rows of the E2 identity
    rows = {}
    def ridx(e):
        v = rows.get(e)
        if v is None:
            v = len(rows); rows[e] = v
        return v
    tri, tci, tvals = [], [], []
    for j, mj in enumerate(Dm):
        mjp = {mj: 1}
        dmj = [pderiv(mjp, i, p) for i in range(NV)]
        for nu, (i, mu) in enumerate(Ncols):
            mup = {mu: 1}
            P = sub(mul(mjp, pderiv(mup, i, p), p), mul(mup, dmj[i], p), p)
            for e, v in P.items():
                tri.append(ridx(e)); tci.append(j * nN + nu); tvals.append(int(v))
    zl = [inst.label(tuple(1 if t == i else 0 for t in range(NV))) for i in range(NV)]
    Ccols = []
    if degC >= 0:
        for i in range(NV):
            for jj in range(i + 1, NV):
                want = inst.addlab(inst.labA, lab, zl[i], zl[jj])
                for mu in monomials(degC):
                    if inst.label(mu) == want:
                        Ccols.append((i, jj, mu))
    nC = len(Ccols)
    Mj_entries = [[] for _ in range(K0)]   # per sector monomial: (row, col, val)
    for j, mj in enumerate(Dm):
        mjp = {mj: 1}
        dmj = [pderiv(mjp, t, p) for t in range(NV)]
        for col, (i, jj, mu) in enumerate(Ccols):
            mup = {mu: 1}
            P = sub(mul(mjp, sub(mul(pderiv(mup, i, p), Fi[jj], p), mul(pderiv(mup, jj, p), Fi[i], p), p), p),
                    mul(mup, sub(mul(Fi[jj], dmj[i], p), mul(Fi[i], dmj[jj], p), p), p), p)
            for e, v in P.items():
                Mj_entries[j].append((ridx(e), col, int(v)))
    labQ = inst.addlab(inst.labA, lab, lab)
    Qcols = [mu for mu in monomials(degQ) if inst.label(mu) == labQ] if degQ >= 0 else []
    nQ = len(Qcols)
    MQ_entries = []
    for jq, mu in enumerate(Qcols):
        for e, v in shift(Fp, mu, p - 1, p).items():
            MQ_entries.append((ridx(e), jq, int(v)))
    R = len(rows)
    info.update({'nN': nN, 'nC': nC, 'nQ': nQ, 'rows': R, 'e1_rows': R1})
    # left null space of MQ
    if nQ:
        MQt = [[0] * R for _ in range(nQ)]
        for r, c, v in MQ_entries:
            MQt[c][r] = (MQt[c][r] + v) % p
        Lvecs = nullspace_cols(nm(MQt, p))          # vectors of length R with MQ^T-nullspace = left null space of MQ
    else:
        Lvecs = [[1 if i == r else 0 for i in range(R)] for r in range(R)]
    L = Lvecs; RL = len(L)
    info['rows_L'] = RL
    # LT = L T  (L residues x T small ints) via numpy int64
    Tsp = sp.csr_matrix((np.array(tvals, dtype=np.int64), (np.array(tri), np.array(tci))), shape=(R, K0 * nN))
    Lnp = np.array(L, dtype=np.int64)                 # RL x R
    LT = (Tsp.T @ Lnp.T).T % p                        # RL x K0*nN  (entries < 2^31 * small * nnz)
    LT = np.ascontiguousarray(LT)
    Bm = nm(B, p); Bt_m = Bm.transpose()
    Nm = nm(Nmat, p, R=nN, C=K)
    # L M^(j) for each j
    LMj = []
    Lm = nm(L, p)
    for j in range(K0):
        Mj = [[0] * nC for _ in range(R)]
        for r, c, v in Mj_entries[j]:
            Mj[r][c] = (Mj[r][c] + v) % p
        LMj.append(table(Lm * nm(Mj, p, R=R, C=nC)) if nC else [[0] * 0 for _ in range(RL)])
    # Koszul gauge dimension at a random a
    rng = random.Random(seed)
    a0 = [rng.randrange(1, p) for _ in range(K)]
    c0 = table(Bm * nm([[x] for x in a0], p))
    c0 = [row[0] for row in c0]
    if nC:
        LM0 = [[0] * nC for _ in range(RL)]
        for j in range(K0):
            if c0[j]:
                for r in range(RL):
                    row = LMj[j][r]
                    LM0[r] = [(LM0[r][c] + c0[j] * row[c]) % p for c in range(nC)]
        rankLM0 = nm(LM0, p, R=RL, C=nC).rank()
        n_g = nC - rankLM0
    else:
        n_g = 0
    info['koszul_gauge_dim_generic'] = n_g
    # equations
    varnames = [f'a{k}' for k in range(K)] + [f'x{c}' for c in range(nC)]
    eqs = []
    for l in range(RL):
        G = nm([list(LT[l, j * nN:(j + 1) * nN]) for j in range(K0)], p, R=K0, C=nN)
        Ql = table(Bt_m * (G * Nm))              # K x K
        terms = []
        for k in range(K):
            for kk in range(k, K):
                cf = (Ql[k][kk] + (Ql[kk][k] if kk != k else 0)) % p if kk != k else Ql[k][k] % p
                if cf:
                    terms.append(f'{cf}*a{k}*a{kk}' if kk != k else f'{cf}*a{k}^2')
        if nC:
            S = [[LMj[j][l][c] for c in range(nC)] for j in range(K0)]   # K0 x nC
            Hl = table(Bt_m * nm(S, p, R=K0, C=nC))                       # K x nC
            for k in range(K):
                for c in range(nC):
                    if Hl[k][c]:
                        terms.append(f'{Hl[k][c]}*a{k}*x{c}')
        if terms:
            eqs.append('+'.join(terms))
    # F-gauge: directions s with B s in span(F g)
    nF = 0
    if d >= 6:
        Gm = inst.sector_monomials(d - 6, lab)
        if Gm:
            FM = [[0] * len(Gm) for _ in range(K0)]
            for jj, g in enumerate(Gm):
                for e, v in shift(Fp, g, 1, p).items():
                    FM[didx[e]][jj] = int(v) % p
            big = [B[j] + [(-FM[j][t]) % p for t in range(len(Gm))] for j in range(K0)]
            Svecs = nullspace_cols(nm(big, p, R=K0, C=K + len(Gm)))
            for s in Svecs:
                sa = s[:K]
                if any(sa):
                    eqs.append('+'.join(f'{v % p}*a{k}' for k, v in enumerate(sa) if v % p)); nF += 1
    info['f_gauge_constraints'] = nF
    lam = [rng.randrange(1, p) for _ in range(K)]
    eqs.append('+'.join(f'{v}*a{k}' for k, v in enumerate(lam)) + f'+{p-1}')
    for g in range(n_g):
        mu = [rng.randrange(1, p) for _ in range(nC)]
        eqs.append('+'.join(f'{v}*x{c}' for c, v in enumerate(mu)))
    info['n_eqs'] = len(eqs); info['n_vars'] = len(varnames); info['build_seconds'] = round(time.time() - t0, 1)
    text = ','.join(varnames) + '\n' + str(p) + '\n' + ',\n'.join(eqs) + '\n'
    if verbose:
        print(f"[msolve-build] d={d} lab={lab} p={p}: K0={K0} K={K} nN={nN} nC={nC} nQ={nQ} rows={R} rows_L={RL} "
              f"koszul-gauge={n_g} F-gauge={nF} -> {len(varnames)} vars, {len(eqs)} eqs ({info['build_seconds']}s)", flush=True)
    return text, info


def run_msolve(text, tag, threads=4, timeout=3600):
    fin = f'{MSDIR}/{tag}.ms'; fout = f'{MSDIR}/{tag}.out'
    open(fin, 'w').write(text)
    t0 = time.time()
    try:
        r = subprocess.run(['msolve', '-f', fin, '-o', fout, '-t', str(threads), '-v', '1'], capture_output=True, text=True, timeout=timeout)
        err = r.stderr[-800:]
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'seconds': round(time.time() - t0, 1)}
    except FileNotFoundError:
        return {'status': 'msolve-not-found', 'seconds': 0}
    out = open(fout).read() if os.path.exists(fout) else ''
    head = out[:200].strip()
    if head.startswith('[-1]'):
        status = 'EMPTY'
    elif head.startswith('[1,') and '-1' in head[:40]:
        status = 'POSITIVE-DIM'
    elif head.startswith('[0,') or head.startswith('[2,'):
        status = 'ZERO-DIM'
    else:
        status = 'unknown'
    # degree (number of solutions over the closure) for zero-dim: third entry
    deg = None
    try:
        parts = head.replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
        if status == 'ZERO-DIM':
            deg = int(parts[2])
    except Exception:
        pass
    return {'status': status, 'degree': deg, 'seconds': round(time.time() - t0, 1), 'head': head[:120], 'stderr': err[-300:]}


if __name__ == '__main__':
    target = sys.argv[1]; d = int(sys.argv[2]); labarg = sys.argv[3]
    p = int(sys.argv[4]) if len(sys.argv) > 4 else 2147483647
    threads = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    if target.startswith('control:'):
        iid = target.split(':', 1)[1]
        inst0 = Instance(iid); dw, Dw, Nw, Qw = core.load_witness(iid)
        subg = symmetry.stabilizer(inst0.G, [inst0.A, Dw]); inst = Instance(iid, gens=subg)
        labs = [label_of_poly(inst, Dw)] if labarg == 'auto' else [tuple(int(x) for x in labarg.split(','))]
        tagbase = f'ctl-{iid}'
    else:
        inst = Instance(target)
        labs = sorted(inst.labels_in_degree(d)) if labarg == 'auto' else [tuple(int(x) for x in labarg.split(','))]
        tagbase = target
    for lab in labs:
        tag = f"{tagbase}-d{d}-{'_'.join(str(x) for x in lab)[:40]}-p{p}"
        text, info = build_system(inst, d, lab, p)
        rec = {'kind': 'msolve', 'instance': inst.iid, 'tag': tag, **info, 'engine': ENGINE}
        if text is None:
            rec['status'] = 'skipped (empty family or E1 lift failure)'
        else:
            rec.update(run_msolve(text, tag, threads=threads))
        log_line(OUT, rec)
        print(json.dumps({k: rec[k] for k in rec if k not in ('kind',)}), flush=True)
