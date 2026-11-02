"""Exact 'Hessian relaxation' of the polar (Jacobian-type) family at degree d.

Every E1-admissible denominator of degree d <= 9 is D = w.grad F (K_d = J_d).  Writing N = A w + c.grad F
(the H z/6 part cancels identically), the certificate identity E = F^2 Q is equivalent to
      Phi(w) + L(c, w) = F Q,   Phi(w) = D (grad A . w + A div w) - A w.grad D,
      L(c, w) = sum_{i<j} F_j[(d_i c_ij) D - c_ij d_i D] - F_i[(d_j c_ij) D - c_ij d_j D],
quadratic in w, bilinear in (c, w), linear in Q.  Relaxation: W_ab := w_a w_b, V_{ga} := c_g w_a are new
independent unknowns; the system is LINEAR in (W, V, Q).  Any true solution (w, c, Q) gives a relaxation solution,
so the set of admissible w is contained in {w : w w^T in pi_W(ker)}.  If pi_W(ker) = 0 over Q, no polar
denominator of this degree (and parity sector) admits a certificate for A -- rigorous (exact integer kernel).

usage: python3 wave2_excl_relax.py <jobname> <d> <parity: 0|1|all> [exact=1] [msolve=1]
"""
import sys, os, json, time, subprocess, traceback
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, q, monomials, mul, add, sub, scale, shift, pderiv, power, logline, var
from fractions import Fraction
from math import gcd
from flint import nmod_mat, fmpz_mat
import random

JOB = sys.argv[1]
D_DEG = int(sys.argv[2])
PAR = sys.argv[3]
EXACT = True
MSOLVE = True
for a in sys.argv[4:]:
    if a.startswith('exact='): EXACT = a[6:] == '1'
    if a.startswith('msolve='): MSOLVE = a[7:] == '1'
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')
P = core.PRIMES[1]
NV = 6
pA = core.parity_of(A)
parz = [0, 0, 0, 1, 1, 1]
Fi = [pderiv(F, i) for i in range(NV)]
Ai = [pderiv(A, i) for i in range(NV)]


def mono(m):
    return {tuple(m): Fraction(1)}


def setup(d, par):
    """Index sets for w, c, Q and the polynomial data."""
    pD = None if par == 'all' else int(par)
    # w basis
    wb = []
    for i in range(NV):
        pw = None if pD is None else (pD - parz[i]) % 2
        for m in monomials(d - 5, pw):
            wb.append((i, m))
    # numerator parities (for restricted sector): par(N_i) = pA - pD - parz[i]
    cb = []
    for i in range(NV):
        for j in range(i + 1, NV):
            pc = None if pD is None else (pA - pD - parz[i] - parz[j]) % 2
            for m in monomials(d - 4, pc):
                cb.append((i, j, m))
    qb = monomials(2 * d - 6, None if pD is None else pA)
    return wb, cb, qb


def D_of(wvec, wb):
    D = {}
    for (i, m), c in zip(wb, wvec):
        if c:
            D = add(D, scale(shift(Fi[i], m, 1), c))
    return D


def build(d, par, p=None):
    """p=None: exact Fractions; p given: all arithmetic mod p (fast, screening / rank only)."""
    wb, cb, qb = setup(d, par)
    nw, nc, nq = len(wb), len(cb), len(qb)
    Fp = core.reduce_mod(F, p) if p else F
    Ap = core.reduce_mod(A, p) if p else A
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    Aip = [pderiv(Ap, i, p) for i in range(NV)]
    mono_ = (lambda m: {tuple(m): 1}) if p else mono
    # polar basis D_a = m_a F_{i_a}, and its partials
    Da = [shift(Fip[i], m, 1, p) for (i, m) in wb]
    Dai = [[pderiv(D, k, p) for k in range(NV)] for D in Da]
    cols = []   # (kind, index tuple, poly)
    # W columns (a <= b): Phi_ab + Phi_ba
    def Phi(a, b):
        ib, mb = wb[b]
        mbp = mono_(mb)
        t = mul(Da[a], add(mul(Aip[ib], mbp, p), mul(Ap, pderiv(mbp, ib, p), p), p), p)    # D_a (dA.e_b + A div e_b)
        t = sub(t, mul(Ap, mul(mbp, Dai[a][ib], p), p), p)                                # - A e_b . grad D_a
        return t
    for a in range(nw):
        for b in range(a, nw):
            t = Phi(a, b) if a == b else add(Phi(a, b), Phi(b, a), p)
            cols.append(('W', (a, b), t))
    # V columns (gamma, a): L(e_gamma, e_a)
    for g, (i, j, m) in enumerate(cb):
        mp_ = mono_(m)
        dmi, dmj = pderiv(mp_, i, p), pderiv(mp_, j, p)
        for a in range(nw):
            t = mul(Fip[j], sub(mul(dmi, Da[a], p), mul(mp_, Dai[a][i], p), p), p)
            t = sub(t, mul(Fip[i], sub(mul(dmj, Da[a], p), mul(mp_, Dai[a][j], p), p), p), p)
            cols.append(('V', (g, a), t))
    for m in qb:
        cols.append(('Q', m, scale(shift(Fp, m, 1, p), -1, p)))
    rows = {}
    for _, _, P_ in cols:
        for e in P_:
            rows.setdefault(e, len(rows))
    return wb, cb, qb, cols, rows


def selftest(d, par, wb, cb):
    """Check E(w, c) = F (Phi(w) + L(c,w)) exactly for random small w, c."""
    random.seed(5)
    wvec = [Fraction(random.randint(-2, 2)) for _ in wb]
    cvec = [Fraction(random.randint(-2, 2)) for _ in cb]
    D = D_of(wvec, wb)
    N = [mul(A, D_of_w(wvec, wb, i)) for i in range(NV)]
    cmat = [[{} for _ in range(NV)] for _ in range(NV)]
    for (i, j, m), c in zip(cb, cvec):
        if c:
            cmat[i][j] = add(cmat[i][j], scale(mono(m), c))
            cmat[j][i] = sub(cmat[j][i], scale(mono(m), c))
    for i in range(NV):
        for j in range(NV):
            if cmat[i][j]:
                N[i] = add(N[i], mul(cmat[i][j], Fi[j]))
    E = core.defect(F, A, D, q, N)
    # Phi + L
    Phi = sub(mul(D, add(sum_(mul(Ai[i], D_of_w(wvec, wb, i)) for i in range(NV)), mul(A, sum_(pderiv(D_of_w(wvec, wb, i), i) for i in range(NV))))),
              mul(A, sum_(mul(D_of_w(wvec, wb, i), pderiv(D, i)) for i in range(NV))))
    L = {}
    for i in range(NV):
        for j in range(NV):
            if cmat[i][j]:
                L = add(L, mul(Fi[j], sub(mul(pderiv(cmat[i][j], i), D), mul(cmat[i][j], pderiv(D, i)))))
    return sub(E, mul(F, add(Phi, L))) == {}


def sum_(it):
    out = {}
    for t in it:
        out = add(out, t)
    return out


def D_of_w(wvec, wb, i):
    out = {}
    for (ii, m), c in zip(wb, wvec):
        if ii == i and c:
            out = add(out, scale(mono(m), c))
    return out


def main():
    t0 = time.time()
    d, par = D_DEG, PAR
    wb, cb, qb, cols, rows = build(d, par, None if EXACT else P)
    nw, nc, nq = len(wb), len(cb), len(qb)
    nW = nw * (nw + 1) // 2
    nV = nc * nw
    R, C = len(rows), len(cols)
    ok = selftest(d, par, wb, cb) if d <= 5 or EXACT else 'skipped'
    tb = time.time() - t0
    rec = dict(exp='relax', d=d, parity=par, nw=nw, nc=nc, nq=nq, nW=nW, nV=nV, rows=R, cols=C, selftest_E_eq_F_PhiL=ok, job=JOB)
    print(json.dumps(rec), flush=True)
    # modular ranks first
    def rp(v):
        return (v.numerator * pow(v.denominator, -1, P)) % P if isinstance(v, Fraction) else int(v) % P
    flat = [0] * (R * C)
    for j, (_, _, P_) in enumerate(cols):
        for e, v in P_.items():
            flat[rows[e] * C + j] = rp(v)
    Mp = nmod_mat(R, C, flat, P)
    rank_full = Mp.rank()
    del Mp
    # VQ block
    flatVQ = [0] * (R * (C - nW))
    for j in range(nW, C):
        for e, v in cols[j][2].items():
            flatVQ[rows[e] * (C - nW) + (j - nW)] = rp(v)
    rank_VQ = nmod_mat(R, C - nW, flatVQ, P).rank()
    del flatVQ
    dimW_p = nW - (rank_full - rank_VQ)
    rec.update(prime=P, rank_full_modp=rank_full, rank_VQ_modp=rank_VQ, dim_piW_ker_modp=dimW_p, seconds_build=round(tb, 1), seconds_modp=round(time.time() - t0, 1))
    print(json.dumps({k: rec[k] for k in ('rank_full_modp', 'rank_VQ_modp', 'dim_piW_ker_modp')}), flush=True)
    logline(LOG, dict(rec, stage='modular'))
    if EXACT:
        t1 = time.time()
        # exact integer matrix: scale each column by lcm of denominators
        flatZ = [0] * (R * C)
        for j, (_, _, P_) in enumerate(cols):
            den = 1
            for v in P_.values():
                den = den * v.denominator // gcd(den, v.denominator)
            for e, v in P_.items():
                flatZ[rows[e] * C + j] = int(v * den)
        MZ = fmpz_mat(R, C, flatZ)
        del flatZ
        X, nul = MZ.nullspace()
        del MZ
        rank_exact = C - nul
        # projection to W coordinates (first nW rows of X), exact rank
        Wrows = [[int(X[a, k]) for k in range(nul)] for a in range(nW)]
        Wm = fmpz_mat(nW, nul, [x for r in Wrows for x in r])
        rW = Wm.rank()
        rec.update(rank_full_exact=rank_exact, nullity_exact=nul, dim_piW_ker_exact=rW, seconds_exact=round(time.time() - t1, 1))
        print(json.dumps({k: rec[k] for k in ('rank_full_exact', 'nullity_exact', 'dim_piW_ker_exact')}), flush=True)
        if rW == 0:
            rec['result'] = 'NO polar denominator of degree %d (sector %s): pi_W(ker) = 0 exactly' % (d, par)
        else:
            # basis of pi_W(ker) and the annihilating quadrics in w
            WT = Wm.transpose()            # nul x nW
            Wr, rk = WT.rref()
            # annihilator: left nullspace of basis -> use nullspace of Wr (rk x nW)^T? we want linear forms lam with lam . v = 0 for all v in rowspace(Wr)
            Bm = fmpz_mat(rk, nW, [int(Wr[i, j]) for i in range(rk) for j in range(nW)])
            Lam, nlam = Bm.nullspace()      # nW x nlam : columns are lam with Bm lam = 0
            quadrics = []
            Wkeys = [(a, b) for a in range(nw) for b in range(a, nw)]
            for k in range(nlam):
                terms = []
                for idx, (a, b) in enumerate(Wkeys):
                    c = int(Lam[idx, k])
                    if c:
                        terms.append(f'{c}*w{a}*w{b}')
                quadrics.append('+'.join(terms).replace('+-', '-'))
            rec.update(dim_piW=rk, n_quadrics=nlam)
            qpath = os.path.join(w2x.LANE, f'wave2-excl-relax-{JOB}-quadrics.ms')
            with open(qpath, 'w') as f:
                f.write(','.join(f'w{a}' for a in range(nw)) + '\n')
                f.write(f'{P}\n' if False else '0\n')
                f.write(',\n'.join(quadrics) + '\n')
            rec['quadrics_file'] = qpath
            rec['w_basis'] = [f'{i}:{" ".join(map(str, m))}' for (i, m) in wb]
            # a basis of pi_W(ker) as symmetric matrices (for inspection)
            rec['piW_basis_sample'] = [{f'{Wkeys[j]}': int(Wr[i, j]) for j in range(nW) if int(Wr[i, j])} for i in range(min(rk, 5))]
            rec['result'] = f'pi_W(ker) has dim {rk}; candidate w must satisfy {nlam} quadrics (file written)'
            if MSOLVE and nlam > 0:
                try:
                    opath = qpath.replace('.ms', '.out')
                    r = subprocess.run(['msolve', '-v', '1', '-t', '8', '-f', qpath, '-o', opath], capture_output=True, text=True, timeout=7200)
                    rec['msolve_stdout_tail'] = r.stdout[-600:]
                    rec['msolve_stderr_tail'] = r.stderr[-400:]
                    out = open(opath).read() if os.path.exists(opath) else ''
                    rec['msolve_out_head'] = out[:600]
                    rec['msolve_dim_positive'] = ('positive' in r.stdout.lower() or 'positive' in r.stderr.lower())
                except Exception:
                    rec['msolve_error'] = traceback.format_exc()[-400:]
    rec['seconds'] = round(time.time() - t0, 1)
    logline(LOG, dict(rec, stage='final'))
    print(json.dumps({k: v for k, v in rec.items() if k not in ('w_basis', 'piW_basis_sample')}), flush=True)


if __name__ == '__main__':
    main()
