"""Rigorous completion of the degree-d polar relaxation (default d = 6) in staged, light steps.

  python3 wave2_excl_relax6.py <par> modp <primeindex>   -> nullspace of the relaxation matrix mod p (saved .npy)
  python3 wave2_excl_relax6.py <par> recon <pi1,pi2,...> -> CRT + rational reconstruction + EXACT verification;
                                                           writes exact pi_W basis and the annihilating quadrics (msolve file)
  python3 wave2_excl_relax6.py <par> msolve              -> candidate enumeration (homogeneous dimension; affine charts)
Rigor: the reconstructed vectors are verified exactly (fmpz arithmetic) to lie in the rational kernel; their number
equals the modular nullity, hence rank_Q = rank_p for both the full matrix and the VQ block, and pi_W(ker_Q) is exactly
the span of the verified W-parts.  Any admissible w satisfies w w^T in pi_W(ker_Q).
"""
import sys, os, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, logline
_argv = sys.argv; sys.argv = ['x', 'imported', '6', sys.argv[1]]
import wave2_excl_relax as R
sys.argv = _argv
from fractions import Fraction
from math import gcd
import numpy as np
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq, fmpz

PAR = sys.argv[1]; MODE = sys.argv[2]
D_DEG = int(os.environ.get('RELAX_D', '6'))
LANE = w2x.LANE
TAG = f'relax{D_DEG}{"e" if PAR == "0" else ("o" if PAR == "1" else "a")}'
LOG = os.path.join(LANE, f'wave2-excl.{TAG}cert.jsonl')


def build_cols(p):
    wb, cb, qb, cols, rows = R.build(D_DEG, PAR, p)
    return wb, cb, qb, cols, rows


def stage_modp(pi):
    p = core.PRIMES[pi] if pi >= 0 else {-1: 1073741789, -2: 1000000007, -3: 536870909}[pi]
    t0 = time.time()
    wb, cb, qb, cols, rows = build_cols(p)
    Rn, C = len(rows), len(cols)
    flat = [0] * (Rn * C)
    for j, (_, _, P_) in enumerate(cols):
        for e, v in P_.items():
            flat[rows[e] * C + j] = int(v) % p
    M = nmod_mat(Rn, C, flat, p)
    del flat
    X, nul = M.nullspace()
    del M
    arr = np.zeros((C, nul), dtype=np.int64)
    for j in range(nul):
        for i in range(C):
            arr[i, j] = int(X[i, j])
    np.save(os.path.join(LANE, f'wave2-excl-{TAG}-null-p{p}.npy'), arr)
    # canonical structure: pivot columns are those rows i of X that are standard basis... record which coordinates are 'free'
    rec = dict(exp='relax6cert', stage='modp', d=D_DEG, parity=PAR, prime=p, rows=Rn, cols=C, nullity=nul, rank=C - nul,
               nW=len(wb) * (len(wb) + 1) // 2, seconds=round(time.time() - t0, 1))
    logline(LOG, rec); print(json.dumps(rec), flush=True)


def stage_recon(pis):
    t0 = time.time()
    primes = [core.PRIMES[i] for i in pis]
    arrs = [np.load(os.path.join(LANE, f'wave2-excl-{TAG}-null-p{p}.npy')) for p in primes]
    C, nul = arrs[0].shape
    assert all(a.shape == (C, nul) for a in arrs), 'nullities differ between primes'
    # identify the free (identity) rows: python-flint nullspace gives X with X[free_j, k] = 1 for k-th free column
    # CRT entrywise
    mod = 1
    acc = [0] * (C * nul)
    for a, p in zip(arrs, primes):
        flat = a.ravel().tolist()
        if mod == 1:
            acc = flat; mod = p
        else:
            acc = [core.crt_pair(x, mod, y, p)[0] for x, y in zip(acc, flat)]
            mod *= p
    rat = [core.ratrecon(x, mod) for x in acc]
    nfail = sum(1 for r in rat if r is None)
    rec = dict(exp='relax6cert', stage='recon', d=D_DEG, parity=PAR, primes=primes, cols=C, nullity=nul, recon_failures=nfail,
               seconds_recon=round(time.time() - t0, 1))
    print(json.dumps(rec), flush=True)
    if nfail:
        rec['result'] = 'reconstruction failed for %d entries: add primes' % nfail
        logline(LOG, rec); return
    maxh = max(max(abs(r.numerator), r.denominator) for r in rat)
    rec['max_height'] = maxh
    # exact matrix (scaled columns to integers)
    t1 = time.time()
    wb, cb, qb, cols, rows = build_cols(None)
    Rn = len(rows)
    assert len(cols) == C
    scal = []
    flatZ = [0] * (Rn * C)
    for j, (_, _, P_) in enumerate(cols):
        den = 1
        for v in P_.values():
            den = den * v.denominator // gcd(den, v.denominator)
        scal.append(den)
        for e, v in P_.items():
            flatZ[rows[e] * C + j] = int(v * den)
    MZ = fmpz_mat(Rn, C, flatZ)
    del flatZ
    # kernel vectors of the UNscaled M: v (rational). For the scaled matrix MZ = M diag(s): u_j = v_j / s_j.
    Vint = []
    for k in range(nul):
        u = [rat[i * nul + k] / scal[i] for i in range(C)]
        L = 1
        for x in u:
            L = L * x.denominator // gcd(L, x.denominator)
        Vint.append([int(x * L) for x in u])
    VZ = fmpz_mat(C, nul, [Vint[k][i] for i in range(C) for k in range(nul)])
    Prod = MZ * VZ
    zero = all(int(Prod[i, j]) == 0 for i in range(Rn) for j in range(nul)) if Rn * nul < 5_000_000 else (Prod == fmpz_mat(Rn, nul))
    rec.update(exact_verified=bool(zero), seconds_verify=round(time.time() - t1, 1))
    if not zero:
        rec['result'] = 'EXACT VERIFICATION FAILED'
        logline(LOG, rec); print(json.dumps(rec)); return
    # pi_W(ker_Q): W-parts of the verified vectors (rational), exact rank and basis
    nw = len(wb); nW = nw * (nw + 1) // 2
    Wq = fmpq_mat(nul, nW)
    for k in range(nul):
        for a in range(nW):
            r = rat[a * nul + k]
            if r:
                Wq[k, a] = fmpq(r.numerator, r.denominator)
    Wr, rW = Wq.rref()
    Wbasis = [[Fraction(int(Wr[i, a].p), int(Wr[i, a].q)) for a in range(nW)] for i in range(rW)]
    # annihilator: lam with Wbasis . lam = 0  (exact)
    Bq = fmpq_mat(rW, nW)
    for i in range(rW):
        for a in range(nW):
            if Wbasis[i][a]:
                Bq[i, a] = fmpq(Wbasis[i][a].numerator, Wbasis[i][a].denominator)
    # nullspace over Q: use fmpz after scaling rows
    rowsZ = []
    for i in range(rW):
        L = 1
        for x in Wbasis[i]:
            L = L * x.denominator // gcd(L, x.denominator)
        rowsZ.append([int(x * L) for x in Wbasis[i]])
    BZ = fmpz_mat(rW, nW, [x for r_ in rowsZ for x in r_])
    Lam, nlam = BZ.nullspace()
    Wkeys = [(a, b) for a in range(nw) for b in range(a, nw)]
    quadrics = []
    for k in range(nlam):
        col = [int(Lam[idx, k]) for idx in range(nW)]
        g = 0
        for c in col: g = gcd(g, abs(c))
        g = g or 1
        terms = []
        for idx, (a, b) in enumerate(Wkeys):
            c = col[idx] // g
            if c:
                terms.append(f'{c}*w{a}*w{b}')
        quadrics.append('+'.join(terms).replace('+-', '-'))
    qpath = os.path.join(LANE, f'wave2-excl-{TAG}-quadrics.ms')
    with open(qpath, 'w') as f:
        f.write(','.join(f'w{a}' for a in range(nw)) + '\n0\n' + ',\n'.join(quadrics) + '\n')
    json.dump({'d': D_DEG, 'parity': PAR, 'w_basis': [f'{i}:{" ".join(map(str, m))}' for (i, m) in wb], 'Wkeys': Wkeys,
               'piW_basis': [[str(x) for x in row] for row in Wbasis], 'quadrics': quadrics},
              open(os.path.join(LANE, f'wave2-excl-{TAG}-piW.json'), 'w'))
    rec.update(dim_piW_exact=rW, n_quadrics=nlam, quadrics_file=qpath, seconds=round(time.time() - t0, 1),
               result=f'rank certified; pi_W(ker_Q) has dimension {rW}; admissible w satisfy {nlam} exact quadrics')
    logline(LOG, rec); print(json.dumps({k: v for k, v in rec.items()}), flush=True)


def stage_msolve():
    t0 = time.time()
    info = json.load(open(os.path.join(LANE, f'wave2-excl-{TAG}-piW.json')))
    nw = len(info['w_basis']); quadrics = info['quadrics']
    base = os.path.join(LANE, f'wave2-excl-{TAG}-quadrics')
    out = {}
    def run(tag, extra_eqs, timeout=900):
        path = base + f'-{tag}.ms'
        with open(path, 'w') as f:
            f.write(','.join(f'w{a}' for a in range(nw)) + '\n0\n' + ',\n'.join(quadrics + extra_eqs) + '\n')
        opath = path.replace('.ms', '.out')
        try:
            r = subprocess.run(['msolve', '-v', '1', '-t', '2', '-P', '1', '-f', path, '-o', opath], capture_output=True, text=True, timeout=timeout)
            o = open(opath).read() if os.path.exists(opath) else ''
            return {'stdout_tail': r.stdout[-500:], 'stderr_tail': r.stderr[-300:], 'out_head': o[:400], 'out_len': len(o), 'rc': r.returncode}
        except subprocess.TimeoutExpired:
            return {'timeout': timeout}
    # (1) homogeneous system: dimension of the affine cone
    out['homogeneous'] = run('hom', [])
    # (2) random affine chart: sum r_k w_k = 1 (projective points off a hyperplane become affine points)
    import random
    random.seed(3)
    r1 = [random.randint(1, 97) for _ in range(nw)]
    out['chart1'] = run('chart1', ['+'.join(f'{c}*w{a}' for a, c in enumerate(r1)) + '-1'])
    random.seed(4)
    r2 = [random.randint(1, 97) for _ in range(nw)]
    out['chart2'] = run('chart2', ['+'.join(f'{c}*w{a}' for a, c in enumerate(r2)) + '-1'])
    rec = dict(exp='relax6cert', stage='msolve', d=D_DEG, parity=PAR, nw=nw, n_quadrics=len(quadrics), msolve=out, seconds=round(time.time() - t0, 1))
    logline(LOG, rec); print(json.dumps(rec)[:3000], flush=True)


def stage_msolvep(pi):
    PRIME_MAP = {-1: 1073741789, -2: 1000000007, -3: 536870909}
    """Characteristic-p screening: pi_W(ker_p), its annihilating quadrics mod p, msolve over F_p (homogeneous + charts)."""
    import random
    t0 = time.time()
    p = core.PRIMES[pi] if pi >= 0 else PRIME_MAP[pi]
    arr = np.load(os.path.join(LANE, f'wave2-excl-{TAG}-null-p{p}.npy'))
    C, nul = arr.shape
    wb, cb, qb = R.setup(D_DEG, PAR)
    nw = len(wb); nW = nw * (nw + 1) // 2
    Wp = nmod_mat(nul, nW, [int(arr[a, k]) for k in range(nul) for a in range(nW)], p)
    Wr, rW = Wp.rref()
    Bm = nmod_mat(rW, nW, [int(Wr[i, a]) for i in range(rW) for a in range(nW)], p)
    Lam, nlam = Bm.nullspace()
    Wkeys = [(a, b) for a in range(nw) for b in range(a, nw)]
    quadrics = []
    for k in range(nlam):
        terms = []
        for idx, (a, b) in enumerate(Wkeys):
            c = int(Lam[idx, k])
            if c:
                terms.append(f'{c}*w{a}*w{b}')
        quadrics.append('+'.join(terms))
    base = os.path.join(LANE, f'wave2-excl-{TAG}-quadrics-p{p}')
    out = {}
    def run(tag, extra_eqs, timeout=600):
        path = base + f'-{tag}.ms'
        with open(path, 'w') as f:
            f.write(','.join(f'w{a}' for a in range(nw)) + f'\n{p}\n' + ',\n'.join(quadrics + extra_eqs) + '\n')
        opath = path.replace('.ms', '.out')
        try:
            r = subprocess.run(['msolve', '-v', '1', '-t', '2', '-P', '1', '-f', path, '-o', opath], capture_output=True, text=True, timeout=timeout)
            o = open(opath).read() if os.path.exists(opath) else ''
            return {'stdout_tail': r.stdout[-700:], 'stderr_tail': r.stderr[-300:], 'out_head': o[:300], 'out_tail': o[-300:], 'out_len': len(o), 'rc': r.returncode}
        except subprocess.TimeoutExpired:
            return {'timeout': timeout}
    out['homogeneous'] = run('hom', [])
    random.seed(3); r1 = [random.randint(1, 97) for _ in range(nw)]
    out['chart1'] = run('chart1', ['+'.join(f'{c}*w{a}' for a, c in enumerate(r1)) + f'+{p-1}'])
    random.seed(4); r2 = [random.randint(1, 97) for _ in range(nw)]
    out['chart2'] = run('chart2', ['+'.join(f'{c}*w{a}' for a, c in enumerate(r2)) + f'+{p-1}'])
    rec = dict(exp='relax6cert', stage='msolve-modp', d=D_DEG, parity=PAR, prime=p, nw=nw, dim_piW_modp=rW, n_quadrics=nlam,
               msolve=out, seconds=round(time.time() - t0, 1), note='screening in characteristic p')
    logline(LOG, rec); print(json.dumps(rec)[:4000], flush=True)


def stage_leftker(pi):
    """Left kernel of the (V,Q)-block mod p: functionals y on S_{2d}^odd with y.L_V = y.L_Q = 0.  Saved as .npy
    (rows x nullity).  Each such y gives the necessary quadric y.Phi(w) = 0 on admissible w."""
    p = core.PRIMES[pi] if pi >= 0 else {-1: 1073741789, -2: 1000000007, -3: 536870909}[pi]
    t0 = time.time()
    wb, cb, qb, cols, rows = build_cols(p)
    Rn, C = len(rows), len(cols)
    nw = len(wb); nW = nw * (nw + 1) // 2
    # transpose of the VQ block: (C - nW) x Rn
    flat = [0] * ((C - nW) * Rn)
    for j in range(nW, C):
        for e, v in cols[j][2].items():
            flat[(j - nW) * Rn + rows[e]] = int(v) % p
    T = nmod_mat(C - nW, Rn, flat, p)
    del flat
    Y, nul = T.nullspace()
    del T
    arr = np.zeros((Rn, nul), dtype=np.int64)
    for k in range(nul):
        for i in range(Rn):
            arr[i, k] = int(Y[i, k])
    np.save(os.path.join(LANE, f'wave2-excl-{TAG}-leftker-p{p}.npy'), arr)
    # row (monomial) order must be reproducible: save the monomial list
    mons = sorted(rows, key=lambda e: rows[e])
    json.dump([' '.join(map(str, e)) for e in mons], open(os.path.join(LANE, f'wave2-excl-{TAG}-rows.json'), 'w'))
    rec = dict(exp='relax6cert', stage='leftker', d=D_DEG, parity=PAR, prime=p, rows=Rn, cols_VQ=C - nW, left_nullity=nul,
               rank_VQ=Rn - nul, seconds=round(time.time() - t0, 1))
    logline(LOG, rec); print(json.dumps(rec), flush=True)


if __name__ == '__main__':
    if MODE == 'modp':
        stage_modp(int(sys.argv[3]))
    elif MODE == 'recon':
        stage_recon([int(x) for x in sys.argv[3].split(',')])
    elif MODE == 'msolve':
        stage_msolve()
    elif MODE == 'msolvep':
        stage_msolvep(int(sys.argv[3]))
    elif MODE == 'leftker':
        stage_leftker(int(sys.argv[3]))
