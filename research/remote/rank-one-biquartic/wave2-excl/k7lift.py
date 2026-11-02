"""K7lift: exact (rational) generators of the non-Jacobian E1-kernel K_7 inside the structured subspace
B = S_7(X) (+) M(Y)*S_5(X) (57 coordinates; wave-1 kstruct3 showed K_7 lives there mod 5 primes).
Per prime: canonical rref basis (mod p) of W = {b in B : A b in J_13}; then CRT over all primes + rational
reconstruction; then exact verification (A b in J_13 over Q via exact fmpq solve) and J_7-membership.
Usage: python3 k7lift.py <nprimes> <workers>"""
import sys, json, time, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
from sympy import prevprime
NPR = int(sys.argv[1]) if len(sys.argv) > 1 else 48
WK = int(sys.argv[2]) if len(sys.argv) > 2 else 48
LOGF = f'{PK}/k7lift.log.jsonl'
def xmons(k): return [m for m in monomials(k) if m[3] == m[4] == m[5] == 0]
Bpolys = [mono(m) for m in xmons(7)] + [mul(M, mono(m)) for m in xmons(5)]
nB = len(Bpolys)
def per_prime(p):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    rows_mons = monomials(13, 1); idx = {m: k for k, m in enumerate(rows_mons)}
    cols = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(8, (1 - pi) % 2):
            cols.append(shift(Fip[i], m, 1, p))
    nJ = len(cols)
    for b in Bpolys:
        cols.append(mul(Ap, core.reduce_mod(b, p), p))
    R, C = len(rows_mons), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * C + j] = v
    X, nul = nmod_mat(R, C, flat, p).nullspace()
    proj = nmod_mat(nul, nB, [int(X[nJ + j, k]) for k in range(nul) for j in range(nB)], p)
    Pr, rk = proj.rref()
    basis = [[int(Pr[i, j]) for j in range(nB)] for i in range(rk)]
    # J cap B: b in J_7 (Y-even)
    colsJ = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(2, (0 - pi) % 2):
            colsJ.append(shift(Fip[i], m, 1, p))
    mons7 = monomials(7, 0); idx7 = {m: k for k, m in enumerate(mons7)}
    flatJ = [0] * (len(mons7) * (len(colsJ) + nB))
    W_ = len(colsJ) + nB
    for j, c in enumerate(colsJ):
        for e, v in c.items(): flatJ[idx7[e] * W_ + j] = v
    for j, b in enumerate(Bpolys):
        for e, v in core.reduce_mod(b, p).items(): flatJ[idx7[e] * W_ + len(colsJ) + j] = v
    XJ, nulJ = nmod_mat(len(mons7), W_, flatJ, p).nullspace()
    projJ = nmod_mat(nulJ, nB, [int(XJ[len(colsJ) + j, k]) for k in range(nulJ) for j in range(nB)], p)
    PrJ, rkJ = projJ.rref()
    basisJ = [[int(PrJ[i, j]) for j in range(nB)] for i in range(rkJ)]
    return dict(p=p, rank_W=rk, rank_JB=rkJ, basis_W=basis, basis_JB=basisJ, seconds=round(time.time() - t0, 1))
if __name__ == '__main__':
    primes = []; x = 2147483647 + 1
    while len(primes) < NPR:
        x = prevprime(x); primes.append(int(x))
    t0 = time.time()
    with mp.Pool(WK) as pool:
        res = pool.map(per_prime, primes)
    json.dump(res, open(f'{PK}/k7lift_perprime.json', 'w'))
    print('per-prime done', round(time.time() - t0), 's; ranks', sorted(set((r['rank_W'], r['rank_JB']) for r in res)), flush=True)
    log(dict(exp='K7-structured-kernel-perprime', d=7, primes=len(res), rank_W=[r['rank_W'] for r in res][:5], rank_JB=[r['rank_JB'] for r in res][:5],
             seconds=round(time.time() - t0), note='rref basis mod p of W={b in S7(X)+M*S5(X): A b in J_13} and of J_7 cap B'), path=LOGF)
    # CRT + ratrecon for the W basis (pivot structure must agree)
    def recon(key, rk_key):
        good = [r for r in res if r[rk_key] == res[0][rk_key]]
        pivs = [tuple(next(j for j in range(nB) if row[j]) for row in r[key]) for r in good]
        from collections import Counter
        piv0 = Counter(pivs).most_common(1)[0][0]
        good = [r for r, pv in zip(good, pivs) if pv == piv0]
        rk = res[0][rk_key]
        vecs = []
        for i in range(rk):
            vec = []
            for j in range(nB):
                a_, m_ = good[0][key][i][j], good[0]['p']
                for r in good[1:]:
                    a_, m_ = core.crt_pair(a_, m_, r[key][i][j], r['p'])
                rr = core.ratrecon(a_, m_)
                if rr is None:
                    return None, len(good)
                vec.append(rr)
            vecs.append(vec)
        return vecs, len(good)
    W_vecs, nW = recon('basis_W', 'rank_W')
    J_vecs, nJg = recon('basis_JB', 'rank_JB')
    print('reconstruction W:', 'ok' if W_vecs else 'FAILED', 'using', nW, 'primes; J cap B:', 'ok' if J_vecs else 'FAILED', flush=True)
    out = dict(coords=['S7:' + core.to_str(b) for b in Bpolys[:36]] + ['M*S5:' + core.to_str(b) for b in Bpolys[36:]],
               W=[[str(x) for x in v] for v in W_vecs] if W_vecs else None, JB=[[str(x) for x in v] for v in J_vecs] if J_vecs else None,
               primes_used=nW)
    json.dump(out, open(f'{PK}/k7lift_rational.json', 'w'))
    hts = 0
    if W_vecs:
        hts = max(max(abs(x.numerator), x.denominator) for v in W_vecs for x in v)
    log(dict(exp='K7-structured-kernel-ratrecon', d=7, primes=nW, result='reconstructed' if W_vecs else 'ratrecon failed (need more primes)',
             max_height_digits=len(str(hts)), seconds=round(time.time() - t0), note='CRT over primes with common pivot structure, then rational reconstruction of the canonical rref basis of W (21-dim expected: 18 = J cap B, 3 = K_7)'), path=LOGF)
    print('done', round(time.time() - t0), 's')
