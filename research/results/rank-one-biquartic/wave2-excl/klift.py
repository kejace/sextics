"""klift: multi-prime rational reconstruction of the structured E1-kernel W_d = {b in B_d : A b in J_{d+6}} for the base-type
subspace B_d = span{ P(X) M(Y)^a g(Y)^b : deg = d } (Y-even).  Usage: python3 klift.py <d> <nprimes> <workers>
Output: klift_d<d>_perprime.json, klift_d<d>_rational.json (same format as k7lift)."""
import sys, json, time, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
from sympy import prevprime
d = int(sys.argv[1]); NPR = int(sys.argv[2]) if len(sys.argv) > 2 else 48; WK = int(sys.argv[3]) if len(sys.argv) > 3 else 48
LOGF = f'{PK}/klift_d{d}.log.jsonl'
def xmons(k): return [m for m in monomials(k) if m[3] == m[4] == m[5] == 0]
Bpolys = []; Bnames = []
for b_ in range(0, d // 6 + 1):
    for a_ in range(0, (d - 6 * b_) // 2 + 1):
        k = d - 2 * a_ - 6 * b_
        MG = mul(core.power(M, a_), core.power(g, b_))
        for m in xmons(k):
            Bpolys.append(mul(MG, mono(m))); Bnames.append(f'S{k}.M{a_}.g{b_}:' + core.to_str(mono(m)))
nB = len(Bpolys)
print(f'd={d}: base-type subspace dim {nB}', flush=True)
def per_prime(p):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    rows_mons = monomials(d + 6, 1); idx = {m: k for k, m in enumerate(rows_mons)}
    cols = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (1 - pi) % 2):
            cols.append(shift(Fip[i], m, 1, p))
    nJ = len(cols)
    for b in Bpolys: cols.append(mul(Ap, core.reduce_mod(b, p), p))
    R, C = len(rows_mons), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[idx[e] * C + j] = v
    X, nul = nmod_mat(R, C, flat, p).nullspace()
    proj = nmod_mat(nul, nB, [int(X[nJ + j, k]) for k in range(nul) for j in range(nB)], p)
    Pr, rk = proj.rref()
    basis = [[int(Pr[i, j]) for j in range(nB)] for i in range(rk)]
    colsJ = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d - 5, (0 - pi) % 2): colsJ.append(shift(Fip[i], m, 1, p))
    monsd = monomials(d, 0); idxd = {m: k for k, m in enumerate(monsd)}
    W_ = len(colsJ) + nB
    flatJ = [0] * (len(monsd) * W_)
    for j, c in enumerate(colsJ):
        for e, v in c.items(): flatJ[idxd[e] * W_ + j] = v
    for j, b in enumerate(Bpolys):
        for e, v in core.reduce_mod(b, p).items(): flatJ[idxd[e] * W_ + len(colsJ) + j] = v
    XJ, nulJ = nmod_mat(len(monsd), W_, flatJ, p).nullspace()
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
    json.dump(res, open(f'{PK}/klift_d{d}_perprime.json', 'w'))
    print('per-prime done', round(time.time() - t0), 's; ranks', sorted(set((r['rank_W'], r['rank_JB']) for r in res)), flush=True)
    log(dict(exp='structured-kernel-perprime', d=d, primes=len(res), dimB=nB, rank_W=res[0]['rank_W'], rank_JB=res[0]['rank_JB'], kernel_mod_J=res[0]['rank_W'] - res[0]['rank_JB'],
             seconds=round(time.time() - t0), note='W={b in base-type B_d: A b in J}, J cap B; per-prime rref bases'), path=LOGF)
    from collections import Counter
    def recon(key, rk_key):
        good = [r for r in res if r[rk_key] == res[0][rk_key]]
        pivs = [tuple(next(j for j in range(nB) if row[j]) for row in r[key]) for r in good]
        piv0 = Counter(pivs).most_common(1)[0][0]
        good = [r for r, pv in zip(good, pivs) if pv == piv0]
        vecs = []
        for i in range(res[0][rk_key]):
            vec = []
            for j in range(nB):
                a_, m_ = good[0][key][i][j], good[0]['p']
                for r in good[1:]:
                    a_, m_ = core.crt_pair(a_, m_, r[key][i][j], r['p'])
                rr = core.ratrecon(a_, m_)
                if rr is None: return None, len(good)
                vec.append(rr)
            vecs.append(vec)
        return vecs, len(good)
    W_vecs, nW = recon('basis_W', 'rank_W'); J_vecs, nJg = recon('basis_JB', 'rank_JB')
    print('reconstruction W:', 'ok' if W_vecs else 'FAILED', 'using', nW, 'primes; J cap B:', 'ok' if J_vecs else 'FAILED', flush=True)
    out = dict(coords=Bnames, W=[[str(x) for x in v] for v in W_vecs] if W_vecs else None, JB=[[str(x) for x in v] for v in J_vecs] if J_vecs else None, primes_used=nW)
    json.dump(out, open(f'{PK}/klift_d{d}_rational.json', 'w'))
    hts = max(max(abs(x.numerator), x.denominator) for v in W_vecs for x in v) if W_vecs else 0
    log(dict(exp='structured-kernel-ratrecon', d=d, primes=nW, result='reconstructed' if W_vecs else 'ratrecon failed (need more primes)', max_height_digits=len(str(hts)), seconds=round(time.time() - t0)), path=LOGF)
    print('done', round(time.time() - t0), 's')
