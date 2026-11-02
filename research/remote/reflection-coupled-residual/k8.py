"""Extract the E1-kernel class at degree d (default 8): Y-even D with A*D in J(F), modulo J_d.
Canonical representative = normal form w.r.t. rref of J_d rows; computed for two primes and rationally reconstructed.
Usage: python3 k8.py <d> [A-expr]"""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv
from flint import nmod_mat
from fractions import Fraction
d = int(sys.argv[1]) if len(sys.argv) > 1 else 8
Aexpr = sys.argv[2] if len(sys.argv) > 2 else 'X0**3*Y0**3'
F, _, q, rec = core.load_instance('reflection-coupled-residual')
A = core.parse(Aexpr)
parA = core.parity_of(A)

def kernel_mod_p(P):
    t0 = time.time()
    Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
    Fi = [pderiv(Fp, i, P) for i in range(NV)]
    parD = 0; parT = (parD + parA) % 2
    monsT = monomials(d + 6, parT); idxT = {m: i for i, m in enumerate(monsT)}
    monsD = monomials(d, parD); idxD = {m: i for i, m in enumerate(monsD)}
    Jrows = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(d + 1, (parT - pi) % 2):
            Jrows.append(shift(Fi[i], m, 1, P))
    Arows = [shift(Ap, m, 1, P) for m in monsD]
    rows = Jrows + Arows
    R, C = len(rows), len(monsT)
    flat = [0] * (R * C)
    for r, poly in enumerate(rows):
        base = r * C
        for e, v in poly.items():
            flat[base + idxT[e]] = v
    M = nmod_mat(R, C, flat, P)
    # left nullspace: nullspace of M^T
    X, nul = M.transpose().nullspace()
    nJ = len(Jrows)
    # project to A-coordinates (D coefficients)
    proj = [[int(X[nJ + j, k]) for k in range(nul)] for j in range(len(monsD))]
    Pm = nmod_mat(len(monsD), nul, [proj[j][k] for j in range(len(monsD)) for k in range(nul)], P).transpose()
    Pr, rk = Pm.rref()
    cand = []
    for i in range(rk):
        poly = {monsD[j]: int(Pr[i, j]) for j in range(len(monsD)) if int(Pr[i, j])}
        cand.append(poly)
    # J_d rows (even part)
    Jd = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        if d - 5 >= 0:
            for m in monomials(d - 5, (parD - pi) % 2):
                Jd.append(shift(Fi[i], m, 1, P))
    RJ, CD = len(Jd), len(monsD)
    flatJ = [0] * (RJ * CD)
    for r, poly in enumerate(Jd):
        for e, v in poly.items(): flatJ[r * CD + idxD[e]] = v
    MJ = nmod_mat(RJ, CD, flatJ, P); MJr, rJ = MJ.rref()
    piv = []
    for i in range(rJ):
        j = 0
        while int(MJr[i, j]) == 0: j += 1
        piv.append(j)
    def nf(poly):
        v = [0] * CD
        for e, c in poly.items(): v[idxD[e]] = c % P
        for i, j in enumerate(piv):
            if v[j]:
                c = v[j]
                for k in range(CD):
                    a = int(MJr[i, k])
                    if a: v[k] = (v[k] - c * a) % P
        return v
    nfs = [nf(c) for c in cand]
    K = nmod_mat(len(nfs), CD, [x for v in nfs for x in v], P)
    Kr, rk2 = K.rref()
    out = []
    for i in range(rk2):
        out.append({monsD[j]: int(Kr[i, j]) for j in range(CD) if int(Kr[i, j])})
    print(f'p={P}: d={d} |S_d^even|={CD} rank J_d={rJ} dim R_d={CD-rJ} | E1-solutions incl J: {rk} | kernel mod J: {rk2} ({time.time()-t0:.0f}s)', flush=True)
    return out

P0, P1 = core.PRIMES[0], core.PRIMES[1]
k0 = kernel_mod_p(P0); k1 = kernel_mod_p(P1)
assert len(k0) == len(k1), (len(k0), len(k1))
res = []
for a, b in zip(k0, k1):
    keys = set(a) | set(b)
    poly = {}
    for e in sorted(keys):
        x, _m = core.crt_pair(a.get(e, 0), P0, b.get(e, 0), P1)
        r = core.ratrecon(x, P0 * P1)
        if r is None:
            print('ratrecon failed at', e); r = None
        poly[e] = r
    res.append(poly)
    print('kernel element:', core.to_str({e: v for e, v in poly.items() if v is not None}))
json.dump([core.to_str({e: v for e, v in poly.items() if v is not None}) for poly in res],
          open(f'/tmp/hodge/results/reflection-coupled-residual/kernel_d{d}_{Aexpr.replace("*","")}.json', 'w'), indent=1)
