"""Diagnostic: dimension of the Koszul span Wtot (per E2 sector) versus pi(J) and the X2<=5 quotient, using random family
members (sum of images over r random theta spans Wtot w.h.p.).  usage: python3 wspan.py d p r threads [sector c0,c1]"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from multiprocessing import Pool
d = int(sys.argv[1]); p = int(sys.argv[2]); r = int(sys.argv[3]); threads = int(sys.argv[4])
sector = tuple(int(x) for x in sys.argv[5].split(',')) if len(sys.argv) > 5 else None
S = Setup()
p7inv = pow(7, -1, p); Fz_rest = {e: v % p for e, v in S.Fz.items() if e[2] < 6}
def reduce_F(P):
    """Reduce P (dict, int coeffs) modulo F over F_p: result has X2-exponent <= 5.  (7 F = 7 X2^6 + rest.)"""
    P = {e: v % p for e, v in P.items() if v % p}
    while True:
        hi = sorted((e for e in P if e[2] >= 6), key=lambda e: -e[2])
        if not hi: return P
        for e in hi:
            v = P.pop(e, None)
            if v is None: continue
            c = (v * p7inv) % p
            m = (e[0], e[1], e[2] - 6, e[3], e[4], e[5])
            for e2, v2 in Fz_rest.items():
                k = (m[0]+e2[0], m[1]+e2[1], m[2]+e2[2], m[3]+e2[3], m[4]+e2[4], m[5]+e2[5])
                nv = (P.get(k, 0) - c * v2) % p
                if nv: P[k] = nv
                else: P.pop(k, None)
# family basis mod p: Jacobian part + kernel elements (mod p, canonical)
Dk = []
if sector is None:
    sectors = sorted({S.chi(m) for m in mons(d)})
else:
    sectors = [sector]
for c in sectors:
    Dm, vecs, rJ = kernel_canonical(S, d, c, p)
    for v in vecs:
        Dk.append({Dm[j]: v[j] for j in range(len(Dm)) if v[j]})
for i in range(6):
    for m in mons(d - 5):
        if sector is not None and S.chi(m) != cadd(sector, S.chiE[i]): continue
        Dk.append(imul({m: 1}, S.dF[i]))
print(f'd={d}: family size {len(Dk)} (sector {sector})', flush=True)
if sector is None:
    Sbasis = [(i, j, m) for i in range(6) for j in range(i + 1, 6) for m in mons(d - 4)]
else:
    target = cadd(S.chiA, sector)
    Sbasis = [(i, j, m) for i in range(6) for j in range(i + 1, 6) for m in mons(d - 4) if S.chi(m) == cadd(target, S.chiE[i], S.chiE[j])]
random.seed(7)
thetas = [[random.randrange(1, p) for _ in Dk] for _ in range(r)]
def Dof(theta):
    D = {}
    for t_, Dk_ in zip(theta, Dk):
        for e, v in Dk_.items(): D[e] = (D.get(e, 0) + t_ * v) % p
    return {e: v for e, v in D.items() if v}
Ds = [Dof(th) for th in thetas]
dDs = [[ideriv(D, i) for i in range(6)] for D in Ds]
def kvec(args):
    (i, j, m), q = args
    D = Ds[q]
    divN = core.sub(imul(ideriv({m: 1}, i), S.dF[j]), imul(ideriv({m: 1}, j), S.dF[i]))
    ND = core.sub(imul(imul({m: 1}, S.dF[j]), dDs[q][i]), imul(imul({m: 1}, S.dF[i]), dDs[q][j]))
    return reduce_F(core.sub(imul(D, divN), ND))
t0 = time.time()
with Pool(threads) as pool:
    W = pool.map(kvec, [(s, q) for s in Sbasis for q in range(r)], chunksize=32)
print(f'computed {len(W)} Koszul images ({time.time()-t0:.1f}s)', flush=True)
bych = {}
for w in W:
    comp = {}
    for e, v in w.items(): comp.setdefault(S.chi(e), {})[e] = v
    for c, wc in comp.items(): bych.setdefault(c, []).append(wc)
out = {}
for c, ws in sorted(bych.items()):
    monsc = [m for m in mons(2 * d) if S.chi(m) == c and m[2] <= 5]; idx = {m: i for i, m in enumerate(monsc)}; C = len(monsc)
    flat = [0] * (len(ws) * C)
    for rr, w in enumerate(ws):
        for e, v in w.items(): flat[rr * C + idx[e]] = v
    rkW = nmod_mat(len(ws), C, flat, p).rank()
    # pi(J) in this sector: J generators reduced mod F
    gens = [reduce_F(g) for g in jac_gens(S, 2 * d, c)]
    flatJ = [0] * (len(gens) * C)
    for rr, g in enumerate(gens):
        for e, v in g.items(): flatJ[rr * C + idx[e]] = v
    rkJ = nmod_mat(len(gens), C, flatJ, p).rank() if gens else 0
    out[str(c)] = {'quotient_dim': C, 'rank_W': rkW, 'rank_piJ': rkJ, 'extra_functionals': rkJ - rkW, 'hessian_functionals': C - rkJ}
    print(f'  sector {c}: quotient {C}, Koszul span {rkW}, pi(J) {rkJ}  => Koszul relaxation functionals {C - rkW} (Hessian: {C - rkJ})', flush=True)
tot = {k: sum(v[k] for v in out.values()) for k in ('quotient_dim', 'rank_W', 'rank_piJ', 'extra_functionals', 'hessian_functionals')}
print('TOTAL', tot, flush=True)
json.dump({'d': d, 'p': p, 'r': r, 'sector': sector, 'family_size': len(Dk), 'n_koszul': len(Sbasis), 'sectors': out, 'total': tot}, open(f'{PK}/wspan_d{d}' + (f'_s{sector[0]}{sector[1]}' if sector else '') + '.json', 'w'), indent=1)
log({'exp': 'koszul-span-dimension', 'family': f'E1-admissible degree {d} ' + (f'sector {sector}' if sector else 'all'), 'd': d, 'prime': p,
     'result': f"Koszul span rank {tot['rank_W']} vs pi(J) rank {tot['rank_piJ']} in quotient dim {tot['quotient_dim']}: Koszul relaxation has {tot['quotient_dim']-tot['rank_W']} functionals vs Hessian {tot['hessian_functionals']}",
     'rank': tot['rank_W'], 'seconds': round(time.time() - t0, 1), 'job': f'excl-ctr-wspan{d}', 'note': f'{r} random family members; family size {len(Dk)}'})
