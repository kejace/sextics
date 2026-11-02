"""Referee replay: (a) exact lifted certificates L*D7^+ in J_8, M*D7^- in J_8 (excl_lib.exact_membership, verified with
Fractions); (b) Hilbert values from the series (1+t+..+t^4)^6 and Theta=(1+t^2+..+t^8)^3; (c) A not in J_6 at a fresh prime."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
from hodge import core
from hodge.core import NV, pderiv, mul, reduce_mod, monomials
OUT = '/tmp/hodge/results/short-power-presentation/wave2-excl/referee/referee.jsonl'
def log(**r):
    r['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); r['lens'] = 'referee'
    open(OUT, 'a').write(json.dumps(r) + '\n'); print(json.dumps(r)[:420], flush=True)
F, A, q, rec = W.load()
ex = json.load(open('/tmp/hodge/results/short-power-presentation/wave2-excl/referee/excl_kernel_D7.json'))
Lp = core.parse('X0+2*X1+3*X2'); Mp = core.parse('Y0+3*Y1+5*Y2')
for par, lin, name in (('0', Lp, 'L*D7+'), ('1', Mp, 'M*D7-')):
    D = core.parse(ex[par]['D']); P = mul(lin, D); t0 = time.time()
    g, info = W.exact_membership(F, P, primes=W.PRIMES31, max_primes=6, verbose=True)
    log(exp='referee-exact-membership-J8', what=name, exact=bool(g is not None), info={k: v for k, v in info.items() if k != 'pivots'}, seconds=round(time.time() - t0, 1))
# Hilbert values
def coeffs(poly, n):
    return poly
import itertools
def series_pow(base, k, nmax):
    out = [0] * (nmax + 1); out[0] = 1
    for _ in range(k):
        new = [0] * (nmax + 1)
        for i, a in enumerate(out):
            if a:
                for j, b in enumerate(base):
                    if i + j <= nmax: new[i + j] += a * b
        out = new
    return out
H = series_pow([1, 1, 1, 1, 1], 6, 24); Th = series_pow([1, 0, 1, 0, 1, 0, 1, 0, 1], 3, 24)
vals = {d: ((H[d] + Th[d]) // 2, (H[d] - Th[d]) // 2) for d in (5, 6, 7, 8, 9, 12, 13, 14)}
log(exp='referee-hilbert-values', R_plus_minus=vals, dimR_lib={d: (W.dimR(d, 0), W.dimR(d, 1)) for d in (6, 7, 8, 12, 13, 14)},
    note='R_d^+/- from the tau-equivariant Koszul series; used as the (R1) prime-goodness reference')
# A not in J_6 at a fresh prime
p = 1073741741; Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
g, info = W.membership_mod_p(Fp, Ap, p)
gens = W.jac_gens(Fp, 6, 1, p); mons = monomials(6, 1); idx = {m: i for i, m in enumerate(mons)}
rJ = W.cols_matrix([x[2] for x in gens], idx, p).rank()
log(exp='referee-A-notin-J', prime=p, A_in_J6_modp=bool(info['consistent']), rankJ6minus_modp=rJ, hilbert=W.dimJ(6, 1),
    note='rank_p[J|A] = rank_p J + 1 with rank_p J at the Hilbert value proves A not in J over Q')
