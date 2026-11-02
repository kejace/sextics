"""Exact (over Q) basis of the degree-9 E1 kernel modulo J(F) (Y-even and Y-odd), by CRT over several primes of the
canonical rref representatives in standard monomials, plus structure probes: (01)-action on the kernel, membership of
representatives in powers of linear ideals (modulo J), and sparsest representatives.  usage: python3 w2_k9_exact.py [d] [pD]
Writes w2-kernel-d<d>-p<pD>-exact.json and appends to w2-k9exact.jsonl."""
import os, sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hodge import core, jacobian
from hodge.core import monomials, NV, pderiv, shift, mul, add, scale, to_str
from flint import nmod_mat
from fractions import Fraction
LANE = os.path.dirname(os.path.abspath(__file__))
d = int(sys.argv[1]) if len(sys.argv) > 1 else 9
pD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
F, A, q, rec = core.load_instance('diagonal-permutation')
t0 = time.time()
def log(e):
    e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(os.path.join(LANE, 'w2-k9exact.jsonl'), 'a') as fh: fh.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:500], flush=True)
primes = core.PRIMES[:4]
bases = {}
std0 = None
for p in primes:
    K = jacobian.mult_kernel(F, A, d, p=p, parity=pD, verbose=True, return_reducer=True)
    std = K['reducer'].standard
    if std0 is None: std0 = std
    assert std == std0, 'standard monomials differ between primes'
    # canonical: basis is rref over standard monomials -> same pivots expected
    bases[p] = K['basis']
    log({'exp': 'kernel-modp', 'd': d, 'pD': pD, 'p': p, 'dim': K['dim'], 'seconds': round(time.time() - t0, 1)})
dim = len(bases[primes[0]])
sidx = {m: i for i, m in enumerate(std0)}
# CRT per entry
exact = []
for k in range(dim):
    vecs = []
    for p in primes:
        v = [0] * len(std0)
        for m, c in bases[p][k].items(): v[sidx[m]] = c
        vecs.append(v)
    acc, mod = vecs[0], primes[0]
    for v, p in zip(vecs[1:], primes[1:]):
        acc = [core.crt_pair(a, mod, b, p)[0] for a, b in zip(acc, v)]
        mod *= p
    rat = [core.ratrecon(a, mod) for a in acc]
    ok = all(r is not None for r in rat)
    exact.append({std0[j]: r for j, r in enumerate(rat) if r} if ok else None)
nok = sum(1 for e in exact if e is not None)
log({'exp': 'kernel-exact', 'd': d, 'pD': pD, 'dim': dim, 'reconstructed': nok, 'primes': primes, 'seconds': round(time.time() - t0, 1)})
# verify E1 exactly mod a fresh prime for reconstructed reps
from hodge import fixed_d
pv = core.PRIMES[5]
ver = []
for e in exact:
    if e is None: ver.append(None); continue
    ok, info = fixed_d.e1_solvable(F, A, q, e, p=pv)
    ver.append(bool(ok))
log({'exp': 'kernel-exact-verify', 'd': d, 'pD': pD, 'verified_e1_mod_fresh_prime': ver, 'max_height_digits': max((max(len(str(v.numerator)) + len(str(v.denominator)) for v in e.values()) for e in exact if e), default=None)})
json.dump({'d': d, 'pD': pD, 'standard_monomials': [','.join(map(str, m)) for m in std0],
           'basis_exact': [to_str(e) if e else None for e in exact]}, open(os.path.join(LANE, f'w2-kernel-d{d}-p{pD}-exact.json'), 'w'))
# ---- structure probes mod p: (01) action and vanishing ideals (D in I^k + J ?) for the kernel space
p = primes[0]
Fp = core.reduce_mod(F, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
mons = monomials(d, pD); idx = {m: j for j, m in enumerate(mons)}; CA = len(mons)
def vec(P):
    v = [0] * CA
    for e, c in core.reduce_mod(P, p).items(): v[idx[e]] = c
    return v
Jrows = []
for i in range(NV):
    for m in monomials(d - 5, (pD + PAR[i]) % 2) if (PAR := [0,0,0,1,1,1]) else []:
        Jrows.append(vec(shift(Fi[i], m, 1, p)))
rkJ = nmod_mat(len(Jrows), CA, sum(Jrows, []), p).rank()
Krows = [vec(b) for b in bases[p]]
rkJK = nmod_mat(len(Jrows) + len(Krows), CA, sum(Jrows, []) + sum(Krows, []), p).rank()
assert rkJK == rkJ + dim
from w2sym import sig_poly
sK = [vec(sig_poly(b)) for b in bases[p]]
rk_s = nmod_mat(len(Jrows) + len(Krows) + len(sK), CA, sum(Jrows, []) + sum(Krows, []) + sum(sK, []), p).rank()
log({'exp': 'kernel-(01)-stable', 'd': d, 'pD': pD, 'rkJ': rkJ, 'dim': dim, 'stable_under_(01)': rk_s == rkJK})
# (01)-eigenspace dimensions inside the kernel: dim of (K ± sigma K) mod J
for sign, nm in ((1, 'invariant'), (-1, 'anti')):
    rows = [vec(add(b, scale(sig_poly(b), sign))) for b in bases[p]]
    rk = nmod_mat(len(Jrows) + len(rows), CA, sum(Jrows, []) + sum(rows, []), p).rank() - rkJ
    log({'exp': f'kernel-(01)-{nm}-dim', 'd': d, 'pD': pD, 'dim': rk})
# vanishing-ideal probes: dim of (I^k + J) ∩ (K + J) modulo J  = dim(K+J) + dim(I^k+J) - dim(K + I^k + J)
X = [core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
lin = {'X0': X[0], 'X1': X[1], 'X2': X[2], 'Y0': X[3], 'Y1': X[4], 'Y2': X[5],
       'r': add(X[0], scale(X[1], -1)), 's': add(X[3], scale(X[4], -1)), 'p': add(X[0], X[1]), 'q': add(X[3], X[4]),
       'r2': add(X[1], scale(X[2], -1)), 's2': add(X[4], scale(X[5], -1))}
def ideal_rows(names, k):
    rows = []
    forms = [lin[n] for n in names]
    for combo in itertools.combinations_with_replacement(range(len(forms)), k):
        P = {(0,)*6: 1}
        for c in combo: P = mul(P, forms[c], p)
        for m in monomials(d - k, None):
            Q = shift(P, m, 1, p)
            if core.parity_of(Q) == pD: rows.append(vec(Q))
    return rows
tests = [(['X0','X1','X2'], 3), (['X0','X1','X2'], 4), (['Y0','Y1','Y2'], 3), (['Y0','Y1','Y2'], 4), (['X0','Y0'], 2), (['X0','Y0'], 3),
         (['X2','Y2'], 2), (['X2','Y2'], 3), (['X0','X1'], 2), (['Y0','Y1'], 2), (['r','s'], 2), (['r','s'], 3), (['p','q'], 2),
         (['X0','X1','Y0','Y1'], 4), (['X0','X1','X2','Y0'], 4), (['r','s','X2'], 3), (['r','s','Y2'], 3), (['r','r2','s','s2'], 4), (['r','r2','s','s2'], 3),
         (['r','r2'], 2), (['s','s2'], 2), (['X0','Y0','X1'], 3), (['X0','Y0','Y1'], 3)]
for names, k in tests:
    t1 = time.time()
    rows = ideal_rows(names, k)
    base = nmod_mat(len(Jrows) + len(rows), CA, sum(Jrows, []) + sum(rows, []), p).rank()
    both = nmod_mat(len(Jrows) + len(rows) + len(Krows), CA, sum(Jrows, []) + sum(rows, []) + sum(Krows, []), p).rank()
    inter = (base - rkJ) + dim - (both - rkJ)
    log({'exp': 'kernel-vanishing-probe', 'd': d, 'pD': pD, 'ideal': f'({",".join(names)})^{k}', 'dim_IkJ_modJ': base - rkJ,
         'dim_kernel_in_Ik_modJ': inter, 'seconds': round(time.time() - t1, 1)})
print('K9 EXACT DONE', flush=True)
