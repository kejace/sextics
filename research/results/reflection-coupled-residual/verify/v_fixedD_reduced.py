"""Reduced fixed-D certificate test (independent code), valid at q=1 for Y-even D with gcd(D,F)=1 and A*D in J:
all E1 solutions are N = N0 + a.F (Koszul, a antisymmetric in S_{d-4}); a certificate exists iff
   rho(N0) + rho(a.F) = F*Q   for some a, Q in S_{2d-6}        (E2)
which is a LINEAR system in (a, Q) over S_{2d}^{cT}.  Reads N0/rho0 from ob_<tag>_p<p>.json (v_ob.py).
Reports rank / augmented rank mod p and the nullity, and attempts rational reconstruction of a kernel basis from a
second prime only if asked.  'inconsistent mod p' transfers to Q iff rank_Q = rank_p (certified when nullity = 0,
or when the mod-p kernel lifts exactly).   Usage: python3 v_fixedD_reduced.py <p> [tag]"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]); tag = sys.argv[2] if len(sys.argv) > 2 else 'K8'
Fs, As, qv = load_rec(); F = modp(parse(Fs), p); A = modp(parse(As), p); cA = char_of(A)
if tag == 'K8':
    DQ = parse(json.load(open('/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json'))[0])
else:
    DQ = parse(tag)
D = modp(DQ, p); d = max(sum(e) for e in D); cD = char_of(DQ); cT = chi_mul(cA, cD)
Fi = [deriv(F, i, p) for i in range(6)]; cFi = [char_of(f) for f in Fi]; Di = [deriv(D, i, p) for i in range(6)]
data = json.load(open(f'ob_{tag}_p{p}.json')); assert data['p'] == p
rho0 = {tuple(int(x) for x in k.split(',')): v for k, v in data['rho0'].items()}
t0 = time.time()
def rho_of(N):
    divN = {}
    for i in range(6): divN = add(divN, deriv(N[i], i, p), p)
    ND = {}
    for i in range(6): ND = add(ND, mul(N[i], Di[i], p), p)
    return add(mul(D, divN, p), scale(ND, p - 1, p), p)
# unknowns: a_ij (i<j) in S_{d-4} with char cT*chi_i*chi_j ; Q in S_{2d-6} with char cT
cols = []
for i in range(6):
    for j in range(i + 1, 6):
        for m in monos_chi(d - 4, chi_mul(cT, chi_mul(cFi[i], cFi[j]))):
            cols.append(('a', i, j, m))
na = len(cols)
for m in monos_chi(2 * d - 6, cT): cols.append(('Q', m))
nQ = len(cols) - na
mons2 = monos_chi(2 * d, cT); idx2 = {m: j for j, m in enumerate(mons2)}
R, C = len(mons2), len(cols)
print(f'[{tag}] reduced E2 system: {R} equations, {na} Koszul unknowns + {nQ} Q unknowns', flush=True)
M = nmod_mat(R, C + 1, p)
for jcol, col in enumerate(cols):
    if col[0] == 'a':
        _, i, j, m = col
        N = [dict() for _ in range(6)]
        N[i] = shift(Fi[j], m); N[j] = scale(shift(Fi[i], m), p - 1, p)   # a_ij F_j e_i - a_ij F_i e_j
        v = rho_of(N)
    else:
        v = scale(shift(F, col[1]), p - 1, p)       # -F*Q
    for e, c in v.items(): M[idx2[e], jcol] = c
for e, c in rho0.items(): M[idx2[e], C] = (p - c) % p    # move rho0 to RHS: sum = -rho0
t1 = time.time()
M0 = nmod_mat(R, C, p)
for r in range(R):
    for c in range(C):
        v = M[r, c]
        if v: M0[r, c] = v
rk0 = M0.rank(); rk1 = M.rank()
cons = (rk0 == rk1); nullity = C - rk0
print(f'[{tag}] p={p}: rank {rk0}, augmented rank {rk1}, CONSISTENT mod p = {cons}, nullity {nullity} '
      f'({"rigorous over Q since nullity 0" if (not cons and nullity == 0) else ("transfer needs rank_Q = rank_p" if not cons else "a certificate exists mod p -> lift!")}) ({time.time()-t1:.0f}s)', flush=True)
out = {'check': 'fixedD-reduced', 'D': tag, 'd': d, 'p': p, 'rows': R, 'na': na, 'nQ': nQ, 'rank': rk0, 'rank_aug': rk1, 'consistent_modp': cons, 'nullity': nullity}
logjson(VLOG, out)
if cons:
    Mr, rk = M.rref()
    sol = {}
    for r in range(rk):
        j = 0
        while Mr[r, j] == 0: j += 1
        if j < C and Mr[r, C]: sol[j] = int(Mr[r, C])
    json.dump({'p': p, 'cols': [list(map(str, c)) for c in cols], 'sol': sol}, open(f'fixedD_sol_{tag}_p{p}.json', 'w'))
    print('particular solution saved', flush=True)
elif nullity:
    # save kernel basis for possible rational reconstruction
    K, nul = M0.nullspace()
    ker = [[int(K[i, k]) for i in range(C)] for k in range(nul)]
    json.dump({'p': p, 'cols': [list(map(str, c)) for c in cols], 'kernel': ker}, open(f'fixedD_ker_{tag}_p{p}.json', 'w'))
    print(f'kernel basis ({nul}) saved', flush=True)
