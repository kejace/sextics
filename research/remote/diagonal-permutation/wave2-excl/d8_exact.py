"""Exact rational certificate that A * D8 in J(F) over Q: find cofactors U_i (deg 9) with A*D8 = sum_i U_i F_i.
Method: dense nmod_mat rref mod p1 to fix a pivot set; solve the pivot subsystem mod several primes; CRT + rational
reconstruction; EXACT verification over Q with Fractions. Receipt: d8_cofactors.json (U_i over Q)."""
import time, json
from fractions import Fraction as Fr
from flint import nmod_mat
from common import *
from hodge.core import monomials, shift, pderiv, mul, reduce_mod
F, A, q, rec = core.load_instance(IID)
D8 = core.parse(json.load(open(os.path.join(LANE, 'D8_e1kernel.json')))['D'])
assert core.degree(D8) == 8 and core.parity_of(D8) == 0
AD = mul(A, D8)                     # degree 14, Y-odd, exact
t0 = time.time()
# columns: m*F_i with parity(m) + parity(F_i) = 1 (odd target), deg m = 9
Fi = [pderiv(F, i) for i in range(6)]
cols = []   # (i, m)
for i in range(6):
    par = (1 - (0 if i < 3 else 1)) % 2
    for m in monomials(9, par): cols.append((i, m))
rows = monomials(14, 1)
ridx = {e: k for k, e in enumerate(rows)}
R, C = len(rows), len(cols)
print('system', R, 'x', C, flush=True)
def build(p):
    Fip = [reduce_mod(f, p) for f in Fi]
    flat = [0] * (R * C)
    for j, (i, m) in enumerate(cols):
        for e, v in shift(Fip[i], m, 1, p).items():
            flat[ridx[e] * C + j] = v
    b = [0] * R
    for e, v in reduce_mod(AD, p).items(): b[ridx[e]] = v
    return nmod_mat(R, C, flat, p), b
primes = core.PRIMES
p1 = primes[0]
M, b = build(p1)
Mt = M.transpose()   # rref of transpose gives column pivots? no: we need pivot COLUMNS of M -> rref(M)
Rr, rk = M.rref()
pivots = []
j = 0
for i in range(rk):
    while j < C and int(Rr[i, j]) == 0: j += 1
    pivots.append(j); j += 1
print('rank mod p1', rk, 'pivots', len(pivots), time.time() - t0, flush=True)
# consistency check mod p1: rank of [M|b]
Mb = nmod_mat(R, C + 1, [0] * (R * (C + 1)), p1)
# cheaper: solve via the pivot subsystem below and check residual
sols = []
mods = []
for p in primes:
    Mp, bp = build(p)
    sub = nmod_mat(R, rk, [int(Mp[r, c]) for r in range(R) for c in pivots], p)
    # solve least-squares-free: pick rk independent rows via rref of transpose
    subT = sub.transpose()
    RrT, rkT = subT.rref()
    if rkT != rk:
        print('pivot set rank-deficient mod', p, 'skip', flush=True); continue
    rpiv = []; jj = 0
    for i in range(rkT):
        while jj < R and int(RrT[i, jj]) == 0: jj += 1
        rpiv.append(jj); jj += 1
    Sq = nmod_mat(rk, rk, [int(Mp[r, c]) for r in rpiv for c in pivots], p)
    bb = nmod_mat(rk, 1, [bp[r] for r in rpiv], p)
    x = Sq.solve(bb)
    # verify the full residual mod p
    xs = [int(x[i, 0]) for i in range(rk)]
    res = [0] * R
    for k_, c in enumerate(pivots):
        if xs[k_]:
            i, m = cols[c]
            for e, v in shift(reduce_mod(Fi[i], p), m, xs[k_], p).items():
                res[ridx[e]] = (res[ridx[e]] + v) % p
    bad = sum(1 for r in range(R) if res[r] != bp[r])
    print('prime', p, 'residual mismatches', bad, time.time() - t0, flush=True)
    if bad: 
        log({'exp': 'd8-cofactor-modp', 'p': p, 'consistent': False}); continue
    sols.append(xs); mods.append(p)
    # try rational reconstruction with the primes so far
    Mprod = 1; cur = [0] * rk
    for xs_, p_ in zip(sols, mods):
        if Mprod == 1: cur = list(xs_); Mprod = p_
        else:
            cur = [core.crt_pair(a, Mprod, c_, p_) for a, c_ in zip(cur, xs_)]; Mprod *= p_
    try:
        rat = [core.ratrecon(a, Mprod) for a in cur]
    except Exception as ex:
        print('ratrecon failed with', len(mods), 'primes', ex, flush=True); continue
    if any(r is None for r in rat):
        print('ratrecon incomplete with', len(mods), 'primes', flush=True); continue
    # exact verification over Q
    U = [dict() for _ in range(6)]
    for k_, c in enumerate(pivots):
        if rat[k_]:
            i, m = cols[c]; U[i][m] = Fr(rat[k_])
    total = {}
    for i in range(6):
        for e, v in mul(U[i], Fi[i]).items(): total[e] = total.get(e, 0) + v
    diff = {e: Fr(total.get(e, 0)) - Fr(AD.get(e, 0)) for e in set(total) | set(AD)}
    ok = all(v == 0 for v in diff.values())
    print('EXACT verification with', len(mods), 'primes:', ok, time.time() - t0, flush=True)
    if ok:
        receipt('d8_cofactors.json', {'statement': 'A*D8 = sum_i U_i * dF/dz_i exactly over Q', 'D8': core.to_str(D8),
                 'U': [core.to_str(u) for u in U], 'primes_used': mods, 'rank_mod_p1': rk})
        log({'exp': 'd8-in-J-exact-over-Q', 'ok': True, 'primes_used': mods, 'rank': rk, 'rows': R, 'cols': C,
             'max_height_digits': max(len(str(abs(Fr(r).numerator))) + len(str(Fr(r).denominator)) for r in rat),
             'seconds': round(time.time() - t0, 1)})
        break
else:
    log({'exp': 'd8-in-J-exact-over-Q', 'ok': False, 'note': 'lift did not verify with all primes', 'seconds': round(time.time() - t0, 1)})
