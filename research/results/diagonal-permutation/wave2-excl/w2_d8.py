"""JOB excl-dp-d8: the unique non-Jacobian degree-8 E1-kernel line D8 (wave 1), made rigorous and tested.

1. exact rational u8 with u8(F) = A*D8 (CRT over primes + rational reconstruction, verified with Fractions)
   -> proves A*D8 in J(F) over Q, hence dim ker[A] on R_8^+ is exactly 1 (wave-1: = 1 mod p, i.e. <= 1 over Q).
2. relaxation: Psi(D8;u8) = D8 div u8 - u8(D8)  in J_16 ?  (mod P1, P2 with the exact u8; good reduction checked)
3. reduced E2' system at P1 and P2 (consistent? ranks); rank certification: lift the mod-p nullspace to Q.
4. engine cross-checks at P1: full fixed-D system (compressed) and served target space.
Outputs: d8_<job>.json, log lines exp in {'d8_e1_exact','d8_relax','d8_reducedE2','d8_engine','d8_rankcert'}."""
import sys, os, time, json
import numpy as np
from fractions import Fraction
import w2lib as W
from hodge import core, fixed_d
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree, shift, NV

job = os.environ.get('W2JOB', 'd8')
OUT = f'{W.PKT}/d8_{job}.json'
res = {}


def save():
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


D8 = core.parse(json.load(open(f'{W.LANE}/D8_e1kernel.json'))['D'])
assert degree(D8) == 8 and core.parity_of(D8) == 0
d = 8
F, A = W.F, W.A
Fi = [pderiv(F, i) for i in range(NV)]

# ---------------------------------------------------------------- 1. exact u8
t0 = time.time()
primes = [W.P1, W.P2, W.P3] + [p for p in core.PRIMES[1:] if p not in (W.P1, W.P2, W.P3)]


def e1_rref(D, p):
    """rref of the E1 system [cols | A D] mod p; returns (cols, red, piv)."""
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Dp = reduce_mod(D, p)
    Fip = [pderiv(Fp, i, p) for i in range(NV)]
    cols, colpolys = [], []
    for i in range(NV):
        par = (1 + 0 - (0 if i < 3 else 1)) % 2
        for m in monomials(d + 1, par):
            cols.append((i, m)); colpolys.append(shift(Fip[i], m, 1, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for P in colpolys + [rhs]:
        for e in P:
            rows.setdefault(e, len(rows))
    C = len(cols)
    M = np.zeros((len(rows), C + 1), dtype=np.int64)
    for j, P in enumerate(colpolys):
        for e, v in P.items():
            M[rows[e], j] = v
    for e, v in rhs.items():
        M[rows[e], C] = v
    red, piv = W.rref_np(M, p)
    return cols, red, piv


sols = []; used = []; piv0 = None
lifted = None; u8 = None
for p in primes[:4]:
    cols, red, piv = e1_rref(D8, p)
    C = len(cols)
    assert piv[-1] != C, f'E1 inconsistent mod {p}?!'
    if piv0 is None:
        piv0 = piv
    elif piv != piv0:
        print('pivot structure differs at', p, '-- skipping prime', flush=True)
        continue
    x = np.zeros(C, dtype=np.int64)
    for r, j in enumerate(piv):
        x[j] = red[r, C]
    sols.append(x); used.append(p)
    lifted = W.lift_vector(sols, used)
    print(f'prime {p}: rank {len(piv)} of {C} cols; lift {"ok" if lifted else "not yet"} ({time.time()-t0:.0f}s)', flush=True)
    if lifted is not None:
        u8 = [{} for _ in range(NV)]
        for (i, m), v in zip(cols, lifted):
            if v:
                u8[i][m] = v
        lhs = W.exact_apply_field(u8, F)
        ok = (sub(lhs, mul(A, D8)) == {})
        print('exact verification u8(F) == A*D8 :', ok, flush=True)
        if ok:
            break
        lifted = None
if lifted is not None:
    maxh = max(max(abs(v.numerator), v.denominator) for u in u8 for v in u.values())
    res['u8'] = {'primes': used, 'terms': sum(len(u) for u in u8), 'max_height_bits': maxh.bit_length(), 'seconds': round(time.time() - t0, 1)}
    json.dump({'u8': [core.to_str(u) for u in u8]}, open(f'{W.PKT}/d8_u8_exact.json', 'w'))
    W.log('d8_e1_exact', D='D8', d=8, prime=used, result='A*D8 in J(F) over Q (exact u8 verified)', seconds=res['u8']['seconds'],
          note=f'u8 terms {res["u8"]["terms"]}, height {maxh.bit_length()} bits; with ker_p = 1 this gives dim ker[A] on R_8^+ = 1 over Q')
else:
    res['u8'] = {'primes': used, 'lift': 'failed', 'seconds': round(time.time() - t0, 1)}
    W.log('d8_e1_exact', D='D8', d=8, prime=used, result='rational lift of u8 FAILED with %d primes (heights too large); continuing mod p' % len(used),
          seconds=res['u8']['seconds'], note='E1 holds mod every prime tried (same pivot structure); exact statement pending')
save()

# ---------------------------------------------------------------- 2. relaxation Psi(D8;u8) in J_16 ?
def u8_mod(p):
    if u8 is not None:
        return [reduce_mod(u, p) for u in u8]
    return W.solve_E1(D8, p)

for p in (W.P1, W.P2):
    t0 = time.time()
    red16 = W.Red(16, p, 1, None, verbose=True)
    good = (red16.ring_dim == W.koszul_hilbert(16, 1))
    up = u8_mod(p)
    Pp = W.psi(reduce_mod(D8, p), up, p)
    assert core.parity_of(Pp) == 1 and degree(Pp) == 16
    inJ = red16.in_J(Pp)
    nfv = red16.nf(Pp)
    res[f'relax_p{p}'] = {'in_J': bool(inJ), 'nf_nonzero_coords': int(np.count_nonzero(nfv)), 'ring_dim_16_odd': red16.ring_dim,
                          'good_reduction': good, 'exact_u8': u8 is not None, 'seconds': round(time.time() - t0, 1)}
    W.log('d8_relax', D='D8', d=8, prime=p, result=('Psi in J' if inJ else 'Psi NOT in J (relaxation excludes D8)'),
          rank=red16.rank, seconds=res[f'relax_p{p}']['seconds'],
          note=f'good reduction deg16 odd: {good}; nf nonzero coords {int(np.count_nonzero(nfv))}; u8 exact: {u8 is not None}')
    print('relaxation at', p, res[f'relax_p{p}'], flush=True)
    save()

# ---------------------------------------------------------------- 3. reduced E2' at P1, P2 + rank certification
null_bases = []; null_primes = []; free0 = None
for p in (W.P1, W.P2, W.P3):
    t0 = time.time()
    up = u8_mod(p)
    S = W.ReducedE2(D8, p, u=up, verbose=True)
    r = S.solve(compressed=True, verbose=True)
    res[f'reducedE2_p{p}'] = r
    W.log('d8_reducedE2', D='D8', d=8, prime=p, result=('consistent' if r['consistent'] else 'inconsistent'), rank=r['rank'],
          nullity=r['nullity'], rows=r['rows'], cols=r['cols'], seconds=r['seconds'], note='reduced E2 system (c_ij, Q)')
    print('reducedE2', p, r, flush=True)
    save()
    if r['consistent']:
        continue
    # nullspace of M (without rhs) mod p, full (uncompressed) for certification
    M0 = np.ascontiguousarray(S.M[:, :S.ncols])
    t1 = time.time()
    X, nullity = W.to_nmod(M0, p).nullspace()
    Xa = W.from_nmod(X)[:, :nullity]
    print(f'full nullspace mod {p}: nullity {nullity} ({time.time()-t1:.0f}s)', flush=True)
    res[f'reducedE2_p{p}']['full_nullity'] = int(nullity)
    # identify free variables (rows of Xa that form an identity block)
    null_bases.append(Xa); null_primes.append(p)
    save()

# exact column polynomials for verification
Di = [pderiv(D8, i) for i in range(NV)]
excols, expolys = [], []
for i in range(NV):
    for j in range(i + 1, NV):
        base = sub(mul(Fi[i], Di[j]), mul(Fi[j], Di[i]))
        DFj = mul(D8, Fi[j]); DFi = mul(D8, Fi[i])
        for m in monomials(d - 4):
            P = shift(base, m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                P = add(P, shift(DFj, tuple(m2), m[i]))
            if m[j]:
                m2 = list(m); m2[j] -= 1
                P = sub(P, shift(DFi, tuple(m2), m[j]))
            if P:
                excols.append(('c', (i, j), m)); expolys.append(P)
for m in monomials(2 * d - 6):
    excols.append(('Q', None, m)); expolys.append(shift(F, m, -1))
eqpar = 1
keep = [k for k, P in enumerate(expolys) if core.parity_of(P) == eqpar]
excols = [excols[k] for k in keep]; expolys = [expolys[k] for k in keep]
assert excols == [c for c in S.cols], 'exact/modular column lists differ'

if len(null_bases) >= 2 and all(b.shape[1] == null_bases[0].shape[1] for b in null_bases):
    t0 = time.time()
    nb = null_bases[0].shape[1]
    certified = 0
    for v in range(nb):
        vecs = [b[:, v] for b in null_bases]
        lifted = W.lift_vector(vecs, null_primes)
        if lifted is None:
            print('nullvector', v, 'lift failed with', len(null_primes), 'primes', flush=True)
            break
        acc = {}
        for k, c in enumerate(lifted):
            if c:
                acc = add(acc, scale(expolys[k], c))
        if acc == {}:
            certified += 1
        else:
            print('nullvector', v, 'does NOT verify exactly', flush=True)
            break
    res['rankcert'] = {'nullity_modp': nb, 'certified_exact_nullvectors': certified, 'primes': null_primes, 'seconds': round(time.time() - t0, 1)}
    ok = (certified == nb)
    W.log('d8_rankcert', D='D8', d=8, prime=null_primes, result=('rank_Q = rank_p certified' if ok else 'certification incomplete'),
          nullity=nb, seconds=res['rankcert']['seconds'],
          note='exact rational nullvectors of the reduced E2 matrix; with inconsistency mod p this proves inconsistency over Q')
    print('rank certification', res['rankcert'], flush=True)
    save()

# ---------------------------------------------------------------- 4. engine cross-checks at P1
t0 = time.time()
out = fixed_d.solve_fixed_D_fast(F, A, 1, D8, p=W.P1, verbose=True, want_exact=False)
res['engine_full_P1'] = {k: out[k] for k in ('consistent', 'rank', 'nullity', 'rows', 'cols')}
W.log('d8_engine', D='D8', d=8, prime=W.P1, result=('consistent' if out['consistent'] else 'inconsistent'), rank=out['rank'],
      nullity=out['nullity'], seconds=round(time.time() - t0, 1), note='engine v2 full (N,Q) system, compressed')
save()
t0 = time.time()
dimA, basis, contains = W.served_space(D8, W.P1, verbose=True)
res['served_P1'] = {'dim': dimA, 'contains_A': contains}
W.log('d8_served', D='D8', d=8, prime=W.P1, result=f'served_dim={dimA} containsA={contains}', seconds=round(time.time() - t0, 1),
      note='baseline 18 = dim J_6^- (zero class)')
save()
print('DONE', res, flush=True)
