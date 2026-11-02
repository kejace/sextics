"""Extract the unique (mod J_8^+) degree-8 Y-even E1-kernel element over Q via CRT + rational reconstruction."""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv
from flint import nmod_mat
from fractions import Fraction
F, A, q, rec = core.load_instance('diagonal-permutation')
d = 8
Dmons = monomials(d, 0)
CA = len(Dmons)

def kernel_rep_modp(p):
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    jcols = []
    for i in range(6):
        par_needed = (1 - (0 if i < 3 else 1)) % 2
        for m in monomials(d + 1, par_needed):
            jcols.append(shift(Fi[i], m, 1, p))
    acols = [shift(Ap, m, 1, p) for m in Dmons]
    rows = {}
    for col in jcols + acols:
        for e in col: rows.setdefault(e, len(rows))
    R = len(rows); CJ = len(jcols); C = CJ + CA
    flat = [0] * (R * C)
    for j, col in enumerate(jcols + acols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    # K basis: rows = projections of nullspace vectors on D coords
    K = nmod_mat(nullity, CA, [int(X[CJ + j, k]) for k in range(nullity) for j in range(CA)], p)
    # J_8^+ basis in the monomial coordinates of Dmons
    idx = {m: j for j, m in enumerate(Dmons)}
    jd = []
    for i in range(6):
        par_needed = (0 - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, par_needed):
            jd.append(shift(Fi[i], m, 1, p))
    J = nmod_mat(len(jd), CA, [0] * (len(jd) * CA), p)
    for r, col in enumerate(jd):
        for e, v in col.items(): J[r, idx[e]] = v
    Jr, rkJ = J.rref()
    # reduce K rows modulo row space of J (eliminate pivot columns of Jr)
    piv = []
    for i in range(rkJ):
        j = 0
        while int(Jr[i, j]) == 0: j += 1
        piv.append(j)
    Kr, rkK = K.rref()
    # Stack Jr pivot rows and Kr rows; rref; find rows not in J-space
    S = nmod_mat(rkJ + rkK, CA, [0] * ((rkJ + rkK) * CA), p)
    for i in range(rkJ):
        for j in range(CA): S[i, j] = Jr[i, j]
    for i in range(rkK):
        for j in range(CA): S[rkJ + i, j] = Kr[i, j]
    Sr, rkS = S.rref()
    # the rref rows whose pivot is not a J-pivot give the complement representatives (canonical, reduced mod J)
    out = []
    for i in range(rkS):
        j = 0
        while int(Sr[i, j]) == 0: j += 1
        if j not in piv:
            vec = [int(Sr[i, c]) for c in range(CA)]
            # fully reduce modulo J rows: Sr is rref so pivot columns of other rows are zero already
            out.append((j, vec))
    return rkJ, rkK, rkS, out

res = {}
primes = core.PRIMES[:3]
for p in primes:
    t0 = time.time()
    rkJ, rkK, rkS, out = kernel_rep_modp(p)
    print(f'p={p}: rkJ={rkJ} rkK={rkK} rkS={rkS} complement reps={len(out)} pivots={[o[0] for o in out]} ({time.time()-t0:.1f}s)', flush=True)
    res[p] = out
# CRT + ratrecon on the first complement rep
vecs = [res[p][0][1] for p in primes]
acc, m = vecs[0], primes[0]
for v, p in zip(vecs[1:], primes[1:]):
    acc = [core.crt_pair(a, m, b, p)[0] for a, b in zip(acc, v)]
    m *= p
rat = [core.ratrecon(a, m) for a in acc]
assert all(r is not None for r in rat), 'ratrecon failed'
D8 = {Dmons[j]: r for j, r in enumerate(rat) if r != 0}
s = core.to_str(D8)
print('D8 terms:', len(D8))
print(s)
# verify exactly over Q: A*D8 in J? use e1_solvable (mod p) and also check with different prime
from hodge import fixed_d
for p in core.PRIMES[3:5]:
    ok, info = fixed_d.e1_solvable(F, A, q, D8, p=p)
    print('e1 check prime', p, ok, info)
json.dump({'d': 8, 'D': s, 'terms': len(D8)}, open('/tmp/hodge/results/diagonal-permutation/D8_e1kernel.json', 'w'), indent=1)
