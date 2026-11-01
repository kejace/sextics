"""E1 kernel sweep for Y-ODD D: K_d = {D in S_d^- : A*D in J(F)} modulo J_d^-."""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, mul, shift, pderiv
from flint import nmod_mat
F, A, q, rec = core.load_instance('diagonal-permutation')
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for d in range(1, dmax + 1):
    t0 = time.time()
    Dmons = monomials(d, 1)          # Y-odd degree d
    jcols = []
    for i in range(6):
        par_needed = (0 - (0 if i < 3 else 1)) % 2   # m*F_i must be Y-even
        for m in monomials(d + 1, par_needed):
            jcols.append(shift(Fi[i], m, 1, p))
    acols = [shift(Ap, m, 1, p) for m in Dmons]
    rows = {}
    for col in jcols + acols:
        for e in col:
            rows.setdefault(e, len(rows))
    R = len(rows); CJ = len(jcols); CA = len(acols); C = CJ + CA
    flat = [0] * (R * C)
    for j, col in enumerate(jcols + acols):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    if nullity:
        proj = nmod_mat(CA, nullity, [int(X[CJ + j, k]) for j in range(CA) for k in range(nullity)], p)
        rk = proj.rank()
    else:
        rk = 0
    if d >= 5:
        jd = []
        for i in range(6):
            par_needed = (1 - (0 if i < 3 else 1)) % 2
            for m in monomials(d - 5, par_needed):
                jd.append(shift(Fi[i], m, 1, p))
        rws = {}
        for col in jd:
            for e in col: rws.setdefault(e, len(rws))
        fl = [0] * (len(rws) * len(jd))
        for j, col in enumerate(jd):
            for e, v in col.items(): fl[rws[e] * len(jd) + j] = v
        dimJ = nmod_mat(len(rws), len(jd), fl, p).rank() if jd else 0
    else:
        dimJ = 0
    print(f'[Y-odd D] d={d}: dim S_d^-={CA}, J_{d+6}^+ cols={CJ}, rows={R}, nullity={nullity}, dim K_d={rk}, dim J_d^-={dimJ}, primitive E1 kernel={rk - dimJ}  (build {t1-t0:.1f}s, null {t2-t1:.1f}s)', flush=True)
