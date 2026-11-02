"""Exact (over Q) bases of structural Jacobian-ideal families:
   W(idx,k,d) = { D in J_d^+ : D in (z_idx)^k }  modulo S_{d-6}F, in the ORIGINAL coordinates (idx refers to X0..Y2),
   using fmpq_mat nullspace (small systems).  Saves to families.json."""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, add, scale
from flint import fmpq_mat, fmpq
from fractions import Fraction
F, A, q, rec = core.load_instance('diagonal-permutation')
Fi = [pderiv(F, i) for i in range(6)]

def family(idx, k, d, verbose=True, Mx=None, My=None):
    """If Mx,My given (Z = M z), the ideal is (z_idx)^k in the new coordinates; results are returned in ORIGINAL coords."""
    from jideal import lin_subst, inv3
    global F, Fi
    F0, Fi0 = F, Fi
    if Mx is not None:
        F = lin_subst(F0, Mx, My); Fi = [pderiv(F, i) for i in range(6)]
    try:
        reps = _family(idx, k, d, verbose)
    finally:
        F, Fi = F0, Fi0
    if Mx is not None:
        Mxi, Myi = inv3(Mx), inv3(My)
        reps = [lin_subst(P, Mxi, Myi) for P in reps]
        reps = [{e: v for e, v in P.items() if v} for P in reps]
    return reps

def _family(idx, k, d, verbose=True):
    gens, cols = [], []
    for i in range(6):
        par_needed = (0 - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, par_needed):
            gens.append(('J', i, m)); cols.append(shift(Fi[i], m, 1))
    nJ = len(cols)
    if d >= 6:
        for m in monomials(d - 6, 0):
            gens.append(('F', -1, m)); cols.append(shift(F, m, 1))
    C = len(cols)
    def in_ideal(e): return sum(e[i] for i in idx) >= k
    bad = {}
    for col in cols:
        for e in col:
            if not in_ideal(e): bad.setdefault(e, len(bad))
    M = fmpq_mat(len(bad), C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            if e in bad: M[bad[e], j] = fmpq(v.numerator, v.denominator)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while Mr[i, j] == 0: j += 1
        piv.append(j)
    free = [j for j in range(C) if j not in set(piv)]
    nullity = len(free)
    # D polys for nullspace vectors: x_f = 1, x_{piv_i} = -Mr[i, f]
    Ds = []
    for f in free:
        Dp = dict(cols[f])
        for i, pj in enumerate(piv):
            v = Mr[i, f]
            if v != 0:
                Dp = add(Dp, scale(cols[pj], -Fraction(int(v.p), int(v.q))))
        Dp = {e: v for e, v in Dp.items() if v}
        if Dp: Ds.append(Dp)
    # reduce modulo F-multiples: rref of [F-mults ; Ds] and keep D rows with new pivots
    Fm = [cols[j] for j in range(nJ, C)]
    allm = {}
    for P in Fm + Ds:
        for e in P: allm.setdefault(e, len(allm))
    def mat(polys):
        Mq = fmpq_mat(len(polys), len(allm))
        for r, P in enumerate(polys):
            for e, v in P.items(): Mq[r, allm[e]] = fmpq(v.numerator, v.denominator)
        return Mq
    fp = set()
    if Fm:
        Fr, rkF = mat(Fm).rref()
        for i in range(rkF):
            j = 0
            while Fr[i, j] == 0: j += 1
            fp.add(j)
    reps = []
    if Ds:
        Sr, rk = mat(Fm + Ds).rref()
        inv = {v: e for e, v in allm.items()}
        for i in range(rk):
            j = 0
            while Sr[i, j] == 0: j += 1
            if j not in fp:
                P = {inv[c]: Fraction(int(Sr[i, c].p), int(Sr[i, c].q)) for c in range(len(allm)) if Sr[i, c] != 0}
                reps.append(P)
    if verbose:
        print(f'family idx={idx} k={k} d={d}: nullity={nullity}, nontrivial reps={len(reps)}', flush=True)
    return reps

if __name__ == '__main__':
    out = {}
    for (idx, k, d) in [([0,1,2], 3, 7), ([3,4,5], 4, 8), ([0,1,2], 4, 8)]:
        t0 = time.time()
        reps = family(idx, k, d)
        key = f'idx{idx}_k{k}_d{d}'
        out[key] = [core.to_str(P) for P in reps]
        print(key, len(reps), f'{time.time()-t0:.1f}s', flush=True)
    json.dump(out, open('/tmp/hodge/results/diagonal-permutation/families.json', 'w'), indent=1)
