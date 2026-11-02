"""Find D in J(F)_d^+ (Y-even Jacobian ideal elements of degree d) lying in a power of a linear ideal
I = (l_1..l_m)^k, modulo the trivial part S_{d-6}*F.  Linear forms are given via a block-diagonal
coordinate change (3x3 on X, 3x3 on Y) and a list of new-coordinate indices.
Usage as library: find(Mx, My, idx, k, d) -> dict with dims and basis of D's (original coords, exact Q)."""
import sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, mul, add, scale
from flint import nmod_mat, fmpq_mat, fmpq
from fractions import Fraction
F, A, q, rec = core.load_instance('diagonal-permutation')

def lin_subst(P, Mx, My):
    """Substitute Z = M * Znew (i.e. express P in new coordinates): X_i = sum_j Mx[i][j] X'_j, etc."""
    lin = []
    for i in range(3):
        lin.append({tuple(1 if t == j else 0 for t in range(6)): Fraction(Mx[i][j]) for j in range(3) if Mx[i][j]})
    for i in range(3):
        lin.append({tuple(1 if t == 3 + j else 0 for t in range(6)): Fraction(My[i][j]) for j in range(3) if My[i][j]})
    # powers cache
    out = {}
    cache = {}
    def pw(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = core.power(lin[i], k) if k else {(0,)*6: Fraction(1)}
        return cache[(i, k)]
    for e, c in P.items():
        term = {(0,)*6: Fraction(c)}
        for i in range(6):
            if e[i]: term = mul(term, pw(i, e[i]))
        out = add(out, term)
    return {e: v for e, v in out.items() if v}

def inv3(M):
    import sympy
    Mi = sympy.Matrix(M).inv()
    return [[Fraction(int(Mi[i, j].p), int(Mi[i, j].q)) for j in range(3)] for i in range(3)]

def find(Mx, My, idx, k, d, p=None, verbose=True, extra_filter=None):
    """D in J_d^+ with D in (z_idx)^k in the new coordinates z = M^{-1} Z.  Returns basis (orig coords)."""
    p = p or core.PRIMES[0]
    Mxi, Myi = inv3(Mx), inv3(My)
    # F in new coordinates: Z = M z  => F(Z) = F(M z)
    Fn = lin_subst(F, Mx, My)
    Fnp = core.reduce_mod(Fn, p)
    Fi = [pderiv(Fnp, i, p) for i in range(6)]
    gens = []   # (i, m) generators of J_d^+ : m * F_i
    cols = []
    for i in range(6):
        par_needed = (0 - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, par_needed):
            gens.append(('J', i, m)); cols.append(shift(Fi[i], m, 1, p))
    nJ = len(cols)
    # trivial part: S_{d-6}^+ * F
    if d >= 6:
        for m in monomials(d - 6, 0):
            gens.append(('F', -1, m)); cols.append(shift(Fnp, m, 1, p))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R = len(rows); C = len(cols)
    def in_ideal(e):
        return sum(e[i] for i in idx) >= k and (extra_filter is None or extra_filter(e))
    bad_rows = [rows[e] for e in rows if not in_ideal(e)]
    # unknown combination c (C-vector) of J-generators only (first nJ) such that sum c_j col_j has zero coefficient on bad rows,
    # modulo adding F-multiples: we solve for combos of all cols (J gens and F-multiples) with zero on bad rows,
    # then quotient by those supported purely on F-multiples... simpler: solve on J-gens + F gens, report dim and dim of F-part.
    flat = [0] * (len(bad_rows) * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            r = rows[e]
            if not in_ideal(e):
                flat[bad_rows.index(r) * C + j] = v
    # faster index map
    M = nmod_mat(len(bad_rows), C, flat, p)
    X, nullity = M.nullspace()
    # the nullspace vectors: combos of (J gens, F gens) producing D in ideal.  Project to D polynomials; mod out those with zero J-part?
    # D = sum_j c_j col_j ; two combos giving same D differ by a syzygy; we want the space of D's (in the ideal) modulo S_{d-6}F.
    # Compute D-space: evaluate each nullspace vector to a polynomial (mod p), then rank modulo the F-multiples subspace.
    Dpolys = []
    for kk in range(nullity):
        Dp = {}
        for j in range(C):
            v = int(X[j, kk])
            if v: Dp = add(Dp, scale(cols[j], v, p), p)
        Dpolys.append(Dp)
    Fcols = cols[nJ:]
    allm = {}
    for P in Dpolys + Fcols:
        for e in P: allm.setdefault(e, len(allm))
    def mat(polys):
        fl = [0] * (len(polys) * len(allm))
        for r, P in enumerate(polys):
            for e, v in P.items(): fl[r * len(allm) + allm[e]] = v
        return nmod_mat(len(polys), len(allm), fl, p)
    rkD = mat(Dpolys).rank() if Dpolys else 0
    rkF = mat([P for P in Fcols if all(in_ideal(e) for e in P)]).rank() if Fcols else 0
    rkDF = mat(Dpolys + Fcols).rank() if (Dpolys or Fcols) else 0
    # dim of nontrivial D's = rkD - dim(Dspace ∩ Fspace) = rkD - (rkD + rkF_all - rkDF) ; use rkF_all = rank of all F cols
    rkFall = mat(Fcols).rank() if Fcols else 0
    nontriv = rkDF - rkFall
    if verbose:
        print(f'[jideal] d={d} ideal=(z{idx})^{k}: nullity={nullity}, dim D-space={rkD}, trivial(F-mult in ideal)={rkF}, nontrivial dim={nontriv}', flush=True)
    # extract nontrivial representatives: rref of [Fcols; Dpolys], take D rows with new pivots
    reps = []
    if nontriv > 0:
        S = mat(Fcols + Dpolys)
        Sr, rk = S.rref()
        # pivots of F-part
        Fr, rkF2 = mat(Fcols).rref()
        fpiv = set()
        for i in range(rkF2):
            j = 0
            while int(Fr[i, j]) == 0: j += 1
            fpiv.add(j)
        inv_allm = {v: e for e, v in allm.items()}
        for i in range(rk):
            j = 0
            while int(Sr[i, j]) == 0: j += 1
            if j not in fpiv:
                P = {inv_allm[c]: int(Sr[i, c]) for c in range(len(allm)) if int(Sr[i, c])}
                reps.append(P)
    # convert reps back to original coordinates (mod p symmetric lift then exact?) -- we return mod-p polys in NEW coords plus the
    # transformation; the caller can test in new coordinates directly with the transformed F (same problem up to linear change).
    return {'nontrivial': nontriv, 'reps_newcoords_modp': reps, 'Fn': Fn, 'Mx': Mx, 'My': My, 'Mxi': Mxi, 'Myi': Myi, 'p': p}

def to_orig(P, Mxi, Myi):
    """P given in new coords z; Z = M z  =>  z = M^{-1} Z ; substitute z_i = sum_j Mi[i][j] Z_j."""
    return lin_subst(P, Mxi, Myi)

if __name__ == '__main__':
    I3 = [[1,0,0],[0,1,0],[0,0,1]]
    # original coordinates: ideal (Y0,Y1,Y2)^4 and (X0,X1,X2)^4, d=7,8,9
    for d in (6, 7, 8, 9):
        for name, idx, k in [('(Y)^4', [3,4,5], 4), ('(X)^4', [0,1,2], 4), ('(Y)^3', [3,4,5], 3), ('(X)^3', [0,1,2], 3)]:
            t0 = time.time()
            r = find(I3, I3, idx, k, d, verbose=False)
            print(f'd={d} {name}: nontrivial dim = {r["nontrivial"]}  ({time.time()-t0:.1f}s)', flush=True)
