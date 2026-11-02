"""Localized higher-degree relaxation on a rational linear subspace L of t-space (basis v_1..v_w), exact over Q.
Given the exact bilinear equations eq_r = sum_{a<=b} e_r[(a,b)] t_a t_b + sum_{a,k} e_r[(a,nt+k)] t_a s_k, restrict t = sum x_i v_i and
find all cofactor vectors mu_r(x) of degree k-2 whose combination sum_r mu_r(x) eq_r has no s-terms; the pure-x results are
degree-k consequences.  Returns the list of consequences (dict monomial-exponent-tuple -> Fraction) for the requested degrees."""
import itertools
from fractions import Fraction
from flint import fmpq_mat, fmpq, fmpz_mat, fmpz

def monomials_w(w, deg):
    out = []
    for comb in itertools.combinations(range(deg + w - 1), w - 1):
        e = []; prev = -1
        for x in comb: e.append(x - prev - 1); prev = x
        e.append(deg + w - 2 - prev)
        out.append(tuple(e))
    return out

def restrict(eqs, nt, ns, basis):
    """Returns per equation (alpha: dict (i,j) i<=j -> Fraction, beta: dict (i,k) -> Fraction) in the x-coordinates."""
    w = len(basis)
    out = []
    for e in eqs:
        alpha = {}; beta = {}
        for (a, b), v in e.items():
            if v == 0: continue
            if b < nt:
                for i in range(w):
                    va = basis[i][a]
                    if va == 0: continue
                    for j in range(w):
                        vb = basis[j][b]
                        if vb == 0: continue
                        if a == b:
                            key = (min(i, j), max(i, j)); alpha[key] = alpha.get(key, 0) + v * va * vb if i == j else alpha.get(key, 0) + v * va * vb
                        else:
                            key = (min(i, j), max(i, j)); alpha[key] = alpha.get(key, 0) + v * va * vb
            else:
                k = b - nt
                for i in range(w):
                    va = basis[i][a]
                    if va: beta[(i, k)] = beta.get((i, k), 0) + v * va
        # note: for a == b the double loop over (i,j) counts (i,j) and (j,i) separately -> coefficient of x_i x_j (i != j) accumulates both orders: correct expansion of t_a^2 = (sum x_i v_i[a])^2.
        # for a != b the term t_a t_b = (sum x_i v_i[a])(sum x_j v_j[b]) also accumulates (i,j) and (j,i) into key (min,max): correct.
        out.append(({k: v for k, v in alpha.items() if v}, {k: v for k, v in beta.items() if v}))
    return out

def consequences(restricted, w, ns, degrees=(2, 3, 4)):
    """Exact pure-x consequences of the given degrees."""
    neq = len(restricted)
    result = {}
    for k in degrees:
        dmons = monomials_w(w, k - 2)      # cofactor monomials
        gmons = monomials_w(w, k - 1)      # monomials of the s-part coefficients
        gidx = {m: i for i, m in enumerate(gmons)}
        unknowns = [(r, d) for r in range(neq) for d in range(len(dmons))]
        uidx = {u: i for i, u in enumerate(unknowns)}
        rows = {}
        for r, (alpha, beta) in enumerate(restricted):
            for (i, kk), v in beta.items():
                for di, dm in enumerate(dmons):
                    g = list(dm); g[i] += 1; g = tuple(g)
                    row = (kk, gidx[g])
                    rows.setdefault(row, {})
                    rows[row][uidx[(r, di)]] = rows[row].get(uidx[(r, di)], 0) + v
        R = len(rows); C = len(unknowns)
        # exact nullspace via fmpq_mat rref (clear denominators per row to integers)
        M = fmpq_mat(R, C)
        for ri, (row, ent) in enumerate(rows.items()):
            for ci, v in ent.items():
                if v: M[ri, ci] = fmpq(v.numerator, v.denominator)
        Mr, rank = M.rref()
        piv = []
        for i in range(rank):
            j = 0
            while j < C and Mr[i, j] == 0: j += 1
            piv.append(j)
        pset = set(piv)
        free = [j for j in range(C) if j not in pset]
        cons = []
        kmons = monomials_w(w, k)
        for f in free:
            mu = {f: Fraction(1)}
            for i, j in enumerate(piv):
                if Mr[i, f] != 0: mu[j] = -Fraction(int(Mr[i, f].p), int(Mr[i, f].q))
            g = {}
            for ui, coef in mu.items():
                r, di = unknowns[ui]
                dm = dmons[di]
                alpha = restricted[r][0]
                for (i, j), v in alpha.items():
                    m = list(dm); m[i] += 1; m[j] += 1; m = tuple(m)
                    g[m] = g.get(m, 0) + coef * v
            g = {m: v for m, v in g.items() if v}
            if g: cons.append(g)
        result[k] = {'unknowns': C, 'equations': R, 'kernel': len(free), 'consequences': cons}
    return result
