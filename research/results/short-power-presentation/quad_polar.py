"""Polar-divisor E2 filter.  For D = w(F) in J_d (w a vector field of degree d-5 with the parity making D
Y-pure), the certificate identity forces  A * Phi(w) in J,  Phi(w) = w^T Hess(F) w  (degree 2d-6).
If [A] is injective on R_{2d-6-?}... (precisely: on R_{2d-6}^{parD}) this gives Phi(w) in J_{2d-6}.
We build the quadratic equations for Phi(w) = 0 in R_{2d-6}^{parD} and hand them to Singular (char p).
Also: a direct mod-p check of the identity  B_D(A w) + A*Phi(w) in J  for a random w."""
import sys, time, json, random, subprocess
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from flint import nmod_mat
p = 536870909  # < 2^29 for Singular
RES = '/tmp/hodge/results/short-power-presentation'

def field_basis(d, parD):
    """List of (j, m): component j, monomial m of degree d-5, parity parD - par(F_j)."""
    out = []
    for j in range(6):
        pj = 0 if j < 3 else 1
        for m in core.monomials(d - 5, (parD - pj) % 2):
            out.append((j, m))
    return out

def jac_rref(Fp, e, par):
    """rref of J_e^{par} in the monomial basis of S_e^{par}: returns (mons, index, Rr, pivots, rank)."""
    mons = core.monomials(e, par); idx = {m: i for i, m in enumerate(mons)}
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    gens = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(e - 5, (par - pi) % 2):
            gens.append(core.shift(Fi[i], m, 1, p))
    K = len(gens); n = len(mons)
    flat = [0] * (K * n)
    for k, g in enumerate(gens):
        for m, v in g.items(): flat[k * n + idx[m]] = v
    M = nmod_mat(K, n, flat, p)
    Rr0, r = M.rref()
    Rr = nmod_mat(r, n, [int(Rr0[i, j]) for i in range(r) for j in range(n)], p)
    pivots = []
    for i in range(r):
        j = 0
        while int(Rr[i, j]) == 0: j += 1
        pivots.append(j)
    return mons, idx, Rr, pivots, r

def build(d, parD, verbose=True):
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p)
    e = 2 * d - 6
    basis = field_basis(d, parD)
    nb = len(basis)
    Hess = [[core.pderiv(core.pderiv(Fp, j, p), l, p) for l in range(6)] for j in range(6)]
    mons, idx, Rr, pivots, r = jac_rref(Fp, e, 0)   # Phi is always Y-even
    n = len(mons)
    nonpiv = [j for j in range(n) if j not in set(pivots)]
    pairs = [(a, b) for a in range(nb) for b in range(a, nb)]
    flat = [0] * (len(pairs) * n)
    for k, (a, b) in enumerate(pairs):
        (j, m), (l, m2) = basis[a], basis[b]
        mm = tuple(x + y for x, y in zip(m, m2))
        P = core.shift(Hess[j][l], mm, 1 if a == b else 2, p)   # a_a a_b counted twice for a != b
        for mon, v in P.items(): flat[k * n + idx[mon]] = v
    U = nmod_mat(len(pairs), n, flat, p)
    # reduce mod J: U_red = U - U[:, pivots] * Rr
    Upiv = nmod_mat(len(pairs), r, [int(U[k, pj]) for k in range(len(pairs)) for pj in pivots], p)
    Ured = U - Upiv * Rr
    eqs = []  # list of dict {(a,b): coef}
    for c in nonpiv:
        eq = {}
        for k, (a, b) in enumerate(pairs):
            v = int(Ured[k, c])
            if v: eq[(a, b)] = v
        if eq: eqs.append(eq)
    info = {'d': d, 'parD': parD, 'unknowns': nb, 'deg Phi': e, 'dim S_e': n, 'dim J_e': r,
            'equations': len(eqs), 'p': p}
    if verbose: print(json.dumps(info), flush=True)
    return basis, eqs, info

def singular_file(d, parD, basis, eqs, path):
    nb = len(basis)
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, (a(1..{nb})), dp;\n')
        f.write('ideal I;\n')
        for eq in eqs:
            terms = [f'{v}*a({a+1})*a({b+1})' for (a, b), v in eq.items()]
            f.write('I[size(I)+1] = ' + '+'.join(terms) + ';\n')
        f.write('int t = timer;\n')
        f.write('ideal G = std(I);\n')
        f.write('print("dim = " + string(dim(G)));\n')
        f.write('print("vdim = " + string(vdim(G)));\n')
        f.write('print("time = " + string(timer - t));\n')
        f.write('quit;\n')

def check_identity(d, parD, trials=2):
    """Direct check: for random w, with N = A w, B_D(N) + A Phi(w) in J_{d+6}?  (mod p)"""
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    basis = field_basis(d, parD)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    for t in range(trials):
        w = [{} for _ in range(6)]
        for (j, m) in basis:
            w[j] = core.add(w[j], {m: random.randrange(1, p)}, p)
        D = {}
        for j in range(6): D = core.add(D, core.mul(w[j], Fi[j], p), p)
        N = [core.mul(Ap, w[j], p) for j in range(6)]
        divN = {}
        for j in range(6): divN = core.add(divN, core.pderiv(N[j], j, p), p)
        ND = {}
        for j in range(6): ND = core.add(ND, core.mul(N[j], core.pderiv(D, j, p), p), p)
        B = core.sub(core.mul(D, divN, p), ND, p)
        Phi = {}
        for j in range(6):
            for l in range(6):
                Phi = core.add(Phi, core.mul(core.mul(w[j], w[l], p), core.pderiv(Fi[j], l, p), p), p)
        T = core.add(B, core.mul(Ap, Phi, p), p)
        # is T in J_{deg}?
        e = core.degree(T); par = core.parity_of(T)
        mons, idx, Rr, pivots, r = jac_rref(Fp, e, par)
        vec = nmod_mat(1, len(mons), [T.get(m, 0) for m in mons], p)
        vp = nmod_mat(1, r, [int(vec[0, pj]) for pj in pivots], p)
        red = vec - vp * Rr
        zero = all(int(red[0, c]) == 0 for c in range(len(mons)))
        print(f'identity check d={d} parD={parD} trial {t}: B_D(Aw) + A Phi(w) in J: {zero}', flush=True)

if __name__ == '__main__':
    mode = sys.argv[1]
    d = int(sys.argv[2]); parD = int(sys.argv[3])
    if mode == 'check':
        check_identity(d, parD)
    else:
        basis, eqs, info = build(d, parD)
        path = f'{RES}/polar_d{d}_p{parD}.sing'
        singular_file(d, parD, basis, eqs, path)
        with open(f'{RES}/polar_d{d}_p{parD}.json', 'w') as f: json.dump(info, f)
        print('wrote', path, flush=True)
