"""Structured candidate families vs the E1-admissible space (exact mod p), degrees 6..12.

Uses out/e1kernel_d{d}_p{par}.json (kernel reps mod p from e1spaces.py) + J(F)_d generators:
  W_d = {D : A D in J} = span(kernel reps) + J_d.   For each structured subspace U:  dim(U cap W_d).
Members of U cap W_d (mod p) are written out; individual candidates get an exact e1 test and, for d <= 9,
the engine's fixed-D solve.
  python3 wave2/structured.py --dmin 6 --dmax 10 [--parities 0,1] [--tag structured]
"""
import sys, os, time, json, argparse, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=6)
ap.add_argument('--dmax', type=int, default=10)
ap.add_argument('--parities', default='0,1')
ap.add_argument('--prime-index', type=int, default=0)
ap.add_argument('--tag', default='structured')
ap.add_argument('--solve-dmax', type=int, default=9, help='run the exact fixed-D solve for passing members up to this degree')
args = ap.parse_args()

F, A, q = load()
st = structure(F)
f, g, L, M, Pc, t = st['f'], st['g'], st['L'], st['M'], st['P'], st['t']
p = core.PRIMES[args.prime_index]
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(F, i) for i in range(NV)]
Fip = [core.reduce_mod(x, p) for x in Fi]
Hess = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
LOGP = f'{OUT}/{args.tag}.jsonl'


def mono(e):
    return {tuple(e): Fraction(1)}


def var(i):
    return mono([1 if j == i else 0 for j in range(6)])


def xmons(k, vars_):
    """monomials of degree k in the given variable subset."""
    out = []
    for m in monomials(k):
        if all(m[j] == 0 for j in range(6) if j not in vars_):
            out.append(m)
    return out


def det3(Mx):
    a = Mx
    t1 = mul(a[0][0], sub(mul(a[1][1], a[2][2]), mul(a[1][2], a[2][1])))
    t2 = mul(a[0][1], sub(mul(a[1][0], a[2][2]), mul(a[1][2], a[2][0])))
    t3 = mul(a[0][2], sub(mul(a[1][0], a[2][1]), mul(a[1][1], a[2][0])))
    return add(sub(t1, t2), t3)


def structured_subspaces(d, par):
    """Yield (name, list of polys) of degree d and Y-parity par (exact rational polys)."""
    S = Space(d, par)
    # blocks and variable subsets
    for r in (3, 4, 5):
        for sub_ in itertools.combinations(range(6), r):
            ms = [m for m in xmons(d, sub_) if core.yparity(m) == par]
            if ms:
                name = 'vars[' + ''.join(core.VARS[j] for j in sub_) + ']'
                yield name, [mono(m) for m in ms]
    # bidegree pieces
    for i in range(d + 1):
        ms = [m for m in S.mons if sum(m[:3]) == i]
        if ms:
            yield f'bideg({i},{d-i})', [mono(m) for m in ms]
    # base-type P(X) M^a g^b and Q(Y) L^a f^b
    polys = []
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0:
                continue
            W = power(M, a) if a else {(0,) * 6: Fraction(1)}
            if b:
                W = mul(W, g)
            for m in xmons(k, (0, 1, 2)):
                P_ = mul(mono(m), W)
                if core.parity_of(P_) == par:
                    polys.append(P_)
    if polys:
        yield 'Ybase P(X)M^a g^b', polys
    polys = []
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0:
                continue
            W = power(L, a) if a else {(0,) * 6: Fraction(1)}
            if b:
                W = mul(W, f)
            for m in xmons(k, (3, 4, 5)):
                P_ = mul(mono(m), W)
                if core.parity_of(P_) == par:
                    polys.append(P_)
    if polys:
        yield 'Xbase Q(Y)L^a f^b', polys
    # monomials in f, g, L, M (all Y-even)
    if par == 0:
        polys = []
        for a in range(d // 6 + 1):
            for b in range(d // 6 + 1):
                rest = d - 6 * a - 6 * b
                if rest < 0 or rest % 2:
                    continue
                for c in range(rest // 2 + 1):
                    e = rest // 2 - c
                    P_ = {(0,) * 6: Fraction(1)}
                    for base_, k in ((f, a), (g, b), (L, c), (M, e)):
                        if k:
                            P_ = mul(P_, power(base_, k))
                    polys.append(P_)
        yield 'monomials(f,g,L,M)', polys
    # Hessian minors
    if d == 8:
        polys = []
        for (i, j) in itertools.combinations(range(6), 2):
            for (k, l) in itertools.combinations(range(6), 2):
                P_ = sub(mul(Hess[i][k], Hess[j][l]), mul(Hess[i][l], Hess[j][k]))
                if P_ and core.parity_of(P_) == par:
                    polys.append(P_)
        if polys:
            yield 'Hessian 2x2 minors', polys
    if d == 12:
        polys = []
        for rows in itertools.combinations(range(6), 3):
            for cols in itertools.combinations(range(6), 3):
                sub3 = [[Hess[r][c] for c in cols] for r in rows]
                P_ = det3(sub3)
                if P_ and core.parity_of(P_) == par:
                    polys.append(P_)
        if polys:
            yield 'Hessian 3x3 minors', polys
        if par == 0:
            hf = det3([[pderiv(pderiv(f, i), j) for j in range(3)] for i in range(3)])
            hg = det3([[pderiv(pderiv(g, i), j) for j in range(3, 6)] for i in range(3, 6)])
            yield 'Hess(f), Hess(g)', [hf, hg]
            yield 'products of (f,g,P)', [mul(a_, b_) for a_, b_ in itertools.combinations_with_replacement([f, g, Pc], 2)]
    # radial derivative products of degree d (in J, trivial e1) are skipped on purpose (sum principle)


def admissible_matrix(d, par):
    """Rows spanning W_d = J_d + kernel reps, as nmod_mat over S_d^par monomials; returns (mat, rank, Space)."""
    S = Space(d, par)
    gens = jacobian.jacobian_generators(F, d, p, par)
    rows = [gpoly for (_, _, gpoly) in gens]
    kpath = f'{OUT}/e1kernel_d{d}_p{par}.json'
    nker = 0
    if os.path.exists(kpath):
        kj = json.load(open(kpath))
        assert kj['p'] == p
        for b in kj['basis_modp']:
            rows.append(core.parse(b, p)); nker += 1
    else:
        print(f'[structured] WARNING no kernel file for d={d} par={par}; assuming kernel 0 (true for d <= 6)', file=sys.stderr)
    flat = [0] * (len(rows) * S.n)
    for r_, P_ in enumerate(rows):
        for e, v in P_.items():
            flat[r_ * S.n + S.idx[e]] = v
    Wm = nmod_mat(len(rows), S.n, flat, p)
    rk = Wm.rank()
    return Wm, rk, S, nker


def intersect(Wm, rkW, S, polys):
    """dim(U cap W) and members (mod p) for U = span(polys)."""
    r = len(polys)
    flat = [0] * (r * S.n)
    for i, P_ in enumerate(polys):
        Pp = core.reduce_mod(P_, p)
        for e, v in Pp.items():
            flat[i * S.n + S.idx[e]] = v
    Um = nmod_mat(r, S.n, flat, p)
    rkU = Um.rank()
    # stack
    both = nmod_mat(Wm.nrows() + r, S.n, [int(Wm[i, j]) for i in range(Wm.nrows()) for j in range(S.n)] + flat, p)
    rkB = both.rank()
    dim_int = rkU + rkW - rkB
    members = []
    if dim_int > 0:
        # find coefficient vectors x (len r) with U^T x in W: nullspace of [W^T | -U^T] restricted to U-coordinates
        WT = Wm.transpose(); UT = Um.transpose()
        Cn = Wm.nrows() + r
        flat2 = [0] * (S.n * Cn)
        for i in range(S.n):
            for j in range(Wm.nrows()):
                v = int(WT[i, j])
                if v: flat2[i * Cn + j] = v
            for j in range(r):
                v = int(UT[i, j])
                if v: flat2[i * Cn + Wm.nrows() + j] = (p - v) % p
        X, nul = nmod_mat(S.n, Cn, flat2, p).nullspace()
        proj = nmod_mat(nul, r, [int(X[Wm.nrows() + j, k]) for k in range(nul) for j in range(r)], p)
        Pr, rk = proj.rref()
        for i in range(rk):
            D_ = {}
            for j in range(r):
                v = int(Pr[i, j])
                if v:
                    D_ = add(D_, scale(core.reduce_mod(polys[j], p), v, p), p)
            if D_:
                members.append(D_)
    return dim_int, rkU, members


if __name__ == '__main__':
    pars = [int(x) for x in args.parities.split(',')]
    for d in range(args.dmin, args.dmax + 1):
        for par in pars:
            t0 = time.time()
            Wm, rkW, S, nker = admissible_matrix(d, par)
            print(f'[structured] d={d} par={par}: W dim {rkW} (J_d + {nker} kernel reps) in S_d dim {S.n} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
            for name, polys in structured_subspaces(d, par):
                t1 = time.time()
                polys = [P_ for P_ in polys if P_]
                if not polys:
                    continue
                dim_int, rkU, members = intersect(Wm, rkW, S, polys)
                # how much of the intersection is already J_d (trivial)?  dim(U cap J_d)
                rec = {'batch': 'structured', 'd': d, 'parity': par, 'family': name, 'dim_U': int(rkU), 'dim_W': int(rkW),
                       'dim_U_cap_W': int(dim_int), 'p': p, 'time_s': round(time.time() - t1, 1)}
                if dim_int > 0:
                    # J-part of the intersection
                    gens = jacobian.jacobian_generators(F, d, p, par)
                    Jrows = [gpoly for (_, _, gpoly) in gens]
                    flatJ = [0] * (len(Jrows) * S.n)
                    for r_, P_ in enumerate(Jrows):
                        for e, v in P_.items():
                            flatJ[r_ * S.n + S.idx[e]] = v
                    Jm = nmod_mat(len(Jrows), S.n, flatJ, p) if Jrows else nmod_mat(0, S.n, [], p)
                    rkJ = Jm.rank() if Jrows else 0
                    dimJint, _, _ = intersect(Jm, rkJ, S, polys) if Jrows else (0, 0, [])
                    rec['dim_U_cap_J'] = int(dimJint)
                    rec['dim_new'] = int(dim_int - dimJint)
                    rec['members_modp'] = [core.to_str(jacobian.symmetric(m, p))[:3000] for m in members[:5]]
                print(json.dumps({k: v for k, v in rec.items() if k != 'members_modp'}), flush=True)
                log_line(LOGP, rec)
