"""E2-mod-J quadratic system for the family D = D0 + w(F), w of degree d-5 (pure parity parD), where D0 is a
fixed E1-solution (A*D0 in J) with a particular N0 (N0(F) = A*D0), all mod p.
Condition (necessary for a certificate):  B_{D0}(N0) + D0 w(A) - A w(D0) - N0^T Hess(F) w - A w^T Hess(F) w  in J_{2d}.
Usage: python3 e2_mod_J.py <json with D0 (first item)> <d> <parD> [--with-D0-scale]  -> Singular file."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from quad_polar import field_basis, jac_rref
from flint import nmod_mat
p = 536870909

def solve_N0(Fp, Ap, D0, d):
    """Particular N0 mod p with sum N0_i F_i = A*D0 (pure-parity D0)."""
    pars = fixed_d.numerator_parities(D0, Ap)
    degN = d + 1
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    cols, meta = [], []
    for i in range(6):
        for m in core.monomials(degN, pars[i]):
            cols.append(core.shift(Fi[i], m, 1, p)); meta.append((i, m))
    rhs = core.mul(Ap, D0, p)
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    N0 = [{} for _ in range(6)]
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        assert j < C, 'E1 inconsistent mod p'
        v = int(Mr[i, C])
        if v:
            comp, m = meta[j]; N0[comp][m] = v
    # verify
    chk = {}
    for i in range(6): chk = core.add(chk, core.mul(N0[i], Fi[i], p), p)
    assert core.sub(chk, rhs, p) == {}, 'N0 verification failed'
    return N0

def main():
    path, d, parD = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    D0 = core.reduce_mod(core.parse(json.load(open(path))[0]['D']), p)
    assert core.degree(D0) == d
    t0 = time.time()
    N0 = solve_N0(Fp, Ap, D0, d)
    print(f'N0 found ({time.time()-t0:.0f}s)', flush=True)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    Hess = [[core.pderiv(Fi[j], l, p) for l in range(6)] for j in range(6)]
    basis = field_basis(d, parD); nb = len(basis)
    e = 2 * d
    parB = (core.parity_of(Ap) + 0) % 2  # parity of B = parity of A D^2 = parity of A (D pure) ... D even/odd: A*D^2 parity = par A
    mons, idx, Rr, pivots, r = jac_rref(Fp, e, parB)
    n = len(mons)
    nonpiv = [j for j in range(n) if j not in set(pivots)]
    print(f'R_{e}^{parB}: dim S={n}, dim J={r}, equations={len(nonpiv)}; unknowns={nb}', flush=True)
    # constant term: B_{D0}(N0)
    divN = {}
    for i in range(6): divN = core.add(divN, core.pderiv(N0[i], i, p), p)
    ND = {}
    for i in range(6): ND = core.add(ND, core.mul(N0[i], core.pderiv(D0, i, p), p), p)
    B0 = core.sub(core.mul(D0, divN, p), ND, p)
    terms = {(): B0}   # key: () constant, (a,) linear, (a,b) quadratic
    Ai = [core.pderiv(Ap, i, p) for i in range(6)]
    D0i = [core.pderiv(D0, i, p) for i in range(6)]
    # N0^T Hess w  = sum_j (sum_i N0_i F_ij) w_j ;  precompute HN_j = sum_i N0_i F_ij
    HN = [{} for _ in range(6)]
    for j in range(6):
        for i in range(6): HN[j] = core.add(HN[j], core.mul(N0[i], Hess[i][j], p), p)
    for a, (j, m) in enumerate(basis):
        # D0 * w(A) term: D0 * m * A_j ;  - A * w(D0): - A * m * D0_j ;  - HN_j * m
        t = core.shift(core.mul(D0, Ai[j], p), m, 1, p)
        t = core.sub(t, core.shift(core.mul(Ap, D0i[j], p), m, 1, p), p)
        t = core.sub(t, core.shift(HN[j], m, 1, p), p)
        if t: terms[(a,)] = t
    for a in range(nb):
        for b in range(a, nb):
            (j, m), (l, m2) = basis[a], basis[b]
            mm = tuple(x + y for x, y in zip(m, m2))
            t = core.shift(core.mul(Ap, Hess[j][l], p), mm, (p - 1) if a == b else (p - 2), p)  # -A w^T H w
            if t: terms[(a, b)] = t
    keys = list(terms)
    flat = []
    for k in keys:
        row = [0] * n
        for mon, v in terms[k].items():
            assert core.yparity(mon) == parB, ('parity mismatch', k)
            row[idx[mon]] = v
        flat.extend(row)
    U = nmod_mat(len(keys), n, flat, p)
    Upiv = nmod_mat(len(keys), r, [int(U[i, pj]) for i in range(len(keys)) for pj in pivots], p)
    Ured = U - Upiv * Rr
    eqs = []
    for c in nonpiv:
        eq = {}
        for i, k in enumerate(keys):
            v = int(Ured[i, c])
            if v: eq[k] = v
        if eq: eqs.append(eq)
    out = f'{RES}/e2modJ_d{d}_p{parD}.sing'
    with open(out, 'w') as f:
        f.write(f'ring R = {p}, (a(1..{nb})), dp;\nideal I;\n')
        for eq in eqs:
            ts = []
            for k, v in eq.items():
                if k == (): ts.append(f'{v}')
                elif len(k) == 1: ts.append(f'{v}*a({k[0]+1})')
                else: ts.append(f'{v}*a({k[0]+1})*a({k[1]+1})')
            f.write('I[size(I)+1] = ' + '+'.join(ts) + ';\n')
        f.write('int t = timer;\nideal G = std(I);\nprint("dim = " + string(dim(G)));\nprint("vdim = " + string(vdim(G)));\n')
        f.write('print("time = " + string(timer - t));\nif (dim(G) >= 0) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); L; }\nquit;\n')
    nconst = sum(1 for eq in eqs if () in eq)
    print(f'wrote {out}: {len(eqs)} equations ({nconst} with constant term), {nb} unknowns ({time.time()-t0:.0f}s)', flush=True)

if __name__ == '__main__':
    main()
