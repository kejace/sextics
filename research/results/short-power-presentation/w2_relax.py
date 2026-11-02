"""E2-mod-J (Hessian) relaxation over the full E1-admissible family K_d^par = J_d^par + <kernel reps>, mod p.
Unknowns u = (w_{j,m} for m in S_{d-5}^{parD+par(j)}, t_1..t_k).  D = w(F) + sum t_l D0_l,  N = A w + sum t_l N0_l,
G = sum t_l G0_l (F G0_l = A D0_l - N0_l(F)).  Condition: R = D G + D div N - N(D) in J_{2d}^{parA}  (homogeneous
quadratic system).  Writes Singular and msolve input files.  Usage: python3 w2_relax.py d parD [p]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

def solve_linear_modp(cols, rhs, p):
    """Find x with sum x_j cols[j] = rhs (polys mod p); returns dict j->x_j or None."""
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
    x = {}
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        if j == C: return None
        v = int(Mr[i, C])
        if v: x[j] = v
    return x

def particular_N0(Fp, Ap, D0, d, p):
    pars = fixed_d.numerator_parities(D0, Ap)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    cols, meta = [], []
    for i in range(6):
        for m in core.monomials(d + 1, pars[i]):
            cols.append(core.shift(Fi[i], m, 1, p)); meta.append((i, m))
    rhs = core.mul(Ap, D0, p)
    x = solve_linear_modp(cols, rhs, p)
    assert x is not None, 'E1 inconsistent'
    N0 = [{} for _ in range(6)]
    for j, v in x.items():
        i, m = meta[j]; N0[i][m] = v
    chk = {}
    for i in range(6): chk = core.add(chk, core.mul(N0[i], Fi[i], p), p)
    assert core.sub(chk, rhs, p) == {}
    return N0

def divide_by_F(Fp, P, deg, par, p):
    """G with F*G = P (deg G = deg)."""
    mons = core.monomials(deg, par)
    cols = [core.shift(Fp, m, 1, p) for m in mons]
    x = solve_linear_modp(cols, P, p)
    assert x is not None, 'not divisible by F'
    return {mons[j]: v for j, v in x.items()}

def build(d, parD, p, reps):
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    parA = core.parity_of(Ap)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    Hess = [[core.pderiv(Fi[j], l, p) for l in range(6)] for j in range(6)]
    Ai = [core.pderiv(Ap, i, p) for i in range(6)]
    # w basis
    wb = []
    for j in range(6):
        pj = 0 if j < 3 else 1
        for m in core.monomials(d - 5, (parD + pj) % 2):
            wb.append((j, m))
    nw = len(wb); k = len(reps)
    t0 = time.time()
    # kernel data
    N0s, G0s, HN = [], [], []
    for D0 in reps:
        N0 = particular_N0(Fp, Ap, D0, d, p)
        NF = {}
        for i in range(6): NF = core.add(NF, core.mul(N0[i], Fi[i], p), p)
        rhs = core.sub(core.mul(Ap, D0, p), NF, p)
        G0 = divide_by_F(Fp, rhs, d, (parA + parD) % 2, p) if rhs else {}
        N0s.append(N0); G0s.append(G0)
        H = [{} for _ in range(6)]
        for j in range(6):
            for i in range(6): H[j] = core.add(H[j], core.mul(N0[i], Hess[i][j], p), p)
        HN.append(H)
    print(f'kernel data for {k} reps ({time.time()-t0:.0f}s)', flush=True)
    e = 2 * d
    red = Reducer(Fp, e, parA, p, verbose=True)
    n = len(red.mons)
    terms = {}
    # w-w block: T(a,b)+T(b,a) = -2 A m m' F_jl  (a<b);  T(a,a) = -A m^2 F_jj
    for a in range(nw):
        j, m = wb[a]
        for b in range(a, nw):
            l, m2 = wb[b]
            mm = tuple(x + y for x, y in zip(m, m2))
            t = core.shift(core.mul(Ap, Hess[j][l], p), mm, (p - 1) if a == b else (p - 2), p)
            if t: terms[(a, b)] = t
    # w-t block: T(w_a,t_l)+T(t_l,w_a) = -m HN_l,j + D0_l (A_j m + A d_j m) - A m d_j D0_l
    for l in range(k):
        D0 = reps[l]; D0i = [core.pderiv(D0, i, p) for i in range(6)]
        for a in range(nw):
            j, m = wb[a]
            t = core.shift(HN[l][j], m, p - 1, p)
            t = core.add(t, core.shift(core.mul(D0, Ai[j], p), m, 1, p), p)
            dm = core.pderiv({m: 1}, j, p)
            if dm:
                t = core.add(t, core.mul(core.mul(Ap, D0, p), dm, p), p)
            t = core.sub(t, core.shift(core.mul(Ap, D0i[j], p), m, 1, p), p)
            if t: terms[(a, nw + l)] = t
    # t-t block
    def T(l, l2):
        D0 = reps[l]; N = N0s[l2]; G = G0s[l2]
        divN = {}
        for i in range(6): divN = core.add(divN, core.pderiv(N[i], i, p), p)
        out = core.mul(D0, core.add(G, divN, p), p)
        for i in range(6): out = core.sub(out, core.mul(N[i], core.pderiv(D0, i, p), p), p)
        return out
    for l in range(k):
        for l2 in range(l, k):
            t = T(l, l2) if l == l2 else core.add(T(l, l2), T(l2, l), p)
            if t: terms[(nw + l, nw + l2)] = t
    print(f'{len(terms)} nonzero quadratic coefficient vectors ({time.time()-t0:.0f}s)', flush=True)
    keys = list(terms)
    U = red.nf_rows([terms[key] for key in keys])
    nonpiv = [j for j in range(n) if j not in set(red.pivots)]
    eqs = []
    for c in nonpiv:
        eq = {}
        for i, key in enumerate(keys):
            v = int(U[i, c])
            if v: eq[key] = v
        if eq: eqs.append(eq)
    nu = nw + k
    base = f'{RES}/w2/relax_d{d}_p{parD}'
    with open(base + '.sing', 'w') as f:
        f.write(f'ring R = {p}, (a(1..{nu})), dp;\nideal I;\n')
        for eq in eqs:
            f.write('I[size(I)+1] = ' + '+'.join(f'{v}*a({a+1})*a({b+1})' for (a, b), v in eq.items()) + ';\n')
        f.write('int t = timer;\nideal G = std(I);\nprint("dim = " + string(dim(G)));\nprint("vdim = " + string(vdim(G)));\n')
        f.write('print("time = " + string(timer - t));\nprint("size GB = " + string(size(G)));\n')
        f.write('if (dim(G) >= 0) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); int ii; for (ii=1; ii<=size(L); ii++) { print("component " + string(ii) + " dim " + string(dim(std(L[ii]))) + " ngens " + string(size(L[ii]))); if (size(L[ii]) < 400) { L[ii]; } } }\nquit;\n')
    with open(base + '.ms', 'w') as f:
        f.write(','.join(f'a{i+1}' for i in range(nu)) + '\n' + f'{p}\n')
        f.write(',\n'.join('+'.join(f'{v}*a{a+1}*a{b+1}' for (a, b), v in eq.items()) for eq in eqs) + '\n')
    info = {'d': d, 'parD': parD, 'p': p, 'n_w': nw, 'n_t': k, 'unknowns': nu, 'equations': len(eqs),
            'dim_R_2d': n - red.rank, 'seconds': round(time.time() - t0, 1)}
    print(info, flush=True)
    json.dump({'info': info, 'wbasis': [(j, list(m)) for j, m in wb]}, open(base + '.json', 'w'))
    wlog(batch='relaxation-build', family='K_d+J_d Hessian/E2-mod-J', degree=d, parity=parD, **{k2: v for k2, v in info.items() if k2 not in ('d', 'parD')})
    return base

if __name__ == '__main__':
    d = int(sys.argv[1]); parD = int(sys.argv[2])
    p = int(sys.argv[3]) if len(sys.argv) > 3 else PRIMES[0]
    path = f'{RES}/w2/kernel_d{d}_p{parD}_mod{p}.json'
    reps = [core.reduce_mod(core.parse(s), p) for s in json.load(open(path))['reps']] if os.path.exists(path) else []
    print(f'd={d} parD={parD} p={p}: {len(reps)} kernel reps', flush=True)
    build(d, parD, p, reps)
