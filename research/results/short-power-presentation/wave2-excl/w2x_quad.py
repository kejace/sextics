"""E2-mod-J quadratic relaxation (Jacobian-type reduction), written independently of the wave-1 scripts.

Family (q = 1, deg D = d, Y-parities `pars` of D):  D = sum_i lam_i D0_i + w(F),  w = sum_a x_a m_a d/dz_{j_a}
with {D0_i} a mod-p basis of K_d^par / J_d^par (E1 kernel) and (j_a, m_a) running over monomials m_a of degree d-5
with parity par - par(z_{j_a}).  Every certificate denominator of degree d with these parities lies in this family
(E1: A D in J).  With N0_i(F) = A D0_i, every admissible numerator is N = sum lam_i N0_i + A w + Z G/6 + K(F)
(Koszul syzygies, F smooth); the G and K terms contribute 0 resp. elements of J to
    B_D(N) := D div N - N(D),
and E2 (B_D(N) - D G = F Q) forces   B_D( sum lam_i N0_i + A w )  in  J_{2d}.
This script builds those quadratic equations in (lam, x) modulo J_{2d} (normal form by rref), checks that they are
invariant under the Euler translations w -> w + l Z (l of degree d-6; those give D -> D + 6 l F, not coprime), restricts
to the complement W' = {x_(0, l*X0) = 0}, and writes msolve + Singular inputs.  If the restricted ideal I' is
zero-dimensional at the origin mod p (Hilbert function of S'/I' vanishes in some degree) then, since
dim_Q (S'/I'_Q)_k <= dim_p (S'/I'_p)_k, the same holds over Q, so V_C(I') = {0}: every denominator of the family is
a multiple of F  =>  no certificate with such D.  Ranks of J are compared with the smooth Hilbert function.
Usage:  python3 w2x_quad.py <d> <e|o|m> [--run-singular] [--run-msolve T]"""
import sys, os, time, json, subprocess
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave2-excl')
from w2x_lib import *
from hodge.core import monomials, shift, mul, add, sub, scale, pderiv, parity_of, degree

p = P

def var_parity(j): return 0 if j < 3 else 1

def field_basis(d, parD):
    return [(j, m) for j in range(6) for m in monomials(d - 5, (parD - var_parity(j)) % 2)]

def rref_J(Fp, k, par):
    """rref of J_k^{par} spanned by m*F_i (rows) in the monomial basis of S_k^{par}. Returns mons, idx, Rr, pivots, rank."""
    mons = monomials(k, par); idx = {m: i for i, m in enumerate(mons)}; n = len(mons)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    gens = []
    for i in range(6):
        for m in monomials(k - 5, (par - var_parity(i)) % 2):
            gens.append(shift(Fi[i], m, 1, p))
    flat = [0] * (len(gens) * n)
    for r, g in enumerate(gens):
        for m, v in g.items(): flat[r * n + idx[m]] = v
    Rr0, rk = nmod_mat(len(gens), n, flat, p).rref()
    Rr = nmod_mat(rk, n, [int(Rr0[i, j]) for i in range(rk) for j in range(n)], p)
    pivots = []
    for i in range(rk):
        j = 0
        while int(Rr[i, j]) == 0: j += 1
        pivots.append(j)
    return mons, idx, Rr, pivots, rk

def normal_forms(polys, mons, idx, Rr, pivots, chunk=1500):
    """Rows = normal forms of polys modulo the row space of Rr (dense, chunked)."""
    n = len(mons); r = len(pivots)
    out_rows = []
    for s in range(0, len(polys), chunk):
        blk = polys[s:s + chunk]
        flat = [0] * (len(blk) * n)
        for i, P_ in enumerate(blk):
            for m, v in P_.items(): flat[i * n + idx[m]] = v
        U = nmod_mat(len(blk), n, flat, p)
        Up = nmod_mat(len(blk), r, [int(U[i, pj]) for i in range(len(blk)) for pj in pivots], p)
        Ured = U - Up * Rr
        out_rows.append(Ured)
    return out_rows

def solve_N0(Fp, Ap, D0, d):
    """Particular N0 (mod p) with sum_j N0_j F_j = A D0."""
    rhs = mul(Ap, D0, p); parR = parity_of(rhs)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    cols, meta = [], []
    for i in range(6):
        for m in monomials(d + 1, (parR - var_parity(i)) % 2):
            cols.append(shift(Fi[i], m, 1, p)); meta.append((i, m))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    Mr, rank = nmod_mat(R, C + 1, flat, p).rref()
    N0 = [{} for _ in range(6)]
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        assert j < C, 'E1 inconsistent for D0 (mod p)'
        v = int(Mr[i, C])
        if v:
            comp, m = meta[j]; N0[comp][m] = v
    chk = {}
    for i in range(6): chk = add(chk, mul(N0[i], Fi[i], p), p)
    assert sub(chk, rhs, p) == {}, 'N0 check failed'
    return N0

def divergence(N): 
    out = {}
    for j in range(6): out = add(out, pderiv(N[j], j, p), p)
    return out
def apply(N, D):
    out = {}
    for j in range(6):
        if N[j]: out = add(out, mul(N[j], pderiv(D, j, p), p), p)
    return out

def build(d, mode, kernel_bases=None, verbose=True):
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    pars = [0] if mode == 'e' else [1] if mode == 'o' else [0, 1]
    # ---- unknowns: lam_i (kernel reps) then x_a (vector-field basis) ----
    D0s = []   # (parity, poly)
    for par in pars:
        fn = f'{PK}/e1kernel_d{d}_p{par}.json'
        if kernel_bases is not None and par in kernel_bases:
            basis = kernel_bases[par]
        elif os.path.exists(fn):
            basis = [core.reduce_mod(core.parse(s), p) for s in json.load(open(fn))['basis_modp']]
            if verbose: print(f'kernel basis d={d} par={par} from {fn}: {len(basis)} elements', flush=True)
        else:
            out = jacobian.mult_kernel(F, A, d, p=p, parity=par, verbose=True)
            basis = out['basis']
            with open(fn, 'w') as f:
                json.dump({'rec': {'exp': 'e1kernel(from quad)', 'd': d, 'parD': par, 'p': p, 'kernel mod J': out['dim'],
                                   'dim J_d (mod p)': out['jacobian_dim'], 'dim J_d (expected smooth)': expected_rank_J(d, par)},
                           'basis_modp': [core.to_str(b) for b in basis]}, f, indent=1)
        for b in basis: D0s.append((par, b))
    wb = []
    for par in pars:
        for (j, m) in field_basis(d, par): wb.append((par, j, m))
    nl, nw = len(D0s), len(wb); n = nl + nw
    names = [f'lam{i+1}' for i in range(nl)] + [f'x{a+1}' for a in range(nw)]
    if verbose: print(f'd={d} mode={mode}: {nl} kernel reps, {nw} field coordinates, {n} unknowns', flush=True)
    # ---- D_u, N_u for each unknown ----
    t0 = time.time()
    Du, Nu, paru = [], [], []
    for (par, D0) in D0s:
        Du.append(D0); Nu.append(solve_N0(Fp, Ap, D0, d)); paru.append(par)
    if verbose and nl: print(f'N0 solved for {nl} kernel reps ({time.time()-t0:.0f}s)', flush=True)
    for (par, j, m) in wb:
        Du.append(shift(Fi[j], m, 1, p))
        N = [{} for _ in range(6)]; N[j] = shift(Ap, m, 1, p)
        Nu.append(N); paru.append(par)
    divN = [divergence(N) for N in Nu]
    # ---- bilinear coefficients T_uv, grouped by parity of T ----
    pairs, T = [], []
    for u in range(n):
        for v in range(u, n):
            if u == v:
                t = sub(mul(Du[u], divN[u], p), apply(Nu[u], Du[u]), p)
            else:
                t = sub(mul(Du[u], divN[v], p), apply(Nu[v], Du[u]), p)
                t = add(t, sub(mul(Du[v], divN[u], p), apply(Nu[u], Du[v]), p), p)
            pairs.append((u, v)); T.append(t)
    if verbose: print(f'{len(pairs)} bilinear coefficient polynomials built ({time.time()-t0:.0f}s)', flush=True)
    parT = [(1 + paru[u] + paru[v]) % 2 for (u, v) in pairs]   # parity of A D^2
    eqs = []   # list of dict {(u,v): coeff}
    jinfo = {}
    for par in sorted(set(parT)):
        mons, idx, Rr, pivots, rk = rref_J(Fp, 2 * d, par)
        jinfo[par] = {'deg': 2 * d, 'dim S': len(mons), 'rank J (mod p)': rk, 'rank J (expected)': expected_rank_J(2 * d, par)}
        assert jinfo[par]['rank J (mod p)'] == jinfo[par]['rank J (expected)'], 'Hilbert check failed'
        sel = [k for k in range(len(pairs)) if parT[k] == par]
        for t_, k in zip(T, range(len(T))):
            if parT[k] == par: assert all(core.yparity(m) == par for m in t_), 'parity mismatch'
        blocks = normal_forms([T[k] for k in sel], mons, idx, Rr, pivots)
        nonpiv = [c for c in range(len(mons)) if c not in set(pivots)]
        # equations: one per non-pivot column
        cols = {c: {} for c in nonpiv}
        off = 0
        for U in blocks:
            nr = U.nrows()
            for i in range(nr):
                k = sel[off + i]
                for c in nonpiv:
                    v = int(U[i, c])
                    if v: cols[c][pairs[k]] = v
            off += nr
        for c in nonpiv:
            if cols[c]: eqs.append(cols[c])
        if verbose: print(f'parity {par}: {len(nonpiv)} coordinate equations ({time.time()-t0:.0f}s)', flush=True)
    # ---- row-reduce the equations in the basis of quadratic monomials ----
    pidx = {pr: k for k, pr in enumerate(pairs)}
    flat = [0] * (len(eqs) * len(pairs))
    for i, e in enumerate(eqs):
        for pr, v in e.items(): flat[i * len(pairs) + pidx[pr]] = v
    Er, erk = nmod_mat(len(eqs), len(pairs), flat, p).rref()
    eqs2 = []
    for i in range(erk):
        e = {}
        for k in range(len(pairs)):
            v = int(Er[i, k])
            if v: e[pairs[k]] = v
        eqs2.append(e)
    if verbose: print(f'{len(eqs)} equations, rank {erk} (independent quadrics) ({time.time()-t0:.0f}s)', flush=True)
    # ---- Euler translation invariance check ----
    widx = {(par, j, m): nl + a for a, (par, j, m) in enumerate(wb)}
    euler = []
    for par in pars:
        for l in monomials(d - 6, par):
            vec = {}
            for j in range(6):
                m = tuple(l[i] + (1 if i == j else 0) for i in range(6))
                vec[widx[(par, j, m)]] = 1
            euler.append((par, l, vec))
    viol = 0
    for e in eqs2:
        for (par, l, vec) in euler:
            # linear form  x -> Q(vec, x)
            lin = {}
            inv2 = pow(2, -1, p)
            for (u, v), c in e.items():
                if u == v:
                    if u in vec: lin[u] = (lin.get(u, 0) + c * vec[u]) % p
                else:
                    if u in vec: lin[v] = (lin.get(v, 0) + c * inv2 * vec[u]) % p
                    if v in vec: lin[u] = (lin.get(u, 0) + c * inv2 * vec[v]) % p
            viol += sum(1 for c in lin.values() if c)
    if verbose: print(f'Euler invariance: {len(euler)} directions, violations = {viol}', flush=True)
    # ---- restrict to complement W' : drop x_(par,0,l*X0) ----
    drop = set()
    for (par, l, vec) in euler:
        m = tuple(l[i] + (1 if i == 0 else 0) for i in range(6))
        drop.add(widx[(par, 0, m)])
    keep = [u for u in range(n) if u not in drop]
    kmap = {u: i for i, u in enumerate(keep)}
    eqs3 = []
    for e in eqs2:
        e3 = {}
        for (u, v), c in e.items():
            if u in drop or v in drop: continue
            e3[(kmap[u], kmap[v])] = c
        if e3: eqs3.append(e3)
    names3 = [names[u] for u in keep]
    info = {'exp': 'quad-relaxation', 'd': d, 'mode': mode, 'p': p, 'n_kernel_reps': nl, 'n_field_coords': nw,
            'unknowns_full': n, 'euler_directions': len(euler), 'euler_invariance_violations': viol,
            'unknowns_restricted': len(keep), 'equations_raw': len(eqs), 'equations_independent': erk,
            'equations_restricted': len(eqs3), 'J_checks': jinfo, 'build_seconds': round(time.time() - t0, 1)}
    return names3, eqs3, info, {'names_full': names, 'drop': sorted(drop), 'eqs_full': [[list(k), v] for e in eqs2 for k, v in e.items()]}

def write_msolve(path, names, eqs):
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n'); f.write(f'{p}\n')
        lines = []
        for e in eqs:
            lines.append('+'.join(f'{v}*{names[u]}*{names[v_]}' if u != v_ else f'{v}*{names[u]}^2' for (u, v_), v in e.items()))
        f.write(',\n'.join(lines) + '\n')

def write_singular(path, names, eqs, extra=''):
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, ({",".join(names)}), dp;\noption(redSB);\nideal I;\n')
        for e in eqs:
            f.write('I[size(I)+1] = ' + '+'.join(f'{v}*{names[u]}*{names[v_]}' for (u, v_), v in e.items()) + ';\n')
        f.write('int t = timer;\nideal G = std(I);\nprint("dim = " + string(dim(G)));\nprint("vdim = " + string(vdim(G)));\n')
        f.write('print("time = " + string(timer - t));\nintvec hv = hilb(G, 2);\nprint("hilb2 = " + string(hv));\n')
        f.write('if (dim(G) > 0) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); L; }\n')
        f.write(extra + 'quit;\n')

if __name__ == '__main__':
    d = int(sys.argv[1]); mode = sys.argv[2]
    names, eqs, info, full = build(d, mode)
    tag = f'quad_d{d}_{mode}'
    write_msolve(f'{PK}/{tag}.ms', names, eqs)
    write_singular(f'{PK}/{tag}.sing', names, eqs)
    with open(f'{PK}/{tag}.eqs.json', 'w') as f: json.dump({'names': names, 'eqs': [[list(k), v] for e in eqs for k, v in e.items()], 'eq_sizes': [len(e) for e in eqs], 'full': full}, f)
    with open(f'{PK}/{tag}.info.json', 'w') as f: json.dump(info, f, indent=1)
    print(json.dumps(info), flush=True)
    if '--run-singular' in sys.argv:
        t0 = time.time()
        r = subprocess.run(['Singular', '-q', f'{PK}/{tag}.sing'], capture_output=True, text=True, timeout=36000)
        out = r.stdout + r.stderr
        with open(f'{PK}/{tag}.singular.out', 'w') as f: f.write(out)
        print(out[-3000:], flush=True)
        info['singular_seconds'] = round(time.time() - t0, 1)
        info['singular_dim'] = next((l.split('=')[1].strip() for l in out.splitlines() if l.startswith('dim =')), None)
        info['singular_vdim'] = next((l.split('=')[1].strip() for l in out.splitlines() if l.startswith('vdim =')), None)
        info['singular_hilb2'] = next((l.split('=',1)[1].strip() for l in out.splitlines() if l.startswith('hilb2 =')), None)
        with open(f'{PK}/{tag}.info.json', 'w') as f: json.dump(info, f, indent=1)
        log(info)
