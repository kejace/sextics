"""Bilinear E1+E2 system for eigenvector denominators of degree d and character c.
Unknowns: t (coefficients of D in the colon space modulo F*S_{d-6}) and s (Koszul syzygy coefficients of N).
Equations: E2 = D div N - N(D) = 0 modulo F (cokernel functionals of multiplication by F)."""
import oplib, time, sys, json, random, subprocess, os
from hodge import core
from hodge.core import mul, scale, pderiv, add, sub, monomials, shift, NV
from flint import nmod_mat
PM = 1073741827   # prime for msolve (< 2^31)

def e1_solutions(F, A, q, d, c, p):
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    cA = oplib.char_of(Ap)
    degN = d + 6*q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    kinds = []; cols = []
    Dmons = oplib.monos_char(d, c)
    for m in Dmons:
        cols.append(scale(shift(Ap, m, 1, p), p-1, p)); kinds.append(('D', -1, m))
    for i in range(NV):
        ci = oplib.cadd(oplib.cadd(c, oplib.VARCHAR[i]), cA)
        for m in oplib.monos_char(degN, ci):
            cols.append(shift(Fi[i], m, q, p)); kinds.append(('N', i, m))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*C+j] = v
    M = nmod_mat(R, C, flat, p)
    X, n = M.nullspace()
    # transpose, rref: D-coordinates come first so pivots land on D first
    T = nmod_mat(n, C, [int(X[j, k]) for k in range(n) for j in range(C)], p)
    Tr, rk = T.rref()
    nD = len(Dmons)
    part = []; syz = []
    for i in range(rk):
        D = {}; N = [{} for _ in range(NV)]
        for j in range(C):
            v = int(Tr[i, j])
            if not v: continue
            kind, comp, m = kinds[j]
            if kind == 'D': D[m] = v
            else: N[comp][m] = v
        if D: part.append((D, N))
        else: syz.append(N)
    # gauge complement: D's independent modulo F*S_{d-6}
    gauge = [mul(Fp, {m: 1}, p) for m in (oplib.monos_char(d-6, c) if d >= 6 else [])]
    idx = {m: i for i, m in enumerate(Dmons)}
    def vec(D): 
        v = [0]*nD
        for e, val in D.items(): v[idx[e]] = val
        return v
    base = [vec(g) for g in gauge]
    chosen = []
    cur_rank = nmod_mat(len(base), nD, [x for r in base for x in r], p).rank() if base else 0
    for (D, N) in part:
        rows_ = base + [vec(D_) for (D_, N_) in chosen] + [vec(D)]
        rkk = nmod_mat(len(rows_), nD, [x for r in rows_ for x in r], p).rank()
        if rkk > cur_rank + len(chosen):
            chosen.append((D, N))
    return chosen, syz, {'colon_dim': len(part), 'gauge_dim': cur_rank, 'syz_dim': len(syz), 'nD': nD}

def B_poly(D, N, p):
    """D div N - N(D) mod p."""
    divN = {}; ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

def cokernel_F(F, deg, sector, p):
    """Left nullspace of multiplication by F : S_{deg-6,sector} -> S_{deg,sector}. Returns (row monomials idx, L matrix rows as dict lists)."""
    Fp = core.reduce_mod(F, p)
    tgt = oplib.monos_char(deg, sector)
    src = oplib.monos_char(deg-6, sector) if deg >= 6 else []
    idx = {m: i for i, m in enumerate(tgt)}
    R, C = len(tgt), len(src)
    if C == 0:
        # identity functionals
        return idx, None, R
    flat = [0]*(R*C)
    for j, m in enumerate(src):
        for e, v in shift(Fp, m, 1, p).items():
            flat[idx[e]*C + j] = v
    M = nmod_mat(R, C, flat, p)
    L, k = M.transpose().nullspace()   # C? no: M^T is C x R, nullspace vectors length R
    # L is R x k ; functionals are columns
    Lrows = [[int(L[i, j]) for i in range(R)] for j in range(k)]
    return idx, Lrows, k

def build_system(F, A, q, d, c, p=PM, verbose=True):
    t0 = time.time()
    chosen, syz, info = e1_solutions(F, A, q, d, c, p)
    cA = oplib.char_of(core.reduce_mod(A, p))
    sector = oplib.cadd(cA, oplib.cadd(c, c))
    nt, ns = len(chosen), len(syz)
    if nt == 0:
        return None, info
    idx, Lrows, k = cokernel_F(F, 2*d, sector, p)
    Ds = [D for D, N in chosen]; Ns = [N for D, N in chosen]
    R = len(idx)
    keys = []; vecs = []
    for a in range(nt):
        for b in range(nt):
            P = B_poly(Ds[a], Ns[b], p)
            v = [0]*R
            for e, val in P.items(): v[idx[e]] = val
            keys.append((('t', a), ('t', b)) if a <= b else (('t', b), ('t', a))); vecs.append(v)
        for j in range(ns):
            P = B_poly(Ds[a], syz[j], p)
            v = [0]*R
            for e, val in P.items(): v[idx[e]] = val
            keys.append((('t', a), ('s', j))); vecs.append(v)
    Mv = nmod_mat(len(vecs), R, [x for v in vecs for x in v], p)
    if Lrows is None:
        Pm = Mv
    else:
        Lt = nmod_mat(R, k, [Lrows[j][i] for i in range(R) for j in range(k)], p)
        Pm = Mv * Lt
    eqs = [{} for _ in range(k)]
    for rowi, key in enumerate(keys):
        for r in range(k):
            val = int(Pm[rowi, r])
            if val:
                eqs[r][key] = (eqs[r].get(key, 0) + val) % p
    info.update({'nt': nt, 'ns': ns, 'neq': k, 'build_s': round(time.time()-t0, 1)})
    if verbose: print(f'[bilinear d={d} c={c}] nt={nt} ns={ns} neq={k} info={info}', file=sys.stderr)
    return {'eqs': eqs, 'nt': nt, 'ns': ns, 'Ds': Ds, 'Ns': Ns, 'syz': syz, 'p': p, 'd': d, 'c': c}, info

def rank_defect(sysd, seed=1):
    """At a random t, number of independent linear conditions on s (cokernel dim of M(t)) and consistency."""
    p = sysd['p']; nt, ns = sysd['nt'], sysd['ns']
    rng = random.Random(seed)
    t = [rng.randrange(1, p) for _ in range(nt)]
    rows = []; rhs = []
    for eq in sysd['eqs']:
        row = [0]*ns; r = 0
        for (va, vb), cf in eq.items():
            if va[0] == 't' and vb[0] == 't':
                r = (r + cf * t[va[1]] * t[vb[1]]) % p
            else:
                row[vb[1]] = (row[vb[1]] + cf * t[va[1]]) % p
        rows.append(row); rhs.append(r)
    k = len(rows)
    if ns == 0:
        rk = 0
    else:
        rk = nmod_mat(k, ns, [x for r in rows for x in r], p).rank()
    rka = nmod_mat(k, ns+1, [x for r, b in zip(rows, rhs) for x in (r + [b])], p).rank()
    return {'rank_M': rk, 'rank_aug': rka, 'coker': k - rk, 'consistent_at_random_t': rka == rk}

def write_msolve(sysd, path, seed=1):
    p = sysd['p']; nt, ns = sysd['nt'], sysd['ns']
    names = [f't{i}' for i in range(nt)] + [f's{j}' for j in range(ns)]
    def vn(v): return f't{v[1]}' if v[0] == 't' else f's{v[1]}'
    lines = [','.join(names), str(p)]
    polys = []
    for eq in sysd['eqs']:
        terms = [f'{cf}*{vn(a)}*{vn(b)}' for (a, b), cf in eq.items() if cf]
        if terms: polys.append('+'.join(terms))
    rng = random.Random(seed)
    polys.append('+'.join(f'{rng.randrange(1,p)}*t{i}' for i in range(nt)) + f'+{p-1}')
    lines.append(',\n'.join(polys))
    with open(path, 'w') as f: f.write('\n'.join(lines) + '\n')
    return len(polys)

def run_msolve(path, out, threads=4, timeout=3600):
    t0 = time.time()
    try:
        r = subprocess.run(['msolve', '-f', path, '-o', out, '-t', str(threads), '-v', '0'], capture_output=True, text=True, timeout=timeout)
        res = open(out).read().strip() if os.path.exists(out) else ''
    except subprocess.TimeoutExpired:
        return 'timeout', time.time()-t0, ''
    dt = time.time() - t0
    if res.startswith('[-1]'): return 'no-solution', dt, res
    if res.startswith('[1,'): return 'infinite', dt, res
    return 'finite-or-other', dt, res[:2000]

if __name__ == '__main__':
    F, A, q = oplib.load()
    d = int(sys.argv[1]); 
    chars = [(int(x[0]), int(x[1])) for x in sys.argv[2:]] if len(sys.argv) > 2 else [(a, b) for a in range(6) for b in range(6)]
    threads = int(os.environ.get('MSOLVE_THREADS', '4'))
    for c in chars:
        sysd, info = build_system(F, A, q, d, c)
        if sysd is None:
            print(f'd={d} c={c}: empty colon space (after gauge) info={info}'); sys.stdout.flush(); continue
        rd = rank_defect(sysd)
        path = f'/tmp/hodge/results/off-paired/msolve_d{d}_c{c[0]}{c[1]}.ms'
        npol = write_msolve(sysd, path)
        st, dt, res = run_msolve(path, path.replace('.ms', '.out'), threads=threads, timeout=int(os.environ.get('MSOLVE_TIMEOUT', '900')))
        rec = {'d': d, 'char': c, **info, **rd, 'msolve': st, 'msolve_seconds': round(dt, 1), 'npolys': npol}
        print(json.dumps(rec)); sys.stdout.flush()
        with open('/tmp/hodge/results/off-paired/bilinear_results.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
        oplib.log(f'[all eigenvector D of degree {d}, character {c}]', d, 'bilinear-sector-sweep', True, (st != 'no-solution'), rd['rank_M'], rd['coker'], dt,
                  f"msolve {st}; nt={info['nt']} ns={info['ns']} neq={info['neq']} coker_at_random_t={rd['coker']}")
        if st not in ('no-solution', 'timeout'):
            print('!!! POSSIBLE SOLUTIONS:', res[:500]); sys.stdout.flush()
