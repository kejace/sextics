"""Job D (remote): (1) structured sigma-eigen candidate families up to degree 12: rigorous E1 screen (A*D not in J mod p
=> no certificate over D, any power) and full restricted E2 for survivors; (2) corrected served-target computation in
the character-3 sector (modulo J) for the key denominators (wave-1 ran it in the wrong sector); (3) exact rational
(LLL-reduced) E1-kernel representatives at d = 9 for every sector (wave-1 rational lift failed in 8 of 11 sectors) and
exact rational dual certificates for each of them; (4) Singular E1 screen for high-degree structured candidates
(P_eps, Hessian determinant, products).  Outputs jobD_*/ ; log lines in wave2-excl.log.jsonl."""
import sys, os, json, time, random, subprocess, itertools
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000003
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 16
F, A, q = load()
OUT = PACKET + f'/jobD_p{P}'; os.makedirs(OUT, exist_ok=True); os.makedirs(OUT + '/duals', exist_ok=True)
def say(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
x, y, z, t, r, s = [core.parse(v) for v in core.VARS]
Fi = [core.pderiv(F, i) for i in range(6)]
Fstr = open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0]
Jhat = core.parse(f'7*X0*X2**2*(9*(Y0**6 - Y1**6 + Y2**6) - 9*({Fstr})) + 15*X1*Y2*Y0*(2*Y0**6 - 4*Y1**6 + 8*Y2**6 - 2*({Fstr}))')
def pp(s_): return core.parse(s_)

# ------------------------------------------------------------------ (1) structured candidates
def structured_candidates():
    C = []
    X_, Y_ = 'X0,X1,X2'.split(','), 'Y0,Y1,Y2'.split(',')
    # pure-Y invariants (char 0)
    for e, name in [('Y0**6+Y1**6+Y2**6', 'G Fermat Y-block'), ('Y0*Y1*Y2', 'trs'), ('Y0**2+Y1**2+Y2**2', 'Y-quadric sum'),
                    ('Y0**3+Y1**3+Y2**3', 'Y-Fermat cubic'), ('Y0**4+Y1**4+Y2**4', 'Y-Fermat quartic'), ('Y0**3*Y1+Y1**3*Y2+Y2**3*Y0', 'Y-Klein quartic'),
                    ('Y0*Y1+Y1*Y2+Y2*Y0', 'Y elementary sym 2'), ('Y0**3+Y1**3+Y2**3-3*Y0*Y1*Y2', 'Y Hesse-type cubic'),
                    ('Y0**6+Y1**6+Y2**6-30*Y0**2*Y1**2*Y2**2', 'Y sextic with (trs)^2'), ('(Y0**2+Y1**2+Y2**2)**3', 'Y quadric cubed'),
                    ('Y0**5*Y1+Y1**5*Y2+Y2**5*Y0', 'Y Klein sextic'), ('Y0*Y1*Y2*(Y0**3+Y1**3+Y2**3)', 'trs * Fermat cubic')]:
        C.append((name, pp(e), 'pure-Y'))
    # pure-X sigma-eigen invariants and Klein-block pieces
    for e, name in [('X0*X1*X2', 'xyz'), ('X0**3*X2', 'x^3 z'), ('X0*X1**3', 'x y^3'), ('X1*X2**3', 'y z^3'),
                    ('X0**5*X1+X1**5*X2+X2**5*X0', 'K Klein sextic block'), ('X0**2*X1**2*X2**2', '(xyz)^2'),
                    ('X0**3*X2+X0*X1**3+X1*X2**3', 'sum of invariant quartics'), ('X0**3*X2*X0*X1**3*X1*X2**3', 'product of invariant quartics (deg 12)'),
                    ('X0*X1*X2*(X0**3*X2+X0*X1**3+X1*X2**3)', 'xyz * quartic sum (deg 7)'), ('X0**7+X1**7+X2**7', 'x^7+y^7+z^7'),
                    ('X0**14*X1**0+X1**14+X2**14', 'x^14+y^14+z^14 (char 0, deg 14)')]:
        C.append((name, pp(e), 'pure-X'))
    # block pieces of F
    M = core.sub(core.sub(F, pp('X0**5*X1+X1**5*X2+X2**5*X0')), pp('Y0**6+Y1**6+Y2**6'))
    K = pp('X0**5*X1+X1**5*X2+X2**5*X0'); G = pp('Y0**6+Y1**6+Y2**6')
    C += [('M mixed block', M, 'block'), ('K+G = F-M', core.add(K, G), 'block'), ('F-K', core.sub(F, K), 'block'), ('F-G', core.sub(F, G), 'block'),
          ('K-G', core.sub(K, G), 'block'), ('K+2G', core.add(K, core.scale(G, 2)), 'block')]
    # radial / Euler combinations
    EX = core.add(core.add(core.mul(x, Fi[0]), core.mul(y, Fi[1])), core.mul(z, Fi[2]))
    EY = core.add(core.add(core.mul(t, Fi[3]), core.mul(r, Fi[4])), core.mul(s, Fi[5]))
    C += [('E_X = x F_x + y F_y + z F_z', EX, 'radial'), ('E_X - E_Y', core.sub(EX, EY), 'radial'), ('E_X - 2 E_Y', core.sub(EX, core.scale(EY, 2)), 'radial')]
    # coupling quartics and their products with Y-forms
    for e, name in [('X0*X1**3*Y0*Y1', 'coupling x y^3 t r'), ('X0**3*X2*Y1*Y2', 'coupling x^3 z r s'), ('X1*X2**3*Y0*Y2', 'coupling y z^3 t s'),
                    ('X0*X1**3*Y0*Y1+2*X0**3*X2*Y1*Y2+3*X1*X2**3*Y0*Y2', 'M (7/5 scaled)'), ('X0*X1**3*Y0*Y1-X0**3*X2*Y1*Y2', 'coupling difference'),
                    ('X0**3*Y0**3', 'A itself'), ('X0**3*Y0**3+X1**3*Y1**3+X2**3*Y2**3', 'A + rotations (not eigen)'), ('X0**3*Y0**2', 'x^3 t^2'), ('X0**3*Y0', 'x^3 t'),
                    ('X0*Y0', 'x t'), ('X0**2*Y0**2', 'x^2 t^2'), ('X0*X2**2*Y0**6-X1*Y0*Y2*X2**2*Y1**4', 'Jhat-type monomial pair')]:
        C.append((name, pp(e), 'mixed'))
    # Jhat-adjacent: pieces of Jhat
    C += [('Jhat', Jhat, 'companion'), ('7 x z^2 * 9(t^6 - r^6 + s^6) part of Jhat', pp('63*X0*X2**2*(Y0**6-Y1**6+Y2**6)'), 'companion-piece'),
          ('15 y s t (2t^6-4r^6+8s^6) part of Jhat', pp('15*X1*Y2*Y0*(2*Y0**6-4*Y1**6+8*Y2**6)'), 'companion-piece'),
          ('Jhat with F-terms dropped', pp('63*X0*X2**2*(Y0**6-Y1**6+Y2**6)+15*X1*Y2*Y0*(2*Y0**6-4*Y1**6+8*Y2**6)'), 'companion-piece'),
          ('x z^2 (char 2) * G', core.mul(pp('X0*X2**2'), G), 'companion-piece'), ('y t s (char 2) * G', core.mul(pp('X1*Y0*Y2'), G), 'companion-piece')]
    # Hessian minors (2x2 of degree 8, 3x3 blocks of degree 12)
    H = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    def det(Mx):
        n = len(Mx)
        if n == 1: return Mx[0][0]
        acc = {}
        for j in range(n):
            minor = [row[:j] + row[j+1:] for row in Mx[1:]]
            term = core.mul(Mx[0][j], det(minor))
            acc = core.add(acc, term) if j % 2 == 0 else core.sub(acc, term)
        return acc
    for I in [(0, 1, 2), (3, 4, 5), (0, 3, 4), (0, 1, 3)]:
        for J_ in [I]:
            C.append((f'Hessian 3x3 minor rows {I} cols {J_}', det([[H[i][j] for j in J_] for i in I]), 'hessian-minor'))
    C.append(('Hessian 3x3 minor rows (0,1,2) cols (3,4,5)', det([[H[i][j] for j in (3, 4, 5)] for i in (0, 1, 2)]), 'hessian-minor'))
    for (i, j) in [(0, 1), (0, 3), (3, 4), (0, 0), (3, 3), (1, 2)]:
        C.append((f'Hessian 2x2 minor ({i},{j})', det([[H[a][b] for b in (i, j)] for a in (i, j)]) if i != j else core.sub(core.mul(H[i][i], H[(i+1)%6][(i+1)%6]), core.mul(H[i][(i+1)%6], H[(i+1)%6][i])), 'hessian-minor'))
    # products of structured pieces with coordinates/forms (irreducible-only is not needed for the E1 screen)
    C += [('K * G', core.mul(K, G), 'product'), ('xyz * G', core.mul(pp('X0*X1*X2'), G), 'product'), ('trs * K', core.mul(pp('Y0*Y1*Y2'), K), 'product'),
          ('x^3 z * t^3', pp('X0**3*X2*Y0**3'), 'product'), ('(x y^3)(y z^3)', pp('X0*X1**4*X2**3'), 'product'),
          ('Jhat * x', core.mul(Jhat, x), 'product'), ('Jhat * t', core.mul(Jhat, t), 'product'), ('Jhat * Y0^2', core.mul(Jhat, pp('Y0**2')), 'product'),
          ('Jhat * xyz', core.mul(Jhat, pp('X0*X1*X2')), 'product')]
    out = []
    for name, D, fam in C:
        if not D or not core.is_homogeneous(D) or core.parity_of(D) is None:
            say('skip (not homogeneous / mixed parity):', name); continue
        out.append((name, D, fam))
    return out

def task_struct(args):
    name, Dser, fam = args
    D = {tuple(json.loads(k)): Fraction(v) for k, v in Dser.items()}
    d = core.degree(D)
    try:
        if d > 12:
            e1, info = fixed_d.e1_solvable(F, A, q, D, P)
            log({'exp': 'structured-E1', 'D': name, 'd': d, 'par': core.parity_of(D), 'k': char_of(D), 'prime': P, 'family': fam,
                 'result': 'E1-ok (not tested further)' if e1 else 'E1-fail (RIGOROUS: A*D not in J mod p => not in J over Q)', 'job': 'jobD'})
            return name, {'d': d, 'e1': e1}
        e = test_D_restricted(F, A, q, D, P, name, fam, note='structured candidate', extra={'job': 'jobD'})
        res = {k: v for k, v in e.items() if not k.startswith('_')}
        if e.get('result') == 'inconsistent' and d <= 10 and char_of(D) is not None:
            r = exact_dual(F, A, q, D, P); yv = r.pop('y', None)
            fn = f'{OUT}/duals/{sha(name)}.json'
            json.dump({'D': name, 'D_poly': core.to_str(D), 'certificate': r, 'y': yv}, open(fn, 'w'))
            log({'exp': 'exact-dual', 'D': name, 'd': d, 'k': char_of(D), 'par': core.parity_of(D), 'prime': P, 'family': fam,
                 'result': 'RIGOROUS-INCONSISTENT (exact rational dual verified over Q)' if r['ok'] else 'dual-failed', 'rank': r.get('rank_mod_p'),
                 'seconds': r.get('t_solve'), 'job': 'jobD', 'note': json.dumps({k: v for k, v in r.items()})[:500]})
            res['dual'] = r
        return name, res
    except Exception as ex:
        log({'exp': 'fixedD', 'D': name, 'family': fam, 'prime': P, 'result': 'error', 'note': repr(ex), 'job': 'jobD'})
        return name, {'error': repr(ex)}

# ------------------------------------------------------------------ (2) served targets in the character-3 sector
def served_sector3(D, p=P):
    """Space of Y-odd degree-6 targets of sigma-character 3 served by D (mod p), and its dimension modulo J_6^-[3]."""
    Dp = core.reduce_mod(D, p); Fp = core.reduce_mod(F, p)
    d = core.degree(D); pD = core.parity_of(D); k = char_of(D)
    assert k is not None
    Fi_ = [pderiv(Fp, i, p) for i in range(NV)]; Di = [pderiv(Dp, i, p) for i in range(NV)]
    FD = mul(Fp, Dp, p); D2 = mul(Dp, Dp, p); F2 = scale(power(Fp, 2, p), p - 1, p)
    G = [scale(add(mul(Dp, Fi_[i], p), mul(Fp, Di[i], p), p), p - 1, p) for i in range(NV)]
    cols, keys = [], []
    for i in range(NV):
        for m in wmonomials(d + 1, (1 - pD - PI[i]) % 2, (3 + k + W[i]) % 7):
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(col); keys.append(('N', i, m))
    for m in wmonomials(2 * d - 6, 1, (3 + 2 * k) % 7):
        cols.append(shift(F2, m, 1, p)); keys.append(('Q', m))
    Am = wmonomials(6, 1, 3)
    for m in Am:
        cols.append(shift(D2, m, 1, p)); keys.append(('A', m))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    M = dense_cols(cols, rows, p)
    X, nul = M.nullspace()
    a_idx = [j for j, kk in enumerate(keys) if kk[0] == 'A']
    served = []
    if nul:
        proj = nmod_mat(len(a_idx), nul, [int(X[j, c]) for j in a_idx for c in range(nul)], p).transpose()
        Pr, rk = proj.rref()
        served = [[int(Pr[i, j]) for j in range(len(a_idx))] for i in range(rk)]
    # J_6^-[3] inside the A-monomials
    ctx = Ctx(F, A, p)
    jk, jp = jac_cols(ctx, 6, 1, 3)
    aidx = {m: i for i, m in enumerate(Am)}
    jrows = [[P_.get(m, 0) % p for m in Am] for P_ in jp]
    red = SubspaceReducer(jrows, Am, p)
    # served modulo J: rank of reduced served vectors
    if served:
        V = nmod_mat(len(served), len(Am), [v for row in served for v in row], p)
        NF = red.reduce_matrix(V)
        dim_mod_J = NF.rank()
    else:
        dim_mod_J = 0
    # does served contain J_6^-[3] (trivially served)?
    contains_J = True
    if jrows:
        base = nmod_mat(len(served) if served else 0, len(Am), [v for row in served for v in row], p) if served else None
        for jr in jrows:
            if base is None: contains_J = False; break
            aug = nmod_mat(len(served) + 1, len(Am), [v for row in served for v in row] + jr, p)
            if aug.rank() != len(served): contains_J = False; break
    # is A served?  A = X0^3 Y0^3 in A-coords
    Avec = [0] * len(Am); Avec[aidx[(3, 0, 0, 3, 0, 0)]] = 1
    A_served = False
    if served:
        aug = nmod_mat(len(served) + 1, len(Am), [v for row in served for v in row] + Avec, p)
        A_served = aug.rank() == len(served)
    return {'d': d, 'k': k, 'par': pD, 'rows': len(rows), 'cols': len(cols), 'n_A_monomials': len(Am), 'dim_J6_sector3': red.rank,
            'dim_served_sector3': len(served), 'dim_served_mod_J': dim_mod_J, 'contains_trivial_J': contains_J, 'A_served': A_served,
            'served_mod_J_polys': [core.to_str({Am[red.std[j]]: int(NF[i, j]) for j in range(len(red.std)) if int(NF[i, j])}) for i in range(NF.nrows())][:8] if served else []}

# ------------------------------------------------------------------ (3) exact rational kernel reps at d = 9
def exact_kernel_reps(d, par, k):
    """Q-basis of (J:A)_d[par,k] modulo J_d[par,k], LLL-reduced integer representatives, via fmpz nullspace."""
    F7 = core.scale(F, 7); A1 = A
    Fi7 = [core.pderiv(F7, i) for i in range(NV)]
    Dm = wmonomials(d, par, k)
    if not Dm: return []
    T = wmonomials(d + 6, (par + 1) % 2, (k + 3) % 7)
    ridx = {m: i for i, m in enumerate(T)}
    cols = [shift(A1, m, 1) for m in Dm]
    jk = []
    for i in range(NV):
        for m in wmonomials(d + 1, (par + 1 + PI[i]) % 2, (k + 3 + W[i]) % 7):
            cols.append(shift(Fi7[i], m, 1)); jk.append((i, m))
    R, C = len(T), len(cols)
    Mz = fmpz_mat(R, C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            v = Fraction(v); assert v.denominator == 1
            Mz[ridx[e], j] = int(v)
    X, nul = Mz.nullspace()   # X: C x C with first nul columns a basis of the (rational) nullspace (integer)
    nD = len(Dm)
    if nul == 0: return []
    # projection to D-coordinates: rows = nullspace vectors restricted to first nD coords
    proj = fmpz_mat(nul, nD, [int(X[j, c]) for c in range(nul) for j in range(nD)])
    # row space over Q = (J:A)_d sector.  Reduce modulo J_d sector: J-gens as integer rows
    jrows = []
    for i in range(NV):
        for m in wmonomials(d - 5, (par + PI[i]) % 2, (k + W[i]) % 7):
            Pm = shift(Fi7[i], m, 1); jrows.append([int(Fraction(Pm.get(mm, 0))) for mm in Dm])
    # rank computations over Q with fmpq_mat rref
    Jq = fmpq_mat(len(jrows), nD, [fmpq(v) for row in jrows for v in row]) if jrows else None
    rJ = Jq.rref()[1] if Jq is not None else 0
    allrows = jrows + [[int(proj[i, j]) for j in range(nD)] for i in range(nul)]
    Aq = fmpq_mat(len(allrows), nD, [fmpq(v) for row in allrows for v in row])
    Ar, rk = Aq.rref()
    kernel = rk - rJ
    # choose kernel reps: pivot rows among the colon vectors beyond the J-span (greedy, exact)
    reps = []
    cur = jrows[:]
    currank = rJ
    for i in range(nul):
        v = [int(proj[i, j]) for j in range(nD)]
        test = fmpq_mat(len(cur) + 1, nD, [fmpq(w) for row in cur for w in row] + [fmpq(w) for w in v])
        if test.rref()[1] > currank:
            cur.append(v); currank += 1; reps.append(v)
        if len(reps) == kernel: break
    # LLL-reduce the integer lattice spanned by reps (plus nothing else) for small representatives
    if reps:
        L = fmpz_mat(len(reps), nD, [w for row in reps for w in row])
        # also reduce modulo J-lattice directions: add small multiples? keep simple: LLL on reps alone
        Lr = L.lll()
        reps = [[int(Lr[i, j]) for j in range(nD)] for i in range(Lr.nrows())]
    polys = []
    for v in reps:
        g = 0
        for w in v: g = __import__('math').gcd(g, abs(w))
        polys.append({Dm[j]: Fraction(v[j] // g) for j in range(nD) if v[j]})
    return polys

def task_rational_kernel(args):
    d, par, k = args
    t = time.time()
    try:
        polys = exact_kernel_reps(d, par, k)
        out = []
        for i, D in enumerate(polys):
            ok, info = fixed_d.e1_solvable(F, A, q, D, P)
            out.append({'i': i, 'terms': len(D), 'max_abs_coeff': max(abs(int(v)) for v in D.values()), 'E1_ok_mod_p': ok, 'D': core.to_str(D)})
        log({'exp': 'rational-kernel-reps', 'd': d, 'par': par, 'k': k, 'result': {'kernel_dim_over_Q': len(polys)}, 'seconds': round(time.time() - t, 1), 'job': 'jobD',
             'note': 'exact Q-nullspace (fmpz_mat) of [A*m | m*F_i]; representatives modulo J_d, LLL-reduced'})
        return (d, par, k), out
    except Exception as ex:
        log({'exp': 'rational-kernel-reps', 'd': d, 'par': par, 'k': k, 'result': 'error', 'note': repr(ex), 'job': 'jobD'})
        return (d, par, k), {'error': repr(ex)}

if __name__ == '__main__':
    # (3) rational kernel reps first (cheap), then their exact duals
    secs = [(9, par, k) for par in (0, 1) for k in range(7)]
    ratreps = {}
    with mp.Pool(min(NW, 14)) as pool:
        for key, out in pool.imap_unordered(task_rational_kernel, secs):
            ratreps[f'{key[0]},{key[1]},{key[2]}'] = out
            say('rational kernel', key, (len(out) if isinstance(out, list) else out))
            json.dump(ratreps, open(f'{OUT}/rational_kernel_reps_d9.json', 'w'), indent=1)
    cands = []
    for key, out in ratreps.items():
        if isinstance(out, list):
            for o in out:
                cands.append((f'rational E1-kernel rep (LLL) d=9 (par,k)=({key.split(",")[1]},{key.split(",")[2]}) #{o["i"]}', core.parse(o['D']), 'e1-kernel-rep-rational-lll'))
    # (1) structured candidates
    cands += structured_candidates()
    say(len(cands), 'candidates (rational kernel reps + structured)')
    ctasks = [(n, {json.dumps(list(m)): str(v) for m, v in D.items()}, fam) for n, D, fam in cands]
    results = {}
    with mp.Pool(NW) as pool:
        for name, res in pool.imap_unordered(task_struct, ctasks):
            results[name] = res
            say('cand', name, res.get('result') or res.get('e1'), 'rank', res.get('rank'), 'dual', (res.get('dual') or {}).get('ok'), res.get('error', ''))
            json.dump(results, open(f'{OUT}/structured_results.json', 'w'), indent=1, default=str)
    say('STRUCTURED_DONE')
    # (2) served targets sector 3
    served = {}
    keyD = [('Jhat', Jhat), ('F_0', Fi[0]), ('F_1', Fi[1]), ('F_2', Fi[2]), ('F_3', Fi[3]), ('F_4', Fi[4]), ('F_5', Fi[5]),
            ('X0*F_0', core.mul(x, Fi[0])), ('Y0*F_3', core.mul(t, Fi[3])), ('X2*F_2', core.mul(z, Fi[2])), ('Y0*F_2', core.mul(t, Fi[2])),
            ('E_X', core.add(core.add(core.mul(x, Fi[0]), core.mul(y, Fi[1])), core.mul(z, Fi[2]))), ('F_0*F_1', core.mul(Fi[0], Fi[1])), ('F_3*F_4', core.mul(Fi[3], Fi[4]))]
    for key, out in ratreps.items():
        if isinstance(out, list) and out:
            keyD.append((f'rational kernel rep d=9 ({key}) #0', core.parse(out[0]['D'])))
    for name, D in keyD:
        t = time.time()
        try:
            r = served_sector3(D)
            served[name] = r
            say('served3', name, r)
            log({'exp': 'served-targets-sector3', 'D': name, 'd': r['d'], 'k': r['k'], 'par': r['par'], 'prime': P, 'result': {kk: r[kk] for kk in ('dim_served_sector3', 'dim_served_mod_J', 'contains_trivial_J', 'A_served', 'dim_J6_sector3')},
                 'seconds': round(time.time() - t, 1), 'job': 'jobD', 'note': 'targets A\' in S_6 Y-odd of sigma-character 3 served by D (mod p); dimension modulo J_6; corrected version of wave-1 served-targets (which ran in character 4)'})
        except Exception as ex:
            say('served3 error', name, ex)
            log({'exp': 'served-targets-sector3', 'D': name, 'result': 'error', 'note': repr(ex), 'job': 'jobD'})
        json.dump(served, open(f'{OUT}/served_sector3.json', 'w'), indent=1, default=str)
    # (4) Singular E1 screen for high-degree candidates (mod p): reduce(A*D, std(J)) == 0 ?
    sing = f'''ring R = {P}, (X0,X1,X2,Y0,Y1,Y2), dp;
poly F = {Fstr.replace('**','^')};
poly A = X0^3*Y0^3;
ideal J = jacob(F);
ideal SJ = std(J);
int vd = vdim(SJ);
"vdim", vd;
poly Peps1 = X0^14*Y0^7 - (X1*X2)^7*(Y0^7+Y1^7+Y2^7) + X0^3*Y0^3*X1^7*X2^8;
poly Peps0 = X0^14*Y0^7 - (X1*X2)^7*(Y0^7+Y1^7+Y2^7);
matrix H[6][6] = jacob(jacob(F));
poly HessF = det(H);
poly K = X0^5*X1+X1^5*X2+X2^5*X0; poly G = Y0^6+Y1^6+Y2^6;
list L = list("P_eps=1 (deg 21)", Peps1, "P_eps=0 (deg 21)", Peps0, "Hess(F) (deg 24)", HessF, "K*G (deg 12)", K*G, "K^2 (deg 12)", K^2, "G^2", G^2, "K*G*X0*X1*X2 (deg 15)", K*G*X0*X1*X2, "F_0*F_1*F_2 (deg 15)", diff(F,X0)*diff(F,X1)*diff(F,X2), "F_3*F_4*F_5 (deg 15)", diff(F,Y0)*diff(F,Y1)*diff(F,Y2), "Peps1*X0 (deg 22)", Peps1*X0, "Peps1*Y0 (deg 22)", Peps1*Y0);
int i;
for (i = 1; i <= size(L); i = i + 2) {{ poly r = reduce(A*L[i+1], SJ); "E1", L[i], "deg", deg(L[i+1]), "A*D in J:", (r == 0); }}
quit;
'''
    open(f'{OUT}/e1_screen.sing', 'w').write(sing)
    t = time.time()
    try:
        rr = subprocess.run(['Singular', '-q', f'{OUT}/e1_screen.sing'], capture_output=True, text=True, timeout=7200)
        open(f'{OUT}/e1_screen.out', 'w').write(rr.stdout + '\n' + rr.stderr)
        say('singular E1 screen:\n' + rr.stdout[-2000:])
        log({'exp': 'singular-E1-screen', 'prime': P, 'result': rr.stdout.strip().splitlines()[-14:], 'seconds': round(time.time() - t, 1), 'job': 'jobD',
             'note': 'reduce(A*D, std(jacob(F))) == 0 mod p; a nonzero remainder is a rigorous exclusion (A*D not in J over Q); vdim must be 15625'})
    except Exception as ex:
        say('singular error', ex)
    say('JOBD_DONE')
