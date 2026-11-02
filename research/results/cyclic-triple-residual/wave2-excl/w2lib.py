"""Wave-2 exclusion library for cyclic-triple-residual (q = 1, A = X0^3 Y0^3).

Jacobian-type reduction.  E1 (sum_i N_i F_i = A D) says [A][D] = 0 in R = S/J(F).  In degrees d where
Ann([A])_d = 0 (d <= 7 here) every denominator is D = G(F) = sum_i G_i F_i with G a vector of forms of degree
d - 5, and then N = A G + S.gradF with S antisymmetric (Koszul).  Reducing E2 (D div N - N(D) = F Q) modulo
J(F) kills every term except  A * G^T Hess(F) G, so a necessary condition ("Hessian relaxation") is

        A * sum_{i,k} G_i G_k F_{ik}  in  J(F)_{2d}.                                   (d <= 7)

For d = 8 the E1-admissible denominators are D = c D8 + G(F) (D8 = the unique sporadic annihilator) with
N = c P + A G + S.gradF where P(F) = A D8, and the relaxation reads (mod J_16)

   c^2 [D8 div P - P(D8)] + c [D8 (A div G + G(A)) - P^T Hess(F) G - A G(D8)] - A G^T Hess(F) G  = 0.

Illegal solutions: D = h F  <=>  G = h E (E = Euler field), c = 0.

All normal forms are computed modulo a prime p (sector by sector under the order-18 diagonal group); the
quadrics are therefore modular objects; exact certificates are produced separately (see square_certificate).
"""
import sys, time, json, itertools, math, random, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from fractions import Fraction
from hodge import core
from sect import elems, semi_invariant_subgroup, generators, make_chi, cadd, chi_of, imul, ideriv, mons, E6, integral
from flint import nmod_mat, fmpq_mat, fmpz_mat, fmpq, fmpz

WD = '/tmp/hodge/results/cyclic-triple-residual'
PK = WD + '/wave2-excl'
LOG = WD + '/wave2-excl.log.jsonl'

def log(entry):
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    entry.setdefault('engine', 2)
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')

class Setup:
    def __init__(self):
        F, A, q, rec = core.load_instance('cyclic-triple-residual')
        if os.environ.get('W2_FERMAT'):      # control: Fermat sextic (t = 0), where join-type denominators exist
            F = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
        self.F, self.A, self.q = F, A, q
        self.Fz, self.Az = integral(F), integral(A)          # Fz = 7 F (integer), Az = A
        self.G = semi_invariant_subgroup(elems(self.Fz), [self.Az])
        self.gens = generators(self.G)
        self.chi = make_chi(self.gens)
        self.chiA = chi_of(self.Az, self.chi)
        self.chiE = [self.chi(E6[i]) for i in range(6)]
        self.dF = [ideriv(self.Fz, i) for i in range(6)]
        self.ddF = [[ideriv(self.dF[i], k) for k in range(6)] for i in range(6)]
    def cneg(self, c):
        return tuple((-x) % 6 for x in c)
    def chi_poly(self, P):
        return chi_of(P, self.chi)

def mono_poly(m, c=1):
    return {tuple(m): c}

def jac_gens(S, n, c):
    """Integer polynomials m * Fz_i spanning J(Fz)_n in character sector c."""
    out = []
    if n < 5: return out
    for i in range(6):
        ci = cadd(c, S.chiE[i])
        for m in mons(n - 5):
            if S.chi(m) == ci: continue_ = False
            if S.chi(m) != ci: continue
            out.append(imul({m: 1}, S.dF[i]))
    return out

class Reducer:
    """Row-reduced J(F)_n in one character sector, mod p.  nf(V): normal forms of rows of an nmod_mat."""
    def __init__(self, S, n, c, p):
        self.p, self.n, self.c = p, n, c
        self.mons = [m for m in mons(n) if S.chi(m) == c]
        self.idx = {m: i for i, m in enumerate(self.mons)}
        gens = jac_gens(S, n, c)
        C = len(self.mons)
        self.ngens = len(gens)
        if gens:
            flat = [0] * (len(gens) * C)
            for r, P in enumerate(gens):
                for e, v in P.items():
                    flat[r * C + self.idx[e]] = v % p
            M = nmod_mat(len(gens), C, flat, p)
            Mr, rk = M.rref()
            self.rank = rk
            piv = []
            for i in range(rk):
                j = 0
                while j < C and int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            self.piv = piv
            self.Rr = nmod_mat(rk, C, [int(Mr[i, j]) for i in range(rk) for j in range(C)], p) if rk else None
        else:
            self.rank = 0; self.piv = []; self.Rr = None
        ps = set(self.piv)
        self.standard = [j for j in range(C) if j not in ps]
        self.dimR = len(self.standard)
    def vec(self, P):
        v = [0] * len(self.mons)
        for e, c in P.items():
            v[self.idx[e]] = (v[self.idx[e]] + c) % self.p
        return v
    def nf_rows(self, rows):
        """rows: list of coefficient vectors (mod p).  Returns list of vectors restricted to standard monomials."""
        if not rows: return []
        k, C, p = len(rows), len(self.mons), self.p
        V = nmod_mat(k, C, [x % p for r in rows for x in r], p)
        if self.rank:
            Vp = nmod_mat(k, self.rank, [int(V[i, j]) for i in range(k) for j in self.piv], p)
            V = V - Vp * self.Rr
        return [[int(V[i, j]) for j in self.standard] for i in range(k)]

def theta_vars(d):
    """theta coordinates of the E1-admissible family in degree d: list of ('c',) or (i, m)."""
    th = []
    if d == 8:
        th.append(('c',))
    for i in range(6):
        for m in mons(d - 5):
            th.append((i, m))
    return th

def illegal_forms(d, th):
    """Linear forms (dict var-index -> coeff) vanishing exactly on the illegal locus {c = 0, G = h E}."""
    idx = {t: k for k, t in enumerate(th)}
    forms = []
    if d == 8:
        forms.append({idx[('c',)]: 1})
    for i in range(6):
        for m in mons(d - 5):
            if m[i] == 0:
                forms.append({idx[(i, m)]: 1})
    if d >= 6:
        for m in mons(d - 6):
            e0 = list(m); e0[0] += 1; e0 = tuple(e0)
            for j in range(1, 6):
                ej = list(m); ej[j] += 1; ej = tuple(ej)
                forms.append({idx[(0, e0)]: 1, idx[(j, ej)]: -1})
    return forms

def relaxation_terms(S, d, P=None, D8=None):
    """Dict (a, b) [a <= b, theta indices] -> integer polynomial coefficient of theta_a theta_b in the relaxation
    polynomial (degree 2d).  For d <= 7: A G^T H G.  For d = 8 the full expression in the docstring (P = E1
    solution with P(F) = A D8, integer polys; note Fz = 7F so P here must satisfy sum P_i dFz_i = Az D8z)."""
    th = theta_vars(d)
    terms = {}
    Az = S.Az
    def addt(key, P_):
        if not P_: return
        if key in terms:
            terms[key] = {e: v for e, v in ((e, terms[key].get(e, 0) + c) for e, c in P_.items()) if v} | {e: v for e, v in terms[key].items() if e not in P_}
        else:
            terms[key] = dict(P_)
    # G^T H G part: sum_{a,b} theta_a theta_b m_a m_b F_{i_a i_b}; sign -1 for d = 8 (see docstring)
    sg = -1 if d == 8 else 1
    gi = [(k, t) for k, t in enumerate(th) if t[0] != 'c']
    for x in range(len(gi)):
        ka, (ia, ma) = gi[x]
        for y in range(x, len(gi)):
            kb, (ib, mb) = gi[y]
            Hik = S.ddF[ia][ib]
            if not Hik: continue
            P_ = imul(imul({ma: 1}, {mb: 1}), Hik)
            P_ = imul(Az, P_)
            mult = sg * (1 if x == y else 2)
            P_ = {e: v * mult for e, v in P_.items()}
            addt((ka, kb), P_)
    if d == 8:
        assert P is not None and D8 is not None
        kc = 0
        D8z = D8
        divP = {}
        for i in range(6):
            divP = core.add(divP, ideriv(P[i], i))
        PD8 = {}
        for i in range(6):
            PD8 = core.add(PD8, imul(P[i], ideriv(D8z, i)))
        addt((kc, kc), core.sub(imul(D8z, divP), PD8))
        dA = [ideriv(Az, i) for i in range(6)]
        dD8 = [ideriv(D8z, i) for i in range(6)]
        for ka, (ia, ma) in gi:
            # D8 (A d_i m + m d_i A) - A m d_i D8 - sum_k P_k F_{k i} m
            t1 = imul(D8z, core.add(imul(Az, ideriv({ma: 1}, ia)), imul({ma: 1}, dA[ia])))
            t2 = imul(Az, imul({ma: 1}, dD8[ia]))
            t3 = {}
            for k in range(6):
                t3 = core.add(t3, imul(P[k], imul(S.ddF[k][ia], {ma: 1})))
            addt((kc, ka), core.sub(core.sub(t1, t2), t3))
    return th, terms

def relaxation_quadrics(S, d, p, P=None, D8=None, verbose=True):
    """Return (th, quadrics) where quadrics is a list of dicts {(a,b): coeff mod p}: the coordinates of the
    normal form (mod J(F)_{2d}, sector by sector) of the relaxation polynomial."""
    t0 = time.time()
    th, terms = relaxation_terms(S, d, P, D8)
    n = 2 * d
    # group term keys by character
    bych = {}
    for key, P_ in terms.items():
        c = chi_of(P_, S.chi)
        bych.setdefault(c, []).append(key)
    quadrics = []
    info = {}
    for c, keys in sorted(bych.items()):
        red = Reducer(S, n, c, p)
        rows = [red.vec(terms[k]) for k in keys]
        nfs = red.nf_rows(rows)
        # transpose: for each standard monomial s, quadric sum_k nfs[k][s] theta_key
        for s in range(red.dimR):
            qd = {}
            for k, key in enumerate(keys):
                v = nfs[k][s]
                if v: qd[key] = v
            if qd: quadrics.append(qd)
        info[str(c)] = {'mons': len(red.mons), 'J_rank': red.rank, 'dimR': red.dimR, 'terms': len(keys)}
    if verbose:
        print(f'[relax d={d} p={p}] theta={len(th)} term-pairs={len(terms)} quadrics={len(quadrics)} '
              f'dimR_total={sum(v["dimR"] for v in info.values())} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
    return th, quadrics, info

def quad_str(qd, names, p=None):
    parts = []
    for (a, b), v in sorted(qd.items()):
        if p is not None:
            v = v % p
            if v > p // 2: v -= p
        if v == 0: continue
        mon = f'{names[a]}^2' if a == b else f'{names[a]}*{names[b]}'
        parts.append(f'{"+" if v > 0 else "-"}{abs(v)}*{mon}')
    s = ''.join(parts)
    return s.lstrip('+') if s else '0'

def lin_str(form, names):
    parts = []
    for k, v in sorted(form.items()):
        parts.append(f'{"+" if v > 0 else "-"}{abs(v)}*{names[k]}')
    return ''.join(parts).lstrip('+')

def write_msolve(path, nvars, p, polys, names=None):
    names = names or [f't{k}' for k in range(nvars)]
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n')
        f.write(f'{p}\n')
        f.write(',\n'.join(polys) + '\n')
    return names

def write_singular(path, nvars, p, polys, names=None, extra=''):
    names = names or [f't{k}' for k in range(nvars)]
    with open(path, 'w') as f:
        f.write(f'ring r = {p}, ({",".join(names)}), dp;\n')
        f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
        f.write(extra)
    return names

def run_msolve(inp, out, threads=1, gb=0, timeout=None, extra=()):
    import subprocess
    cmd = ['msolve', '-f', inp, '-o', out, '-t', str(threads), '-g', str(gb), *extra]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return r.returncode, time.time() - t0, r.stdout[-2000:], r.stderr[-2000:]

def parse_msolve_gb(path):
    """Return list of polynomial strings from a -g 2 output file."""
    txt = '\n'.join(l for l in open(path).read().splitlines() if not l.startswith('#'))
    txt = txt.strip()
    if txt.startswith('['): txt = txt[1:]
    if txt.endswith(']:'): txt = txt[:-2]
    elif txt.endswith(']'): txt = txt[:-1]
    polys = [s.strip() for s in txt.split(',\n') if s.strip()]
    return polys

# ---------------------------------------------------------------------------------------------
# Sealed fixed-D exclusion: modular ranks + exact rational null vectors (modular-rank argument)
# ---------------------------------------------------------------------------------------------
def sealed_exclusion(F, A, q, D, primes=(1000003, 2147483587), max_primes=4, verbose=False, exact_fallback=True):
    """Sector-reduced fixed-D system (sect.Sys).  Returns a dict with modular ranks at two primes, the number of
    exactly verified rational null vectors (lower bound on nullity_Q => upper bound on rank_Q), and the conclusion:
      sealed inconsistent over Q  iff  verified_null == nullity_p  and  rank_aug_p > rank_p.
    (rank_p <= rank_Q <= cols - verified_null; rank_aug_Q >= rank_aug_p.)"""
    import sect
    t0 = time.time()
    Sy = sect.Sys(F, A, q, D, verbose=verbose)
    p1, p2 = primes[0], primes[1]
    def null_mod(p):
        M = Sy.dense(p, aug=False)
        X, nullity = M.nullspace()
        vs = [[int(X[j, k]) for j in range(Sy.ncols)] for k in range(nullity)]
        if vs:
            Mv = nmod_mat(len(vs), Sy.ncols, [x for v in vs for x in v], p)
            R, rk = Mv.rref()
            vs = [[int(R[i, j]) for j in range(Sy.ncols)] for i in range(rk)]
        return vs, nullity
    v1, n1 = null_mod(p1); v2, n2 = null_mod(p2)
    ra1 = Sy.dense(p1, aug=True).rank(); ra2 = Sy.dense(p2, aug=True).rank()
    out = {'rows': Sy.nrows, 'cols': Sy.ncols, 'nN': Sy.nN, 'nQ': Sy.nQ, 'G_order': Sy.G_order, 'd': Sy.d,
           'primes': [p1, p2], 'nullity': [n1, n2], 'rank': [Sy.ncols - n1, Sy.ncols - n2], 'rank_aug': [ra1, ra2]}
    sealed = 0
    if n1 == n2 and len(v1) == len(v2):
        mod = p1 * p2
        for a, b in zip(v1, v2):
            xs = [core.crt_pair(x, p1, y, p2)[0] for x, y in zip(a, b)]
            rat = [core.ratrecon(x, mod) for x in xs]
            if any(r is None for r in rat): continue
            L = 1
            for r in rat: L = L * r.denominator // math.gcd(L, r.denominator)
            xi = [int(r * L) for r in rat]
            acc = {}
            for j, col in enumerate(Sy.cols):
                if xi[j] == 0: continue
                for r_, v in col.items(): acc[r_] = acc.get(r_, 0) + v * xi[j]
            if all(v == 0 for v in acc.values()): sealed += 1
    out['verified_null'] = sealed
    out['consistent_mod_p'] = (ra1 == Sy.ncols - n1) or (ra2 == Sy.ncols - n2)
    out['sealed_inconsistent'] = (sealed == n1) and (ra1 > Sy.ncols - n1) and (ra2 > Sy.ncols - n2)
    out['method'] = 'modular-rank+exact-null-vectors'
    maxD = max((abs(Fraction(v).numerator) for v in D.values()), default=0)
    out['max_abs_D_coeff_digits'] = len(str(maxD))
    if not out['sealed_inconsistent'] and not out['consistent_mod_p'] and exact_fallback and Sy.nrows * Sy.ncols <= 12_000_000 and maxD < 10**15:
        # exact fallback: fmpz ranks over Q of the system and of the augmented system
        t1 = time.time()
        rQ = Sy.dense(None, aug=False).rank(); raQ = Sy.dense(None, aug=True).rank()
        out.update(rank_exact=rQ, rank_aug_exact=raQ, t_exact=round(time.time() - t1, 1), method='exact-fmpz-rank')
        out['sealed_inconsistent'] = raQ > rQ
        out['consistent_exact'] = (raQ == rQ)
    out['seconds'] = round(time.time() - t0, 1)
    return out

# ---------------------------------------------------------------------------------------------
# Exact E1-kernel elements by sector (multi-prime CRT + rational reconstruction + exact E1 solve)
# ---------------------------------------------------------------------------------------------
def kernel_canonical(S, d, c, p):
    """Canonical (rref) basis mod p of {D in S_d^c : A D in J_{d+6}} / J_d^c, as coefficient rows over the sector
    monomials (non-pivot-of-J coordinates only are nonzero at pivot positions... we return full rows)."""
    Dm = [m for m in mons(d) if S.chi(m) == c]; idx = {m: i for i, m in enumerate(Dm)}
    target = cadd(S.chiA, c)
    cols = []
    for m in Dm:
        cols.append({tuple(a + b for a, b in zip(m, e)): v % p for e, v in S.Az.items()})
    nA = len(cols)
    for i in range(6):
        ci = cadd(target, S.chiE[i])
        for m in mons(d + 1):
            if S.chi(m) != ci: continue
            col = {}
            for e, v in S.dF[i].items():
                k = tuple(a + b for a, b in zip(m, e)); col[k] = (col.get(k, 0) + v) % p
            cols.append(col)
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    if nullity == 0: return Dm, [], 0
    proj = nmod_mat(nullity, nA, [int(X[j, k]) for k in range(nullity) for j in range(nA)], p)
    Pr, kA = proj.rref()
    vecs = [[int(Pr[i, j]) for j in range(nA)] for i in range(kA)]
    # reduce modulo J_d^c
    red = Reducer(S, d, c, p)
    assert red.mons == Dm
    nfs = red.nf_rows(vecs) if vecs else []
    # nfs are coordinates on standard monomials; embed back and rref
    std = red.standard
    full = []
    for v in nfs:
        row = [0] * nA
        for s, j in enumerate(std): row[j] = v[s]
        if any(row): full.append(row)
    if not full: return Dm, [], red.rank
    Mf = nmod_mat(len(full), nA, [x for r in full for x in r], p)
    Rf, rk = Mf.rref()
    out = [[int(Rf[i, j]) for j in range(nA)] for i in range(rk)]
    return Dm, out, red.rank

def kernel_exact(S, d, c, primes, verbose=False):
    """Exact rational kernel basis in sector c (CRT over primes + ratrecon + exact E1 solve).
    Returns list of dicts {Kz, Pz, L, verified} (integer polys; sum Pz_i dFz_i = Az Kz) or None entries on failure."""
    t0 = time.time()
    bases = []
    for p in primes:
        Dm, vecs, rJ = kernel_canonical(S, d, c, p)
        bases.append((p, vecs))
    k = len(bases[0][1])
    if any(len(b[1]) != k for b in bases):
        return None, {'error': 'kernel dimension differs between primes', 'dims': [len(b[1]) for b in bases]}
    if k == 0: return [], {'dim': 0}
    nA = len(Dm)
    # CRT across primes
    elems_ = []
    for r in range(k):
        xs = list(bases[0][1][r]); mod = bases[0][0]
        for p, vecs in bases[1:]:
            xs = [core.crt_pair(x, mod, y, p)[0] for x, y in zip(xs, vecs[r])]; mod *= p
        rat = [core.ratrecon(x, mod) for x in xs]
        if any(x is None for x in rat):
            elems_.append(None); continue
        L = 1
        for x in rat: L = L * x.denominator // math.gcd(L, x.denominator)
        Kz = {Dm[j]: int(rat[j] * L) for j in range(nA) if rat[j]}
        elems_.append(Kz)
    # exact E1 solve for all elements at once: sum P_i dFz_i = Az Kz  (rows: S_{d+6}^{target})
    target = cadd(S.chiA, c)
    mT = [m for m in mons(d + 6) if S.chi(m) == target]; iT = {m: i for i, m in enumerate(mT)}
    cols = []; meta = []
    for i in range(6):
        ci = cadd(target, S.chiE[i])
        for m in mons(d + 1):
            if S.chi(m) != ci: continue
            cols.append(imul({m: 1}, S.dF[i])); meta.append((i, m))
    good = [e for e in elems_ if e is not None]
    R, C, K = len(mT), len(cols), len(good)
    flat = [0] * (R * (C + K))
    for j, P_ in enumerate(cols):
        for e, v in P_.items(): flat[iT[e] * (C + K) + j] += v
    for kk, Kz in enumerate(good):
        for e, v in imul(S.Az, Kz).items(): flat[iT[e] * (C + K) + C + kk] += v
    Maug = fmpq_mat(fmpz_mat(R, C + K, flat))
    Mr, rk = Maug.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < C + K and Mr[i, j] == 0: j += 1
        piv.append(j)
    rankM = sum(1 for j in piv if j < C)
    results = []
    gi = 0
    for Kz in elems_:
        if Kz is None:
            results.append(None); continue
        kk = gi; gi += 1
        # consistent iff no pivot in column C + kk ... (pivots among rhs columns indicate inconsistency of that rhs
        # only if the pivot row has zeros in the other rhs columns; simpler: solve individually using the rref of M)
        results.append({'Kz': Kz})
    # individual exact solves (robust): reuse rref of [M | all rhs]: a rhs column is consistent iff it is not a pivot
    # column AND it has no nonzero entry in rows beyond rankM where all M-columns are zero. Do the simple thing:
    for kk, res in enumerate([r for r in results if r is not None]):
        colj = C + kk
        # rows i >= rankM of rref have zero M-part; consistency iff Mr[i, colj] == 0 for all such i
        ok = all(Mr[i, colj] == 0 for i in range(rankM, rk))
        if not ok:
            res['verified'] = False; res['Pz'] = None; continue
        xs = [Fraction(0)] * C
        for i, j in enumerate(piv[:rankM]):
            v = Mr[i, colj]; xs[j] = Fraction(int(v.p), int(v.q))
        L = 1
        for x in xs: L = L * x.denominator // math.gcd(L, x.denominator)
        Pz = [dict() for _ in range(6)]
        for (i, m), x in zip(meta, xs):
            if x: Pz[i][m] = int(x * L)
        Kzz = {e: v * L for e, v in res['Kz'].items()}
        lhs = {}
        for i in range(6): lhs = core.add(lhs, imul(Pz[i], S.dF[i]))
        res['verified'] = (lhs == imul(S.Az, Kzz))
        res['Kz'] = Kzz; res['Pz'] = Pz; res['L'] = L
    info = {'dim': k, 'reconstructed': sum(1 for r in results if r is not None),
            'verified': sum(1 for r in results if r is not None and r.get('verified')), 'seconds': round(time.time() - t0, 1),
            'e1_rows': R, 'e1_cols': C}
    return results, info

# ---------------------------------------------------------------------------------------------
# Exact E1-kernel elements via exact integer nullspace + LLL (small-height representatives)
# ---------------------------------------------------------------------------------------------
def integer_kernel_lll(Cmat_rows, n):
    """Saturated, LLL-reduced Z-basis of {x in Z^n : C x = 0} for an integer matrix C (list of rows of length n),
    via LLL of the lattice spanned by the rows of [I_n | N * C^T]."""
    m = len(Cmat_rows)
    if m == 0:
        return [[int(i == j) for j in range(n)] for i in range(n)]
    maxc = max((abs(v) for r in Cmat_rows for v in r), default=1)
    N = 2 ** (maxc.bit_length() + n.bit_length() + 40)
    rows = []
    for i in range(n):
        rows.append([int(i == j) for j in range(n)] + [N * Cmat_rows[r][i] for r in range(m)])
    B = fmpz_mat(rows).lll()
    out = []
    for i in range(B.nrows()):
        row = [int(B[i, j]) for j in range(n + m)]
        if all(v == 0 for v in row[n:]) and any(row[:n]):
            out.append(row[:n])
    return out

def kernel_exact_lll(S, d, c, verbose=False):
    """Exact basis (small integer coefficients) of {D in S_d^c : A D in J_{d+6}} modulo J_d^c, with exact E1 lifts.
    Steps: exact rational nullspace of [A-columns | J_{d+6}-columns] (fmpz_mat), projection to the A-coordinates
    (rational kernel K_Q of [A] on S_d^c), saturation + LLL of K_Q cap Z^{nA} via the LLL kernel trick, choice of an
    independent subset modulo J_d^c (exact), exact E1 solve for the lifts.  Returns (list of {'Kz','Pz','L','verified'}, info)."""
    t0 = time.time()
    Dm = [m for m in mons(d) if S.chi(m) == c]; nA = len(Dm)
    target = cadd(S.chiA, c)
    cols = [imul(S.Az, {m: 1}) for m in Dm]
    for i in range(6):
        ci = cadd(target, S.chiE[i])
        for m in mons(d + 1):
            if S.chi(m) != ci: continue
            cols.append(imul({m: 1}, S.dF[i]))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] += v
    X, nullity = fmpz_mat(R, C, flat).nullspace()
    t1 = time.time()
    proj = [[int(X[j, k]) for j in range(nA)] for k in range(nullity)]
    proj = [r for r in proj if any(r)]
    if not proj:
        return [], {'dim': 0, 'seconds': round(time.time() - t0, 1)}
    # orthogonal complement of the projected kernel in Q^nA, as integer rows
    Y, ny = fmpz_mat(proj).nullspace()          # columns of Y: vectors y with proj * y = 0
    Crows = [[int(Y[i, k]) for i in range(nA)] for k in range(ny)]
    Ksat = integer_kernel_lll(Crows, nA)         # saturated LLL-short basis of the rational kernel K_Q cap Z^nA
    kA = len(Ksat)
    # exact reduction modulo J_d^c: pick an independent subset modulo J
    gensJ = jac_gens(S, d, c)
    idx = {m: i for i, m in enumerate(Dm)}
    Jrows = []
    for P_ in gensJ:
        r_ = [0] * nA
        for e, v in P_.items(): r_[idx[e]] += v
        Jrows.append(r_)
    rJ = fmpz_mat(Jrows).rank() if Jrows else 0
    chosen = []
    cur = list(Jrows); currank = rJ
    for v in Ksat:
        if currank - rJ >= kA - rJ: break
        rk = fmpz_mat(cur + [v]).rank()
        if rk > currank:
            chosen.append(v); cur.append(v); currank = rk
    elems_ = []
    for row in chosen:
        g = 0
        for x in row: g = math.gcd(g, abs(x))
        if g > 1: row = [x // g for x in row]
        elems_.append({Dm[j]: row[j] for j in range(nA) if row[j]})
    t2 = time.time()
    # exact E1 lifts (same as kernel_exact)
    mT = [m for m in mons(d + 6) if S.chi(m) == target]; iT = {m: i for i, m in enumerate(mT)}
    colsE = []; metaE = []
    for i in range(6):
        ci = cadd(target, S.chiE[i])
        for m in mons(d + 1):
            if S.chi(m) != ci: continue
            colsE.append(imul({m: 1}, S.dF[i])); metaE.append((i, m))
    RE, CE, K = len(mT), len(colsE), len(elems_)
    flatE = [0] * (RE * (CE + K))
    for j, P_ in enumerate(colsE):
        for e, v in P_.items(): flatE[iT[e] * (CE + K) + j] += v
    for kk, Kz in enumerate(elems_):
        for e, v in imul(S.Az, Kz).items(): flatE[iT[e] * (CE + K) + CE + kk] += v
    Mr, rkE = fmpq_mat(fmpz_mat(RE, CE + K, flatE)).rref()
    piv = []
    for i in range(rkE):
        j = 0
        while j < CE + K and Mr[i, j] == 0: j += 1
        piv.append(j)
    rankM = sum(1 for j in piv if j < CE)
    results = []
    for kk, Kz in enumerate(elems_):
        colj = CE + kk
        ok = all(Mr[i, colj] == 0 for i in range(rankM, rkE))
        if not ok:
            results.append({'Kz': Kz, 'Pz': None, 'verified': False}); continue
        xs = [Fraction(0)] * CE
        for i, j in enumerate(piv[:rankM]):
            v = Mr[i, colj]; xs[j] = Fraction(int(v.p), int(v.q))
        L = 1
        for x in xs: L = L * x.denominator // math.gcd(L, x.denominator)
        Pz = [dict() for _ in range(6)]
        for (i, m), x in zip(metaE, xs):
            if x: Pz[i][m] = int(x * L)
        Kzz = {e: v * L for e, v in Kz.items()}
        lhs = {}
        for i in range(6): lhs = core.add(lhs, imul(Pz[i], S.dF[i]))
        results.append({'Kz': Kzz, 'Pz': Pz, 'L': L, 'verified': (lhs == imul(S.Az, Kzz))})
    info = {'dim': len(elems_), 'kerA_in_S': kA, 'verified': sum(1 for r in results if r['verified']), 'rank_J_d': rJ,
            'max_abs_coeff': max((abs(v) for r in results for v in r['Kz'].values()), default=0),
            't_nullspace': round(t1 - t0, 1), 't_lll': round(t2 - t1, 1), 'seconds': round(time.time() - t0, 1), 'e1_rows': R, 'e1_cols': C}
    if verbose: print(f'[kernel_exact_lll d={d} c={c}] {info}', file=sys.stderr, flush=True)
    return results, info
