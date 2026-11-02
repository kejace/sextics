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
    txt = open(path).read()
    txt = txt.strip()
    if txt.startswith('['): txt = txt[1:]
    if txt.endswith(']:'): txt = txt[:-2]
    elif txt.endswith(']'): txt = txt[:-1]
    polys = [s.strip() for s in txt.split(',\n') if s.strip()]
    return polys
