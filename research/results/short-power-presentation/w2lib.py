"""Wave-2 library for the short-power-presentation lane (q = 1, A = X0^3 Y0^3).

Exact modular linear algebra (python-flint) for:
  * JRef: row-reduced J(F)_e^{par} (generators m*F_j as rows over the monomial basis), normal forms,
  * e1_kernel: K_d^{par} = {D in S_d^{par} : A*D in J(F)_{d+6}} via the pivot rows (A is a monomial),
    reduced modulo J(F)_d, with multi-prime rational reconstruction,
  * particular numerators P (P(F) = A*G) for E1-admissible generators,
  * the E2-mod-J quadratic relaxation  sum_{k,l} a_k a_l T_kl in J_{2d},  T_kl = G_k div P_l - P_l(G_k),
    for a family of generators G_k (kernel representatives + Jacobian generators m*F_j), written for msolve/Singular,
and logging helpers.  Everything mod p unless stated.
"""
import sys, os, json, time, math, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, reduce_mod, yparity
import flint
from flint import nmod_mat

IID = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
W2LOG = os.path.join(RES, 'wave2-search.log.jsonl')
PGB = [536870909, 536870879, 536870869]   # primes < 2^29 (Singular-safe); checked with sympy.isprime

def w2log(rec):
    rec = dict(rec); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); rec.setdefault('lane', IID)
    with open(W2LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec

def load():
    return core.load_instance(IID)

def set_threads(n):
    try:
        flint.ctx.threads = n
    except Exception:
        pass

def mon_add(a, b):
    return tuple(x + y for x, y in zip(a, b))

def mon_sub(a, b):
    return tuple(x - y for x, y in zip(a, b))

def poly_str(P):
    return core.to_str(P)

# --------------------------------------------------------------------------------------
# row-reduced Jacobian ideal in one degree / parity
# --------------------------------------------------------------------------------------
class JRef:
    """rref of J(F)_e^{par} over GF(p); rows = generators, columns = monomials of degree e and parity par."""
    def __init__(self, F, e, par, p, verbose=True):
        t0 = time.time()
        self.p, self.e, self.par = p, e, par
        mons = monomials(e, par); idx = {m: i for i, m in enumerate(mons)}
        self.mons, self.idx = mons, idx
        gens = jacobian.jacobian_generators(F, e, p, par)
        R, C = len(gens), len(mons)
        flat = [0] * (R * C)
        for r, (_, _, P) in enumerate(gens):
            base = r * C
            for ee, v in P.items():
                flat[base + idx[ee]] = v
        M = nmod_mat(R, C, flat, p)
        del flat
        t1 = time.time()
        Mr, rank = M.rref()
        del M
        t2 = time.time()
        piv = []; j = 0
        for i in range(rank):
            while int(Mr[i, j]) == 0:
                j += 1
            piv.append(j); j += 1
        self.Mr, self.rank, self.piv = Mr, rank, piv
        self.ngens = R
        self.pivset = set(piv)
        self.std = [j for j in range(C) if j not in self.pivset]
        self.rowof = {j: i for i, j in enumerate(piv)}
        self.ncols = C
        self.t_build, self.t_rref = t1 - t0, t2 - t1
        if verbose:
            print(f'[JRef] e={e} par={par} p={p}: gens={R} mons={C} rank={rank} dimR={C-rank} '
                  f'(build {t1-t0:.0f}s, rref {t2-t1:.0f}s)', file=sys.stderr, flush=True)

    def row(self, i):
        """Row i of the rref as a python list (length ncols)."""
        Mr = self.Mr; C = self.ncols
        return [int(Mr[i, j]) for j in range(C)]

    def reduce_rows(self, U):
        """U: nmod_mat (n x ncols).  Returns U - U[:, piv] * Rr  (normal forms; zero on pivot columns)."""
        n = U.nrows(); Rrows = self.Mr.nrows(); p = self.p
        flat = [0] * (n * Rrows)
        piv = self.piv
        for i in range(n):
            base = i * Rrows
            for k, j in enumerate(piv):
                v = int(U[i, j])
                if v:
                    flat[base + k] = v
        Upiv = nmod_mat(n, Rrows, flat, p)
        del flat
        return U - Upiv * self.Mr

    def nf_poly(self, P):
        """Normal form of one poly (dict mod p, degree e, parity par) as dict on standard monomials."""
        C = self.ncols
        flat = [0] * C
        for ee, v in P.items():
            flat[self.idx[ee]] = v % self.p
        U = nmod_mat(1, C, flat, self.p)
        Ur = self.reduce_rows(U)
        return {self.mons[j]: int(Ur[0, j]) for j in range(C) if int(Ur[0, j])}

    def nf_polys(self, polys):
        C = self.ncols; n = len(polys)
        flat = [0] * (n * C)
        for i, P in enumerate(polys):
            base = i * C
            for ee, v in P.items():
                flat[base + self.idx[ee]] = v % self.p
        U = nmod_mat(n, C, flat, self.p)
        del flat
        Ur = self.reduce_rows(U)
        out = []
        for i in range(n):
            out.append({self.mons[j]: int(Ur[i, j]) for j in self.std if int(Ur[i, j])})
        return out

# --------------------------------------------------------------------------------------
# E1 kernel in degree d
# --------------------------------------------------------------------------------------
def e1_kernel(F, A, d, par, p, jref=None, jref_d=None, verbose=True):
    """K_d^{par} = {D in S_d^{par}: A*D in J_{d+6}}, A a monomial.  Returns dict with
    'vectors' (list of dict polys mod p spanning the kernel incl. J_d), 'basis_modJ' (normal forms mod J_d,
    rref-canonical), dims, and the JRef objects used."""
    t0 = time.time()
    assert len(A) == 1
    Aexp, Acoef = next(iter(A.items()))
    pA = core.parity_of(A)
    if jref is None:
        jref = JRef(F, d + 6, (par + pA) % 2, p, verbose=verbose)
    monsD = monomials(d, par)
    Am = [mon_add(m, Aexp) for m in monsD]
    Pk = [k for k, mu in enumerate(Am) if jref.idx[mu] in jref.pivset]
    Tk = [k for k, mu in enumerate(Am) if jref.idx[mu] not in jref.pivset]
    StdA = {jref.idx[Am[k]]: k for k in Tk}
    std_other = [j for j in jref.std if j not in StdA]
    pos_other = {j: t for t, j in enumerate(std_other)}
    # Rsub: |P| x |std| ; equations on std_other, definitions on StdA
    nP = len(Pk)
    rows = [jref.row(jref.rowof[jref.idx[Am[k]]]) for k in Pk]
    t1 = time.time()
    # matrix E (|std_other| x nP): E[t, a] = rows[a][std_other[t]]  -> nullspace gives x_P
    flat = [0] * (len(std_other) * nP)
    W = nP
    for a, r in enumerate(rows):
        for t, j in enumerate(std_other):
            v = r[j]
            if v:
                flat[t * W + a] = v
    E = nmod_mat(len(std_other), nP, flat, p) if nP else None
    del flat
    vectors = []
    if nP:
        X, nul = E.nullspace()
        for c in range(nul):
            xP = [int(X[a, c]) for a in range(nP)]
            D = {}
            for a, k in enumerate(Pk):
                if xP[a]:
                    D[monsD[k]] = xP[a]
            # x_T
            for j, k in StdA.items():
                v = sum(xP[a] * rows[a][j] for a in range(nP)) % p
                if v:
                    D[monsD[k]] = v
            vectors.append(D)
    t2 = time.time()
    # rref of the vectors (canonical) then quotient by J_d
    if jref_d is None:
        jref_d = JRef(F, d, par, p, verbose=verbose)
    basis_modJ = []
    if vectors:
        nfs = jref_d.nf_polys(vectors)
        nfs = [v for v in nfs if v]
        if nfs:
            sm = jref_d.std; sidx = {jref_d.mons[j]: t for t, j in enumerate(sm)}
            flat = [0] * (len(nfs) * len(sm))
            for i, v in enumerate(nfs):
                for ee, c in v.items():
                    flat[i * len(sm) + sidx[ee]] = c
            Q = nmod_mat(len(nfs), len(sm), flat, p)
            Qr, rk = Q.rref()
            for i in range(rk):
                basis_modJ.append({jref_d.mons[sm[j]]: int(Qr[i, j]) for j in range(len(sm)) if int(Qr[i, j])})
    out = {'d': d, 'par': par, 'p': p, 'dim_S_d': len(monsD), 'dim_J_d': jref_d.rank, 'dim_R_d': len(jref_d.std),
           'dim_J_d6': jref.rank, 'dim_R_d6': len(jref.std), 'kernel_dim_incl_J': len(vectors),
           'kernel_mod_J': len(basis_modJ), 'nP': nP, 'nT': len(Tk), 'vectors': vectors, 'basis_modJ': basis_modJ,
           'jref': jref, 'jref_d': jref_d, 't_rows': round(t1 - t0, 1), 't_null': round(t2 - t1, 1),
           't_total': round(time.time() - t0, 1)}
    if verbose:
        print(f'[e1_kernel] d={d} par={par} p={p}: kernel incl J = {len(vectors)}, J_d = {jref_d.rank}, '
              f'mod J = {len(basis_modJ)} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    return out

# --------------------------------------------------------------------------------------
# rational reconstruction across primes
# --------------------------------------------------------------------------------------
def reconstruct_polys(lists_by_prime):
    """lists_by_prime: dict p -> list of polys (same length, canonical).  Returns (list of rational polys or None, info)."""
    primes = list(lists_by_prime)
    n = len(lists_by_prime[primes[0]])
    out = []
    for i in range(n):
        keys = set()
        for p in primes:
            if i >= len(lists_by_prime[p]):
                return None, {'reason': 'length mismatch'}
            keys |= set(lists_by_prime[p][i])
        D = {}
        for m in keys:
            x, mod = lists_by_prime[primes[0]][i].get(m, 0), primes[0]
            for p in primes[1:]:
                x, mod = core.crt_pair(x, mod, lists_by_prime[p][i].get(m, 0), p)
            v = core.ratrecon(x, mod)
            if v is None:
                return None, {'reason': f'ratrecon failed at vector {i} monomial {m}', 'modulus_bits': mod.bit_length()}
            if v:
                D[m] = v
        out.append(D)
    return out, {'modulus_bits': int(math.prod(primes)).bit_length()}

def poly_to_json(P):
    return core.to_str(P)

def poly_from_json(s):
    return core.parse(s)

def modp_poly_to_json(P):
    return {' '.join(map(str, e)): v for e, v in P.items()}

def modp_poly_from_json(o):
    return {tuple(int(x) for x in k.split()): int(v) for k, v in o.items()}

# --------------------------------------------------------------------------------------
# particular numerators: P(F) = A*G  (mod p, dense rref of the augmented system; several G at once)
# --------------------------------------------------------------------------------------
def particular_numerators(F, A, Gs, p, verbose=True):
    """For each G in Gs (dict polys mod p, degree d, pure parity) find P = (P_0..P_5) mod p with sum P_i F_i = A*G.
    Returns list of 6-lists (None where inconsistent)."""
    t0 = time.time()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    d = core.degree(Gs[0]); parG = core.parity_of(Gs[0])
    pars = fixed_d.numerator_parities(Gs[0], Ap)
    cols, meta = [], []
    for i in range(NV):
        for m in monomials(d + 1, pars[i]):
            cols.append(shift(Fi[i], m, 1, p)); meta.append((i, m))
    rhss = [mul(Ap, reduce_mod(G, p), p) for G in Gs]
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for r in rhss:
        for e in r: rows.setdefault(e, len(rows))
    R, C, K = len(rows), len(cols), len(rhss)
    W = C + K
    flat = [0] * (R * W)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * W + j] = v
    for k, r in enumerate(rhss):
        for e, v in r.items(): flat[rows[e] * W + C + k] = v
    M = nmod_mat(R, W, flat, p); del flat
    Mr, rank = M.rref(); del M
    piv = []; j = 0
    for i in range(rank):
        while int(Mr[i, j]) == 0: j += 1
        piv.append(j); j += 1
    out = []
    for k in range(K):
        if (C + k) in set(piv):
            out.append(None); continue
        P = [{} for _ in range(NV)]
        for i, j in enumerate(piv):
            if j >= C: continue
            v = int(Mr[i, C + k])
            if v:
                comp, m = meta[j]; P[comp][m] = v
        # verify
        chk = {}
        for i in range(NV): chk = add(chk, mul(P[i], Fi[i], p), p)
        assert sub(chk, rhss[k], p) == {}, 'numerator verification failed'
        out.append(P)
    if verbose:
        print(f'[numerators] d={d} p={p}: {K} systems, rows={R} cols={C} rank={rank}, '
              f'{sum(1 for o in out if o is not None)} consistent ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    return out

def numerators_rational(F, A, Gs_rat, primes=None, verbose=True):
    """Rational particular numerators for rational E1-admissible G's (CRT over primes + ratrecon + exact check)."""
    primes = primes or core.PRIMES[:4]
    acc = {}
    for k, p in enumerate(primes):
        Ps = particular_numerators(F, A, [reduce_mod(G, p) for G in Gs_rat], p, verbose=verbose)
        acc[p] = Ps
        if k == 0:
            continue
        # try reconstruction of each numerator vector
        ok_all = True; out = []
        for gi in range(len(Gs_rat)):
            if any(acc[q][gi] is None for q in acc):
                out.append(None); continue
            comps = []
            for i in range(NV):
                lists = {q: [acc[q][gi][i]] for q in acc}
                rat, info = reconstruct_polys(lists)
                if rat is None:
                    ok_all = False; break
                comps.append(rat[0])
            if not ok_all:
                break
            # exact check
            lhs = {}
            for i in range(NV): lhs = add(lhs, mul(comps[i], pderiv(F, i)))
            if sub(lhs, mul(A, Gs_rat[gi])) != {}:
                ok_all = False; break
            out.append(comps)
        if ok_all:
            return out, {'primes_used': k + 1}
    return None, {'primes_used': len(primes), 'reason': 'reconstruction failed'}

# --------------------------------------------------------------------------------------
# family of E1-admissible generators and the E2-mod-J quadratic relaxation
# --------------------------------------------------------------------------------------
def jacobian_family(F, d, parD):
    """Generators m*F_j of J_d with the parity parD, as (tag, G, P) with P = A m e_j filled in later."""
    out = []
    for j in range(NV):
        pj = 0 if j < 3 else 1
        for m in monomials(d - 5, (parD - pj) % 2):
            out.append(('J', j, m))
    return out

def spurious_subspace(d, parD):
    """Vector fields w = (h/6) E with h in S_{d-6}^{parD}: coordinates in the Jacobian family basis (j, m) ->
    list of dicts {(j, m): coeff}."""
    out = []
    if d < 6:
        return out
    for h in monomials(d - 6, parD):
        w = {}
        for j in range(NV):
            ej = [0] * NV; ej[j] = 1
            w[(j, mon_add(h, tuple(ej)))] = Fraction(1, 6)
        out.append(w)
    return out

class Family:
    """E1-admissible family: kernel reps (rational G_i with rational numerators N_i) + Jacobian generators.
    Members: list of dicts {tag, G (rational poly), P (6 rational polys), par, label}."""
    def __init__(self, F, A, d, kernel_reps=(), kernel_nums=(), parities=(0,), include_jacobian=True):
        self.F, self.A, self.d = F, A, d
        self.members = []
        for G, P in zip(kernel_reps, kernel_nums):
            self.members.append({'tag': 'K', 'G': G, 'P': P, 'par': core.parity_of(G), 'label': f'K{len(self.members)}'})
        self.nK = len(self.members)
        if include_jacobian:
            Fi = [pderiv(F, i) for i in range(NV)]
            for parD in parities:
                for (_, j, m) in jacobian_family(F, d, parD):
                    G = shift(Fi[j], m, 1)
                    P = [{} for _ in range(NV)]
                    P[j] = shift(A, m, 1)
                    self.members.append({'tag': 'J', 'j': j, 'm': m, 'G': G, 'P': P, 'par': parD,
                                         'label': f'J{j}:{"".join(map(str, m))}'})
        self.K = len(self.members)
        self.parities = sorted({mem['par'] for mem in self.members})

    def spurious_vectors(self):
        """Coordinate vectors (length K) spanning the F-multiple directions."""
        pos = {(mem['j'], mem['m']): k for k, mem in enumerate(self.members) if mem['tag'] == 'J'}
        vecs = []
        for parD in self.parities:
            for w in spurious_subspace(self.d, parD):
                v = [Fraction(0)] * self.K
                ok = True
                for key, c in w.items():
                    if key not in pos:
                        ok = False; break
                    v[pos[key]] = c
                if ok:
                    vecs.append(v)
        return vecs

def T_kl(Gk, Pl, p=None):
    """G_k div P_l - P_l(G_k)."""
    divP = {}
    for i in range(NV):
        divP = add(divP, pderiv(Pl[i], i, p), p)
    out = mul(Gk, divP, p)
    for i in range(NV):
        if Pl[i]:
            out = sub(out, mul(Pl[i], pderiv(Gk, i, p), p), p)
    return out

def relaxation_terms(fam, p, verbose=True):
    """All symmetrized T_kl (k <= l) mod p, using the structure of Jacobian generators for speed.
    Returns dict (k, l) -> poly mod p (coefficient of a_k a_l in B(a), k <= l)."""
    t0 = time.time()
    F, A = fam.F, fam.A
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    Ai = [pderiv(Ap, i, p) for i in range(NV)]
    # precomputed products for Jacobian pairs: FA[j][j'] = F_j * A_{j'}, AF[j] = A F_j, AH[j][j'] = A F_{jj'}
    FA = [[mul(Fi[j], Ai[jp], p) for jp in range(NV)] for j in range(NV)]
    AF = [mul(Ap, Fi[j], p) for j in range(NV)]
    AH = [[mul(Ap, Fij[j][jp], p) for jp in range(NV)] for j in range(NV)]
    mem = fam.members
    Gp = [reduce_mod(m['G'], p) for m in mem]
    Pp = [[reduce_mod(c, p) for c in m['P']] for m in mem]
    terms = {}
    def jac_T(k, l):
        # G_k = m F_j, P_l = A m' e_{j'}:  T_kl = m m' F_j A_{j'} + m (d_{j'} m') A F_j - m' (d_{j'} m) A F_j - m m' A F_{jj'}
        j, m = mem[k]['j'], mem[k]['m']; jp, mp = mem[l]['j'], mem[l]['m']
        mm = mon_add(m, mp)
        out = shift(FA[j][jp], mm, 1, p)
        if mp[jp]:
            e = list(mp); e[jp] -= 1
            out = add(out, shift(AF[j], mon_add(m, tuple(e)), mp[jp], p), p)
        if m[jp]:
            e = list(m); e[jp] -= 1
            out = sub(out, shift(AF[j], mon_add(mp, tuple(e)), m[jp], p), p)
        out = sub(out, shift(AH[j][jp], mm, 1, p), p)
        return out
    K = fam.K
    for k in range(K):
        for l in range(k, K):
            tk, tl = mem[k]['tag'], mem[l]['tag']
            if tk == 'J' and tl == 'J':
                t = jac_T(k, l)
                if l != k:
                    t = add(t, jac_T(l, k), p)
            else:
                t = T_kl(Gp[k], Pp[l], p)
                if l != k:
                    t = add(t, T_kl(Gp[l], Pp[k], p), p)
            if t:
                terms[(k, l)] = t
    if verbose:
        print(f'[relaxation_terms] K={K}: {len(terms)} nonzero symmetrized terms ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    return terms

def reduce_terms(fam, terms, p, jrefs=None, verbose=True):
    """Reduce every term modulo J_{2d} (by parity).  Returns (equations, jrefs) where equations is a list of
    dicts {(k,l): coeff} (one per standard monomial with a nonzero column), and jrefs the reducers used."""
    t0 = time.time()
    F = fam.F; d = fam.d; e = 2 * d
    keys = sorted(terms)
    bypar = {}
    for key in keys:
        par = core.parity_of(terms[key])
        if par is None:
            # mixed term: split
            for pr in (0, 1):
                part = {m: v for m, v in terms[key].items() if yparity(m) == pr}
                if part:
                    bypar.setdefault(pr, {})[key] = part
        else:
            bypar.setdefault(par, {})[key] = terms[key]
    jrefs = jrefs or {}
    equations = []
    for par, tdict in bypar.items():
        if par not in jrefs:
            jrefs[par] = JRef(F, e, par, p, verbose=verbose)
        jr = jrefs[par]
        ks = sorted(tdict)
        nfs = jr.nf_polys([tdict[k] for k in ks])
        # transpose: one equation per standard monomial
        cols = {}
        for key, nf in zip(ks, nfs):
            for m, v in nf.items():
                cols.setdefault(m, {})[key] = v
        for m, eq in cols.items():
            equations.append(eq)
        if verbose:
            print(f'[reduce_terms] par={par}: {len(ks)} terms -> {len(cols)} equations (dim R_{e}^{par} = {len(jr.std)})',
                  file=sys.stderr, flush=True)
    if verbose:
        print(f'[reduce_terms] total {len(equations)} equations ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    return equations, jrefs

def change_basis_quadratic(equations, K, basis_cols, p):
    """Rewrite quadratic forms sum_{k<=l} c_kl a_k a_l with a = B u (B: K x n, columns basis_cols as lists mod p)
    as quadratic forms in u.  Returns list of dicts {(i<=j): coeff}."""
    n = len(basis_cols)
    B = [[basis_cols[c][k] % p for c in range(n)] for k in range(K)]  # K x n
    out = []
    for eq in equations:
        # symmetric matrix C (K x K) from eq
        # compute Q = B^T C B  via flint
        flatC = [0] * (K * K)
        for (k, l), v in eq.items():
            if k == l:
                flatC[k * K + k] = (flatC[k * K + k] + v) % p
            else:
                h = (v * pow(2, -1, p)) % p
                flatC[k * K + l] = (flatC[k * K + l] + h) % p
                flatC[l * K + k] = (flatC[l * K + k] + h) % p
        C = nmod_mat(K, K, flatC, p)
        Bm = nmod_mat(K, n, [B[k][c] for k in range(K) for c in range(n)], p)
        Q = Bm.transpose() * C * Bm
        q = {}
        for i in range(n):
            v = int(Q[i, i])
            if v: q[(i, i)] = v
            for j in range(i + 1, n):
                v = (int(Q[i, j]) + int(Q[j, i])) % p
                if v: q[(i, j)] = v
        if q:
            out.append(q)
    return out

def write_msolve(path, equations, nvars, p, varname='a'):
    with open(path, 'w') as f:
        f.write(','.join(f'{varname}{i+1}' for i in range(nvars)) + '\n')
        f.write(f'{p}\n')
        lines = []
        for eq in equations:
            ts = []
            for key, v in eq.items():
                if key == ():
                    ts.append(f'{v}')
                elif len(key) == 1:
                    ts.append(f'{v}*{varname}{key[0]+1}')
                else:
                    i, j = key
                    ts.append(f'{v}*{varname}{i+1}*{varname}{j+1}' if i != j else f'{v}*{varname}{i+1}^2')
            lines.append('+'.join(ts))
        f.write(',\n'.join(lines) + '\n')

def write_singular(path, equations, nvars, p, varname='a', extra=''):
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, ({varname}(1..{nvars})), dp;\nideal I;\n')
        for eq in equations:
            ts = []
            for key, v in eq.items():
                if key == ():
                    ts.append(f'{v}')
                elif len(key) == 1:
                    ts.append(f'{v}*{varname}({key[0]+1})')
                else:
                    i, j = key
                    ts.append(f'{v}*{varname}({i+1})*{varname}({j+1})')
            f.write('I[size(I)+1] = ' + '+'.join(ts) + ';\n')
        f.write('int t = timer;\noption(redSB);\nideal G = std(I);\nprint("dim = " + string(dim(G)));\n'
                'print("vdim = " + string(vdim(G)));\nprint("size GB = " + string(size(G)));\n'
                'print("time = " + string(timer - t));\n')
        f.write(extra)
        f.write('quit;\n')

# --------------------------------------------------------------------------------------
# full exact fixed-D test helpers
# --------------------------------------------------------------------------------------
def full_test(F, A, q, D, family, note='', want_exact=True, image=False, verbose=True):
    """E1 + compressed consistency (+ exact solve) for one rational pure-parity D.  Logs to wave2 log."""
    t0 = time.time()
    d = core.degree(D)
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    rec = {'kind': 'fixed-D', 'family': family, 'd': d, 'parity': core.parity_of(D), 'e1': bool(ok),
           'e1_info': {k: info[k] for k in ('rank', 'rank_aug', 'rows', 'cols')}, 'note': note,
           'D': core.to_str(D)[:400]}
    if not ok:
        rec.update({'consistent': False, 'seconds': round(time.time() - t0, 1)})
        w2log(rec); return rec
    out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=want_exact)
    rec.update({'consistent': bool(out['consistent']), 'rank': out['rank'], 'nullity': out['nullity'],
                'rows': out['rows'], 'cols': out['cols']})
    if out.get('exact'):
        okv, vinfo = core.verify_witness(F, A, q, D, out['N'], out['Q'])
        rec['exact_verified'] = bool(okv)
        if okv:
            with open(os.path.join(RES, 'witness-candidate.json'), 'w') as f:
                json.dump({'instance': IID, 'd': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']],
                           'Q': core.to_str(out['Q']), 'family': family}, f, indent=1)
            rec['WITNESS'] = True
            print('!!! WITNESS FOUND — saved witness-candidate.json', file=sys.stderr, flush=True)
    if image:
        dim, basis, S = fixed_d.target_image(F, q, D)
        rec['served_dim'] = dim
    rec['seconds'] = round(time.time() - t0, 1)
    w2log(rec)
    return rec
