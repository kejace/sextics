"""Wave-2 EXCLUSION/RIGOR-lens library for the lane `short-power-presentation` (file prefix excl_ — the
lane directory is shared with other lens agents).

All modular linear algebra uses python-flint nmod_mat.  Every rigorous claim produced with this
library is backed by one of two arguments, which the functions record explicitly:

  (R1) Hilbert-series check.  F is smooth over Q, so F_0..F_5 is a regular sequence and the
       tau-equivariant Hilbert series of R = S/J(F) is known exactly:
          H_R(t) = (1+t+t^2+t^3+t^4)^6,   Theta_R(t) = (1+t^2+t^4+t^6+t^8)^3  (trace of tau = Y-sign).
       rank_p of the generator matrix of J_e^{par} is always <= dim_Q J_e^{par}; when the two agree the
       prime p is "good in degree e" and the mod-p quotient S_e/J_e is the reduction of the Q-quotient.
  (R2) Modular rank inequality.  rank_p M <= rank_Q M for every integer matrix M.  An inconsistency
       mod p of  M x = b  proves inconsistency over Q as soon as rank_p M is known to equal rank_Q M,
       e.g. because rank_p M = #cols - (dimension of an explicitly known Q-nullspace).

Conventions: polys are dicts {exponent 6-tuple: coeff}; parity = Y-degree mod 2; A = X0^3 Y0^3 (odd).
"""
import sys, time, json, os, itertools, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, reduce_mod
from fractions import Fraction
from flint import nmod_mat

IID = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
LOG = os.path.join(RES, 'wave2-excl.log.jsonl')

# primes below 2^31 (Singular / msolve friendly), none of them 2,3,5,7
PRIMES31 = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549]
PRIMES30 = [1073741789, 1073741783, 1073741741, 1073741723, 1073741719, 1073741717]
PRIMES29 = [536870909, 536870879, 536870869, 536870849, 536870839, 536870837]


def wlog(**rec):
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('lens', 'excl')
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print('LOG ' + json.dumps(rec), flush=True)


def load():
    return core.load_instance(IID)


# --------------------------------------------------------------------------------------------
# (R1) Hilbert-series bookkeeping
# --------------------------------------------------------------------------------------------
def _pmul(a, b, nmax):
    out = [0] * (nmax + 1)
    for i, x in enumerate(a):
        if x == 0 or i > nmax:
            continue
        for j, y in enumerate(b):
            if i + j > nmax:
                break
            out[i + j] += x * y
    return out


def _ppow(a, k, nmax):
    out = [1] + [0] * nmax
    for _ in range(k):
        out = _pmul(out, a, nmax)
    return out


_NMAX = 60
_H = _ppow([1, 1, 1, 1, 1], 6, _NMAX)                  # dim R_e
_TH = _ppow([1, 0, 1, 0, 1, 0, 1, 0, 1], 3, _NMAX)     # trace of tau on R_e


def dimS(e, par=None):
    """dim S_e^{par} over Q (par None = all)."""
    if e < 0:
        return 0
    tot = math.comb(e + 5, 5)
    if par is None:
        return tot
    tr = math.comb(e // 2 + 2, 2) if e % 2 == 0 else 0   # trace of tau on S_e = coeff of 1/(1-t^2)^3
    return (tot + tr) // 2 if par == 0 else (tot - tr) // 2


def dimR(e, par=None):
    if e < 0 or e > _NMAX:
        return 0
    if par is None:
        return _H[e]
    return (_H[e] + _TH[e]) // 2 if par == 0 else (_H[e] - _TH[e]) // 2


def dimJ(e, par=None):
    return dimS(e, par) - dimR(e, par)


# --------------------------------------------------------------------------------------------
# basic modular linear algebra helpers
# --------------------------------------------------------------------------------------------
def jac_gens(Fp, deg, par, p):
    """[(i, m, m*F_i)] spanning J(F)_deg^{par} mod p."""
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    out = []
    if deg < 5:
        return out
    for i in range(NV):
        pi = 0 if i < 3 else 1
        mpar = None if par is None else (par - pi) % 2
        for m in monomials(deg - 5, mpar):
            out.append((i, m, shift(Fi[i], m, 1, p)))
    return out


def cols_matrix(polys, idx, p):
    """nmod_mat with the polys as COLUMNS over the row index idx (monomial -> row)."""
    R, C = len(idx), len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items():
            flat[idx[e] * C + j] = v
    return nmod_mat(R, C, flat, p)


def rows_matrix(polys, idx, p):
    """nmod_mat with the polys as ROWS over the column index idx (monomial -> column)."""
    R, C = len(polys), len(idx)
    flat = [0] * (R * C)
    for i, P in enumerate(polys):
        base = i * C
        for e, v in P.items():
            flat[base + idx[e]] = v
    return nmod_mat(R, C, flat, p)


class JReducer:
    """Row-reduced J(F)_e^{par} mod p in the monomial basis; normal forms and the (R1) check."""

    def __init__(self, Fp, e, par, p):
        self.p, self.e, self.par = p, e, par
        self.mons = monomials(e, par)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        gens = jac_gens(Fp, e, par, p)
        n = len(self.mons)
        if not gens:
            self.rank, self.pivots, self.Rr = 0, [], None
        else:
            M = rows_matrix([g[2] for g in gens], self.idx, p)
            Rr0, r = M.rref()
            self.rank = r
            self.Rr = nmod_mat(r, n, [int(Rr0[i, j]) for i in range(r) for j in range(n)], p) if r else None
            piv = []
            for i in range(r):
                j = 0
                while int(Rr0[i, j]) == 0:
                    j += 1
                piv.append(j)
            self.pivots = piv
        self.pivset = set(self.pivots)
        self.nonpiv = [j for j in range(n) if j not in self.pivset]
        self.expected = dimJ(e, par)
        self.good = (self.rank == self.expected)

    def reduce_matrix(self, U):
        """U: nmod_mat (k x n) of polys as rows -> normal forms (k x n) = U - U[:,piv] * Rr."""
        if self.rank == 0:
            return U
        k = U.nrows()
        Up = nmod_mat(k, self.rank, [int(U[i, pj]) for i in range(k) for pj in self.pivots], self.p)
        return U - Up * self.Rr

    def normal_form_vec(self, P):
        """normal form of a single poly (dict) -> dict over non-pivot monomials."""
        U = rows_matrix([P], self.idx, self.p)
        Ur = self.reduce_matrix(U)
        return {self.mons[j]: int(Ur[0, j]) for j in self.nonpiv if int(Ur[0, j])}


# --------------------------------------------------------------------------------------------
# E1 kernel  K_d^{par} = {D in S_d^{par} : A D in J}  and its dimension modulo J_d^{par}
# --------------------------------------------------------------------------------------------
def e1_kernel(F, A, d, parD, p, want_basis=False):
    """Returns (info, basis) where basis (if requested) is a list of normal-form representatives mod p
    (dicts) of a basis of K_d/J_d, in rref form w.r.t. the monomial order of core.monomials.

    Rigor (R1+R2): dim_Q(K_d/J_d) <= value reported here provided info['hilbert_ok'] is True."""
    t0 = time.time()
    Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
    degT, parT = d + 6, (parD + 1) % 2
    monsT = monomials(degT, parT)
    idxT = {m: i for i, m in enumerate(monsT)}
    gens = jac_gens(Fp, degT, parT, p)
    monsD = monomials(d, parD)
    AD = [shift(Ap, m, 1, p) for m in monsD]
    MJ = cols_matrix([g[2] for g in gens], idxT, p)
    rJ = MJ.rank()
    Mall = cols_matrix([g[2] for g in gens] + AD, idxT, p)
    rall = Mall.rank()
    kdim = len(monsD) - (rall - rJ)
    red = JReducer(Fp, d, parD, p)
    info = {'d': d, 'parD': parD, 'p': p, 'dim S_d': len(monsD), 'dim J_d': red.rank, 'dim J_d expected': red.expected,
            'rank J_{d+6}': rJ, 'rank J_{d+6} expected': dimJ(degT, parT), 'rank [J|AD]': rall,
            'kernel dim (incl J_d)': kdim, 'kernel mod J': kdim - red.rank,
            'hilbert_ok': bool(red.good and rJ == dimJ(degT, parT)), 'rows': len(monsT), 'cols': len(gens) + len(monsD)}
    basis = None
    if want_basis and kdim > red.rank:
        X, nul = Mall.nullspace()
        nJ = len(gens)
        nD = len(monsD)
        Dv = nmod_mat(nul, nD, [int(X[nJ + j, k]) for k in range(nul) for j in range(nD)], p)   # rows = D parts
        Dr0, rk = Dv.rref()
        Dr = nmod_mat(rk, nD, [int(Dr0[i, j]) for i in range(rk) for j in range(nD)], p)
        # normal forms modulo J_d then rref -> canonical representatives
        NF = red.reduce_matrix(Dr)
        NFr, rk2 = NF.rref()
        basis = []
        for i in range(rk2):
            D = {monsD[j]: int(NFr[i, j]) for j in range(nD) if int(NFr[i, j])}
            if D:
                basis.append(D)
        info['basis_count'] = len(basis)
    info['seconds'] = round(time.time() - t0, 1)
    return info, basis


def lift_polys(polys_by_prime):
    """polys_by_prime: list of (p, [poly mod p ...]) with identical shapes -> list of exact polys or None."""
    k = len(polys_by_prime[0][1])
    out = []
    for i in range(k):
        keys = set()
        for p, polys in polys_by_prime:
            keys |= set(polys[i])
        D = {}
        for m in keys:
            x, mod = 0, 1
            for p, polys in polys_by_prime:
                x, mod = core.crt_pair(x, mod, polys[i].get(m, 0) % p, p)
            v = core.ratrecon(x, mod)
            if v is None:
                return None
            if v:
                D[m] = v
        out.append(D)
    return out


# --------------------------------------------------------------------------------------------
# exact Jacobian membership certificates:  P = sum_i g_i F_i  over Q
# --------------------------------------------------------------------------------------------
def membership_mod_p(Fp, P, p):
    """Particular solution (g_i) mod p of sum g_i F_i = P (free variables 0), plus pivot structure."""
    deg = core.degree(P)
    par = core.parity_of(P)
    gens = jac_gens(Fp, deg, par, p)
    mons = monomials(deg, par)
    idx = {m: i for i, m in enumerate(mons)}
    M = cols_matrix([g[2] for g in gens] + [P], idx, p)
    Mr, rank = M.rref()
    C = len(gens)
    piv = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        piv.append(j)
    if piv and piv[-1] == C:
        return None, {'rank': rank - 1, 'rank_aug': rank, 'consistent': False, 'cols': C, 'rows': len(mons)}
    g = [{} for _ in range(NV)]
    for i, j in enumerate(piv):
        v = int(Mr[i, C])
        if v:
            gi, m, _ = gens[j]
            g[gi][m] = v
    return g, {'rank': rank, 'consistent': True, 'pivots': piv, 'cols': C, 'rows': len(mons)}


def exact_membership(F, P, primes=None, max_primes=6, verbose=True):
    """Exact (g_i) with sum g_i F_i = P over Q, or None.  Verified exactly with Fractions."""
    primes = primes or PRIMES31
    acc = []
    piv0 = None
    Fi = [pderiv(F, i) for i in range(NV)]
    for p in primes[:max_primes]:
        Fp = reduce_mod(F, p)
        Pp = reduce_mod(P, p)
        g, info = membership_mod_p(Fp, Pp, p)
        if g is None:
            return None, {'inconsistent_mod': p, **info}
        if piv0 is None:
            piv0 = info['pivots']
        elif piv0 != info['pivots']:
            if verbose:
                print('[exact_membership] pivot structure differs between primes; restarting', flush=True)
            acc = []
            piv0 = info['pivots']
        acc.append((p, g))
        lifted = lift_polys(acc)
        if lifted is None:
            continue
        chk = {}
        for i in range(NV):
            chk = add(chk, mul(lifted[i], Fi[i]))
        if sub(chk, P) == {}:
            return lifted, {'primes': len(acc), 'rank': info['rank']}
        elif verbose:
            print(f'[exact_membership] lift after {len(acc)} primes did not verify; continuing', flush=True)
    return None, {'lift_failed': True}


# --------------------------------------------------------------------------------------------
# rigorous fixed-D test (compressed modular consistency + (R2) trivial-nullity accounting)
# --------------------------------------------------------------------------------------------
def trivial_nullity(d, q, pA, pD):
    """dim of the explicitly known Q-nullspace {N = F u + G E, Q determined} of the fixed-D system."""
    degN = d + 6 * q - 5
    parsN = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
    nu = sum(len(monomials(degN - 6, parsN[i])) for i in range(NV))
    parG = (pA - pD) % 2
    return nu + len(monomials(degN - 1, parG)) - len(monomials(degN - 7, parG))


def rigorous_fixed_D(F, A, q, D, p, k_extra=60, verbose=True):
    """Compressed consistency test of the full certificate system for a fixed pure-parity D.

    Returns dict with 'consistent' (mod p, compressed), 'rank' (of the compressed system, a LOWER bound
    for rank_p and rank_Q), 'trivial_nullity', and 'rigorous_exclusion' = (not consistent) and
    (cols - rank == trivial_nullity).  Compressed inconsistency alone already proves inconsistency mod p."""
    S = fixed_d.FixedDSystem(F, A, q, D, p, verbose=verbose)
    r = fixed_d.solve_compressed(S, k_extra=k_extra, verbose=verbose)
    pA, pD = core.parity_of(A), core.parity_of(D)
    tn = trivial_nullity(S.d, q, pA, pD)
    out = {'d': S.d, 'parD': pD, 'p': p, 'rows': S.nrows, 'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ,
           'consistent_mod_p': bool(r['consistent']), 'rank_compressed': r['rank'], 'nullity_compressed': r['nullity'],
           'trivial_nullity': tn, 'k_compressed': r['k']}
    out['rigorous_exclusion'] = (not r['consistent']) and (r['nullity'] == tn)
    out['x'] = r['x'] if r['consistent'] else None
    out['pivots'] = r['pivots']
    return out


# --------------------------------------------------------------------------------------------
# The quadratic (Hessian-type) relaxation of E2 for the family  D = sum_a lam_a D_a + w(F)
# --------------------------------------------------------------------------------------------
def field_basis(d, parD):
    """Unknown coordinates of a polynomial vector field w of degree d-5 with w(F) of parity parD."""
    out = []
    for j in range(NV):
        pj = 0 if j < 3 else 1
        for m in monomials(d - 5, (parD - pj) % 2):
            out.append((j, m))
    return out


def build_relaxation(F, A, d, kern, pars, p, slice_euler=True, verbose=True):
    """Quadratic system (mod p) in the unknowns (lam_a) + (w_{j,m}) expressing

        Xi(lam, w) := sum_{a,b} lam_a lam_b B_{D_a}(N_b)
                      + sum_a lam_a ( D_a w(A) - A w(D_a) - N_a^T Hess(F) w )
                      - A * w^T Hess(F) w                               in  J(F)_{2d}

    which is a necessary condition for a certificate with denominator D = sum lam_a D_a + w(F)
    (derivation in the wave-2 report).  kern = [(D_a, N_a)] with N_a(F) = A D_a exactly (mod p);
    pars = list of parities for w (e.g. [0], [1], [0,1]).  With slice_euler the coordinates w_{0,m},
    X0 | m, are dropped (they parametrise the trivial solutions w = G*E, which act on the solution set).

    Returns dict(unknowns, eqs (list of {(i,j): c} with i<=j), info)."""
    t0 = time.time()
    Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    H = [[pderiv(Fi[j], l, p) for l in range(NV)] for j in range(NV)]
    Ai = [pderiv(Ap, i, p) for i in range(NV)]
    kern = [(reduce_mod(Da, p), [reduce_mod(n, p) for n in Na]) for Da, Na in kern]
    for Da, Na in kern:   # sanity: N_a(F) = A D_a
        chk = {}
        for i in range(NV):
            chk = add(chk, mul(Na[i], Fi[i], p), p)
        assert sub(chk, mul(Ap, Da, p), p) == {}, 'kernel particular solution does not satisfy N(F) = A D'
    unknowns = [('lam', a) for a in range(len(kern))]
    for parD in pars:
        for (j, m) in field_basis(d, parD):
            if slice_euler and j == 0 and m[0] > 0:
                continue
            unknowns.append(('w', j, m))
    nu = len(unknowns)
    # coefficient polynomials of each unordered pair of unknowns
    terms = {}   # (i, j) i<=j -> poly

    def addterm(i, j, P):
        if not P:
            return
        if i > j:
            i, j = j, i
        terms[(i, j)] = add(terms.get((i, j), {}), P, p)

    # lam-lam terms: B_{D_a}(N_b)
    divN = []
    for Da, Na in kern:
        dv = {}
        for i in range(NV):
            dv = add(dv, pderiv(Na[i], i, p), p)
        divN.append(dv)
    Dgrad = [[pderiv(Da, i, p) for i in range(NV)] for Da, _ in kern]
    for a, (Da, Na) in enumerate(kern):
        for b, (Db, Nb) in enumerate(kern):
            B = mul(Da, divN[b], p)
            for i in range(NV):
                B = sub(B, mul(Nb[i], Dgrad[a][i], p), p)
            addterm(a, b, B)
    # lam-w terms: D_a * Z^m A_j - A Z^m D_{a,j} - (sum_i N_{a,i} F_{ij}) Z^m
    HN = []
    for Da, Na in kern:
        hn = []
        for j in range(NV):
            s = {}
            for i in range(NV):
                s = add(s, mul(Na[i], H[i][j], p), p)
            hn.append(s)
        HN.append(hn)
    ADaj = [[sub(mul(Da, Ai[j], p), mul(Ap, Dgrad[a][j], p), p) for j in range(NV)] for a, (Da, _) in enumerate(kern)]
    widx = [(k, u) for k, u in enumerate(unknowns) if u[0] == 'w']
    for a in range(len(kern)):
        for k, (_, j, m) in widx:
            P = shift(sub(ADaj[a][j], HN[a][j], p), m, 1, p)
            addterm(a, k, P)
    # w-w terms: -A w^T H w  -> coefficient of w_{j,m} w_{l,m'}:  -(1 or 2) A Z^{m+m'} F_{jl}
    AH = [[mul(Ap, H[j][l], p) for l in range(NV)] for j in range(NV)]
    for x in range(len(widx)):
        k1, (_, j, m) = widx[x]
        for y in range(x, len(widx)):
            k2, (_, l, m2) = widx[y]
            mm = tuple(u + v for u, v in zip(m, m2))
            c = (p - 1) if x == y else (p - 2)
            addterm(k1, k2, shift(AH[j][l], mm, c, p))
    t1 = time.time()
    # reduce every coefficient polynomial modulo J_{2d} (both parities if needed) and read off equations
    e = 2 * d
    keys = list(terms)
    parities_needed = sorted({core.yparity(mon) for P in terms.values() for mon in P})
    eqs = []
    parts = {}
    for par in parities_needed:
        red = JReducer(Fp, e, par, p)
        assert red.good, f'prime {p} not good in degree {e} parity {par}: rank {red.rank} != {red.expected}'
        n = len(red.mons)
        polys = [{mon: v for mon, v in terms[k].items() if core.yparity(mon) == par} for k in keys]
        U = rows_matrix(polys, red.idx, p)
        Ur = red.reduce_matrix(U)
        cnt = 0
        # equations: one per non-pivot monomial column
        for c in red.nonpiv:
            eq = {}
            for i, k in enumerate(keys):
                v = int(Ur[i, c])
                if v:
                    eq[k] = v
            if eq:
                eqs.append(eq)
                cnt += 1
        parts[par] = {'dim S_e': n, 'dim J_e': red.rank, 'dim R_e': n - red.rank, 'equations': cnt}
    info = {'d': d, 'pars': pars, 'p': p, 'n_kernel': len(kern), 'unknowns': nu, 'slice_euler': slice_euler,
            'pairs': len(keys), 'equations': len(eqs), 'parts': parts, 'build_seconds': round(t1 - t0, 1),
            'reduce_seconds': round(time.time() - t1, 1)}
    if verbose:
        print(json.dumps(info), flush=True)
    return {'unknowns': unknowns, 'eqs': eqs, 'info': info}


def write_singular(sysd, path, p, extra=''):
    """Singular script: ring in p with variables u(1..n); ideal I; std, dimension, degree, Hilbert polynomial."""
    nu = len(sysd['unknowns'])
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, (u(1..{nu})), dp;\nideal I;\n')
        for eq in sysd['eqs']:
            ts = []
            for (i, j), v in eq.items():
                ts.append(f'{v}*u({i+1})*u({j+1})')
            f.write('I[size(I)+1] = ' + '+'.join(ts) + ';\n')
        f.write('int t0 = timer;\noption(redSB);\nideal G = std(I);\n'
                'print("std_time = " + string(timer - t0));\nprint("dim = " + string(dim(G)));\n'
                'print("vdim = " + string(vdim(G)));\nprint("degree:"); degree(G);\n'
                'intvec hv = hilb(G, 1);\nprint("hilb1 = " + string(hv));\n')
        f.write(extra)
        f.write('quit;\n')


def write_msolve(sysd, path, p):
    nu = len(sysd['unknowns'])
    with open(path, 'w') as f:
        f.write(','.join(f'u{i+1}' for i in range(nu)) + '\n')
        f.write(f'{p}\n')
        lines = []
        for eq in sysd['eqs']:
            ts = []
            for (i, j), v in eq.items():
                ts.append(f'{v}*u{i+1}^2' if i == j else f'{v}*u{i+1}*u{j+1}')
            lines.append('+'.join(ts))
        f.write(',\n'.join(lines) + '\n')


def save_system(sysd, path):
    with open(path, 'w') as f:
        json.dump({'unknowns': [list(u) if u[0] == 'lam' else [u[0], u[1], list(u[2])] for u in sysd['unknowns']],
                   'eqs': [[[list(k), v] for k, v in eq.items()] for eq in sysd['eqs']],
                   'info': sysd['info']}, f)


# --------------------------------------------------------------------------------------------
# exact dependency certificates:  dim_Q (span{P_k} + J_e)/J_e  <=  #P - #verified dependencies
# --------------------------------------------------------------------------------------------
def _mpoly_ctx():
    from flint import fmpq_mpoly_ctx, Ordering
    return fmpq_mpoly_ctx.get(('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2'), Ordering.degrevlex)


def _to_mpoly(P, ctx):
    from flint import fmpq
    return ctx.from_dict({e: fmpq(Fraction(v).numerator, Fraction(v).denominator) for e, v in P.items()})


def exact_dependencies(F, polys, e, par, primes=None, max_primes=4, verbose=True):
    """polys: exact rational homogeneous polys of degree e and parity par.  Finds, mod p, a basis of the space of
    dependencies  sum_k c_k polys_k  in  J(F)_e  (c-part of the nullspace of [J gens | polys]), lifts each vector
    (c, h) with sum_k c_k polys_k = sum_i h_i F_i to Q by CRT/rational reconstruction over several primes and verifies
    the identity exactly with fmpq_mpoly.  Returns dict(n_polys, rank_modp (= dim of span mod J mod p), n_deps_verified,
    dim_span_mod_J_exact_upper = n_polys - n_deps_verified, ...).  If every dependency verifies, dim_Q = rank_modp."""
    primes = primes or PRIMES31
    t0 = time.time()
    K = len(polys)
    mons = monomials(e, par); idx = {m: i for i, m in enumerate(mons)}
    Fi = [pderiv(F, i) for i in range(NV)]
    # generator bookkeeping (same order every prime)
    gen_meta = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        for m in monomials(e - 5, (par - pi) % 2):
            gen_meta.append((i, m))
    nG = len(gen_meta)
    acc = []
    piv0 = None
    deps_vecs = None
    for p in primes[:max_primes]:
        Fp = reduce_mod(F, p)
        gens = [shift(pderiv(Fp, i, p), m, 1, p) for (i, m) in gen_meta]
        pp = [reduce_mod(P, p) for P in polys]
        M = cols_matrix(gens + pp, idx, p)
        X, nul = M.nullspace()
        # restrict to vectors with nonzero c-part: rref of the c-part, then re-express
        Xm = nmod_mat(nG + K, nul, [int(X[i, k]) for i in range(nG + K) for k in range(nul)], p)
        C = nmod_mat(nul, K, [int(X[nG + k, j]) for j in range(nul) for k in range(K)], p)   # rows = nullvectors, cols = c-part
        # choose rref of C to get independent c-parts, keeping track of the combination to rebuild the h-part
        # simpler: rref the whole transposed nullspace with c-columns FIRST
        T = nmod_mat(nul, K + nG, [int(X[nG + k, j]) for j in range(nul) for k in range(K)] and
                     [ (int(X[nG + c, j]) if c < K else int(X[c - K, j])) for j in range(nul) for c in range(K + nG)], p)
        Tr, r = T.rref()
        rows = []
        for i in range(r):
            j = 0
            while j < K + nG and int(Tr[i, j]) == 0:
                j += 1
            if j >= K:
                break   # remaining rows have zero c-part (pure Koszul syzygies)
            rows.append({c: int(Tr[i, c]) for c in range(K + nG) if int(Tr[i, c])})
        ndep = len(rows)
        piv = [min(rw) for rw in rows]
        if piv0 is None:
            piv0 = piv
        elif piv != piv0:
            if verbose:
                print('[exact_dependencies] pivot structure differs between primes; restarting', flush=True)
            acc = []; piv0 = piv
        acc.append((p, rows))
        if len(acc) < 2:
            continue
        # lift
        lifted = []
        ok = True
        for i in range(ndep):
            keys = set()
            for p_, rws in acc:
                keys |= set(rws[i])
            vec = {}
            for c in keys:
                x, mod = 0, 1
                for p_, rws in acc:
                    x, mod = core.crt_pair(x, mod, rws[i].get(c, 0) % p_, p_)
                v = core.ratrecon(x, mod)
                if v is None:
                    ok = False; break
                if v:
                    vec[c] = v
            if not ok:
                break
            lifted.append(vec)
        if not ok:
            if verbose:
                print(f'[exact_dependencies] lift failed with {len(acc)} primes; continuing', flush=True)
            continue
        # exact verification with fmpq_mpoly
        ctx = _mpoly_ctx()
        Fm = [_to_mpoly(f, ctx) for f in Fi]
        Pm = [_to_mpoly(P, ctx) for P in polys]
        Zm = {}
        nver = 0
        bad = 0
        for vec in lifted:
            lhs = ctx.from_dict({})
            hpoly = [{} for _ in range(NV)]
            for c, v in vec.items():
                if c < K:
                    lhs += _to_mpoly({(0,) * NV: v}, ctx) * Pm[c]
                else:
                    i, m = gen_meta[c - K]
                    hpoly[i][m] = hpoly[i].get(m, 0) + v
            for i in range(NV):
                if hpoly[i]:
                    lhs += _to_mpoly(hpoly[i], ctx) * Fm[i]
            if lhs == ctx.from_dict({}):
                nver += 1
            else:
                bad += 1
        out = {'n_polys': K, 'n_gens': nG, 'rank_modp_span_mod_J': K - ndep, 'n_deps_modp': ndep, 'n_deps_verified': nver,
               'n_deps_failed': bad, 'primes': [p_ for p_, _ in acc], 'dim_span_mod_J_exact_upper': K - nver,
               'exact_equals_modp': (nver == ndep), 'seconds': round(time.time() - t0, 1)}
        if verbose:
            print('[exact_dependencies]', json.dumps(out), flush=True)
        return out
    return {'n_polys': K, 'lift_failed': True, 'seconds': round(time.time() - t0, 1)}
