"""Exact (over Q) emptiness certificates for bilinear charts.

For a sector system built mod p (w2_bilinear.build/build_fast) the chart a (t_b = 0 for b < a, t_a = 1) is EMPTY over C if
there is a rational linear functional ell on S_{2d}[sector] with ell(F * S_{2d-6}[sector]) = 0 and
      ell( B_{D(t)}(N(t,s)) ) |_chart  ==  1          (constant-cofactor Nullstellensatz certificate),
where D_a, N_a (rational, reconstructed from the modular basis and verified: N_a(F) = A D_a exactly) and S_k (rational Koszul
syzygies, verified S_k(F) = 0) are the rational data.  Since every certificate (D, N) with D = sum t_a D_a in the chart and
N = sum t_a N_a + sum s_k S_k satisfies B_D(N) = F * (something), ell(B_D(N)) = 0, contradiction.
Cascade: if instead ell(B)|_chart = v^2 (or a linear form) for a chart variable v, then v = 0 on the solution set; substitute and
repeat.  Everything modular is only used to FIND ell; the identity is checked exactly with flint fmpq_mpoly.
Output per chart: {'cert': 'constant'|'cascade'|None, 'steps': [...], 'verified_exact': bool}.
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from w2lib import char_of, cadd, VARCHAR, monos_char
from hodge import core
from hodge.core import NV, VARS, monomials, mul, add, sub, pderiv, shift, reduce_mod
from fractions import Fraction
from flint import nmod_mat, fmpq_mpoly_ctx, fmpq, Ordering

def rat_poly(P, p):
    out = {}
    for e, v in P.items():
        r = core.ratrecon(v, p)
        if r is None: return None
        if r: out[e] = r
    return out

def rational_data(F, A, q, Ds, Ns, syz, p):
    """Reconstruct D_a, N_a, S_k over Q and verify N_a(F) = A D_a, S_k(F) = 0 exactly.  Returns (Dq, Nq, Sq) or raises."""
    Fi = [pderiv(F, i) for i in range(NV)]
    Dq = []; Nq = []; Sq = []
    for D, N in zip(Ds, Ns):
        d_ = rat_poly(D, p); n_ = [rat_poly(N[i], p) for i in range(NV)]
        if d_ is None or any(x is None for x in n_): raise ValueError('lift reconstruction failed')
        lhs = {}
        for i in range(NV): lhs = add(lhs, mul(n_[i], Fi[i]))
        if lhs != mul(A, d_): raise ValueError('lift verification failed')
        Dq.append(d_); Nq.append(n_)
    for S in syz:
        s_ = [rat_poly(S[i], p) for i in range(NV)]
        if any(x is None for x in s_): raise ValueError('syzygy reconstruction failed')
        lhs = {}
        for i in range(NV): lhs = add(lhs, mul(s_[i], Fi[i]))
        if lhs: raise ValueError('syzygy verification failed')
        Sq.append(s_)
    return Dq, Nq, Sq

class ExactB:
    """Exact rational B_{D_a}(N_b), B_{D_a}(S_k) as dict polys (via fmpq_mpoly), cached."""
    def __init__(self, Dq, Nq, Sq):
        self.ctx = fmpq_mpoly_ctx.get(VARS, ordering=Ordering.lex)
        ctx = self.ctx
        def mp(P): return ctx.from_dict({tuple(int(x) for x in e): fmpq(v.numerator, v.denominator) for e, v in P.items()}) if P else ctx.from_dict({})
        self.D = [mp(D) for D in Dq]; self.dD = [[D.derivative(i) for i in range(NV)] for D in self.D]
        def prep(N):
            Nm = [mp(N[i]) for i in range(NV)]
            return Nm, sum((Nm[i].derivative(i) for i in range(NV)), ctx.from_dict({}))
        self.N = [prep(N) for N in Nq]; self.S = [prep(S) for S in Sq]
        self.cache = {}
    def B(self, a, kind, b):
        key = (a, kind, b)
        if key not in self.cache:
            Nm, div = (self.N if kind == 't' else self.S)[b]
            P = self.D[a] * div
            for i in range(NV):
                if Nm[i]: P = P - Nm[i] * self.dD[a][i]
            self.cache[key] = {e: Fraction(int(v.p), int(v.q)) for e, v in P.to_dict().items()}
        return self.cache[key]

def chart_variables(nt, ns, a, zero=()):
    return [('t', i) for i in range(a + 1, nt) if ('t', i) not in zero] + [('s', k) for k in range(ns) if ('s', k) not in zero]

def chart_equations_modp(sysd, a, zero):
    """Equations on chart a with extra variables in `zero` set to 0, as dicts monomial(tuple of var ids)->coeff mod p,
    keyed by equation index."""
    p = sysd['p']; nt, ns = sysd['nt'], sysd['ns']
    out = []
    for r, eq in enumerate(sysd['eqs']):
        P = {}
        for (k1, i1, k2, i2), cf in eq.items():
            vs = []; dead = False
            for kk, ii in ((k1, i1), (k2, i2)):
                if kk == 't' and ii < a: dead = True; break
                if kk == 't' and ii == a: continue
                if (kk, ii) in zero: dead = True; break
                vs.append((kk, ii))
            if dead: continue
            mon = tuple(sorted(vs))
            P[mon] = (P.get(mon, 0) + cf) % p
        P = {m: v for m, v in P.items() if v}
        out.append((r, P))
    return out

def find_combination(eqs, target, p):
    """lambda mod p with sum_r lambda_r eq_r == target (as polynomials).  Returns dict r->lambda or None."""
    mons = {}
    for r, P in eqs:
        for m in P: mons.setdefault(m, len(mons))
    for m in target: mons.setdefault(m, len(mons))
    R, C = len(mons), len(eqs)
    if C == 0: return None
    flat = [0] * (R * (C + 1))
    for j, (r, P) in enumerate(eqs):
        for m, v in P.items(): flat[mons[m] * (C + 1) + j] = v
    for m, v in target.items(): flat[mons[m] * (C + 1) + C] = v % p
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    if piv and piv[-1] == C: return None
    lam = {}
    for i, j in enumerate(piv):
        v = int(Mr[i, C])
        if v: lam[eqs[j][0]] = v
    return lam

def functional_from_lambda(lam, red, p):
    """Rational functional ell on the sector monomials from modular lambda over standard coordinates (via the reducer)."""
    # ell[std_r] = lam_r ; ell[piv_i] = - sum_r R[i, std_r] lam_r
    M = len(red.mons)
    ell = [0] * M
    for r, v in lam.items():
        ell[red.std[r]] = v % p
    if red.rank:
        lamvec = nmod_mat(red.nstd, 1, [lam.get(r, 0) % p for r in range(red.nstd)], p)
        prod = red.Rstd * lamvec   # rank x 1
        for i, j in enumerate(red.pivots):
            ell[j] = (-int(prod[i, 0])) % p
    out = {}
    for j, v in enumerate(ell):
        if v:
            r_ = core.ratrecon(v, p)
            if r_ is None: return None
            out[red.mons[j]] = r_
    return out

def apply_functional(ell, P):
    return sum((ell[e] * v for e, v in P.items() if e in ell), Fraction(0))

def verify_certificate(ell, F, mons_low, red, EB, nt, ns, a, zero, target_poly):
    """Exact checks: ell(F*m) = 0 for all m in mons_low; sum over chart pairs of ell(B) * monomial == target_poly (dict var-monomial -> Fraction)."""
    for m in mons_low:
        if apply_functional(ell, shift(F, m, 1)) != 0:
            return False, 'ell(F m) != 0'
    acc = {}
    for b in range(nt):
        for kind, n in (('t', nt), ('s', ns)):
            for c in range(n):
                if kind == 't' and c < b: continue  # symmetric pairs handled once below
                pass
    # explicit enumeration of chart monomials
    def var_ok(kind, i):
        if kind == 't' and i < a: return False
        return (kind, i) not in zero
    for b in range(a, nt):
        if not var_ok('t', b): continue
        # t-t pairs (b <= c): coefficient of t_b t_c is ell(B(b,t,c)) + ell(B(c,t,b)) for b != c, ell(B(b,t,b)) for b == c
        for c in range(b, nt):
            if not var_ok('t', c): continue
            val = apply_functional(ell, EB.B(b, 't', c)) + (apply_functional(ell, EB.B(c, 't', b)) if c != b else 0)
            if val == 0: continue
            vs = [x for x in (('t', b), ('t', c)) if x[1] != a]
            mon = tuple(sorted(vs))
            acc[mon] = acc.get(mon, 0) + val
        for k in range(ns):
            if not var_ok('s', k): continue
            val = apply_functional(ell, EB.B(b, 's', k))
            if val == 0: continue
            vs = ([('t', b)] if b != a else []) + [('s', k)]
            mon = tuple(sorted(vs))
            acc[mon] = acc.get(mon, 0) + val
    acc = {m: v for m, v in acc.items() if v}
    return acc == {m: v for m, v in target_poly.items() if v}, acc

def certify_chart(sysd, red, EB, F, mons_low, a, verbose=False, max_steps=60):
    """Try constant certificate; else cascade of square/linear implications.  Returns record."""
    p = sysd['p']; nt, ns = sysd['nt'], sysd['ns']
    zero = set(); steps = []
    for step in range(max_steps):
        eqs = chart_equations_modp(sysd, a, zero)
        eqs = [(r, P) for r, P in eqs if P]
        # constant certificate
        lam = find_combination(eqs, {(): 1}, p)
        if lam is not None:
            ell = functional_from_lambda(lam, red, p)
            if ell is None: return {'cert': None, 'steps': steps, 'fail': 'ell reconstruction', 'verified_exact': False}
            ok, info = verify_certificate(ell, F, mons_low, red, EB, nt, ns, a, zero, {(): Fraction(1)})
            if ok:
                steps.append({'type': 'constant', 'support': len(lam), 'ell_terms': len(ell)})
                return {'cert': 'constant' if len(steps) == 1 else 'cascade', 'steps': steps, 'verified_exact': True, 'zero_vars': sorted(zero)}
            return {'cert': None, 'steps': steps, 'fail': 'constant certificate failed exact check', 'verified_exact': False}
        # square implications v^2 in span, then linear forms
        progressed = False
        for v in chart_variables(nt, ns, a, zero):
            lam = find_combination(eqs, {(v, v): 1}, p)
            if lam is None: continue
            ell = functional_from_lambda(lam, red, p)
            if ell is None: continue
            ok, info = verify_certificate(ell, F, mons_low, red, EB, nt, ns, a, zero, {(v, v): Fraction(1)})
            if ok:
                zero.add(v); steps.append({'type': 'square', 'var': v, 'support': len(lam)}); progressed = True
                break
        if progressed: continue
        for v in chart_variables(nt, ns, a, zero):
            lam = find_combination(eqs, {(v,): 1}, p)
            if lam is None: continue
            ell = functional_from_lambda(lam, red, p)
            if ell is None: continue
            ok, info = verify_certificate(ell, F, mons_low, red, EB, nt, ns, a, zero, {(v,): Fraction(1)})
            if ok:
                zero.add(v); steps.append({'type': 'linear', 'var': v, 'support': len(lam)}); progressed = True
                break
        if not progressed:
            return {'cert': None, 'steps': steps, 'fail': 'no constant/square/linear certificate', 'verified_exact': False, 'zero_vars': sorted(zero)}
    return {'cert': None, 'steps': steps, 'fail': 'max steps', 'verified_exact': False}
