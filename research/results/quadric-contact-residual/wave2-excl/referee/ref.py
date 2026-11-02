"""Referee replay (own implementation) for the quadric-contact-residual wave-2 exclusion packet.

Prime: REF_P (default 1000000021, not used by any earlier run).  Exact rational parts use python-flint fmpq_mat/fmpz_mat.
Polynomials: dict {exponent6: coeff}; variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5.  Only hodge.core (polynomial
arithmetic + instance loader) is shared with the packet; sectors, families, relaxations and solvers are re-implemented.
Usage: python3 ref.py <cmd> [args]   (cmds: hilbert, e1, d5exact, d8star, d6, d7pair s t, d8sector chi, d8inv,
                                       fixedD-exact, served-exact NAME [NAME...], served-exact-full NAME)
"""
import sys, os, json, time, subprocess, re, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
os.environ.setdefault('HODGE_QUIET', '1')
from fractions import Fraction
from hodge import core
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, power, shift
import flint
from flint import nmod_mat, fmpz_mat, fmpq_mat

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.environ.get('REF_LOG', os.path.join(HERE, 'referee.log.jsonl'))
P = int(os.environ.get('REF_P', '1000000021'))
F, A, Q_ORDER, _rec = core.load_instance('quadric-contact-residual')
assert Q_ORDER == 1
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]
FiQ = [pderiv(F, i) for i in range(NV)]
HessQ = [[pderiv(FiQ[i], j) for j in range(NV)] for i in range(NV)]
Hess = [[pderiv(Fi[i], j, P) for j in range(NV)] for i in range(NV)]

def unit(i): return tuple(1 if j == i else 0 for j in range(NV))
def chi(e): return ((e[0] + e[3]) % 2, (e[2] + e[5]) % 2, (e[1] - e[4]) % 6, (e[3] + e[4] + e[5]) % 2)
def cadd(a, b): return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2, (a[2] + b[2]) % 6, (a[3] + b[3]) % 2)
def cneg(a): return ((-a[0]) % 2, (-a[1]) % 2, (-a[2]) % 6, (-a[3]) % 2)
def chi_of(Pl):
    cs = {chi(e) for e in Pl}
    return cs.pop() if len(cs) == 1 else None
CHI_Z = [chi(unit(i)) for i in range(NV)]
CHI_A = chi_of(A); CHI0 = (0, 0, 0, 0)
assert chi_of(F) == CHI0 and CHI_A == (0, 0, 0, 1)
ALL_CHARS = sorted({chi(e) for e in monomials(12)})

def log(exp, **kw):
    r = {'exp': exp, 'prime': P, 'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'who': 'referee'}
    r.update(kw)
    with open(LOG, 'a') as f: f.write(json.dumps(r) + '\n')
    print('[LOG]', json.dumps(r)[:600], file=sys.stderr, flush=True)
    return r

def ci_hilbert(m):
    """Hilbert function of a complete intersection of six quintics in six variables."""
    poly = [1]
    for _ in range(6):
        new = [0] * (len(poly) + 4)
        for i, c in enumerate(poly):
            for k in range(5): new[i + k] += c
        poly = new
    return poly[m] if m < len(poly) else 0

# ------------------------------------------------------------------ modular sectors
class JSec:
    """degree-m monomials of G_aff character c; the subspace J_m[c] spanned by m'*F_i; normal forms mod p."""
    def __init__(self, m, c, p=P):
        self.m, self.c, self.p = m, c, p
        self.mons = [e for e in monomials(m) if chi(e) == c]
        self.idx = {e: i for i, e in enumerate(self.mons)}
        self.n = len(self.mons)
        rows = []
        Fip = Fi if p == P else [pderiv(core.reduce_mod(F, p), i, p) for i in range(NV)]
        if m >= 5:
            for mm in monomials(m - 5):
                cm = chi(mm)
                for i in range(NV):
                    if cadd(cm, cneg(CHI_Z[i])) == c:
                        rows.append(shift(Fip[i], mm, 1, p))
        self.nrows = len(rows)
        if rows:
            flat = [0] * (len(rows) * self.n)
            for r, Pl in enumerate(rows):
                for e, v in Pl.items(): flat[r * self.n + self.idx[e]] = v
            M = nmod_mat(len(rows), self.n, flat, p)
            R, rk = M.rref()
            self.rank = rk
            piv = []
            for i in range(rk):
                j = 0
                while int(R[i, j]) == 0: j += 1
                piv.append(j)
            self.piv = piv
            self.R = nmod_mat(rk, self.n, [int(R[i, j]) for i in range(rk) for j in range(self.n)], p) if rk else None
        else:
            self.rank, self.piv, self.R = 0, [], None
        ps = set(self.piv)
        self.nonpiv = [j for j in range(self.n) if j not in ps]
        self.dimR = self.n - self.rank
    def vec(self, Pl):
        v = [0] * self.n
        for e, cv in Pl.items(): v[self.idx[e]] = cv % self.p
        return v
    def nf(self, vecs):
        """normal forms (restricted to non-pivot coordinates) of a list of coefficient vectors: k x dimR"""
        k = len(vecs)
        if k == 0: return nmod_mat(0, len(self.nonpiv), [], self.p)
        V = nmod_mat(k, self.n, [x for v in vecs for x in v], self.p)
        if self.rank:
            Pm = nmod_mat(k, self.rank, [v[j] for v in vecs for j in self.piv], self.p)
            V = V - Pm * self.R
        return nmod_mat(k, len(self.nonpiv), [int(V[i, j]) for i in range(k) for j in self.nonpiv], self.p)

_SEC = {}
def sec(m, c, p=P):
    if (m, c, p) not in _SEC: _SEC[(m, c, p)] = JSec(m, c, p)
    return _SEC[(m, c, p)]

# ------------------------------------------------------------------ vector fields
def vfield(N, Pl, p=None):
    out = {}
    for i in range(NV):
        if N[i]: out = add(out, mul(N[i], pderiv(Pl, i, p), p), p)
    return out
def divg(N, p=None):
    out = {}
    for i in range(NV): out = add(out, pderiv(N[i], i, p), p)
    return out
def BD(D, N, p=None):
    return sub(mul(D, divg(N, p), p), vfield(N, D, p), p)

def relaxation_quadrics(Ds, Ns, d, p=P):
    """Hessian relaxation: D = sum c_k D_k, N = sum c_k N_k; B_D(N) in J_{2d} gives quadrics in c (one per coordinate
    of R_{2d}[char]).  Returns list of {(k,l): coeff} with k <= l."""
    n = len(Ds)
    cD = [chi_of(D) for D in Ds]; assert all(c is not None for c in cD)
    divN = [divg(N, p) for N in Ns]
    groups = {}
    for k in range(n):
        for l in range(k, n):
            groups.setdefault(cadd(CHI_A, cadd(cD[k], cD[l])), []).append((k, l))
    quads = []
    for cc, pairs in groups.items():
        S = sec(2 * d, cc, p)
        vecs = []
        for (k, l) in pairs:
            B = sub(mul(Ds[k], divN[l], p), vfield(Ns[l], Ds[k], p), p)
            if l != k:
                B = add(B, sub(mul(Ds[l], divN[k], p), vfield(Ns[k], Ds[l], p), p), p)
            vecs.append(S.vec(B))
        NF = S.nf(vecs)
        for j in range(NF.ncols()):
            cf = {pairs[r]: int(NF[r, j]) for r in range(len(pairs)) if int(NF[r, j])}
            if cf: quads.append(cf)
    return quads

def quad_rank(quads, n, p=P):
    pairs = [(k, l) for k in range(n) for l in range(k, n)]
    pidx = {pr: i for i, pr in enumerate(pairs)}
    if not quads: return 0
    flat = [0] * (len(quads) * len(pairs))
    for r, cf in enumerate(quads):
        for pr, v in cf.items(): flat[r * len(pairs) + pidx[pr]] = v
    return nmod_mat(len(quads), len(pairs), flat, p).rank()

def write_ms_sing(quads, keep, prefix, p=P):
    """Write msolve (c^2 syntax!) and Singular inputs for the quadrics restricted to the unknowns in `keep`."""
    ks = set(keep); names = [f'c{k}' for k in keep]
    eqs = []
    for cf in quads:
        terms = [(f'{v}*c{k}^2' if k == l else f'{v}*c{k}*c{l}') for (k, l), v in cf.items() if k in ks and l in ks]
        if terms: eqs.append('+'.join(terms))
    with open(prefix + '.ms', 'w') as f:
        f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(eqs) + '\n')
    with open(prefix + '.sing', 'w') as f:
        f.write(f'ring r = {p}, ({",".join(names)}), dp;\nideal I = {",".join(eqs) if eqs else "0"};\n')
        f.write('ideal G = std(I);\nprint("SDIM " + string(dim(G)));\nif (dim(G) == 0) { print("SVDIM " + string(vdim(G))); }\nquit;\n')
    return len(eqs)

def msolve_lead(prefix, threads=2, timeout=600):
    out = prefix + '.lead'
    t0 = time.time()
    r = subprocess.run(['timeout', str(timeout), 'msolve', '-v', '0', '-t', str(threads), '-g', '1', '-f', prefix + '.ms', '-o', out], capture_output=True, text=True)
    res = {'rc': r.returncode, 'seconds': round(time.time() - t0, 1)}
    if r.returncode != 0 or not os.path.exists(out):
        res['verdict'] = 'msolve-failed-or-timeout'; return res
    txt = open(out).read()
    vars_ = open(prefix + '.ms').readline().strip().split(',')
    m = re.search(r'\[(.*)\]:', txt, re.S)
    lead = [s.strip() for s in m.group(1).replace('\n', '').split(',') if s.strip()]
    powers = {}
    for mono in lead:
        mm = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', mono)
        if mm: powers[mm.group(1)] = min(powers.get(mm.group(1), 10**9), int(mm.group(2) or 1))
    missing = [v for v in vars_ if v not in powers]
    res.update({'n_lead': len(lead), 'missing_variables': missing, 'max_pure_power': max(powers.values()) if powers else None,
                'verdict': 'ONLY-ORIGIN (m-primary leading ideal)' if not missing else 'POSITIVE-DIMENSIONAL'})
    return res

def singular_dim(prefix, timeout=600):
    t0 = time.time()
    r = subprocess.run(['timeout', str(timeout), 'Singular', '-q', prefix + '.sing'], capture_output=True, text=True)
    sd = re.search(r'SDIM (-?\d+)', r.stdout); sv = re.search(r'SVDIM (\d+)', r.stdout)
    return {'sdim': int(sd.group(1)) if sd else None, 'svdim': int(sv.group(1)) if sv else None, 'rc': r.returncode, 'seconds': round(time.time() - t0, 1)}

# ------------------------------------------------------------------ exact building blocks (Fractions)
def Xv(i): return {unit(i): Fraction(1)}
X0, X1, X2, Y0, Y1, Y2 = [Xv(i) for i in range(NV)]
S0 = add(mul(X0, X0), mul(Y0, Y0)); P0 = mul(X0, Y0); DELTA = sub(mul(X0, X0), mul(Y0, Y0))
S1 = add(mul(X1, X1), mul(Y1, Y1)); P1 = mul(X1, Y1); S2 = add(mul(X2, X2), mul(Y2, Y2)); P2 = mul(X2, Y2)
Q1 = add(add(mul(P1, P1), mul(P1, P2)), mul(P2, P2))
Q2 = mul(P1, add(mul(X2, X2), scale(mul(Y2, Y2), 2)))
GS = add(add(scale(mul(S0, S0), 3), scale(mul(P0, P0), 2)), scale(Q1, 5))
GP = add(scale(mul(P0, S0), 4), scale(Q2, Fraction(5, 7)))
D8 = sub(scale(mul(GS, GS), 4), mul(GP, GP))
NSTAR = [{} for _ in range(NV)]
_h = scale(mul(mul(X0, X0), mul(Y0, Y0)), Fraction(1, 2))
NSTAR[0] = mul(_h, sub(scale(mul(Y0, GS), 2), mul(X0, GP)))
NSTAR[3] = mul(_h, sub(scale(mul(X0, GS), 2), mul(Y0, GP)))
f0 = add(power(S0, 3), scale(mul(mul(P0, P0), S0), 2))
rad = [add(mul(Xv(i), FiQ[i]), mul(Xv(i + 3), FiQ[i + 3])) for i in range(3)]
BLOCKS = {'X0': X0, 'Y0': Y0, 'X1': X1, 'Y1': Y1, 'X2': X2, 'Y2': Y2,
          'S0': S0, 'delta': DELTA, 'S1': S1, 'S2': S2, 'S0+S2': add(S0, S2), 'S0-S2': sub(S0, S2), 'X0^2+X2^2': add(mul(X0, X0), mul(X2, X2)),
          'S0+S1': add(S0, S1), 'X0^2-X1^2': sub(mul(X0, X0), mul(X1, X1)), 'X1^2-Y1^2': sub(mul(X1, X1), mul(Y1, Y1)), 'X2^2-Y2^2': sub(mul(X2, X2), mul(Y2, Y2)),
          'GS': GS, 'GP': GP, 'Q1': Q1, 'T0': add(power(X0, 4), power(Y0, 4)), 'S0^2+2P0^2': add(mul(S0, S0), scale(mul(P0, P0), 2)), 'Q2': Q2,
          'X1^4+Y1^4': add(power(X1, 4), power(Y1, 4)),
          'F_X0': FiQ[0], 'F_Y0': FiQ[3], 'F_X1': FiQ[1], 'F_Y1': FiQ[4], 'F_X2': FiQ[2], 'F_Y2': FiQ[5],
          'f0': f0, 'rad0': rad[0], 'rad1': rad[1], 'rad2': rad[2], 'X1^6+Y1^6': add(power(X1, 6), power(Y1, 6)), 'X2^6+Y2^6': add(power(X2, 6), power(Y2, 6)),
          'D8*': D8}

def fmpz_rows(rows, ncols, idx):
    """rows: list of dict{key: Fraction}; clear denominators per row; return fmpz_mat"""
    flat = [0] * (len(rows) * ncols)
    for r, row in enumerate(rows):
        den = 1
        for v in row.values(): den = den * v.denominator // __import__('math').gcd(den, v.denominator)
        for k, v in row.items(): flat[r * ncols + idx[k]] = int(v * den)
    return fmpz_mat(len(rows), ncols, flat)
def exact_rank(M):
    return fmpq_mat(M).rref()[1]

# ------------------------------------------------------------------ commands
def cmd_hilbert():
    """dim J_m mod p for m = 5..18 vs the complete-intersection Hilbert function (all 48 characters), for several primes."""
    primes = [int(x) for x in os.environ.get('REF_PRIMES', str(P)).split(',') if x]
    for p in primes:
        t0 = time.time(); bad = []; dims = {}
        for m in range(5, 19):
            chars_m = sorted({chi(e) for e in monomials(m)})
            tot = sum(sec(m, c, p).rank for c in chars_m)
            expect = len(monomials(m)) - ci_hilbert(m)
            dims[m] = tot
            if tot != expect: bad.append((m, tot, expect))
        log('hilbert-J', prime_checked=p, dimJ=dims, bad=bad, ok=(not bad), seconds=round(time.time() - t0, 1),
            note='dim J_m mod p equals dim S_m - CI Hilbert function => J_m mod p is the reduction of the saturated rational lattice (used by every mod-p transfer)')
        _SEC.clear()

def cmd_e1():
    known = [0, 0, 0, 0, 0, 6, 36, 126, 337, 760]
    out = {}
    for d in range(0, 10):
        t0 = time.time(); colon = 0; dimJ = 0; extras = {}
        for c in sorted({chi(e) for e in monomials(d)}):
            mons = [e for e in monomials(d) if chi(e) == c]
            S6 = sec(d + 6, cadd(c, CHI_A))
            NF = S6.nf([S6.vec(shift(Ap, e, 1, P)) for e in mons])
            nul = len(mons) - NF.rank()
            dj = sec(d, c).rank
            colon += nul; dimJ += dj
            if nul > dj: extras[str(c)] = nul - dj
        out[d] = colon
        log('E1-colon', d=d, dim_colon=colon, dim_J=dimJ, extras_per_sector=extras, agrees=(colon == known[d]), seconds=round(time.time() - t0, 1))
    return out

def cmd_d5exact():
    """exact over Q: the 21 columns A*d_i d_j F are independent modulo J_10 (sector-wise fmpq rank)."""
    t0 = time.time()
    cols = {}
    for i in range(NV):
        for j in range(i, NV):
            cols.setdefault(cadd(CHI_A, cadd(cneg(CHI_Z[i]), cneg(CHI_Z[j]))), []).append(mul(A, HessQ[i][j]))
    total = 0; detail = {}
    for c, cl in cols.items():
        mons = [e for e in monomials(10) if chi(e) == c]; idx = {e: i for i, e in enumerate(mons)}
        jrows = [shift(FiQ[i], mm, 1) for mm in monomials(5) for i in range(NV) if cadd(chi(mm), cneg(CHI_Z[i])) == c]
        rJ = exact_rank(fmpz_rows(jrows, len(mons), idx)); rAll = exact_rank(fmpz_rows(jrows + cl, len(mons), idx))
        detail[str(c)] = (len(cl), rAll - rJ); total += rAll - rJ
    log('d5-hessian-exact', n_columns=21, independent_mod_J10_over_Q=total, per_sector=detail, result='EXCLUDED over C' if total == 21 else 'NOT CLOSED',
        seconds=round(time.time() - t0, 1), note='exact fmpq ranks; with (J:A)_5 = J_5 this proves: no degree-5 denominator (c_i c_j = 0 for all i,j)')

def cmd_d8star():
    """exact identities: N*(F) = A*D8*, D8*(X0^2-Y0^2) = F_X0^2 - F_Y0^2, and D8* not in J_8 (mod p + Hilbert)."""
    lhs = vfield(NSTAR, F); ok1 = (sub(lhs, mul(A, D8)) == {})
    ok2 = (sub(mul(D8, DELTA), sub(mul(FiQ[0], FiQ[0]), mul(FiQ[3], FiQ[3]))) == {})
    S = sec(8, CHI0); nf = S.nf([S.vec(core.reduce_mod(D8, P))]); notinJ = any(int(nf[0, j]) for j in range(nf.ncols()))
    log('D8star-identities', Nstar_F_equals_A_D8=ok1, D8_times_delta_equals_FX0sq_minus_FY0sq=ok2, D8_not_in_J8_mod_p=notinJ, n_terms=len(D8))

def jac_family(d, char_filter=None, gauge_drop=None):
    """Jacobian family D = m F_j (m in S_{d-5}); C = m e_j; returns (Ds, Ns, names, chars)."""
    Ds, Ns, names = [], [], []
    for m in monomials(d - 5):
        for j in range(NV):
            c = cadd(chi(m), cneg(CHI_Z[j]))
            if char_filter and not char_filter(c): continue
            if gauge_drop and gauge_drop(m, j): continue
            Ds.append(shift(Fi[j], m, 1, P)); N = [{} for _ in range(NV)]; N[j] = {m: Ap[(3, 0, 0, 3, 0, 0)]}
            N[j] = shift(Ap, m, 1, P); Ns.append(N); names.append((m, j))
    return Ds, Ns, names

def cmd_d6():
    t0 = time.time()
    Ds, Ns, names = jac_family(6, gauge_drop=lambda m, j: (m == unit(5) and j == 5))
    assert len(Ds) == 35
    quads = relaxation_quadrics(Ds, Ns, 6)
    rk = quad_rank(quads, 35)
    prefix = os.path.join(HERE, 'd6')
    neq = write_ms_sing(quads, list(range(35)), prefix)
    v = msolve_lead(prefix); s = singular_dim(prefix)
    log('closure-d6', n_unknowns=35, n_quadrics=neq, quadric_rank=rk, msolve=v, singular=s,
        result='EXCLUDED over C' if (v.get('verdict', '').startswith('ONLY') and s['sdim'] == 0) else 'OPEN', seconds=round(time.time() - t0, 1),
        note='gauge: coefficient of Y2*F_Y2 dropped (F = (1/6) sum Z_i F_i has coefficient 1/6 there)')

def d7_family():
    Ds, Ns, names = jac_family(7, gauge_drop=lambda m, j: (j == 5 and m[5] >= 1))
    assert len(Ds) == 120
    chars = [cadd(chi(m), cneg(CHI_Z[j])) for m, j in names]
    return Ds, Ns, names, chars

_D7 = {}
def d7_quads():
    if 'q' not in _D7:
        Ds, Ns, names, chars = d7_family()
        quads = relaxation_quadrics(Ds, Ns, 7)
        _D7['q'] = (quads, names, chars, quad_rank(quads, 120))
        json.dump({'p': P, 'names': [[list(m), j] for m, j in names], 'chi': [list(c) for c in chars], 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads]},
                  open(os.path.join(HERE, 'd7mixed.quads.json'), 'w'))
    return _D7['q']

def cmd_d7pair(s, t, threads=2, timeout=900):
    t0 = time.time()
    quads, names, chars, rk = d7_quads()
    par = lambda c: (c[0], c[1], c[3])
    S = {tuple(int(x) for x in s), tuple(int(x) for x in t)}
    keep = [k for k, c in enumerate(chars) if par(c) in S]
    prefix = os.path.join(HERE, f'd7_pair_{s}_{t}')
    neq = write_ms_sing(quads, keep, prefix)
    sub_q = [{(k, l): v for (k, l), v in cf.items() if k in set(keep) and l in set(keep)} for cf in quads]
    sub_q = [cf for cf in sub_q if cf]
    kidx = {k: i for i, k in enumerate(keep)}
    srk = quad_rank([{(kidx[k], kidx[l]): v for (k, l), v in cf.items()} for cf in sub_q], len(keep))
    v = msolve_lead(prefix, threads, timeout)
    sg = singular_dim(prefix, timeout) if len(keep) <= 28 else None
    diff = sum(1 for a, b in zip(sorted(S)[0], sorted(S)[1]) if a != b) if len(S) == 2 else 0
    kind = {0: 'vertex', 1: 'edge', 2: 'face diagonal', 3: 'space diagonal'}[diff]
    log('d7-parity-pair', pair=[s, t], kind=kind, n_unknowns=len(keep), n_quadrics=neq, independent_quadrics=srk, full_system_rank=rk, msolve=v, singular=sg,
        result='EXCLUDED over C' if v.get('verdict', '').startswith('ONLY') else v.get('verdict'), seconds=round(time.time() - t0, 1))

def cmd_d8sector(tag, threads=2, timeout=900):
    """pure G_aff sector chi at d=8: family = gauge (s F, s in S_2[chi]) + J_8[chi] (+ extras if any); gauge coefficients set to 0."""
    t0 = time.time()
    c = tuple(int(x) for x in tag)
    Sd = sec(8, c)
    Ds, Ns, vecs, kinds = [], [], [], []
    def try_add(D, N, kind):
        v = Sd.vec(D)
        M = nmod_mat(len(vecs) + 1, Sd.n, [x for u in vecs + [v] for x in u], P)
        if M.rank() == len(vecs) + 1:
            vecs.append(v); Ds.append(D); Ns.append(N); kinds.append(kind); return True
        return False
    inv6 = pow(6, P - 2, P)
    for s in monomials(2):
        if chi(s) != c: continue
        N = [shift(Ap, tuple(s[j] + (1 if j == i else 0) for j in range(NV)), inv6, P) for i in range(NV)]
        try_add(shift(Fp, s, 1, P), N, 'gauge')
    ng = len(Ds)
    for m in monomials(3):
        for j in range(NV):
            if cadd(chi(m), cneg(CHI_Z[j])) != c: continue
            N = [{} for _ in range(NV)]; N[j] = shift(Ap, m, 1, P)
            try_add(shift(Fi[j], m, 1, P), N, 'jac')
    assert len(Ds) == Sd.rank, (len(Ds), Sd.rank)
    # E1 kernel in this sector
    S14 = sec(14, cadd(c, CHI_A))
    NF = S14.nf([S14.vec(shift(Ap, e, 1, P)) for e in Sd.mons])
    e1 = Sd.n - NF.rank()
    extra = e1 - Sd.rank
    if extra:
        assert c == CHI0, 'extras only expected in the trivial sector'
        try_add(core.reduce_mod(D8, P), [core.reduce_mod(x, P) for x in NSTAR], 'extra')
        assert len(Ds) == e1
    n = len(Ds)
    quads = relaxation_quadrics(Ds, Ns, 8)
    keep = list(range(ng, n))
    prefix = os.path.join(HERE, f'd8_chi{tag}')
    neq = write_ms_sing(quads, keep, prefix)
    kidx = {k: i for i, k in enumerate(keep)}
    sub_q = [{(kidx[k], kidx[l]): v for (k, l), v in cf.items() if k in kidx and l in kidx} for cf in quads]
    sub_q = [cf for cf in sub_q if cf]
    rk = quad_rank(sub_q, len(keep))
    v = msolve_lead(prefix, threads, timeout) if neq else {'verdict': 'NO EQUATIONS'}
    log('closure-pure-sector', d=8, chi=str(c), e1_kernel=e1, n_gauge=ng, n_params=n, n_unknowns=len(keep), n_quadrics=neq, quadric_rank=rk, msolve=v,
        result='EXCLUDED over C' if v.get('verdict', '').startswith('ONLY') else v.get('verdict'), seconds=round(time.time() - t0, 1))

def perm_poly(Pl, g, p=P):
    return {tuple(e[g[j]] for j in range(NV)): c for e, c in Pl.items()}
def perm_field(N, g, p=P):
    return [perm_poly(N[g[i]], g, p) for i in range(NV)]
GRP = [(0, 1, 2, 3, 4, 5), (3, 1, 2, 0, 4, 5), (0, 4, 2, 3, 1, 5), (3, 4, 2, 0, 1, 5)]   # id, s0, r1, s0 r1

def cmd_d8inv(threads=2, timeout=900):
    """Gamma-invariant family at d=8 (G_aff-trivial character, s0- and r1-invariant): gauge (3) + Jacobian + D8*."""
    t0 = time.time()
    inv4 = pow(4, P - 2, P); inv6 = pow(6, P - 2, P)
    def symP(Pl):
        out = {}
        for g in GRP: out = add(out, perm_poly(Pl, g), P)
        return scale(out, inv4, P)
    def symN(N):
        out = [{} for _ in range(NV)]
        for g in GRP:
            gN = perm_field(N, g); out = [add(out[i], gN[i], P) for i in range(NV)]
        return [scale(x, inv4, P) for x in out]
    assert symP(Fp) == Fp
    S8 = sec(8, CHI0)
    Ds, Ns, vecs, kinds = [], [], [], []
    def try_add(D, N, kind):
        v = S8.vec(D)
        M = nmod_mat(len(vecs) + 1, S8.n, [x for u in vecs + [v] for x in u], P)
        if M.rank() == len(vecs) + 1:
            vecs.append(v); Ds.append(D); Ns.append(N); kinds.append(kind); return True
        return False
    for s in monomials(2):
        if chi(s) != CHI0: continue
        sP = symP({s: 1})
        N = [shift(mul(Ap, sP, P), unit(i), inv6, P) for i in range(NV)]
        try_add(mul(sP, Fp, P), N, 'gauge')
    ng = len(Ds)
    for m in monomials(3):
        for j in range(NV):
            if cadd(chi(m), cneg(CHI_Z[j])) != CHI0: continue
            C = [{} for _ in range(NV)]; C[j] = {m: 1}
            Cs = symN(C)
            D = {}
            for i in range(NV): D = add(D, mul(Cs[i], Fi[i], P), P)
            if D: try_add(D, [mul(Ap, x, P) for x in Cs], 'jac')
    nj = len(Ds) - ng
    assert symP(core.reduce_mod(D8, P)) == core.reduce_mod(D8, P)
    try_add(core.reduce_mod(D8, P), [core.reduce_mod(x, P) for x in NSTAR], 'extra')
    n = len(Ds)
    for D, N in zip(Ds, Ns):
        assert sub(vfield(N, Fp, P), mul(Ap, D, P), P) == {}, 'E1 lift check'
        assert symP(D) == D
    quads = relaxation_quadrics(Ds, Ns, 8)
    keep = list(range(ng, n)); kidx = {k: i for i, k in enumerate(keep)}
    prefix = os.path.join(HERE, 'd8_inv')
    neq = write_ms_sing(quads, keep, prefix)
    sub_q = [{(kidx[k], kidx[l]): v for (k, l), v in cf.items() if k in kidx and l in kidx} for cf in quads]
    sub_q = [cf for cf in sub_q if cf]
    rk = quad_rank(sub_q, len(keep))
    v = msolve_lead(prefix, threads, timeout)
    log('closure-Gamma-invariant', d=8, n_gauge=ng, n_jac=nj, n_extra=n - ng - nj, n_params=n, n_unknowns=len(keep), n_quadrics=neq, quadric_rank=rk, msolve=v,
        result='EXCLUDED over C' if v.get('verdict', '').startswith('ONLY') else v.get('verdict'), seconds=round(time.time() - t0, 1))

# ------------------------------------------------------------------ exact fixed-D / served (sector-wise, over Q)
def fixedD_exact_sector(D, tchar, with_targets=False):
    """Exact (fmpq) sector system for denominator D (G_aff eigenvector) and target character tchar:
    unknowns N_i in S_{d+1}[tchar + chi_D + chi(Z_i)], Q in S_{2d-6}[tchar + 2 chi_D]; equations in S_{2d+6}[tchar + 2 chi_D].
    Returns (rank M, rank [M|B] or rank [M|b]) where B = target columns (all monomials of S_6[tchar], Y-odd part) or b = A."""
    d = core.degree(D); cD = chi_of(D); assert cD is not None
    eqc = cadd(tchar, cadd(cD, cD))
    rows_mons = [e for e in monomials(2 * d + 6) if chi(e) == eqc]; idx = {e: i for i, e in enumerate(rows_mons)}
    cols = []
    Di = [pderiv(D, i) for i in range(NV)]; FD = mul(F, D); F2 = mul(F, F)
    for i in range(NV):
        ci = cadd(tchar, cadd(cD, CHI_Z[i]))
        G = scale(add(mul(D, FiQ[i]), mul(F, Di[i])), -1)          # - D F_i m - F D_i m  (+ F D d_i m)
        for m in monomials(d + 1):
            if chi(m) != ci: continue
            col = shift(G, m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i]))
            cols.append(col)
    if 2 * d - 6 >= 0:
        for m in monomials(2 * d - 6):
            if chi(m) == eqc: cols.append(shift(scale(F2, -1), m, 1))
    D2 = mul(D, D)
    if with_targets:
        tcols = [shift(D2, e, 1) for e in monomials(6, 1) if chi(e) == tchar]
    else:
        assert tchar == CHI_A
        tcols = [mul(A, D2)]
    for col in cols + tcols:
        assert all(chi(e) == eqc for e in col)
    # exact ranks via columns-as-rows (rank is transpose-invariant)
    rM = exact_rank(fmpz_rows(cols, len(rows_mons), idx)) if cols else 0
    rMB = exact_rank(fmpz_rows(cols + tcols, len(rows_mons), idx))
    return rM, rMB, len(cols), len(tcols), len(rows_mons)

def cmd_fixedD_exact(name='D8*'):
    t0 = time.time(); D = BLOCKS[name]
    rM, rMb, ncols, _, nrows = fixedD_exact_sector(D, CHI_A, with_targets=False)
    log('fixedD-exact-Q', D=name, d=core.degree(D), sector=str(CHI_A), rows=nrows, cols=ncols, rank_M=rM, rank_augmented=rMb,
        consistent_over_Q=(rM == rMb), nullity=ncols - rM, seconds=round(time.time() - t0, 1),
        note='exact rational rank (fmpq_mat.rref) of the G_aff-sector-restricted fixed-D system for target A; sector restriction is lossless (averaging over G_aff)')

def cmd_served_exact(names):
    """served_Q(D) for all 24 Y-odd target characters, exactly: dim served[chi'] = dim S_6^odd[chi'] - (rank[M|B] - rank M)."""
    odd_chars = sorted({chi(e) for e in monomials(6, 1)})
    for name in names:
        t0 = time.time(); D = BLOCKS[name]
        if chi_of(D) is None:
            log('served-exact-Q', D=name, result='skipped (not a G_aff eigenvector)'); continue
        per = {}; tot = 0; totJ = 0
        for tc in odd_chars:
            rM, rMB, ncols, ntc, nrows = fixedD_exact_sector(D, tc, with_targets=True)
            dim_served = ntc - (rMB - rM)
            dimJ = sec(6, tc).rank if tc in {chi(e) for e in monomials(6, 1)} else 0
            # J_6^odd[tc]: Jacobian part of character tc (all of J_6[tc] is Y-odd iff tc[3]==1 since F is Y-even)
            per[str(tc)] = [dim_served, dimJ]; tot += dim_served; totJ += dimJ
        log('served-exact-Q', D=name, d=core.degree(D), dim_served_Q=tot, dim_J6odd=totJ, per_char=per,
            result='served(D) = J_6^odd (trivial only)' if tot == totJ else f'NONTRIVIAL served targets: {tot - totJ}', seconds=round(time.time() - t0, 1))

def cmd_served_exact_full(name):
    """full (non-sector) exact served space for a (possibly mixed-character) block"""
    t0 = time.time(); D = BLOCKS[name]; d = core.degree(D)
    rows_mons = [e for e in monomials(2 * d + 6) if core.yparity(e) == 1]; idx = {e: i for i, e in enumerate(rows_mons)}
    Di = [pderiv(D, i) for i in range(NV)]; FD = mul(F, D); F2 = mul(F, F)
    cols = []
    for i in range(NV):
        G = scale(add(mul(D, FiQ[i]), mul(F, Di[i])), -1)
        par_i = 1 if i >= 3 else 0   # N_i must have Y-parity so that N_i F_i is Y-odd: F_i is Y-odd for i>=3
        for m in monomials(d + 1):
            if core.yparity(m) != (1 - par_i): continue
            col = shift(G, m, 1)
            if m[i]:
                m2 = list(m); m2[i] -= 1; col = add(col, shift(FD, tuple(m2), m[i]))
            cols.append(col)
    for m in monomials(2 * d - 6, 1): cols.append(shift(scale(F2, -1), m, 1))
    D2 = mul(D, D); tcols = [shift(D2, e, 1) for e in monomials(6, 1)]
    for col in cols + tcols: assert all(core.yparity(e) == 1 for e in col)
    rM = exact_rank(fmpz_rows(cols, len(rows_mons), idx)); rMB = exact_rank(fmpz_rows(cols + tcols, len(rows_mons), idx))
    dim_served = len(tcols) - (rMB - rM)
    log('served-exact-Q-full', D=name, d=d, rows=len(rows_mons), cols=len(cols), n_targets=len(tcols), rank_M=rM, rank_MB=rMB, dim_served_Q=dim_served,
        result='served(D) = J_6^odd (trivial only)' if dim_served == 18 else f'dim served = {dim_served} (J_6^odd has 18)', seconds=round(time.time() - t0, 1))

if __name__ == '__main__':
    cmd = sys.argv[1]; args = sys.argv[2:]
    if cmd == 'hilbert': cmd_hilbert()
    elif cmd == 'e1': cmd_e1()
    elif cmd == 'd5exact': cmd_d5exact()
    elif cmd == 'd8star': cmd_d8star()
    elif cmd == 'd6': cmd_d6()
    elif cmd == 'd7pair': cmd_d7pair(args[0], args[1], *(int(x) for x in args[2:]))
    elif cmd == 'd8sector':
        for tag in args[0].split(','): cmd_d8sector(tag, *(int(x) for x in args[1:]))
    elif cmd == 'd8inv': cmd_d8inv(*(int(x) for x in args))
    elif cmd == 'fixedD-exact': cmd_fixedD_exact(*args)
    elif cmd == 'served-exact': cmd_served_exact(args)
    elif cmd == 'served-exact-full': cmd_served_exact_full(args[0])
    else: raise SystemExit('unknown command')
