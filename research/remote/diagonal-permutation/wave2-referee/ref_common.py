"""Referee replay library for lane diagonal-permutation (independent of the `hodge` engine and of w2lib).
Polynomials: dict {exponent-tuple (6): coefficient}.  Exact coefficients are Fractions; modular ones are ints.
Variables X0,X1,X2,Y0,Y1,Y2 (indices 0..5); Y-parity of a monomial = (e3+e4+e5) mod 2.
"""
import re, itertools, json, time
from fractions import Fraction
import sympy

NV = 6
VARS = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
REC = '/home/kejace/devel/sextics/data/diagonal-permutation/instance.rec'


def load_instance():
    txt = open(REC).read()
    lines = txt.splitlines()
    fields = {}
    key = None
    for ln in lines:
        if ln.startswith('#'):
            continue
        m = re.match(r'^([A-Za-z]+):\s*(.*)$', ln)
        if m and not ln.startswith(' '):
            key = m.group(1); fields[key] = m.group(2)
        elif key is not None:
            fields[key] += ' ' + ln.strip()
    F = sympy.Poly(sympy.sympify(fields['F']), *VARS)
    A = sympy.Poly(sympy.sympify(fields['A']), *VARS)
    q = int(fields['q'])
    def todict(P):
        return {tuple(int(e) for e in mon): Fraction(int(c.p), int(c.q)) for mon, c in zip(P.monoms(), P.coeffs())}
    return todict(F), todict(A), q


def monomials(n, par=None):
    out = []
    for e in itertools.product(range(n + 1), repeat=NV):
        if sum(e) == n and (par is None or (e[3] + e[4] + e[5]) % 2 == par):
            out.append(e)
    return out


def pderiv(P, i):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            out[tuple(e2)] = out.get(tuple(e2), 0) + c * e[i]
    return {e: c for e, c in out.items() if c != 0}


def mul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            out[e] = out.get(e, 0) + c1 * c2
    return {e: c for e, c in out.items() if c != 0}


def add(P, Q, s=1):
    out = dict(P)
    for e, c in Q.items():
        out[e] = out.get(e, 0) + s * c
    return {e: c for e, c in out.items() if c != 0}


def scale(P, c):
    return {e: c * v for e, v in P.items() if c * v != 0}


def degree(P):
    return max(sum(e) for e in P) if P else -1


def modp(c, p):
    c = Fraction(c)
    return (c.numerator % p) * pow(c.denominator % p, -1, p) % p


def koszul_series(nmax=30):
    """dim R_n^{+}, dim R_n^{-} for R = S/J(F), F smooth sextic in 6 vars (Koszul Hilbert function),
    refined by the Y-sign involution tau = diag(1,1,1,-1,-1,-1):  chi_g(t) = det(1 - t^5 g)/det(1 - t g)."""
    t = sympy.symbols('t')
    chi1 = ((1 - t**5) / (1 - t))**6
    chitau = ((1 - t**5) / (1 - t))**3 * ((1 + t**5) / (1 + t))**3
    s1 = sympy.series(chi1, t, 0, nmax + 1).removeO()
    s2 = sympy.series(chitau, t, 0, nmax + 1).removeO()
    plus = [int((sympy.Poly(s1, t).coeff_monomial(t**n) + sympy.Poly(s2, t).coeff_monomial(t**n)) / 2) for n in range(nmax + 1)]
    minus = [int((sympy.Poly(s1, t).coeff_monomial(t**n) - sympy.Poly(s2, t).coeff_monomial(t**n)) / 2) for n in range(nmax + 1)]
    return plus, minus


def sector_series(nmax=30):
    """dim R_n^{par,sigma} for the group <tau> x <(01)>, (01) = simultaneous swap of (X0,X1) and (Y0,Y1).
    Eigenvalues on variables: tau (1,1,1,-1,-1,-1); (01): (1,-1,1,1,-1,1); (01)tau: (1,-1,1,-1,1,-1)."""
    t = sympy.symbols('t')
    def chi(eigs):
        num = 1; den = 1
        for lam in eigs:
            num *= (1 - lam * t**5); den *= (1 - lam * t)
        return sympy.series(num / den, t, 0, nmax + 1).removeO()
    chis = {'1': chi([1]*6), 'tau': chi([1, 1, 1, -1, -1, -1]), 's': chi([1, -1, 1, 1, -1, 1]), 'stau': chi([1, -1, 1, -1, 1, -1])}
    co = {k: sympy.Poly(v, t) for k, v in chis.items()}
    def c(k, n):
        return int(co[k].coeff_monomial(t**n))
    table = {}
    for n in range(nmax + 1):
        for par in (0, 1):          # 0 = Y-even (tau = +1), 1 = Y-odd (tau = -1)
            for sig in (1, -1):     # (01)-eigenvalue
                tp = 1 if par == 0 else -1
                val = (c('1', n) + tp * c('tau', n) + sig * c('s', n) + tp * sig * c('stau', n))
                assert val % 4 == 0
                table[(n, par, sig)] = val // 4
    return table


def nmod_rank_from_columns(cols, row_index, p):
    """cols: list of polynomials (dicts) in the same degree; rows indexed by row_index (dict monomial -> row).
    Returns (rank, nrows, ncols) using python-flint nmod_mat (dense, built from nonzeros)."""
    from flint import nmod_mat
    M = nmod_mat(len(row_index), len(cols), p)
    for j, P in enumerate(cols):
        for e, c in P.items():
            M[row_index[e], j] = modp(c, p) if not isinstance(c, int) else c % p
    return M.rank(), len(row_index), len(cols)
