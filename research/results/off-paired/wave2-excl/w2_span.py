"""Fallback rigor for non-m-primary relaxation ideals: certify the linear span of V(J2) over Q.

Given exact rational quadrics quads (dict (i,j)->Fraction, i<=j) in nt variables:
 1. Singular mod p: minimal associated primes of J2 mod p; the linear forms common to all components (degree-1 generators present in
    every component, intersected as subspaces) define the span L_p of V(J2^p).
 2. Reconstruct the linear forms rationally (small heights expected), verify they reduce to the modular ones.
 3. For each form ell: exact rational cofactor solve  ell^k = sum_r c_r(t) g_r  (Macaulay system in degree k, fmpq_mat), k from the
    modular normal form.  Success => ell in rad(J2^Q) => V(J2^Q) subset of L_Q = V(ell_1..ell_c).
Returns (basis of L_Q as Fraction vectors, certificate info) or (None, info)."""
import sys, re, subprocess, itertools, time
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
from fractions import Fraction
from flint import fmpq_mat, fmpq, nmod_mat
from hodge import core

def monomials_n(n, deg):
    out = []
    for comb in itertools.combinations(range(deg + n - 1), n - 1):
        e = []; prev = -1
        for x in comb: e.append(x - prev - 1); prev = x
        e.append(deg + n - 2 - prev)
        out.append(tuple(e))
    return out

def frac_str(v):
    return f'({v.numerator}/{v.denominator})' if v.denominator != 1 else f'({v.numerator})'

def quads_modp_str(quads, p):
    out = []
    for Qd in quads:
        terms = []
        for (i, j), v in Qd.items():
            c = (v.numerator * pow(v.denominator, -1, p)) % p
            if c: terms.append(f'{c}*t{i}*t{j}' if i != j else f'{c}*t{i}^2')
        if terms: out.append('+'.join(terms))
    return out

def parse_linear(P, nt, p=None):
    """Parse a Singular polynomial string; return coefficient vector if it is a linear form, else None."""
    s = P.replace(' ', '')
    row = [0] * nt
    for term in re.findall(r'[+-]?[^+-]+', s):
        if not term: continue
        m = re.match(r'([+-]?)(\d+(?:/\d+)?)?\*?t(\d+)$', term)
        if not m: return None
        sign = -1 if m.group(1) == '-' else 1
        coef = Fraction(m.group(2)) if m.group(2) else Fraction(1)
        row[int(m.group(3))] += sign * coef
    return row

def modp_components_linear_forms(quads, nt, p, workdir, timeout=600):
    """Singular mod p: minAssGTZ(J2); returns list of components, each a list of linear-form rows (mod p ints) plus flag nonlinear-gens."""
    polys = quads_modp_str(quads, p)
    script = f'''LIB "primdec.lib";
ring r = {p}, ({",".join(f"t{i}" for i in range(nt))}), dp;
ideal J2 = {",\n".join(polys)};
list L = minAssGTZ(J2);
int i; int j;
for (i = 1; i <= size(L); i++) {{
  string s = "COMP " + string(i) + ":";
  for (j = 1; j <= size(L[i]); j++) {{ s = s + " " + string(L[i][j]) + " ;"; }}
  print(s);
}}
quit;
'''
    path = f'{workdir}/J2_modp_minass.sing'; open(path, 'w').write(script)
    try:
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    comps = []
    for line in r.stdout.splitlines():
        if line.startswith('COMP '):
            gens = [g.strip() for g in line.split(':', 1)[1].split(';') if g.strip()]
            lin = []; nonlin = 0
            for g in gens:
                row = parse_linear(g, nt)
                if row is None: nonlin += 1
                else: lin.append([int(x) % p if x.denominator == 1 else (x.numerator * pow(x.denominator, -1, p)) % p for x in row])
            comps.append({'linear': lin, 'nonlinear': nonlin})
    return comps

def span_forms_modp(comps, nt, p):
    """Linear forms vanishing on the union of the components = intersection of the row spaces of the components' linear parts."""
    # intersection of subspaces U_i (row spaces): forms vanishing on V(P_i) for all i.  Compute as kernel of the stacked annihilators.
    # annihilator of U_i: nullspace of the matrix of U_i's rows; forms in all U_i = forms orthogonal to all annihilators.
    anns = []
    for comp in comps:
        rows = comp['linear']
        if not rows:
            return []   # a component with no linear forms spans everything
        M = nmod_mat(len(rows), nt, [x for r in rows for x in r], p)
        X, nul = M.nullspace()
        for k in range(nul):
            anns.append([int(X[j, k]) for j in range(nt)])
    if not anns:
        return [[1 if i == j else 0 for i in range(nt)] for j in range(nt)]
    Mann = nmod_mat(len(anns), nt, [x for r in anns for x in r], p)
    X, nul = Mann.nullspace()
    forms = []
    T = nmod_mat(nul, nt, [int(X[j, k]) for k in range(nul) for j in range(nt)], p)
    Tr, rk = T.rref()
    for i in range(rk):
        forms.append([int(Tr[i, j]) for j in range(nt)])
    return forms

def reconstruct_forms(forms, p):
    out = []
    for f in forms:
        row = []
        for x in f:
            r = core.ratrecon(x, p) if x else Fraction(0)
            if r is None: return None
            row.append(r)
        out.append(row)
    return out

def nilpotent_certificate(quads, nt, ell, k):
    """Exact: is ell^k in the ideal generated by quads (in degree k)?  fmpq_mat consistency of the Macaulay system."""
    cmons = monomials_n(nt, k - 2); kmons = monomials_n(nt, k); kidx = {m: i for i, m in enumerate(kmons)}
    cols = []
    for Qd in quads:
        for dm in cmons:
            col = {}
            for (i, j), v in Qd.items():
                m = list(dm); m[i] += 1; m[j] += 1; m = tuple(m)
                col[kidx[m]] = col.get(kidx[m], 0) + v
            cols.append(col)
    # ell^k expanded
    powk = {tuple(0 for _ in range(nt)): Fraction(1)}
    for _ in range(k):
        new = {}
        for m, v in powk.items():
            for i, c in enumerate(ell):
                if c == 0: continue
                mm = list(m); mm[i] += 1; mm = tuple(mm)
                new[mm] = new.get(mm, 0) + v * c
        powk = {m: v for m, v in new.items() if v}
    R, C = len(kmons), len(cols)
    M = fmpq_mat(R, C + 1)
    for j, col in enumerate(cols):
        for r, v in col.items(): M[r, j] = fmpq(v.numerator, v.denominator)
    for m, v in powk.items(): M[kidx[m], C] = fmpq(v.numerator, v.denominator)
    M0 = fmpq_mat(R, C)
    for j, col in enumerate(cols):
        for r, v in col.items(): M0[r, j] = fmpq(v.numerator, v.denominator)
    return M.rank() == M0.rank(), {'rows': R, 'cols': C}

def certify_span(quads, nt, p, workdir, kmax=5, timeout=600, verbose=True):
    t0 = time.time()
    comps = modp_components_linear_forms(quads, nt, p, workdir, timeout)
    if comps is None: return None, {'fail': 'singular-modp-timeout'}
    forms_p = span_forms_modp(comps, nt, p)
    forms = reconstruct_forms(forms_p, p)
    if forms is None: return None, {'fail': 'form reconstruction', 'ncomps': len(comps)}
    certs = []
    for ell in forms:
        done = False
        for k in range(2, kmax + 1):
            ok, info = nilpotent_certificate(quads, nt, ell, k)
            if ok:
                certs.append({'k': k, **info}); done = True; break
        if not done:
            return None, {'fail': 'nilpotency certificate failed', 'form': [str(x) for x in ell], 'ncomps': len(comps)}
    # basis of L_Q = common kernel of the forms
    M = fmpq_mat(len(forms), nt)
    for i, row in enumerate(forms):
        for j, v in enumerate(row): M[i, j] = fmpq(v.numerator, v.denominator)
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j < nt and Mr[i, j] == 0: j += 1
        piv.append(j)
    free = [j for j in range(nt) if j not in piv]
    basis = []
    for f in free:
        v = [Fraction(0)] * nt; v[f] = Fraction(1)
        for i, j in enumerate(piv): v[j] = -Fraction(int(Mr[i, f].p), int(Mr[i, f].q))
        basis.append(v)
    info = {'ncomps': len(comps), 'comps_nonlinear_gens': [c['nonlinear'] for c in comps], 'n_forms': len(forms), 'certs': certs,
            'span_dim': len(basis), 'seconds': round(time.time() - t0, 1)}
    if verbose: print('[span]', info, file=sys.stderr, flush=True)
    return basis, info
