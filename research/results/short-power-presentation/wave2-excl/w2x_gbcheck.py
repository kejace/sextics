"""Zero-dimensionality check of a homogeneous ideal from a msolve Groebner basis (-g 2 output, grevlex).
A reduced GB mod p is an exact object; the ideal I' is zero-dimensional (affine) iff for every variable some
leading monomial is a pure power.  Then I' contains every monomial of degree >= sum_i (a_i - 1) + 1, so the
Hilbert function of S'/I'_p vanishes there, and by rank semicontinuity so does that of S'/I'_Q: V_C(I') = {0}.
Usage: python3 w2x_gbcheck.py <msolve input .ms> <msolve output>  -> prints JSON"""
import sys, re, json

def parse_vars(ms_path):
    return open(ms_path).readline().strip().split(',')

def grevlex_key(exp):
    # larger = bigger in grevlex: compare total degree, then the LAST nonzero exponent difference reversed
    return (sum(exp), tuple(-e for e in reversed(exp)))

def leading_monomials(out_path, names):
    txt = open(out_path).read()
    # msolve -g 2 prints: #...\n[poly1,\npoly2,\n...]:  strip brackets
    body = txt[txt.find('['):txt.rfind(']')]
    body = body.strip('[]')
    polys = [s.strip() for s in body.split(',\n') if s.strip()]
    idx = {n: i for i, n in enumerate(names)}
    lms = []
    for poly in polys:
        best = None
        for term in re.split(r'(?=[+-])', poly.replace(' ', '')):
            if not term or term in '+-': continue
            exp = [0] * len(names)
            for fac in term.split('*'):
                m = re.fullmatch(r'[+-]?(\d*)?(?:([A-Za-z_]\w*)(?:\^(\d+))?)?', fac)
                if m and m.group(2):
                    exp[idx[m.group(2)]] += int(m.group(3) or 1)
            key = grevlex_key(exp)
            if best is None or key > best[0]: best = (key, tuple(exp))
        if best: lms.append(best[1])
    return polys, lms

if __name__ == '__main__':
    ms, out = sys.argv[1], sys.argv[2]
    names = parse_vars(ms)
    polys, lms = leading_monomials(out, names)
    n = len(names)
    pure = {}
    for lm in lms:
        nz = [i for i, e in enumerate(lm) if e]
        if len(nz) == 1:
            i = nz[0]; pure[i] = min(pure.get(i, 10 ** 9), lm[i])
    zero_dim = all(i in pure for i in range(n))
    bound = sum(a - 1 for a in pure.values()) + 1 if zero_dim else None
    lin = sum(1 for lm in lms if sum(lm) == 1)
    res = {'gb_size': len(polys), 'n_vars': n, 'zero_dimensional': zero_dim, 'pure_power_vars': len(pure),
           'missing_vars': [names[i] for i in range(n) if i not in pure][:20], 'linear_leading_terms': lin,
           'max_lead_degree': max((sum(lm) for lm in lms), default=0), 'degree_bound_all_monomials_in_I': bound}
    print(json.dumps(res))
