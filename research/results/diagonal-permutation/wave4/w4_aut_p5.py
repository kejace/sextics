"""Monomial linear automorphisms of F (exact, over Q(zeta_6) and over F_p-bar).  GIVEN that the rank-1 Hessian locus of F
is exactly the six coordinate points, every linear automorphism permutes them, hence is monomial: z_i -> lambda_i z_{pi(i)}.
F o g = F forces lambda_i^6 = 1 (Fermat part) and c_m lambda^m = c_{pi*m} for every monomial m; since the coefficients are
rational, lambda^m in {+1,-1}.  We enumerate all (pi, lambda) exactly as a Z/6 system.  Output: the group, its order,
the list of elements with eigenvalue pattern (1^5,-1) (reflections) -- expected none -- and a check of the known S3 x <tau>.
Usage: python3 w4_aut_p5.py [prime]   (prime only used to double-check the coefficient-ratio argument mod p)"""
import sys, itertools, json
import sympy as sp
import w4lib as L
p = int(sys.argv[1]) if len(sys.argv)>1 else L.P1
P = sp.Poly(L.Fs, *L.Z6)
terms = dict(P.terms())           # monomial exponent tuple -> coefficient (Rational)
supp = list(terms.keys())
n = 6
elements = []
for pi in itertools.permutations(range(n)):
    # g: z_i -> lam_i z_{pi(i)}; monomial m = prod z_i^{m_i} -> lam^m prod z_{pi(i)}^{m_i}: new exponent e with e[pi(i)] = m_i
    def img(mm):
        e = [0]*n
        for i in range(n): e[pi[i]] = mm[i]
        return tuple(e)
    ok = True; signs = {}
    for mm in supp:
        e = img(mm)
        if e not in terms: ok = False; break
        ratio = terms[e]/terms[mm]
        if ratio == 1: signs[mm] = 0
        elif ratio == -1: signs[mm] = 3
        else: ok = False; break   # ratio must be a 6th root of unity in Q, i.e. +-1
    if not ok: continue
    # solve lam^m = zeta^{signs} : with lam_i = zeta^{k_i}, sum_i m_i k_i = signs[m] mod 6, for all m in supp
    sols = []
    for ks in itertools.product(range(6), repeat=n):
        if all(sum(mm[i]*ks[i] for i in range(n)) % 6 == signs[mm] for mm in supp):
            sols.append(ks)
    for ks in sols:
        elements.append({'pi': pi, 'k': ks})
# projective classes: (pi, k) ~ (pi, k + c(1,...,1))
classes = {}
for el in elements:
    k = el['k']; c = k[0]; kn = tuple((ki - c) % 6 for ki in k)
    classes[(el['pi'], kn)] = el
order = len(classes)
# eigenvalue patterns: matrix g with g[pi(i), i] = zeta^{k_i} (column i has lambda_i at row pi(i))
zeta = sp.exp(2*sp.pi*sp.I/6)
refl = []; summary = []
for (pi, kn) in classes:
    g = sp.zeros(n, n)
    for i in range(n): g[pi[i], i] = sp.nsimplify(zeta**kn[i])
    # projective eigenvalue pattern: normalise so that g^ord = 1; compute char poly roots numerically
    ev = [complex(sp.N(x)) for x in g.eigenvals(multiple=True)]
    # pattern as multiset of arguments (in units of 1/6 turn) relative to the most frequent eigenvalue
    import cmath
    args = [round((cmath.phase(x)/(2*cmath.pi))*6) % 6 for x in ev]
    from collections import Counter
    cnt = Counter(args)
    base = cnt.most_common(1)[0][0]
    rel = sorted(((a - base) % 6) for a in args)
    summary.append({'pi': pi, 'k': kn, 'rel_eigen_sixths': rel})
    if rel == [0,0,0,0,0,3]: refl.append((pi, kn))
res = {'order_projective': order, 'n_linear_solutions': len(elements), 'reflections': refl, 'elements': summary}
print(json.dumps(res, default=str), flush=True)
L.log({'what':'monomial automorphism group of F (exact over Q(zeta6), given rank-1 locus = coordinate points)','prime':None,
       'result':{'order_projective':order, 'reflections':refl, 'elements':summary}})
# modular double check of the coefficient-ratio argument: the ratios of distinct |coefficients| are not 6th roots of unity mod p
vals = sorted(set(abs(c) for c in terms.values()))
bad = []
for c1 in vals:
    for c2 in vals:
        if c1 != c2:
            r = (sp.Rational(c1)/sp.Rational(c2)); rp = int(r.p) * pow(int(r.q), -1, p) % p
            if pow(rp, 6, p) == 1: bad.append((str(c1), str(c2)))
print('coefficient magnitudes', [str(v) for v in vals], 'ratios that are 6th roots of unity mod p:', bad)
L.log({'what':'mod-p check: coefficient ratios that are 6th roots of unity (would allow extra monomial automorphisms over F_p-bar)','prime':p,'result':bad})
