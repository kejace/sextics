"""Candidate denominators (5-variable, Y2-free, rho-even) organised by family."""
import sys, itertools, pickle, random
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from fractions import Fraction
random.seed(11)
P = core.parse
inv = {  # skew-adapted invariants (name -> (poly, degree))
 'G0': P('X0**6+X1**6+X2**6'), 'c': P('X0**3*X1+X1**3*X2'), 'e': P('X0*X1'), 'p': P('Y0*Y1'),
 't': P('Y0**2-Y1**2'), 'm': P('Y0**2+Y1**2'), 's': P('Y0**3-Y1**3'), 'S': P('Y0**3+Y1**3'),
 'Y6': P('Y0**6+Y1**6'), 'f': f, 'f0': sub(f, P('X0*X1*Y0*Y1*(Y0**2-Y1**2)/7')),
 'skew': P('X0*X1*Y0*Y1*(Y0**2-Y1**2)/7'), 'cp': P('(X0**3*X1+X1**3*X2)*Y0*Y1'),
 'pure': P('X0**6+X1**6+X2**6+Y0**6+Y1**6'), 'c1': P('X0**3*X1*Y0*Y1'), 'c2': P('X1**3*X2*Y0*Y1'),
 'X0': X[0], 'X1': X[1], 'X2': X[2], 'Y0': X[3], 'Y1': X[4],
}
def cands():
    out = []
    # (3) structural pieces of f and their simple products
    pieces = ['G0', 'Y6', 'cp', 'skew', 'f0', 'pure', 'c1', 'c2']
    for nm in pieces:
        out.append((nm, inv[nm], 'struct:piece'))
    out.append(('f-G0', sub(f, inv['G0']), 'struct:piece'))
    out.append(('f-Y6', sub(f, inv['Y6']), 'struct:piece'))
    out.append(('f-c1', sub(f, inv['c1']), 'struct:piece'))
    out.append(('f-c2', sub(f, inv['c2']), 'struct:piece'))
    # partials, products of partials, radial block derivatives, Hessian minors
    for i in range(5):
        out.append((f'f_{core.VARS[i]}', fi(i), 'struct:partial'))
    for i in range(5):
        for j in range(i, 5):
            out.append((f'f_{core.VARS[i]}*f_{core.VARS[j]}', mul(fi(i), fi(j)), 'struct:partial-product'))
    for i in range(3):
        out.append((f'radial{i}', add(mul(X[i], fi(i)), mul(X[3 + min(i, 1)], fi(3 + min(i, 1)))), 'struct:radial'))
    H = [[pderiv(fi(i), j) for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(i + 1, 5):
            out.append((f'Hess[{i},{j}]', sub(mul(H[i][i], H[j][j]), mul(H[i][j], H[i][j])), 'struct:hessian-2x2'))
    # (2)+(E6) products of invariants (degrees 7..10, rho-even), incl. variables
    names = ['G0', 'c', 'e', 'p', 't', 'm', 's', 'S', 'Y6', 'f', 'X0', 'X1', 'X2', 'Y0', 'Y1', 'cp', 'skew']
    seen = set()
    for k in (2, 3):
        for combo in itertools.combinations_with_replacement(names, k):
            Pp = {(0,)*6: Fraction(1)}
            for nm in combo: Pp = mul(Pp, inv[nm])
            d = core.degree(Pp)
            if d < 7 or d > 10 or rparity_of(Pp) != 0: continue
            key = tuple(sorted(Pp.items()))
            if key in seen: continue
            seen.add(key)
            out.append(('*'.join(combo), Pp, f'invariant-product:k{k}'))
    # sums of invariants (Fermat-type combos) at degree 8: G0*q + Y6*q', m*c*e ... a few hand-picked
    hand = [
      ('G0*p+Y6*e', add(mul(inv['G0'], inv['p']), mul(inv['Y6'], inv['e'])), 'invariant-sum'),
      ('G0*m-Y6*e', sub(mul(inv['G0'], inv['m']), mul(inv['Y6'], inv['e'])), 'invariant-sum'),
      ('c*c+G0*p', add(mul(inv['c'], inv['c']), mul(inv['G0'], inv['p'])), 'invariant-sum'),
      ('c*c-3*G0*p', sub(mul(inv['c'], inv['c']), scale(mul(inv['G0'], inv['p']), 3)), 'invariant-sum'),
      ('c*c*t+G0*e*p', add(mul(mul(inv['c'], inv['c']), inv['t']), mul(mul(inv['G0'], inv['e']), inv['p'])), 'invariant-sum'),
      ('f*p+c*e*t', add(mul(f, inv['p']), mul(mul(inv['c'], inv['e']), inv['t'])), 'invariant-sum'),
      ('f*t+c*e*p', add(mul(f, inv['t']), mul(mul(inv['c'], inv['e']), inv['p'])), 'invariant-sum'),
      ('G0*t+c*e*p', add(mul(inv['G0'], inv['t']), mul(mul(inv['c'], inv['e']), inv['p'])), 'invariant-sum'),
      ('s*S*e+G0*p', add(mul(mul(inv['s'], inv['S']), inv['e']), mul(inv['G0'], inv['p'])), 'invariant-sum'),
      ('c*s*Y0+G0*p', add(mul(mul(inv['c'], inv['s']), X[3]), mul(inv['G0'], inv['p'])), 'invariant-sum'),
      # binary-cubic Hessian-type covariants of P(U) = G0 U^3 + cp U^2 + (e p t/7) U + Y6, "homogenised" by X0^2 <-> U
      ('hessP:b2-3ac', sub(mul(inv['cp'], inv['cp']), scale(mul(inv['G0'], inv['skew']), 3)), 'binary-cubic-covariant'),
      ('hessP:bc-9ad', sub(mul(inv['cp'], inv['skew']), scale(mul(inv['G0'], inv['Y6']), 9)), 'binary-cubic-covariant'),
      ('hessP:c2-3bd', sub(mul(inv['skew'], inv['skew']), scale(mul(inv['cp'], inv['Y6']), 3)), 'binary-cubic-covariant'),
    ]
    out += hand
    # cones: 4-variable and 3-variable sparse forms of degree 7, 8
    cone = []
    for drop in range(5):
        keep = [i for i in range(5) if i != drop]
        for d in (7, 8):
            # Fermat-type sum over kept variables with even exponent on Y's as needed
            terms = {}
            for i in keep:
                e = [0]*6; e[i] = d
                if rpar(tuple(e)) == 0: terms[tuple(e)] = Fraction(1)
            if terms: cone.append((f'cone-drop{core.VARS[drop]}-fermat{d}', terms, 'cone:4var-fermat'))
            # f restricted to kept vars times a kept variable / quadric
            fk = {e: c for e, c in f.items() if e[drop] == 0}
            for i in keep:
                Pp = mul(fk, X[i]) if d == 7 else mul(fk, mul(X[i], X[i]))
                if rparity_of(Pp) == 0: cone.append((f'cone-drop{core.VARS[drop]}-f|*{core.VARS[i]}^{d-6}', Pp, 'cone:4var-f-restricted'))
    for keep in itertools.combinations(range(5), 3):
        for d in (7, 8):
            terms = {}
            for i in keep:
                e = [0]*6; e[i] = d
                if rpar(tuple(e)) == 0: terms[tuple(e)] = Fraction(1)
            # add a mixed term to avoid pure Fermat degeneracy
            if len(terms) >= 2:
                cone.append((f'cone3-{"".join(core.VARS[i] for i in keep)}-fermat{d}', terms, 'cone:3var-fermat'))
    out += cone
    # random dense D of degree 7 and 8 (generic calibration), 2 each
    for d in (7, 8):
        for k in range(2):
            Pp = {m: Fraction(random.randint(-3, 3)) for m in mons5(d, 0)}
            Pp = {e: c for e, c in Pp.items() if c}
            out.append((f'random-dense-d{d}-#{k}', Pp, 'random-generic'))
    # filter: rho-even, nonzero, degree range
    final = []
    for name, Pp, fam in out:
        if not Pp or rparity_of(Pp) != 0: continue
        final.append((name, Pp, fam))
    return final
if __name__ == '__main__':
    C = cands()
    from collections import Counter
    print(len(C), Counter((fam, core.degree(Pp)) for _, Pp, fam in C))
