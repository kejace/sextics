#!/usr/bin/env python3
"""qcr: the reflection towers X -> X/s -> P^4 for the real reflections s in Gamma^- (anti-invariant linear form v):
F is even in v: F = g(v^2) with g a cubic over C(w, other coords).  Compute, at random rational specialisations of the
five invariant coordinates: irreducibility of g, its Galois group, whether disc(g) is a square (A3 vs S3), whether
-g(0) is a square (product z1 z2 z3 rational => Galois group of degree-6 poly drops from B3 (48) to S4 (24)),
and the Galois group of the degree-6 polynomial in v (sympy).  Also the three maximal 'reaching' quotients' adjunction.
"""
import sys, json, random, re
import sympy as sp
from sympy.polys.numberfields.galoisgroups import galois_group
X = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
txt = open('/home/kejace/devel/sextics/data/quadric-contact-residual/instance.rec').read()
Fs = re.search(r'^F: (.*(?:\n  .*)*)', txt, re.M).group(1).replace('\n',' ')
F = sp.expand(7*sp.sympify(Fs))
v, w = sp.symbols('v w')
towers = {
 's0 (X0<->Y0), v=X0-Y0':   {X[0]: (w+v)/2, X[3]: (w-v)/2},
 's0*sigma0, v=X0+Y0':      {X[0]: (w+v)/2, X[3]: (v-w)/2},
 'r1 (X1<->Y1), v=X1-Y1':   {X[1]: (w+v)/2, X[4]: (w-v)/2},
 'r1*h1^3, v=X1+Y1':        {X[1]: (w+v)/2, X[4]: (v-w)/2},
}
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
out = {}
for name, sub in towers.items():
    G = sp.expand(F.subs(sub, simultaneous=True))
    Pv = sp.Poly(G, v)
    assert all(m[0] % 2 == 0 for m in Pv.monoms()), 'not even in v'
    rec = {'even_in_v': True, 'deg_v': Pv.degree(), 'specialisations': []}
    others = [s for s in X if s not in sub] + [w]
    for trial in range(3):
        vals = {s: sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for s in others}
        g6 = sp.Poly(G.subs(vals), v)
        u = sp.symbols('u')
        g3 = sp.Poly(sum(c*u**(m[0]//2) for m, c in zip(g6.monoms(), g6.coeffs())), u)
        g3m = sp.Poly(g3.as_expr()/g3.LC(), u)
        irr3 = g3m.is_irreducible
        disc3 = sp.discriminant(g3m.as_expr(), u)
        c0 = g3m.eval(0)
        g6m = sp.Poly(g6.as_expr()/g6.LC(), v)
        try:
            gg = galois_group(g6m, by_name=True)
            gname = str(gg[0]) if isinstance(gg, tuple) else str(gg)
        except Exception as ex:
            gname = 'ERR '+str(ex)[:80]
        try:
            gg3 = galois_group(g3m, by_name=True); g3name = str(gg3[0]) if isinstance(gg3, tuple) else str(gg3)
        except Exception as ex:
            g3name = 'ERR '+str(ex)[:60]
        rec['specialisations'].append({'cubic_irreducible': bool(irr3), 'cubic_galois': g3name,
            'disc_cubic_is_square': bool(sp.sqrt(disc3).is_rational), 'minus_c0_is_square': bool(sp.sqrt(-c0).is_rational),
            'sextic_irreducible': bool(g6m.is_irreducible), 'sextic_galois': gname})
    out[name] = rec
    print(name, json.dumps(rec), flush=True)
json.dump(out, open('/tmp/hodge/results/_cross/trio-qcr-galois.json','w'), indent=1)
