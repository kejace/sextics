import sympy as sp
x,y,z,t,r,s = sp.symbols('x y z t r s')
U,V,W,T = sp.symbols('U V W T')
F = x**5*y + y**5*z + z**5*x + sp.Rational(5,7)*(x*y**3*t*r + 2*x**3*z*r*s + 3*y*z**3*t*s) + t**6 + r**6 + s**6
# (1) multiplicities of tangent hyperplane sections
# at e_x: tangent hyperplane y=0, chart x=1
f_ex = F.subs({x:1, y:0})
print('e_x tangent section (x=1,y=0):', sp.expand(f_ex))
# lowest total degree in (z,t,r,s)
P = sp.Poly(sp.expand(f_ex), z,t,r,s); print('  min degree =', min(sum(m) for m in P.monoms()))
# at Fermat point (0,0,0,1,r0,0) with r0^6=-1: tangent hyperplane r = r0 (chart t=1)
r0 = sp.Symbol('r0')
f_fp = sp.expand(F.subs({t:1, r:r0}).subs(r0**6, -1))
f_fp = sp.expand(f_fp.subs(r0**6, -1))
P2 = sp.Poly(f_fp, x,y,z,s); print('Fermat point tangent section: min degree in (x,y,z,s) =', min(sum(m) for m in P2.monoms()), ' constant term:', P2.coeff_monomial(1))
# (2) Kummer exponent: x^14 = U^5 W^3 / V with U=x^2y/z, V=y^2z/x, W=z^2x/y
Uv, Vv, Wv = x**2*y/z, y**2*z/x, z**2*x/y
print('U^5 W^3 / V =', sp.simplify(Uv**5*Wv**3/Vv))
# (3) Geiser involution from e_U on the weighted W: P = U^2W+V^2U+W^2V + (5/7)(UV tr + 2UW rs + 3VW ts) + t^6+r^6+s^6
Pw = U**2*W + V**2*U + W**2*V + sp.Rational(5,7)*(U*V*t*r + 2*U*W*r*s + 3*V*W*t*s) + t**6 + r**6 + s**6
grad = [sp.diff(Pw, v).subs({U:1,V:0,W:0,t:0,r:0,s:0}) for v in (U,V,W,t,r,s)]
print('grad P at e_U:', grad)
Ug = -U - (V**2 + sp.Rational(5,7)*V*t*r + sp.Rational(10,7)*W*r*s)/W
print('P(gamma) - P == 0 ?', sp.simplify(Pw.subs(U, Ug) - Pw) == 0)
Ugg = Ug.subs(U, Ug)
print('gamma^2 = id ?', sp.simplify(Ugg - U) == 0)
# (4) pullback of gamma U to X equals phi/(x^3 z) modulo F, phi = z^5 x + (15/7) y z^3 t s + t^6+r^6+s^6
phi = z**5*x + sp.Rational(15,7)*y*z**3*t*s + t**6 + r**6 + s**6
pull = Ug.subs({U:Uv, V:Vv, W:Wv})
diff = sp.simplify(pull - phi/(x**3*z))
print('pi^*(gamma U) - phi/(x^3 z) =', sp.factor(diff), '  (should be -F/(x^3 z))')
# (5) Kummer generator of the twisted cover: (gamma U)^5 W^3 / V pulled back
gen = sp.simplify((phi/(x**3*z))**5 * Wv**3 / Vv)
print('pulled-back gamma^* h =', gen)
