#!/usr/bin/env python3
"""Explicit triple-cover elliptic surface S3 over the t = v^2 line (book §2.9 / R14 normalisation):
cubic U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f z^2 x with (a,b,c) = (1+v^6, 1+c1^6 v^6, 1+c2^6 v^6),
(d,e,f) = v^2/2 (c1, c1 c2, c2); s4 = a e^2 f + b d f^2 + c d^2 e, s6 = 27a^2b^2c^2 + 4a^2 c e^3 + 4 a b^2 f^3
+ 6 a b c d e f + 4 b c^2 d^3 - d^2 e^2 f^2; Jacobian eta^2 = xi^3 + s4 xi - s6/4.  Substituting v^2 = t:
s4, s6 are polynomials in t; the RES over T = t^3 needs them to be polynomials in t^3 -- check which
monomials in t occur, the degrees, chi = ceil(max(deg s4/4, deg s6/6)), p_g = chi - 1, the 2-forms
t^j dt dxi/eta (j <= chi-2) and their characters under t -> omega t."""
import sympy as sp
t,c1,c2,xi=sp.symbols('t c1 c2 xi')
a,b,c=1+t**3,1+c1**6*t**3,1+c2**6*t**3
d,e,f=t*c1/2,t*c1*c2/2,t*c2/2
s4=sp.expand(a*e**2*f+b*d*f**2+c*d**2*e)
s6=sp.expand(27*a**2*b**2*c**2+4*a**2*c*e**3+4*a*b**2*f**3+6*a*b*c*d*e*f+4*b*c**2*d**3-d**2*e**2*f**2)
P4=sp.Poly(s4,t); P6=sp.Poly(s6,t)
print("deg_t s4 =",P4.degree()," t-exponents:",sorted(m[0] for m in P4.monoms()))
print("deg_t s6 =",P6.degree()," t-exponents:",sorted(m[0] for m in P6.monoms()))
print("s4 =",sp.factor(s4))
print("t-exponents not divisible by 3 in s6:",sorted(m[0] for m in P6.monoms() if m[0]%3))
disc=sp.expand(4*s4**3+27*(s6/4)**2)   # for y^2 = x^3 + s4 x - s6/4 : 4 s4^3 + 27 (s6/4)^2
Pd=sp.Poly(disc,t); print("deg_t discriminant =",Pd.degree(),"; t-adic valuation at t=0:",min(m[0] for m in Pd.monoms()))
# RES check at generic (c1,c2): is s4 a polynomial in t^3 ? (the RES over T=t^3 exists iff yes)
print("s4 polynomial in t^3:",all(m[0]%3==0 for m in P4.monoms()),"; s6 polynomial in t^3:",all(m[0]%3==0 for m in P6.monoms()))
# chi of the elliptic surface over P^1_t with the (t -> T=t^3) base change understood: degrees 4n, 6n
import math
n=max(math.ceil(P4.degree()/4),math.ceil(P6.degree()/6)); print("n = chi =",n,", p_g =",n-1,", 2-forms t^j dt dxi/eta, j=0..",n-2)
w=sp.exp(2*sp.pi*sp.I/3)
for j in range(n-1):
    print(f"  omega_{j} = t^{j} dt dxi/eta  ->  under t->w t: factor w^{j+1}  (character omega^{(j+1)%3})")
