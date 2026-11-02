import sympy as sp
r,w,t,U,beta,a,b,c,d = sp.symbols('r w t U beta a b c d')
A = 1 + r**6 + w**6; B = t*(r + r**3*w); C = r*t*(1 - t**2)/7; Dd = 1 + t**6
R = 3*a**2*beta**4 - 4*a*b*beta**3 + 6*a*c*beta**2 - 12*a*d*beta + 4*b*d - c**2
disc = sp.discriminant(R, beta)
print("disc_beta(R) factored over Q[a,b,c,d]:", sp.factor(disc))
delta = sp.discriminant(a*U**3 + b*U**2 + c*U + d, U)
print("disc_U(P) =", sp.factor(delta))
sub = {a:A, b:B, c:C, d:Dd}
for name, expr in [("a",A),("b",B),("c",C),("d",Dd),("3ac-b^2",3*a*c-b**2),("delta",delta)]:
    e = sp.expand(expr.subs(sub)); P = sp.Poly(e, r, w, t)
    print(f"total degree of {name} in (r,w,t):", P.total_degree(), "; factored:", sp.factor(e) if name in ("3ac-b^2",) else "")
# norm of M = b - 3 a beta over the quartic: resultant_beta(R, b - 3a beta)
normM = sp.resultant(R, b - 3*a*beta, beta)
print("N(M) = Res_beta(R, b-3a beta) factored:", sp.factor(normM))
