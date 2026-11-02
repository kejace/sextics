import sympy as sp, json, time
r,w,t,U,beta,a,b,c,d = sp.symbols('r w t U beta a b c d')
A = 1 + r**6 + w**6; B = t*(r + r**3*w); C = r*t*(1 - t**2)/7; Dd = 1 + t**6
R = 3*a**2*beta**4 - 4*a*b*beta**3 + 6*a*c*beta**2 - 12*a*d*beta + 4*b*d - c**2
P = a*U**3 + b*U**2 + c*U + d
out = {}
disc = sp.factor(sp.discriminant(R, beta)); out['disc_beta_R'] = str(disc)
delta = sp.discriminant(P, U); out['disc_U_P'] = str(sp.factor(delta))
M = b - 3*a*beta; N = (c - 3*a*beta**2)/2; k = N - M*beta
out['k_simplified'] = str(sp.factor(k))
out['Pprime_at_minus_beta'] = str(sp.factor(sp.diff(P,U).subs(U,-beta)))
out['N_M'] = str(sp.factor(sp.resultant(R, M, beta)))
out['N_k'] = str(sp.factor(sp.resultant(R, sp.numer(sp.together(k)), beta)))
Cnum = M**2*(3*a*beta**2*k - 2*N**2)
out['N_Cnum'] = str(sp.factor(sp.resultant(R, sp.numer(sp.together(Cnum)), beta)))
sub = {a:A, b:B, c:C, d:Dd}
degs = {}
for name, expr in [("a",A),("b",B),("c",C),("d",Dd),("3ac-b^2",3*a*c-b**2),("delta",delta)]:
    e = sp.expand(expr.subs(sub)); Pl = sp.Poly(e, r, w, t)
    degs[name] = {'total_degree': Pl.total_degree(), 'degs_rwt': list(Pl.degree_list())}
    if name == "3ac-b^2": degs[name]['factored'] = str(sp.factor(e))
out['degrees_in_rwt'] = degs
# irreducibility of delta(r,w,t) over Q
e = sp.expand(delta.subs(sub)); fl = sp.factor_list(e)
out['delta_rwt_factor_degrees'] = [(sp.Poly(f, r, w, t).total_degree(), m) for f, m in fl[1]]
json.dump(out, open('hurwitz_rederive.json','w'), indent=1)
for kk,v in out.items(): print(kk, ':', v if len(str(v))<400 else str(v)[:400]+'...')
