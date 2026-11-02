"""Book 2.7 linear-exchange pattern test (any linear form, any complement).
Pattern: F(u + s p) = c s^6 + s B5(u) + C(u) on the hyperplane H = ker(d_p^5 F), i.e. the polars
d_p^k F restricted to H vanish for k = 2,3,4.  Equivalently with l = d_p^5 F (linear form):
  l | d_p^4 F,  l | d_p^3 F,  l | d_p^2 F.
We solve for p in P^5 over Q with sympy Groebner (small system) via the 'remainder' formulation:
  d_p^k F mod l == 0 for k=2,3,4, computed by substituting the variable in which l has the
  nonzero coefficient.  To keep it polynomial we use the resultant-free trick: l divides Q iff
  Q vanishes on {l=0}; parametrize {l=0} generically is messy, so instead we use the condition
  Q(v) * l(w) - Q(w) * l(v) ... too big.  Simpler: brute force structure:
  we classify p by the support pattern using the characteristic 'c s^6' requirement: F(p) != 0.
Here we do the exact algebra with sympy on the 2nd polar only (quartic) which already pins p.
"""
import sympy as sp, itertools, json, time
x,y,z,t,r,s = X = sp.symbols('x y z t r s')
F = x**5*y + y**5*z + z**5*x + sp.Rational(5,7)*(x*y**3*t*r + 2*x**3*z*r*s + 3*y*z**3*t*s) + t**6 + r**6 + s**6
P = sp.symbols('p0:6')
def dpk(k):
    G = F
    for _ in range(k):
        G = sum(P[i]*sp.diff(G, X[i]) for i in range(6))
    return sp.expand(G)
l = dpk(5)            # linear form in X with coefficients quintic in p
d4 = dpk(4); d3 = dpk(3); d2 = dpk(2)
# condition "l divides d_k": for each variable X_j: eliminate using the coefficient of X_j in l.
# Use: l | Q  <=>  Res-type:  Q * (coefficient) ... we use the substitution X_j = -(l - c_j X_j)/c_j
# multiplied through by c_j^deg to stay polynomial; valid where c_j != 0; cover all j.
def div_conditions(Q, j):
    cj = l.coeff(X[j])
    rest = sp.expand(l - cj*X[j])
    deg = sp.Poly(Q, *X).total_degree()
    # Q(X with X_j -> -rest/cj) * cj^deg
    Qs = sp.expand(sp.Poly(Q, X[j]).as_expr().subs(X[j], -rest/cj) * cj**deg)
    Qs = sp.expand(sp.cancel(Qs))
    return [c for c in sp.Poly(Qs, *X).coeffs()]
t0=time.time()
# Work chart by chart: p_j = 1 and use that same j for the division (c_j = coefficient of X_j in l).
results = {}
for j in range(6):
    subs = {P[j]: 1}
    eqs = []
    for Q in (d4, d3, d2):
        for c in div_conditions(Q, j):
            e = sp.expand(c.subs(subs))
            if e != 0: eqs.append(e)
    eqs = list(set(eqs))
    vars_ = [P[i] for i in range(6) if i != j]
    G = sp.groebner(eqs, *vars_, order='grevlex', domain='QQ')
    results[j] = [str(g) for g in G.exprs]
    print('chart', j, 'n_eqs', len(eqs), 'GB', results[j][:12], 'time', round(time.time()-t0,1), flush=True)
json.dump(results, open('/tmp/hodge/results/klein-coupled-character/w2c/exchange_dirs.json','w'), indent=1)
