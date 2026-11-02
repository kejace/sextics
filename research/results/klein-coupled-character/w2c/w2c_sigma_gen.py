"""Generate msolve inputs for the locus  Sigma_all = { p : Hess_F(p) restricted to ker(grad F(p)) == 0 }
formulated with auxiliary vector m:  H_ij(p) - F_i(p) m_j - m_i F_j(p) = 0  (i<=j), in affine charts p_j = 1.
Variants: 'X' (the fourfold), 'S' (surface fibre S_Y: x,y,z,lam with Y=(t0,r0,s0) fixed),
'V' (threefold fibre over [r:s]=[alpha:beta]: x,y,z,t,T).
Optional extra equations: F=0 (triple-point locus on the variety) or 1-u*F (points off the variety).
"""
import sympy as sp, sys, os, itertools
import os
PR = int(os.environ.get('W2C_PR','1000003'))
def polys_for(kind, Yvals=(1,2,3), ab=(1,2)):
    if kind=='X':
        x,y,z,t,r,s = V = sp.symbols('x y z t r s')
        F = x**5*y + y**5*z + z**5*x + sp.Rational(5,7)*(x*y**3*t*r + 2*x**3*z*r*s + 3*y*z**3*t*s) + t**6 + r**6 + s**6
    elif kind=='S':
        x,y,z,l = V = sp.symbols('x y z l')
        t0,r0,s0 = Yvals
        F = x**5*y + y**5*z + z**5*x + l**2*sp.Rational(5,7)*(x*y**3*t0*r0 + 2*x**3*z*r0*s0 + 3*y*z**3*t0*s0) + l**6*(t0**6+r0**6+s0**6)
    elif kind=='V':
        x,y,z,t,T = V = sp.symbols('x y z t T')
        a,b = ab
        F = x**5*y + y**5*z + z**5*x + sp.Rational(5,7)*(x*y**3*t*a*T + 2*x**3*z*a*b*T**2 + 3*y*z**3*t*b*T) + t**6 + (a**6+b**6)*T**6
    F = sp.expand(F*7)   # clear the 7
    return V, F
def msolve_file(kind, chart, extra, path, Yvals=(1,2,3), ab=(1,2)):
    V, F = polys_for(kind, Yvals, ab)
    n = len(V)
    m = sp.symbols('m0:%d'%n)
    grad = [sp.diff(F, v) for v in V]
    eqs = []
    for i in range(n):
        for j in range(i, n):
            eqs.append(sp.expand(sp.diff(F, V[i], V[j]) - grad[i]*m[j] - m[i]*grad[j]))
    u = sp.Symbol('u')
    if extra=='on': eqs.append(F)
    elif extra=='off': eqs.append(sp.expand(1 - u*F))
    subs = {V[chart]: 1}
    eqs = [sp.expand(e.subs(subs)) for e in eqs]
    vars_ = [v for k,v in enumerate(V) if k!=chart] + list(m) + ([u] if extra=='off' else [])
    def fmt(e):
        P = sp.Poly(e, *vars_)
        terms = []
        for mon, c in P.terms():
            c = int(c) % PR if PR else int(c)
            if c==0: continue
            mstr = '*'.join('%s^%d'%(str(v),k) if k>1 else str(v) for v,k in zip(vars_, mon) if k>0)
            terms.append(str(c) + ('*'+mstr if mstr else ''))
        return '+'.join(terms) if terms else '0'
    eqs = [e for e in eqs if e!=0]
    with open(path,'w') as f:
        f.write(','.join(str(v) for v in vars_)+'\n%d\n'%PR)
        f.write(',\n'.join(fmt(e) for e in eqs)+'\n')
    return len(vars_), len(eqs)
if __name__=='__main__':
    kind, chart, extra, path = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    print(msolve_file(kind, chart, extra, path))
