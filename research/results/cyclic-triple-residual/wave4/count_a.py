"""a = number of lines of P^4 through a general point that are tritangent to B' (= lines of W through a general point / 2).
System: l(s,t) = s*P + t*v, v=(v0,v1,v3,v4,v5) with two random linear normalisations; c = cubic (4 unknowns)."""
import sys, random, sympy as sp
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave4')
from lines import *
p = int(sys.argv[1]); seed = int(sys.argv[2]); mode = sys.argv[3]
random.seed(seed)
R = lambda: random.randrange(1, 1000)
v = sp.symbols('v0 v1 v3 v4 v5'); cc = sp.symbols('c0 c1 c2 c3')
cub = cc[0]*s**3 + cc[1]*s**2*t + cc[2]*s*t**2 + cc[3]*t**3
if mode == 'a':      # lines through a general point P
    P = [R() for _ in range(5)]
    lin = [s*P[i] + t*v[i] for i in range(5)]
    eqs = incidence(lin, cub)
    eqs.append(sum(R()*vi for vi in v) - 1)
    eqs.append(sum(R()*vi for vi in v))
    vars_ = list(v) + list(cc)
elif mode == 'b':    # lines in a general hyperplane meeting a general line m of it (sigma_{2,1})
    u = sp.Symbol('u')
    P = [R() for _ in range(5)]; Q = [R() for _ in range(5)]
    lin = [s*(P[i] + u*Q[i]) + t*v[i] for i in range(5)]
    eqs = incidence(lin, cub)
    h = [R() for _ in range(5)]
    # hyperplane through m: h(P)=h(Q)=0 -> choose h in the kernel: build random h with h.P = h.Q = 0
    import itertools
    # solve for h: random h3,h4,h5 then h0,h1 from the two equations
    h3, h4, h5 = R(), R(), R()
    M = sp.Matrix([[P[0], P[1]], [Q[0], Q[1]]]); rhs = sp.Matrix([-(P[2]*h3+P[3]*h4+P[4]*h5), -(Q[2]*h3+Q[3]*h4+Q[4]*h5)])
    sol = M.solve(rhs); h = [sol[0], sol[1], h3, h4, h5]
    eqs.append(sum(h[i]*v[i] for i in range(5)))
    eqs.append(sum(R()*vi for vi in v) - 1)
    eqs.append(sum(R()*vi for vi in v))
    vars_ = [u] + list(v) + list(cc)
elif mode == 'ab':   # lines meeting a general line m and a general plane Pi (sigma_2 sigma_1): count = a+b
    u = sp.Symbol('u')
    P = [R() for _ in range(5)]; Q = [R() for _ in range(5)]
    lin = [s*(P[i] + u*Q[i]) + t*v[i] for i in range(5)]
    eqs = incidence(lin, cub)
    # plane Pi = {g1 = g2 = 0}; line meets Pi iff some point s0*(P+uQ)+t0*v on it: g_j(P+uQ) + w*g_j(v) = 0, w = t0/s0
    w = sp.Symbol('w')
    g1 = [R() for _ in range(5)]; g2 = [R() for _ in range(5)]
    eqs.append(sum(g1[i]*(P[i]+u*Q[i]) for i in range(5)) + w*sum(g1[i]*v[i] for i in range(5)))
    eqs.append(sum(g2[i]*(P[i]+u*Q[i]) for i in range(5)) + w*sum(g2[i]*v[i] for i in range(5)))
    eqs.append(sum(R()*vi for vi in v) - 1)
    eqs.append(sum(R()*vi for vi in v))
    vars_ = [u, w] + list(v) + list(cc)
fn = f'count_{mode}_{p}_{seed}'
msolve_file(fn + '.ms', eqs, vars_, p)
singular_file(fn + '.sing', eqs, vars_, p, 'ideal J = std(I); print("DIM"); print(dim(J)); print("VDIM"); print(vdim(J)); quit;\n')
print('wrote', fn, 'unknowns', len(vars_), 'eqs', len(eqs))
