import sys, random, sympy as sp
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave4')
from lines import *
p = int(sys.argv[1]); seed = int(sys.argv[2]); mode = sys.argv[3]; random.seed(seed)
R = lambda: random.randrange(1, 1000)
a2,a3,a4,b2,b3,b4 = sp.symbols('a2 a3 a4 b2 b3 b4'); cc = sp.symbols('c0 c1 c2 c3')
cub = cc[0]*s**3 + cc[1]*s**2*t + cc[2]*s*t**2 + cc[3]*t**3
while True:
    T = sp.Matrix(5, 5, lambda i, j: random.randrange(-5, 6))
    if T.det() != 0: break
w = sp.Matrix([s, t, a2*s+b2*t, a3*s+b3*t, a4*s+b4*t]); lin = list(T*w)
def coeffs(expr):
    P = sp.Poly(sp.expand(expr), s, t); return [sp.expand(P.coeff_monomial(s**k*t**(6-k))) for k in range(7)]
sub = {x0: lin[0], x1: lin[1], y0: lin[2], y1: lin[3], y2: lin[4]}
if mode == 'linesB':        # lines in B = {M=0} ⊂ P^4
    eqs = coeffs(Msex.subs(sub, simultaneous=True)); vars_ = [a2,a3,a4,b2,b3,b4]
elif mode == 'linesBp':     # lines in B' = {L^2-4M=0} (double lines of W)
    eqs = coeffs((Lcub**2 - 4*Msex).subs(sub, simultaneous=True)); vars_ = [a2,a3,a4,b2,b3,b4]
elif mode == 'randb':       # control: sigma_21 count for a RANDOM sextic 3-fold (W_rand = {z^2 = H_rand})
    import itertools
    vv = [x0, x1, y0, y1, y2]
    H = sum(random.randrange(-9, 10)*sp.Mul(*[v**e for v, e in zip(vv, ex)]) for ex in itertools.product(range(7), repeat=5) if sum(ex) == 6)
    u = sp.Symbol('u'); v = sp.symbols('v0 v1 v3 v4 v5')
    P = [R() for _ in range(5)]; Q = [R() for _ in range(5)]
    lin2 = [s*(P[i] + u*Q[i]) + t*v[i] for i in range(5)]
    e = sp.expand((cub**2 - H.subs({x0: lin2[0], x1: lin2[1], y0: lin2[2], y1: lin2[3], y2: lin2[4]}, simultaneous=True)))
    Pz = sp.Poly(e, s, t); eqs = [sp.expand(Pz.coeff_monomial(s**k*t**(6-k))) for k in range(7)]
    h3, h4, h5 = R(), R(), R()
    M = sp.Matrix([[P[0], P[1]], [Q[0], Q[1]]]); rhs = sp.Matrix([-(P[2]*h3+P[3]*h4+P[4]*h5), -(Q[2]*h3+Q[3]*h4+Q[4]*h5)])
    sol = M.solve(rhs); h = [sol[0], sol[1], h3, h4, h5]
    eqs.append(sum(h[i]*v[i] for i in range(5))); eqs.append(sum(R()*vi for vi in v) - 1); eqs.append(sum(R()*vi for vi in v))
    vars_ = [u] + list(v) + list(cc)
elif mode == 'sweepB':      # degree of the divisor of X swept by the cuspidal cubics E_L, L tangent to B (disc c = 0)
    eqs = incidence(lin, cub)
    c0, c1, c2, c3 = cc
    eqs.append(c1**2*c2**2 - 4*c0*c2**3 - 4*c1**3*c3 - 27*c0**2*c3**2 + 18*c0*c1*c2*c3)
    t0, xx = sp.symbols('t0 xx')
    eqs.append(xx**3 - cub.subs({s: 1, t: t0}))
    pt = [lin[0], lin[1], xx, lin[2], lin[3], lin[4]]       # (X0,X1,X2,Y0,Y1,Y2) in original coordinates
    pt = [sp.expand(q.subs({s: 1, t: t0})) for q in pt]
    for _ in range(3):
        eqs.append(sum(random.randrange(-99, 100)*q for q in pt))
    vars_ = [a2,a3,a4,b2,b3,b4] + list(cc) + [t0, xx]
fn = f'ctrl_{mode}_{p}_{seed}'
msolve_file(fn + '.ms', eqs, vars_, p)
singular_file(fn + '.sing', eqs, vars_, p, 'ideal J = std(I); print("DIM"); print(dim(J)); if (dim(J)==0) { print("VDIM"); print(vdim(J)); } else { print("MULT"); print(mult(J)); }\nquit;\n')
print('wrote', fn, 'unknowns', len(vars_), 'eqs', len(eqs))
